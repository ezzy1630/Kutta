// Design studies: sweep values, which settings may vary, variant bookkeeping,
// storage, and a study run back to back against the fixture CLI.

#include "core/Session.h"
#include "core/Study.h"

#include <QSignalSpy>
#include <QTemporaryDir>
#include <QtTest>

using namespace kapp;

class TestStudy : public QObject {
    Q_OBJECT
    QTemporaryDir m_dir;
    QString m_car;
    const QString m_inlet = QStringLiteral("boundary/inlet/velocity-magnitude");

    std::unique_ptr<Session> open(const QString& data)
    {
        CliInfo cli;
        cli.program = QStringLiteral(KUTTA_FIXTURE_CLI);
        cli.fixture = true;
        auto s = std::make_unique<Session>(cli, m_dir.filePath(data));
        QSignalSpy opened(s.get(), &Session::caseOpened);
        s->openCase(m_car);
        if (!opened.wait(20000)) return nullptr;
        for (const QString& a : s->caseSummary().approximations()) s->setApproximationApproved(a, true);
        return s;
    }

private slots:
    void initTestCase()
    {
        QCoreApplication::setOrganizationName("KuttaTest");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, m_dir.filePath("settings"));
        qputenv("KUTTA_FIXTURE_SPEED", "400");
        m_car = m_dir.filePath("fsae_car.cas.h5");
        QFile f(m_car);
        QVERIFY(f.open(QIODevice::WriteOnly));
    }

    void sweepValuesParse()
    {
        QString why;
        QCOMPARE(parseSweepValues("16, 20, 24", &why), (std::vector<double>{16, 20, 24}));
        QCOMPARE(parseSweepValues("12 to 28 step 4", &why), (std::vector<double>{12, 16, 20, 24, 28}));
        QCOMPARE(parseSweepValues("0.1 to 0.3 step 0.1", &why).size(), std::size_t(3));  // float steps keep the end
        QVERIFY(parseSweepValues("16", &why).empty() && why.contains("two"));
        QVERIFY(parseSweepValues("16, 16", &why).empty() && why.contains("twice"));
        QVERIFY(parseSweepValues("16, fast", &why).empty() && why.contains("fast"));
        QVERIFY(parseSweepValues("28 to 12 step 4", &why).empty());
        QVERIFY(parseSweepValues("1 to 100 step 1", &why).empty() && why.contains("50"));
        QVERIFY(parseSweepValues("1 to 2 step 0", &why).empty() && why.contains("positive"));
    }

    void bookkeepingFollowsRuns()
    {
        Study s;
        s.id = "study-0001";
        s.parameter = m_inlet;
        s.base = {{m_inlet, 99.0}, {"reference/velocity", 99.0}, {"reference/area", 1.1}};
        s.followers = {"reference/velocity"};
        s.variants = {{12, {{m_inlet, 12.0}, {"reference/velocity", 12.0}}, {}, VariantStatus::Queued},
                      {16, {}, {}, VariantStatus::Queued},
                      {20, {{m_inlet, 20.0}}, {}, VariantStatus::Queued}};
        // The varied setting and its followers come from the variant, never the base.
        QCOMPARE(s.overridesFor(s.variants[1]), (QJsonObject{{"reference/area", 1.1}}));
        QCOMPARE(s.overridesFor(s.variants[0]).value(m_inlet).toDouble(), 12.0);

        QCOMPARE(s.nextQueued(), 0);
        markStarted(s, 0, "run-0004");
        QCOMPARE(s.variants[0].status, VariantStatus::Running);
        QCOMPARE(s.nextQueued(), 1);
        markFinished(s, "run-0004", VariantStatus::Done);
        markStarted(s, 1, "run-0005");
        markFinished(s, "run-0005", VariantStatus::Failed);
        cancelQueued(s);
        QCOMPARE(s.variants[2].status, VariantStatus::Cancelled);
        QVERIFY(s.finished());
        // Running again repeats only what did not complete.
        requeue(s);
        QCOMPARE(s.variants[0].status, VariantStatus::Done);
        QCOMPARE(s.variants[1].status, VariantStatus::Queued);
        QVERIFY(s.variants[1].runId.isEmpty());
        QCOMPARE(s.variants[2].status, VariantStatus::Queued);

        // A study saved while a variant ran comes back stopped, not running.
        markStarted(s, 1, "run-0006");
        const auto back = Study::fromJson(s.toJson());
        QVERIFY(back);
        QCOMPARE(back->variants[1].status, VariantStatus::Stopped);
        QCOMPARE(back->variants[1].runId, QString("run-0006"));
        QCOMPARE(back->base, s.base);
        QCOMPARE(back->followers, s.followers);

        StudyStore store;
        store.setDir(m_dir.filePath("store"));
        QCOMPARE(store.nextNumber(), 1);
        s.number = 1;
        QVERIFY(store.save(s));
        QCOMPARE(store.nextNumber(), 2);
        QCOMPARE(int(store.load().size()), 1);
        QVERIFY(store.remove(s.id));
        QVERIFY(store.load().empty());
    }

    // Only editable numbers that change the case may vary; a sweep of the
    // inlet speed carries the reference velocity when the case ties them.
    void sweepIsHonestAboutWhatCanVary()
    {
        auto s = open("data-sweep");
        QVERIFY(s);
        const SettingsModel& m = *s->settings();
        QVERIFY(studyParameterAllowed(*m.find(m_inlet)));
        QVERIFY(!studyParameterAllowed(*m.find("run/iterations")));
        QVERIFY(!studyParameterAllowed(*m.find("monitors/report/drag/stop-criterion")));
        for (const Setting& x : m.settings())
            if (!x.editable) QVERIFY2(!studyParameterAllowed(x), qPrintable(x.path));

        const double caseSpeed = m.find(m_inlet)->value.toDouble();
        QCOMPARE(referenceFollowers(m, m_inlet), QStringList{"reference/velocity"});
        const auto study = makeSweep(m, m_inlet, {12, caseSpeed}, referenceFollowers(m, m_inlet));
        QVERIFY2(study.ok(), qPrintable(study ? QString() : study.error()));
        QCOMPARE(study.value().variants[0].overrides, (QJsonObject{{m_inlet, 12.0}, {"reference/velocity", 12.0}}));
        // The case's own value changes nothing, so that run is not MODIFIED.
        QVERIFY(study.value().variants[1].overrides.isEmpty());
        // Values the inventory rejects are refused up front.
        QVERIFY(!makeSweep(m, m_inlet, {12, -5}, {}).ok());
        QVERIFY(!makeSweep(m, "run/iterations", {100, 200}, {}).ok());
    }

    // A study runs its variants one after another; each is an ordinary run
    // labelled with its value, and only changed values are MODIFIED.
    void studyRunsBackToBack()
    {
        auto s = open("data-run");
        QVERIFY(s);
        SettingsModel& m = *s->settings();
        QCOMPARE(m.setValue("run/iterations", 200), QString());
        const double caseSpeed = m.find(m_inlet)->value.toDouble();
        auto study = makeSweep(m, m_inlet, {caseSpeed, 24}, referenceFollowers(m, m_inlet));
        QVERIFY(study.ok());
        study.value().name = "Speed sweep";
        QString why;
        const QString id = s->addStudy(study.value(), &why);
        QVERIFY2(!id.isEmpty(), qPrintable(why));
        QVERIFY2(s->startStudy(id, &why), qPrintable(why));
        QVERIFY(!s->startStudy(id, &why));  // one at a time
        QTRY_VERIFY_WITH_TIMEOUT(s->activeStudy().isEmpty() && !s->running(), 120000);

        const Study* done = s->study(id);
        QVERIFY(done);
        for (const StudyVariant& v : done->variants) QCOMPARE(v.status, VariantStatus::Done);
        const RunRecord* a = s->history().find(done->variants[0].runId);
        const RunRecord* b = s->history().find(done->variants[1].runId);
        QVERIFY(a && b);
        QCOMPARE(a->studyId, id);
        QCOMPARE(a->iterationsRequested, 200);
        QVERIFY(!a->modified());
        QVERIFY(b->modified());
        QCOMPARE(b->note, QString("Speed sweep · 24"));

        // Coefficients stay referenced to the speed run: the faster variant
        // has more downforce but a similar CL.
        std::optional<Result<AeroResult>> ra, rb;
        s->aero(a->id, [&](Result<AeroResult> r) { ra = std::move(r); });
        s->aero(b->id, [&](Result<AeroResult> r) { rb = std::move(r); });
        QTRY_VERIFY_WITH_TIMEOUT(ra && rb, 20000);
        QVERIFY(ra->ok() && rb->ok());
        const AeroRow* ta = ra->value().total();
        const AeroRow* tb = rb->value().total();
        QVERIFY(tb->downforce > 1.5 * ta->downforce);
        QVERIFY(qAbs(tb->coefficients->cl / ta->coefficients->cl - 1) < 0.05);

        // The study is kept with the case and reloads.
        auto again = open("data-run");
        QVERIFY(again && again->study(id));
        QCOMPARE(again->study(id)->variants[1].runId, b->id);
    }

    void stopCancelsTheRest()
    {
        auto s = open("data-stop");
        QVERIFY(s);
        QCOMPARE(s->settings()->setValue("run/iterations", 2000), QString());
        auto study = makeSweep(*s->settings(), m_inlet, {14, 18, 22}, {});
        QVERIFY(study.ok());
        const QString id = s->addStudy(study.value());
        QVERIFY(s->startStudy(id));
        QTRY_VERIFY_WITH_TIMEOUT(s->live().last && s->live().last->iter >= 20, 20000);
        s->stopStudy();
        QTRY_VERIFY_WITH_TIMEOUT(s->activeStudy().isEmpty() && !s->running(), 20000);
        const Study* st = s->study(id);
        QCOMPARE(st->variants[0].status, VariantStatus::Stopped);
        QCOMPARE(st->variants[1].status, VariantStatus::Cancelled);
        QCOMPARE(st->variants[2].status, VariantStatus::Cancelled);
    }
};

QTEST_MAIN(TestStudy)
#include "tst_study.moc"
