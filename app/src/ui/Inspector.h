#pragma once

// The task page for the selected outline node, like Fluent's task panel.
// Setup and solution pages edit the case's settings (as overrides); results
// pages edit the active viewport and request post-processing.

#include "core/Session.h"

#include <QPointer>
#include <QScrollArea>

class QLabel;

namespace kapp {

class ChartsPanel;
class GraphicsArea;

struct InspectorContext {
    Session* session;
    GraphicsArea* graphics;
    ChartsPanel* charts;
    std::function<void(const QString& id)> navigate;
    std::function<void(const QString& filter)> showSettings;  // All Case Settings, filtered
    std::function<void()> openDialog;                          // File ▸ Open Case…
    std::function<void(const QString& casePath)> openPath;
    std::function<void()> newCaseDialog;                       // File ▸ New Case from Mesh…
    std::function<void()> createCaseDialog;                    // where to write the new case
    std::function<void()> importGeometryDialog;                // File ▸ New Case from Geometry…
    std::function<void()> addGeometryDialog;                   // more STL files for the project
    std::function<void(int fileIndex)> replaceGeometryDialog;  // swap one STL file
    std::function<void()> sampleGeometry;                      // the fixture CLI's synthetic car as STL
    std::function<void()> duplicateVariant;                    // Duplicate as Design Variant
    QStringList samples;                                       // sample cases the CLI can serve
    std::function<void()> refresh;                             // rebuild the current page
    std::function<QString()> current;                          // the current page's id
};

class Inspector : public QScrollArea {
    Q_OBJECT
public:
    Inspector(Session* session, GraphicsArea* graphics, ChartsPanel* charts, QWidget* parent = nullptr);
    InspectorContext& context() { return m_ctx; }
    void show(const QString& id);
    QString current() const { return m_id; }
    void showSettingsFiltered(const QString& filter);
    void refresh() { show(m_id); }
    // The current page laid out at each page width in turn, as the panel would
    // lay it out, and its ScreenChecks layout problems; the page is then put
    // back at the panel's width. For scripted sessions and their screenshots.
    QStringList layoutProblemsAtWidths(const QList<int>& widths);

signals:
    void navigate(const QString& id);

protected:
    void resizeEvent(QResizeEvent* e) override;

private:
    void fitPageWidth();
    QWidget* build(const QString& id);
    InspectorContext m_ctx;
    QString m_id;
    QString m_settingsFilter;
    bool m_meshBusy = false;
};

// Page builders (Inspector.cpp: setup and solution; ResultsPages.cpp: the rest).
namespace pages {
QWidget* header(const QString& title, const QString& crumb, QWidget* parent);
QWidget* general(const InspectorContext& c);
QWidget* models(const InspectorContext& c);
QWidget* materials(const InspectorContext& c, const QString& only);
QWidget* cellZones(const InspectorContext& c);
QWidget* boundaries(const InspectorContext& c);
QWidget* boundary(const InspectorContext& c, const QString& zone);
QWidget* reference(const InspectorContext& c);
QWidget* methods(const InspectorContext& c);
QWidget* controls(const InspectorContext& c);
QWidget* reportDefinitions(const InspectorContext& c);
QWidget* monitors(const InspectorContext& c);
QWidget* initialization(const InspectorContext& c);
QWidget* run(const InspectorContext& c);
QWidget* contours(const InspectorContext& c);
QWidget* vectors(const InspectorContext& c);
QWidget* pathlines(const InspectorContext& c);
QWidget* isoSurfaces(const InspectorContext& c);
QWidget* planes(const InspectorContext& c);
QWidget* mesh(const InspectorContext& c);
QWidget* xyPlot(const InspectorContext& c);
QWidget* residuals(const InspectorContext& c);
QWidget* reportPlots(const InspectorContext& c);
QWidget* forces(const InspectorContext& c);
QWidget* integrals(const InspectorContext& c);
QWidget* aero(const InspectorContext& c);
QWidget* runHistory(const InspectorContext& c);
QWidget* studies(const InspectorContext& c);
QWidget* compare(const InspectorContext& c);
QWidget* allSettings(const InspectorContext& c, const QString& filter);
QWidget* noCase(const InspectorContext& c, const QString& title);
// Meshing: New Case from Geometry (MeshPages.cpp).
QWidget* geometry(const InspectorContext& c);
QWidget* tunnel(const InspectorContext& c);
QWidget* refinement(const InspectorContext& c);
QWidget* sizing(const InspectorContext& c);
QWidget* generateMesh(const InspectorContext& c);
QWidget* meshReport(const InspectorContext& c);
QWidget* meshSensitivity(const InspectorContext& c);
QWidget* variants(const InspectorContext& c);
// File ▸ New Case from Mesh (NewCasePage.cpp).
QWidget* newCase(const InspectorContext& c);
// Where the live or last run stands against its convergence criteria.
QWidget* convergenceStatus(const InspectorContext& c);
} // namespace pages

} // namespace kapp
