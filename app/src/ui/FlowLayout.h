#pragma once

// A row of buttons or badges that wraps onto more lines when the panel is too
// narrow for it, rather than squeezing them below their size (the inspector
// is as narrow as 320 px). Items keep their size hints; stretches are ignored.

#include <QLayout>

namespace kapp {

class FlowLayout : public QLayout {
public:
    explicit FlowLayout(QWidget* parent = nullptr, int hSpacing = 6, int vSpacing = 6);
    ~FlowLayout() override;

    void addItem(QLayoutItem* item) override;
    int count() const override { return int(m_items.size()); }
    QLayoutItem* itemAt(int index) const override { return m_items.value(index); }
    QLayoutItem* takeAt(int index) override;
    Qt::Orientations expandingDirections() const override { return {}; }
    bool hasHeightForWidth() const override { return true; }
    int heightForWidth(int width) const override { return arrange(QRect(0, 0, width, 0), false); }
    QSize minimumSize() const override;
    QSize sizeHint() const override;
    void setGeometry(const QRect& rect) override;

private:
    int arrange(const QRect& rect, bool apply) const;
    QList<QLayoutItem*> m_items;
    int m_h, m_v;
};

} // namespace kapp
