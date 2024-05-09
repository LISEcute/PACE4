#ifndef TFUSION3_H
#define TFUSION3_H

#include <QWidget>
#include <QRadioButton>
#include <QLineEdit>
#include <QComboBox>
#include <QLabel>
#include <QGridLayout>
#include <QStackedWidget>
#include "pace.h"
class tfusion3 : public QWidget
{
    Q_OBJECT
public:
    explicit tfusion3(PACE *p, QWidget *parent = 0);
    int CurrentPage;
    void setPage(int, PACE *p);
    void setVariables(PACE *p);
    void readPage(int, PACE *p);
    QStackedWidget *pages;
    QWidget *w[4];
   // PACE* pace; // = &p;
signals:

public slots:
    void opt_changed();
private:
    QGridLayout *layout_3[4];
    QRadioButton *opt_sys[4];
    QRadioButton *opt_man[4];
    QLineEdit *edits[15][4];
   // void setPage(int);
    QComboBox *imag_box[4];
    void ShowHideEdit(int, int page);
    QLabel *name[4];
    PACE *pace;
};

#endif // TFUSION3_H
