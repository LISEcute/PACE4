#include "p_about.h"
#include "ui_p_about.h"
#include "g_PACE/ftype.h"

#include <QDesktopServices>
#include <QUrl>

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
About::About(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::About)
{
    ui->setupUi(this);
    ui->label_Version->setText(PACE_version);
    ui->label_Date->setText(PACE_date);

    connect(ui->label_LISE, SIGNAL(clicked()), this, SLOT(CmLISE()));
    connect(ui->label_Charge, SIGNAL(clicked()), this, SLOT(CmCharge()));
    connect(ui->label_mail, SIGNAL(clicked()), this, SLOT(CmMail()));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
About::~About() {    delete ui;}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::CmLISE() {    QDesktopServices::openUrl(QUrl(ui->label_LISE->text()));}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::CmCharge() {    QDesktopServices::openUrl(QUrl(ui->label_Charge->text()));}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::CmMail()
{
QString ss("mailto:tarasov@frib.msu.edu?subject=PACE4 ");
ss.append(ui->label_Version->text());
QDesktopServices::openUrl(QUrl(ss));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
