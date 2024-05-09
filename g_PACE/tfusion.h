#ifndef TFUSION_H
#define TFUSION_H

#include <QString>
#include <QLineEdit>
#include <QMainWindow>
#include <QFile>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QStackedWidget>
#include <QSqlDatabase>

#include "tfusion_2.h"
#include "tfusion_2a.h"
#include "tfusion_2a2.h"
#include "tfusion3.h"
#include "particles.h"
#include "results.h"
#include "pace.h"

#ifndef ParticleStatesH
  #include "ParticleStates.h"
#endif



//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
namespace Ui {
class tfusion;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class tfusion : public QMainWindow
{
    Q_OBJECT

public:
    explicit tfusion(QWidget *parent = 0);
    ~tfusion();

    QString FFileName;
    QString FFileNameOut;
    QString FFileNameEve;
    QString FFileNameParticles;
    QString FFileNameCS;
    QString FFileNameHtml;
    QString Buffer;

    particles *part_page;
    results *res_page;
    PACE *p;


protected:
    void keyPressEvent(QKeyEvent *e);


    //bool db_found;
    //bool db_open;

private:
    bool   Modified;
    QLineEdit* lineEdit[5];
    QLineEdit* limits[2];
    QStackedWidget *stackedw;
    tfusion_2 *page_2;
    tfusion_2a *page_2a;
    tfusion_2a2 *page_2a2;
    tfusion3 *page3;
   // QString html_filename;
   // QRadioButton *inp[5];
    QComboBox *inp_combo;
    QRadioButton *idist[3];
    QRadioButton *mdir[2];
    QRadioButton *itrac[2];
    QRadioButton *noshl[2];
    QCheckBox *partOutput[5];
    QCheckBox *use_gate;
    QLineEdit *gateAZ[2];
    QPushButton *q;

    QLabel *Zl;
    QLabel *Al;

    void setFileName(const QString& FileName);
    void RunPACE3();
    void RunPACE3local();
    void format_error();
    void setVariables();

    bool init_base();
    void readFile(const QString& file);

private slots:

    void getVariables();

    void createFile_clicked();
    void idist_clicked();
    void use_gate_clicked();
    void q_clicked();

    void on_actionExecute_triggered();
    void on_actionNext_triggered();
    void on_action_Previous_triggered();
    void on_actionWeb_Documentation_triggered();
    void on_actionOpen_triggered();
    void on_actionSave_triggered();
    void on_actionExit_triggered();
//    void on_actionPrint_Results_File_triggered();
    void on_actionSave_As_triggered();
    void on_actionAbout_triggered();

private:
    Ui::tfusion *ui;
};

#endif // TFUSION_H
