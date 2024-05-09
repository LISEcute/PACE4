#ifndef TFUSION_2_H
#define TFUSION_2_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QCheckBox>
#include <QRadioButton>
class tfusion_2 : public QWidget
{
    Q_OBJECT
public:
    explicit tfusion_2(QWidget *parent = 0);

    bool FlagPermit;
    QLineEdit *fields[12];
    void setVariables(bool setAll=true);

private:	// User declarations

    QLabel *labels[8][3];
    QString field_labels[3][4];
    QLineEdit *qcn_edit;
    QLineEdit *beamE_edit;
    QLineEdit *Elab_max;
    QLineEdit *E_Nsteps;
    QLineEdit *calc_edit[3];
    QCheckBox *batch_mode;
    QLabel *elab_label[3];
    QLineEdit *n_steps;
    QLineEdit *edits[5];
    QRadioButton *radio[2];
    QLineEdit *h_edit;
    QWidget *h_w;


signals:
    //void changed(int field);

public slots:
    void getVariables();
    void changedVariables(int);
    void batch_clicked(bool);
    void trans_prob_changed();
    void about_trans_clicked();
};

#endif // TFUSION_2_H
