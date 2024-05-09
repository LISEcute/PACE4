#ifndef TFUSION_2A_H
#define TFUSION_2A_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QCheckBox>

class tfusion_2a : public QWidget
{
    Q_OBJECT
public:
    explicit tfusion_2a(QWidget *parent = 0);
    void setVariables();
signals:

public slots:
    void getVariables();
    void setVariablesa(int);
    void batch_clicked(bool);

private:
    QString field_labels[4];
    QLabel *labels[8];
    QLineEdit *fields[12];
    QLabel *elab_label[3];
    QLineEdit *n_steps;
    QLineEdit *beamE_edit;
    QLineEdit *Elab_max;
    QLineEdit *E_Nsteps;
    QCheckBox *batch_mode;
    QLineEdit *line[3];

};

#endif // TFUSION_2A_H
