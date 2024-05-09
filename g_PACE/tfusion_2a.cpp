#include "tfusion_2a.h"
#include "tfusion.h"

#include <QGridLayout>
#include <QGroupBox>
#include <QSignalMapper>
#include <QMessageBox>
#include <QSpacerItem>
#include <QPainter>

extern int _IZC, _IAC, _LMINN, _INPUT;
extern double _EEXCN, _EREC, _AJNUC,_EEXCN_MAX;
//extern char* ElementName(int IZ) ;
extern QString ElementName(int IZ);
extern int _BatchMode;
extern int _BatchNsteps;

extern double MASSES(int IZ, int IN, int &opt);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

tfusion_2a::tfusion_2a(QWidget *parent) :
    QWidget(parent)
{

    QString isoBox_labels[8] = {"A = ", "N = ", "Z = ", "<b>Ca</b>","Spin (gs) = ", "ME (MeV) = ", "-44.214", "<font color=\"green\"><b>DB</b></font>"};

    int NC=_IAC-_IZC;
    double MC;
    QGridLayout *gridLayout2a = new QGridLayout;
    gridLayout2a->setSpacing(25);
//    gridLayout2a->setContentsMargins(50,50,50,50);
    QSpacerItem *vspace = new QSpacerItem(100,1500);
    gridLayout2a->addItem(vspace,0,0);

    QGroupBox *isoBox;
    QGridLayout *iBox;
    field_labels[0] = QString::number(_IAC);
    field_labels[2] = QString::number(_IZC);
    field_labels[1] = QString::number(_IAC - _IZC);

    QSignalMapper *signalMapper = new QSignalMapper(this);
    connect(signalMapper,SIGNAL(mappedInt(int)),this, SLOT(setVariablesa(int)));
    isoBox = new QGroupBox("Compound");
    iBox= new QGridLayout;

    for(int j=0;j<4;j++)
        {
        fields[j] = new QLineEdit(field_labels[j]);
        signalMapper->setMapping(fields[j],j);
        connect(fields[j],SIGNAL(textEdited(QString)),signalMapper, SLOT(map()));
        }

    for(int j=0;j<8;j++)
        {
        if(j == 3){
                QString elc = ElementName(_IZC);
                labels[j] = new QLabel("<b>" +elc+"</b>");
                break;
                }
        else if (j == 6){
                int optC;
                MC=MASSES(_IZC, NC, optC);
                labels[j] = new QLabel(QString::number(MC));
                break;
                }
        else    {
            if(!(j == 4)){
                labels[j] = new QLabel(isoBox_labels[j],0);
                }
            else {
                labels[j] = new QLabel(" ",0);
                }
            }
        }
    fields[1]->setDisabled(true);

    iBox->addWidget(labels[0],0,0);
    iBox->addWidget(fields[0],0,1);
    iBox->addWidget(labels[1],0,2);
    iBox->addWidget(fields[1],0,3);
    iBox->addWidget(labels[2],1,0);
    iBox->addWidget(fields[2],1,1);
    iBox->addWidget(labels[3],1,2,1,2,Qt::AlignHCenter);

    isoBox->setLayout(iBox);
    isoBox->setMaximumSize(200,200);

    gridLayout2a->addWidget(isoBox,1,0,1,3, Qt::AlignHCenter);

    QGroupBox *beamE_box  = new QGroupBox("Excitation energy of compound nucleus (0 < Ex < 2000 MeV)");
    QGridLayout *beamE_layout = new QGridLayout;
    elab_label[0] = new QLabel("EEXCN = ", 0);
    elab_label[1] = new QLabel("EEXCN_max = ", 0);
    elab_label[2] = new QLabel("N of Steps = ", 0);
    beamE_layout->addWidget(elab_label[0],0,0);
    beamE_layout->addWidget(elab_label[1],1,0);
    beamE_layout->addWidget(elab_label[2],1,2);
    n_steps = new QLineEdit(QString::number(_BatchNsteps));

    beamE_edit = new QLineEdit(QString::number(_EEXCN));
    Elab_max = new QLineEdit(QString::number(_EEXCN_MAX));

    beamE_layout->addWidget(beamE_edit,0,1);
    beamE_layout->addWidget(Elab_max,1,1);
    beamE_layout->addWidget(n_steps,1,3);
    elab_label[1]->hide();
    elab_label[2]->hide();
    Elab_max->hide();
    n_steps->hide();
    batch_mode = new QCheckBox("Batch Mode");

    connect(batch_mode, SIGNAL(clicked(bool)), this, SLOT(batch_clicked(bool)));

    beamE_layout->addWidget(batch_mode,0,2);
    beamE_box->setLayout(beamE_layout);
    beamE_box->setMaximumSize(400,200);
    gridLayout2a->addWidget(beamE_box,2,0,1,3,Qt::AlignHCenter);
    QWidget *sec = new QWidget;
    QGridLayout *secLayout = new QGridLayout;
    QLabel *l[3][2];
    QString l_vars[3] = {"EREC = ", "AJNUC = ", "LMINN = "};
    QString l_desc[3] = {"recoil energy of compound nucleus","maximum spin of compound nucleus","Lowest partial wave L in calculation."
                          "Partial waves from L=0 to LMINN <br> excluded,"
                          "enabling low-L non-fusion window in reaction<br>"
                          "calculation. (AJNUC > LMINN)"};

    for(int i=0;i<3;i++)
        {
        l[i][0] = new QLabel(l_vars[i]);        secLayout->addWidget(l[i][0],i,0);
        line[i] = new QLineEdit;                secLayout->addWidget(line[i],i,1);
        l[i][1] = new QLabel(l_desc[i]);        secLayout->addWidget(l[i][1],i,2);
        }

    setVariables();

    sec->setLayout(secLayout);
    sec->setFixedHeight(160);
    gridLayout2a->addWidget(sec,3,0,1,3,Qt::AlignHCenter);
    gridLayout2a->setAlignment(this, Qt::AlignCenter);
    this->setFixedHeight(480);

   // gridLayout2a->setMargin(50);

    setLayout(gridLayout2a);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2a::setVariablesa(int index)
{
    QString val = fields[index]->text();
    int i = index;
    if(i == 0){
        _IAC = val.toInt();
    } else if (i == 2){
        _IZC = val.toInt();
    }
    setVariables();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2a::setVariables(){

    field_labels[0] = QString::number(_IAC);
    field_labels[2] = QString::number(_IZC);
    field_labels[1] = QString::number(_IAC - _IZC);

        for(int j=0;j<4;j++){
            fields[j]->setText(field_labels[j]);
        }
        int j = 3;
        //char* elc;
        QString elc;
        //char *el = new char[2];
        elc = ElementName(_IZC);
       // el[0] = elc[0];
       // el[1] = elc[1];
        labels[j]->setText("<b>" +elc+"</b>");//] = new QLabel("<b>" +QString::fromUtf8(el,2)+"</b>");
        line[0]->setText(QString::number(_EREC));
        line[1]->setText(QString::number(_AJNUC));
        line[2]->setText(QString::number(_LMINN));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2a::getVariables()
{
    QString val;
    for(int i=0;i<4;i++){
        val = fields[i]->text();
        if(i == 0){
            _IAC = val.toInt();
        } else if (i == 2){
            _IZC = val.toInt();
        }
    }

    if(_IZC >130 || _IZC <=0) {
        QMessageBox errorMessage;
        errorMessage.critical(0,"Error","Compound Element number is wrong.");//,MB_OK | MB_ICONERROR);
         //  return false;
    }
    _EEXCN = beamE_edit->text().toDouble();
    if(Elab_max->isVisible()){
        _EEXCN_MAX = Elab_max->text().toDouble();
        _BatchNsteps = n_steps->text().toDouble();
    }
    for(int i=0;i<3;i++){
        val = line[i]->text();
        if(i == 0) _EREC = val.toDouble();
        else if(i == 1) _AJNUC = val.toDouble();
        else if(i == 2) _LMINN = val.toDouble();
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2a::batch_clicked(bool b){
    if (b == true){
        _BatchMode = 1;
        elab_label[0]->setText("EEXCN_min = ");
        elab_label[1]->show();
        elab_label[2]->show();
        Elab_max->show();
        n_steps->show();
    } else {
        _BatchMode = 0;
        elab_label[0]->setText("EEXCN = ");
        elab_label[1]->hide();
        elab_label[2]->hide();
        Elab_max->hide();
        n_steps->hide();
    }
//getVariables();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

