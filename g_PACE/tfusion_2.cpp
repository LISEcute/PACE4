#include "tfusion.h"
#include "tfusion_2.h"

#include <QMessageBox>
#include <QSignalMapper>
#include <QDesktopServices>

double _QcLoc2;

extern double _ELAB,_EXPSIGr,_AGRAZr,_ELOSS,_ELAB_MAX;
extern int _LMINN, _JCMAXr, _NOSHL;

extern int _IZP,_IAP,_IZT,_IAT;
extern int _IZC, _IAC;
extern double _SPr, _ST, _QCNr;
//extern char* ElementName(int IZ) ;
extern QString ElementName(int IZ) ;
extern double MASSES(int IZ, int IN, int &opt);
extern double _EEXCN;
extern int _BarrierMode;
extern int _BatchMode;
extern int _BatchNsteps;
extern double _h_omega;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

tfusion_2::tfusion_2(QWidget *parent) :
    QWidget(parent)
{

        QGridLayout *gridLayout2 = new QGridLayout;

        _IAC = _IAP + _IAT;
        _IZC = _IZP + _IZT;

        FlagPermit=false;

        QGroupBox *isoBox[3];
        QString isoBox_names[3] = {"Projectile","Target","Compound"};

        QString isoBox_labels[8] = {"A = ", "N = ", "Z = ", "<b>Ca</b>","Spin (gs) = ", "ME (MeV) = ", "-44.214",
                                    "<font color=\"green\"><b>DB</b></font>"};

        QGridLayout *iBox[3];
        field_labels[0][0] = QString::number(_IAP);
        field_labels[0][2] = QString::number(_IZP);
        field_labels[0][1] = QString::number(_IAP - _IZP);
        field_labels[0][3] = QString::number(_SPr);
        field_labels[1][0] = QString::number(_IAT);
        field_labels[1][2] = QString::number(_IZT);
        field_labels[1][3] = QString::number(_ST);
        field_labels[1][1] = QString::number(_IAT - _IZT);
        field_labels[2][0] = QString::number(_IAC);
        field_labels[2][2] = QString::number(_IZC);
        field_labels[2][1] = QString::number(_IAC - _IZC);

        QSignalMapper *signalMapper = new QSignalMapper(this);

        int NP=_IAP-_IZP;
        int NT=_IAT-_IZT;
        int NC=_IAC-_IZC;
        double MT=0, MP=0, MC=0;

        for(int i=0; i<3; i++)
            {
            isoBox[i] = new QGroupBox(isoBox_names[i]);
            iBox[i]= new QGridLayout;

            for(int j=0; j<4; j++)
                {
                fields[4*i+j] = new QLineEdit(field_labels[i][j]);
                signalMapper->setMapping(fields[4*i+j],4*i+j);
                connect(fields[4*i+j],SIGNAL(textEdited(QString)),signalMapper, SLOT(map()));
                }

            for(int j=0; j<8; j++)
                {
           //     labels[j][i] = new QLabel(isoBox_labels[j],0,0);
                if(j == 3)
                    {
                    QString elc;
                    switch (i) {
                        case 0:
                            elc = ElementName(_IZP);
                            labels[j][i] = new QLabel("<b>" +elc+ "</b>");
                            break;
                        case 1:
                            elc = ElementName(_IZT);
                            labels[j][i] = new QLabel("<b>" +elc+ "</b>");
                            break;
                        case 2:
                            elc = ElementName(_IZC);
                            labels[j][i] = new QLabel("<b>" +elc+ "</b>");
                            break;
                            }
                        }
                else if (j == 6)
                    {
                    int optP, optT, optC;
                    switch(i) {
                        case 0:
                            MP=MASSES(_IZP, NP, optP);
                            labels[j][i] = new QLabel(QString::number(MP));
                            break;
                        case 1:
                            MT=MASSES(_IZT, NT, optT);
                            labels[j][i] = new QLabel(QString::number(MT));
                            break;
                        case 2:
                            MC=MASSES(_IZC, NC, optC);
                            labels[j][i] = new QLabel(QString::number(MC));
                            break;
                        }
                     }
                else {
                     if(! (i == 2 && j == 4))
                            {
                            labels[j][i] = new QLabel(isoBox_labels[j],0);
                            if(j ==7)
                                labels[j][i]->setFixedWidth(22);
                            labels[j][i]->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
                            }
                     else   {
                            labels[j][i] = new QLabel(" ",0);
                            }
                    }


                if (j == 1){
                        fields[4*i+j]->setDisabled(true);
                    }

            }


            iBox[i]->addWidget(labels[0][i],0,0);
            iBox[i]->addWidget(fields[4*i],0,1);
            iBox[i]->addWidget(labels[1][i],0,2);
            iBox[i]->addWidget(fields[1+4*i],0,3);
            iBox[i]->addWidget(labels[2][i],1,0);
            iBox[i]->addWidget(fields[2+4*i],1,1);
            iBox[i]->addWidget(labels[3][i],1,2,1,2,Qt::AlignHCenter);
            iBox[i]->addWidget(labels[4][i],2,0,1,2,Qt::AlignHCenter);

            if(i != 2){
               // fields[3+4*i]->setFixedWidth(70);
                iBox[i]->addWidget(fields[3+4*i],2,2,1,2);
                }
            iBox[i]->addWidget(labels[5][i],3,0,1,2, Qt::AlignHCenter);
            iBox[i]->addWidget(labels[6][i],3,2,1,2,Qt::AlignLeft);
            iBox[i]->addWidget(labels[7][i],3,4, Qt::AlignHCenter);
            if(i == 2)
                {
                for(int j=0;j<4;j++)
                fields[i*4+j]->setDisabled(true);
                }
        }

        connect(signalMapper,SIGNAL(mappedInt(int)),this, SLOT(changedVariables(int)));  // moved here in v.4.33.1

        _QcLoc2=MP+MT-MC;
        if(_QCNr !=0 ) _QcLoc2 = _QCNr;

        double ECM=_ELAB*double(_IAT)/double(_IAC);

        _EEXCN=ECM+_QcLoc2;
        for(int k=0; k<3; k++)
            {
            isoBox[k]->setLayout(iBox[k]);
            gridLayout2->addWidget(isoBox[k],0,k);
            }

        QWidget *qcn_widget = new QWidget;
        QLabel *qcn_label = new QLabel("Q<sub>CN</sub> = ");
        qcn_edit = new QLineEdit(QString::number(_QCNr));
        signalMapper->setMapping(qcn_edit,222);
        connect(qcn_edit,SIGNAL(textEdited(QString)),signalMapper, SLOT(map()));
        QLabel *qcn_desc = new QLabel("Q-value of reaction (MeV) <br> if == 0, calculated from <br> mass tables.");

        QGridLayout *qcn_layout = new QGridLayout;
        qcn_layout->addWidget(qcn_label,0,0);
        qcn_layout->addWidget(qcn_edit,0,1);
        qcn_layout->addWidget(qcn_desc,0,2);
        qcn_widget->setLayout(qcn_layout);
        QGroupBox *calcBox = new QGroupBox("Calculation");
        QString calcLabels[3] = {"Q<sub>CN</sub> = ","E<sub>CM</sub> = ","E<sub>X</sub> = "};
        calc_edit[0]= new QLineEdit(QString::number(_QcLoc2));
        calc_edit[1]= new QLineEdit(QString::number(ECM));
        calc_edit[2]= new QLineEdit(QString::number(_EEXCN));
        QFormLayout *calc_layout = new QFormLayout;
        calc_layout->setLabelAlignment(Qt::AlignRight);

        for (int i=0;i<3;i++)
            {
            QLabel *lab = new QLabel(calcLabels[i]);
            lab->setAlignment(Qt::AlignRight);
            calc_edit[i]->setEnabled(false);
            calc_edit[i]->setAlignment(Qt::AlignTop);
           // calc_layout->addWidget(lab,i,0);
           // calc_layout->addWidget(calc_edit[i],i,1);
            calc_layout->addRow(calcLabels[i],calc_edit[i]);
            }

        calcBox->setLayout(calc_layout);
        QGroupBox *beamE_box  = new QGroupBox("Beam Energy (MeV)  [Lab]");
        beamE_box->setMinimumWidth(210);
        beamE_box->setFixedHeight(100);
        QGridLayout *beamE_layout = new QGridLayout;
        elab_label[0] = new QLabel("E = ", 0);
        elab_label[1] = new QLabel("E <sub>MAX</sub> = ", 0);
        elab_label[2] = new QLabel("N <sub>steps</sub> = ", 0);
        beamE_layout->addWidget(elab_label[0],0,0);
        beamE_layout->addWidget(elab_label[1],1,0);
        beamE_layout->addWidget(elab_label[2],1,2);
        n_steps = new QLineEdit(QString::number(_BatchNsteps));
      //  gridLayout2->addWidget(qcn_widget,1,0,1,2,Qt::AlignLeft);


        beamE_edit = new QLineEdit(QString::number(_ELAB, 'f', 2));
        Elab_max   = new QLineEdit(QString::number(_ELAB_MAX, 'f', 2));

        signalMapper->setMapping(beamE_edit,333);
        signalMapper->setMapping(Elab_max,  333);
        connect(beamE_edit,SIGNAL(textEdited(QString)),signalMapper, SLOT(map()));
        connect(Elab_max,SIGNAL(textEdited(QString)),signalMapper, SLOT(map()));

        beamE_layout->addWidget(beamE_edit,0,1);
        beamE_layout->addWidget(Elab_max,1,1);
        beamE_layout->addWidget(n_steps,1,3);

        elab_label[1]->hide();
        elab_label[2]->hide();
        Elab_max->hide();
        n_steps->hide();
        batch_mode = new QCheckBox("Batch Mode");
        connect(batch_mode, SIGNAL(clicked(bool)), this, SLOT(batch_clicked(bool)));
        beamE_layout->addWidget(batch_mode,0,2,1,2);
        beamE_box->setLayout(beamE_layout);
        gridLayout2->addWidget(calcBox,1,2,2,1);
        gridLayout2->addWidget(beamE_box,1,0,2,1);
        gridLayout2->addWidget(qcn_widget,1,1,2,1);//,Qt::AlignLeft);

        QGridLayout *gl = new QGridLayout;
        QWidget *box = new QWidget;
        QString label_array[10] = {"EXPSIG","JCMAX","AGRAZ","ELOSS","LMINN",
                                   "experimental fusion cross section if known. <br>TL-S from optical model shifted to "
                                   "reproduce this value if inputted, <br>preserving the L-diffuseness. if == 0  Bass "
                                   "model (PRL 1977) fusion  cross section being used.","experimental fusion cross section if known.<br>"
                                   "TL-S from optical model shifted to reproduce this value if inputted,<br> preserving the L-diffuseness.<br> "
                                   "if == 0  Bass model (PRL 1977) fusion  cross section being used.", "To bypass input channel optical model"
                                   "routine (TLOM) specify L-diffuseness of "
                                  "fusion cross section. <br> If == 0 diffuseness will be set to 0.5  which is essentially "
                                  "sharp cutoff. ","energy loss of beam thru full target width. (total dE)"
                                   " energies will be distributed between <br>Ebeam & Ebeam-Elos","Lowest partial wave L in calculation. "
                                   "Partial waves from L=0 to LMINN "
                                  "excluded,<br> enabling low-L non-fusion window in reaction calculation." };

        edits[0] = new QLineEdit(QString::number(_EXPSIGr));
        edits[1] = new QLineEdit(QString::number(_JCMAXr));
        edits[2] = new QLineEdit(QString::number(_AGRAZr));
        edits[3] = new QLineEdit(QString::number(_ELOSS));
        edits[4] = new QLineEdit(QString::number(_LMINN));

        for(int i=0;i<5;i++)
            {
            QLabel *l = new QLabel(label_array[i]);
            l->setStyleSheet("max-width:60px;");
            edits[i]->setStyleSheet("max-width:60px;");
            gl->addWidget(l,i,0);
            gl->addWidget(new QLabel(label_array[i+5]),i,2);
            gl->addWidget(edits[i],i,1);
            }

        box->setLayout(gl);

        //------------------------------------------------------
        gridLayout2->addWidget(box,3,0,1,4);
        QGroupBox *transBox = new QGroupBox("Transmission probability for a one-dimensional barrier (O.T.)");
        radio[0] = new QRadioButton(" Classical (use above the barrier)");
        connect(radio[0],SIGNAL(clicked()),this,SLOT(trans_prob_changed()));
        radio[1] = new QRadioButton(" Quantum-mechanical\n[D.Hill & J.Wheeler, PhysRev 89(1953) 1105]");
        connect(radio[1],SIGNAL(clicked()),this,SLOT(trans_prob_changed()));

        QPushButton *q = new QPushButton("?");
        q->setFixedWidth(25);
        connect(q,SIGNAL(clicked()),this,SLOT(about_trans_clicked()));
        radio[_BarrierMode]->setChecked(true);
        QGridLayout *vbox = new QGridLayout;
        vbox->addWidget(radio[0],0,0);
        vbox->addWidget(radio[1],1,0);
        vbox->addWidget(q,0,1,2,1);
        vbox->setColumnStretch(0,4);
        transBox->setLayout(vbox);
        h_w = new QWidget;
        QGridLayout *h_layout = new QGridLayout;

        QLabel *hw = new QLabel("&#8463;&omega;  = ");
        hw->setTextFormat(Qt::RichText);
        h_layout->addWidget(hw,0,0);
        h_edit = new QLineEdit(QString::number(_h_omega));
        //h_edit->setStyleSheet("max-width:45px;");
        h_layout->addWidget(h_edit,0,1);
        h_layout->addWidget(new QLabel("Curvature parameter of the parabolic potential<br> describing the  barrier (default "
                                       "value 3 MeV) <br>(only for Quantum-Mechanical mode)"),0,2);
        h_w->setLayout(h_layout);
        gridLayout2->addWidget(transBox,4,0,1,2);
        gridLayout2->addWidget(h_w,5,0,1,2);
        gridLayout2->addWidget(new QLabel(" Note: If you are running at high bombarding<br> energies for which the grazing angular <br> momentum "
                                          "is above 75 hbar, it is recomended <br> to input AGRAZ > 0, and to specify an<br> arbitrary value for "
                                          "EXPSIG (or 0 = Bass) <br>which corresponds to afusion cross section  with <br>a limiting  L-value  around "
                                          "80. This will give you <br> all the evaporation residue data and the fission <br>probabilities you need. For "
                                          "J>80 all nuclei will <br> fission anyway and you will run out of <br> dimensions if you try."),
                               4,2,2,1,Qt::AlignTop);

        if(radio[0]->isChecked()) h_w->hide();
        setLayout(gridLayout2);


setVariables();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void tfusion_2::setVariables(bool setAll)
{
    FlagPermit=false;

    _IAC = _IAP + _IAT;
    _IZC = _IZP + _IZT;
    int NP=_IAP-_IZP;
    int NT=_IAT-_IZT;
    int NC=_IAC-_IZC;
    //double MT, MP, MC;
    int optP, optT, optC;
    double MP=MASSES(_IZP, NP, optP);
    double MT=MASSES(_IZT, NT, optT);
    double MC=MASSES(_IZC, NC, optC);

    field_labels[0][0] = QString::number(_IAP);
    field_labels[0][2] = QString::number(_IZP);
    field_labels[0][1] = QString::number(_IAP - _IZP);
    field_labels[0][3] = QString::number(_SPr);
    field_labels[1][0] = QString::number(_IAT);
    field_labels[1][2] = QString::number(_IZT);
    field_labels[1][3] = QString::number(_ST);
    field_labels[1][1] = QString::number(_IAT - _IZT);
    field_labels[2][0] = QString::number(_IAC);
    field_labels[2][2] = QString::number(_IZC);
    field_labels[2][1] = QString::number(_IAC - _IZC);
    field_labels[2][3] = QString::number(0);

    for(int i=0;i<3;i++)
        {
        for(int j=0;j<4;j++){fields[4*i+j]->setText(field_labels[i][j]); }
        int j = 3;

        QString elc;
        switch (i)
            {
            case 0:
                elc = ElementName(_IZP);
                labels[j][i]->setText("<b>" +elc+"</b>");
                break;
            case 1:
                elc = ElementName(_IZT);
                labels[j][i]->setText("<b>" +elc+"</b>");
                break;
            case 2:
                elc = ElementName(_IZC);
                labels[j][i]->setText("<b>" +elc+"</b>");//] = new QLabel("<b>" +QString::fromUtf8(el,2)+"</b>");
                break;
            }

        j = 6;
        int optP, optT, optC;
        switch(i)
            {
            case 0:
                MP=MASSES(_IZP, NP, optP);
                labels[j][i]->setText(QString::number(MP));
                labels[j+1][i]->setText(optP==0 ? "<font color=\"green\"><b>DB</b></font>" : "<font color=\"red\"><b>calc</b></font>");
                break;
            case 1:
                MT=MASSES(_IZT, NT, optT);
                labels[j][i]->setText(QString::number(MT));
                labels[j+1][i]->setText(optT==0 ? "<font color=\"green\"><b>DB</b></font>" : "<font color=\"red\"><b>calc</b></font>");
                break;
            case 2:
                MC=MASSES(_IZC, NC, optC);
                labels[j][i]->setText(QString::number(MC));
                labels[j+1][i]->setText(optC==0 ? "<font color=\"green\"><b>DB</b></font>" : "<font color=\"red\"><b>calc</b></font>");
                break;
            }
        }

if(setAll)
    {
    qcn_edit->setText(QString::number(_QCNr, 'f', 2));
    beamE_edit->setText(QString::number(_ELAB, 'f', 2));

    if(Elab_max->isVisible()){
        Elab_max->setText(QString::number(_ELAB_MAX, 'f', 2));
        n_steps->setText(QString::number(_BatchNsteps));
        }
    }

    edits[0]->setText(QString::number(_EXPSIGr));
    edits[1]->setText(QString::number(_JCMAXr));
    edits[2]->setText(QString::number(_AGRAZr));
    edits[3]->setText(QString::number(_ELOSS));
    edits[4]->setText(QString::number(_LMINN));
    _QcLoc2=MP+MT-MC;

    if(_QCNr !=0 ) _QcLoc2 = _QCNr;

    double ECM=_ELAB*double(_IAT)/double(_IAC);

    _EEXCN=ECM+_QcLoc2;
    calc_edit[0]->setText(QString::number(_QcLoc2));
    calc_edit[1]->setText(QString::number(ECM));
    calc_edit[2]->setText(QString::number(_EEXCN));

    FlagPermit=true;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/** Called before evaluation or page change **/
void tfusion_2::getVariables()
{
    QString val;
    for(int i=0;i<4;i++){
        for(int j=0;j<3;j++){
            val = fields[4*j+i]->text();
            if(i == 0){
                if(j==0) _IAP = val.toInt();
                else if (j==1) _IAT = val.toInt();
                else if (j==3) _IAC = val.toInt();
            } else if (i == 2){
                if(j==0) _IZP = val.toInt();
                else if (j==1) _IZT = val.toInt();
                else if (j==3) _IZC = val.toInt();
            } else if(i == 3) {
                if(j==0) _SPr = val.toDouble();
                else if (j==1) _ST = val.toDouble();
            }
        }
    }
    QString v = qcn_edit->text();
    _QCNr = v.toDouble();
    if(_IAP - _IZP < 0 ) {
        QMessageBox errorMessage;
        errorMessage.critical(0, "Error","Projectile Element number is wrong.");//,MB_OK | MB_ICONERROR);
           //return false;
        }
    if(_IAT - _IZT < 0) {
        QMessageBox errorMessage;
        errorMessage.critical(0,"Error","Target Element number is wrong.");//,MB_OK | MB_ICONERROR);
        // return false;
        }
    if(_IZC >130 || _IZC <=0) {
        QMessageBox errorMessage;
        errorMessage.critical(0,"Error","Compound Element number is wrong.");//,MB_OK | MB_ICONERROR);
        //  return false;
        }

    if(_EEXCN <= 0) {
        QMessageBox errorMessage;
        errorMessage.critical(0,"Error", "Excitation energy is negative.");//,MB_OK | MB_ICONERROR);
        //  return false;
        }

    if(_BatchMode)
    {
        _ELAB_MAX = Elab_max->text().toDouble();
        _BatchNsteps = n_steps->text().toDouble();
        if(_BatchNsteps < 2 || _BatchNsteps > 1000 || _ELAB_MAX <= _ELAB)
        {
            QMessageBox errorMessage;
            errorMessage.critical(0,"Error","Batch mode settings are wrong!");
            //return false;
        }
    }
    for(int i=0;i<5;i++){
        val = edits[i]->text();
        if(i==0) _EXPSIGr = val.toDouble();
        else if(i ==1) _JCMAXr = val.toInt();
        else if (i == 2) _AGRAZr = val.toDouble();
        else if(i == 3) _ELOSS = val.toDouble();
        else if(i == 4) _LMINN = val.toInt();
    }
    if(_BarrierMode == 1) {
        _h_omega = h_edit->text().toDouble();
    }

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2::changedVariables(int index)
{
if(!FlagPermit) return;

    QString val;
    if(index != 222 && index != 333)
        {
        val = fields[index]->text();

        int i = index%4;
        int j = (int)index/4;
        if(i == 0){
                 if(j==0) {_IAP = val.toInt();}
            else if (j==1) _IAT = val.toInt();
            else if (j==3) _IAC = val.toInt();
            setVariables();
            }
        else if (i == 2)
            {
                 if (j==0) _IZP = val.toInt();
            else if (j==1) _IZT = val.toInt();
            else if (j==3) _IZC = val.toInt();
            setVariables();
            }
        else if(i == 3)
            {
            if(j==0) _SPr = val.toDouble();
            else if (j==1) _ST = val.toDouble();
            }
        }
    else if(index == 333)
        {
        _ELAB =  beamE_edit->text().toDouble();
        setVariables(false);
        }
    else if(index == 222)
        {
        _QCNr = qcn_edit->text().toDouble();
        setVariables(false);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2::batch_clicked(bool b){
    if (b)
        {
        _BatchMode = 1;
        elab_label[0]->setText("E <sub>MIN</sub> = ");
        elab_label[1]->show();
        elab_label[2]->show();
        Elab_max->show();
        n_steps->show();
        }
   else {
        _BatchMode = 0;
        elab_label[0]->setText("E = ");
        elab_label[1]->hide();
        elab_label[2]->hide();
        Elab_max->hide();
        n_steps->hide();
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2::trans_prob_changed(){
    if(radio[0]->isChecked()){
        h_w->hide();
        _BarrierMode = 0;
    } else if(radio[1]->isChecked()){
        h_w->show();
        _BarrierMode = 1;
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void tfusion_2::about_trans_clicked()
{
    QDesktopServices::openUrl(QUrl("http://lise.nscl.msu.edu/6_2/lise++_6_2.pdf#page=22"));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
