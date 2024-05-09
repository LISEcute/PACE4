#include "tfusion3.h"
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>

#include "pace.h"
extern  int _INPUT;
extern double _V0D[5],_R0RD[5],_ARD[5],_R0CD[5],_RCCC[4],_VQS[4],_VQ0[4],_W0D[5],_R0ID[5],_AID[5],
       _RMCHD[5],_V01D[5],_V02D[5],_W01D[5],_W02D[5];
extern int _NPD[5],_IMAG[5],_IRAD[5];//, _ITAKE[5];
extern int _IZC, _IAC;

//extern void OPTPOT(int I,int IN, int IZ);


const char *NameForm[4]={"<h3 align=center>neutron</h3>","<h3 align=center>proton</h3>","<h3 align=center>alpha</h3>","<h3 align=center>incoming channel</h3>"};
char NameLabelCard[8];


tfusion3::tfusion3(PACE *p, QWidget *parent) :
    QWidget(parent)
{
    pace = p;
    pages = new QStackedWidget;

    QHBoxLayout *lay = new QHBoxLayout;

    QGroupBox *optBox[4];
    QVBoxLayout *opt_layout[4];
    QLabel *vars[15][4];
    QString var_names[15] = {"V0D(I)","V01D(I)","V02D(I)","R0RD(I)","ARD(I)","R0CD(I)","W0D(I)","W1D(I)","W2D(I)","R0ID(I)","AID(I)","RMCHS(I)",
                             "NPD(I)","IMAG(I)","IRAD(I)"};
    QLabel *text[4];
    for(int i=0;i<4;i++){
        p->_ITAKE[i+1] = false;
        w[i] = new QWidget;
        layout_3[i]= new QGridLayout;

        name[i] = new QLabel(NameForm[i]);
        name[i]->setStyleSheet("background-color: #0f67e0; border-radius: 10px; margin: 30px; color: #fff;");
        layout_3[i]->addWidget(name[i],0,0,1,2);
        optBox[i]= new QGroupBox("Optical model potential");
        opt_layout[i] = new QVBoxLayout;
        opt_sys[i] = new QRadioButton("from systematics");
        connect(opt_sys[i],SIGNAL(clicked()),this,SLOT(opt_changed()));
        opt_layout[i]->addWidget(opt_sys[i]);
        opt_man[i] = new QRadioButton("manual");
        connect(opt_man[i],SIGNAL(clicked()),this,SLOT(opt_changed()));
        opt_layout[i]->addWidget(opt_man[i]);
        optBox[i]->setLayout(opt_layout[i]);
        optBox[i]->setFixedHeight(100);
        layout_3[i]->addWidget(optBox[i],0,2,1,2);
        for(int j=0;j<15;j++){
            vars[j][i] = new QLabel(var_names[j]);
            vars[j][i]->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
            layout_3[i]->addWidget(vars[j][i],(int)j%6 + 1,2*(int)(j/6));

            if(j!=13){
                edits[j][i] = new QLineEdit;
                layout_3[i]->addWidget(edits[j][i],(int)j%6 + 1,1+2*(int)(j/6));
            }else {
                imag_box[i] = new QComboBox;
                imag_box[i]->insertItem(1,"Volume");
                imag_box[i]->insertItem(2,"Surface");
                layout_3[i]->addWidget(imag_box[i],(int)j%6 + 1,1+2*(int)(j/6));
            }
        }

        text[i] = new QLabel("<div style=\"margin: 0 auto;\"><table cellpadding=5><tr><td colspan=2> Depth of real nuclear well</td><td> V0=V0D(I)+V01D(I)*ECM+V02D(I)*ECM<sup>2</sup></td></tr>"
                             "<tr><td colspan=2>Depth of imaginary nuclear well</td><td>  W0=W0D(I)+W01D(I)*ECM+W02D(I)*ECM**2</td></tr></table>  "
                             "<p style=\"margin-left: 5px;\">IMAG(I) imaginary well is potential surface or volume</p>"
                             "<p style=\"margin-left: 5px;\">Radial extension of nuclear well </p>"
                             "<table cellpadding=5><tr><td>IF IRAD(I)=0</td><td>  Radius of real nuclear well</td><td> R0RD(I)</td></tr>"
                             "<tr><td></td><td>Radius of imaginary nuclear well</td><td>  R0ID(I)</td></tr>"
                             "<tr><td>IF IRAD(I)=1 </td><td> Radius of real well </td><td> R0RD(I)* AMT**(1/3)</td></tr>"
                             "<tr><td></td><td>  Radius of imaginary well</td><td> R0ID(I)* AMT**(1/3)</td></tr>"
                             "<tr><td>IF IRAD(I)=2</td><td>  Radius of real well</td>    "
                             "<td>  R0RD(I)*(AMT**(1/3)+AMP**(1/3))</td></tr>"
                             "<tr><td></td><td>             Radius of imaginary well</td>    "
                             "<td> R0ID(I)*(AMT**(1/3)+AMP**(1/3))</td></tr></table>"
                             "<p style=\"margin-left: 5px;\">  Radius of Coulomb potential is R0CD(I) and is controlled by IRAD(I)</p></div>");
        layout_3[i]->addWidget(text[i],7,0,1,6);
        w[i]->setLayout(layout_3[i]);
        w[i]->setFixedWidth(600);

        pages->addWidget(w[i]);
        setPage(i+1,p);
        p->OPTPOT(i+1,_IAC-_IZC,_IZC);
        //ShowHideEdit(_ITAKE[i+1],i+1);
    }
    setVariables(p);

    lay->addWidget(pages);
    CurrentPage=1;
    setLayout(lay);
  //  setPage(CurrentPage);
}
void tfusion3::setVariables(PACE *p){
    for(int page=1;page<=4;page++){
        if( p->_ITAKE[page] == true){
            opt_man[page-1]->setChecked(true);
           // opt_changed();
        } else {
            opt_sys[page-1]->setChecked(true);
            //opt_changed();
        }

        if(!p->_ITAKE[page] && page!=3) p->OPTPOT(page,_IAC-_IZC,_IZC);

        imag_box[page-1]->setCurrentIndex(_IMAG[page]);
        edits[0][page-1]->setText(QString::number(_V0D[page]));
        edits[1][page-1]->setText(QString::number(_V01D[page]));
        edits[2][page-1]->setText(QString::number(_V02D[page]));
        edits[3][page-1]->setText(QString::number(_R0RD[page]));
        edits[4][page-1]->setText(QString::number(_ARD[page]));
        edits[5][page-1]->setText(QString::number(_R0CD[page]));
        edits[6][page-1]->setText(QString::number(_W0D[page]));
        edits[7][page-1]->setText(QString::number(_W01D[page]));
        edits[8][page-1]->setText(QString::number(_W02D[page]));
        edits[9][page-1]->setText(QString::number(_R0ID[page]));
        edits[10][page-1]->setText(QString::number(_AID[page]));
        edits[11][page-1]->setText(QString::number(_RMCHD[page]));
        edits[12][page-1]->setText(QString::number(_NPD[page]));
        edits[14][page-1]->setText(QString::number(_IRAD[page]));
        ShowHideEdit( p->_ITAKE[page],page);
    }
}
void tfusion3::setPage(int page, PACE* p){


    if(! p->_ITAKE[page] && page!=3) p->OPTPOT(page,_IAC-_IZC,_IZC);
    //setVariables();
    //ShowHideEdit(_ITAKE[page],page);

    if(page <= 4){
        imag_box[page-1]->setCurrentIndex(_IMAG[page]);
        edits[0][page-1]->setText(QString::number(_V0D[page]));
        edits[1][page-1]->setText(QString::number(_V01D[page]));
        edits[2][page-1]->setText(QString::number(_V02D[page]));
        edits[3][page-1]->setText(QString::number(_R0RD[page]));
        edits[4][page-1]->setText(QString::number(_ARD[page]));
        edits[5][page-1]->setText(QString::number(_R0CD[page]));
        edits[6][page-1]->setText(QString::number(_W0D[page]));
        edits[7][page-1]->setText(QString::number(_W01D[page]));
        edits[8][page-1]->setText(QString::number(_W02D[page]));
        edits[9][page-1]->setText(QString::number(_R0ID[page]));
        edits[10][page-1]->setText(QString::number(_AID[page]));
        edits[11][page-1]->setText(QString::number(_RMCHD[page]));
        edits[12][page-1]->setText(QString::number(_NPD[page]));
        edits[14][page-1]->setText(QString::number(_IRAD[page]));
        ShowHideEdit( p->_ITAKE[page],page);
        pages->setCurrentWidget(w[page-1]);
    }
}

void tfusion3::ShowHideEdit(int flag, int page){
    for(int i=0;i<15;i++){
        if(i !=13){
            edits[i][page-1]->setEnabled(flag);
        }else {
            imag_box[page-1]->setEnabled(flag);
        }
    }
}
void tfusion3::opt_changed(){
    int i = CurrentPage-1;
    if(opt_sys[i]->isChecked()){
         pace->_ITAKE[CurrentPage] = false;
         pace->OPTPOT(i+1,_IAC-_IZC,_IZC);
         setVariables(pace);
    } else {
         pace->_ITAKE[CurrentPage] = true;
    }
    ShowHideEdit(pace->_ITAKE[CurrentPage],CurrentPage);
}
void tfusion3::readPage(int page, PACE* p) {
    if(opt_sys[page-1]->isChecked()){
         //p->_ITAKE[CurrentPage] = false;
         p->_ITAKE[page] = false;
    } else {
        // p->_ITAKE[CurrentPage] = true;
        p->_ITAKE[page] = true;
        _V0D[page] = edits[0][page-1]->text().toDouble();
        _V01D[page] = edits[1][page-1]->text().toDouble();
        _V02D[page] = edits[2][page-1]->text().toDouble();
        _R0RD[page] = edits[3][page-1]->text().toDouble();
        _ARD[page] = edits[4][page-1]->text().toDouble();
        _R0CD[page] = edits[5][page-1]->text().toDouble();
        _W0D[page] = edits[6][page-1]->text().toDouble();
        _W01D[page] = edits[7][page-1]->text().toDouble();
        _W02D[page] =  edits[8][page-1]->text().toDouble();
        _R0ID[page] = edits[9][page-1]->text().toDouble();
        _AID[page] = edits[10][page-1]->text().toDouble();
        _RMCHD[page] = edits[11][page-1]->text().toDouble();
        _NPD[page] = edits[12][page-1]->text().toDouble();
        _IRAD[page] = edits[14][page-1]->text().toDouble();
    }
}
