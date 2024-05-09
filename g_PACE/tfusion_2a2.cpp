#include "tfusion_2a2.h"
#include "ftype.h"

#include <QLabel>
#include <QGridLayout>
#include <QGroupBox>

extern double _ALEVB[Max_MOM];
extern double _AJNUC;
extern int _LMINN, _IAC;

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
tfusion_2a2::tfusion_2a2(QWidget *parent) :
    QWidget(parent)
{
    Current = 1;
    Current=max(1,min(_LMINN, (int)_AJNUC));
    for(int i=1; i<=min(_LMINN, (int)_AJNUC); i++)  _ALEVB[i]=0;
    QGridLayout *layout2a2 = new QGridLayout;
    QLabel *title = new QLabel("Partial Cross Sections");
    layout2a2->addWidget(title);
    QGroupBox *data_box = new QGroupBox("Data");
    QGridLayout *data_layout = new QGridLayout;
    QLabel *head[2];
    head[0] = new QLabel("J, hbar");
    head[1] = new QLabel("CS, mb");
    data_layout->addWidget(head[0],0,1);
    data_layout->addWidget(head[1],0,2);
    up = new QPushButton("Up");
    connect(up,SIGNAL(clicked()),this,SLOT(up_clicked()));
    down = new QPushButton("Down");
    connect(down,SIGNAL(clicked()),this,SLOT(down_clicked()));
    data_layout->addWidget(up,1,0,2,1);
    data_layout->addWidget(down,3,0,2,1);
    signalMapper = new QSignalMapper(this);

    for(int i=0;i<4;i++)
        {
        for(int j=0;j<2;j++)
            {
            edit[i][j] = new QLineEdit;
            if(j == 0){
                edit[i][j]->setDisabled(true);
                }
            else {
                connect(edit[i][j], SIGNAL(editingFinished()), signalMapper, SLOT(map()));
                signalMapper->setMapping(edit[i][j],i);
                }
            data_layout->addWidget(edit[i][j],i+1,j+1);
            }
        }
    connect(signalMapper, SIGNAL(mappedInt(int)),this,SLOT(edited(int)));
    data_box->setLayout(data_layout);
    layout2a2->addWidget(data_box);

    Current=max(1,min(_LMINN, (int)_AJNUC));
    for(int i=1; i<=min(_LMINN, (int)_AJNUC); i++)  _ALEVB[i]=0;
    setPage();
    //layout2a2->setMargin(200);
    layout2a2->setContentsMargins(20,0,20,40);  //Qt6
    setLayout(layout2a2);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2a2::setPage()
{

    double r=(_IAC%2)*0.5;

    for(int i=0;i<4;i++)
      {
        edit[i][0]->setText(QString::number(Current+(i-1)+r,'f',1));
        edit[i][1]->setText(QString::number(_ALEVB[Current+i],'f',3));
        edit[i][1]->setEnabled(Current+i > _LMINN);
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2a2::up_clicked()
{
    if(Current<=1) return;
    Current--;
    setPage();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2a2::down_clicked()
{
    if(Current>=_AJNUC-2-(_IAC%2)) return;
    Current++;
    setPage();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion_2a2::edited(int k)
{
    _ALEVB[Current+k] = edit[k][1]->text().toDouble();

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

