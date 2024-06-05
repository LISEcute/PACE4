#include <QFormLayout>
#include <QFileDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QDebug>
#include <QStackedWidget>
#include <QDesktopServices>
#include <QUrl>
#include <QMessageBox>
#include <QKeyEvent>
#include <QToolButton>

#include "p_about.h"
#include "tfusion.h"
#include "ui_tfusion.h"
#include "ftype.h"
#include "ParticleStates.h"


#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))

extern int _NCASC, _IDIST, _INPUT, _MDIR, _ITRAC, _NOSHL;
extern double _FYRST, _ARATIO, _FACLA, _BARFACr;


extern double _ELAB,_EXPSIGr,_AGRAZr,_ELOSS,_ELAB_MAX;
extern int _LMINN, _JCMAXr;
extern double _LowLimit, _HighLimit;


extern int _IZP,_IAP,_IZT,_IAT;
extern double _SPr, _ST, _QCNr;

extern int _IZC, _IAC;
extern double _EEXCN, _EREC, _AJNUC, _EEXCN_MAX;

extern double  _FGE1r,_FGM1r,_FGE2r,_FGM2r;
extern double _ALEVB[Max_MOM+1];

extern double _V0D[5],_R0RD[5],_ARD[5],_R0CD[5],_RCCC[4],_VQS[4],_VQ0[4],_W0D[5],_R0ID[5],_AID[5],
_RMCHD[5],_V01D[5],_V02D[5],_W01D[5],_W02D[5];
extern int _NPD[5],_IMAG[5],_IRAD[5];//, _ITAKE[5];

extern int _BarrierMode;
extern double _h_omega;

extern int _BatchMode; // 0-No, 1-Yes
extern int _CurrentBatchStep;
extern int _BatchNsteps;
extern const char *FileNameAbsent;

extern QString LISErootPATH;
extern QString localPATH;
extern const char *FileNameAbsent;
extern QString FileArg;
extern FILE *mfopen(const QString& filename, const char* operand);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
extern bool init_base(QString name);
extern void InitRandom();
extern FILE *f10;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
extern char *GetNextSymbol(char *s);
extern char *GetNextDelimeter(char *s);
extern char *GetIntFromString(int &V, char *s);
extern char *GetDoubleFromString(double &V, char *s);

QSqlDatabase db;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

tfusion::tfusion(QWidget *parent):
  QMainWindow(parent),
  ui(new Ui::tfusion)
{
  ui->setupUi(this);
  setWindowFlags( Qt::Window | Qt::CustomizeWindowHint |
                  Qt::WindowTitleHint | Qt::WindowSystemMenuHint |
                  Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint );
  p = new PACE();

  InitRandom();
  init_base();
  FFileNameEve =  localPATH + "/PACE4evt.dat";
  //qDebug() << FFileNameEve;

  //--------------------------------------------------------

  _NCASC=1000;
  _IDIST=1;
  _INPUT=1;
  _MDIR=0;
  _ITRAC=0;
  _NOSHL=0;
  _FYRST=0;
  _ARATIO=1;
  _FACLA=10;
  _BARFACr=0;
  _LowLimit = 2;
  _HighLimit = 100;

  _FGE1r=_FGM1r=_FGE2r=_FGM2r=0;

  for(int i=0; i<=Max_MOM; i++) _ALEVB[i]=0;

  /************* Added from TMFusion2 constructor **********************/
  _ELAB  = 300;
  _ELAB_MAX  = 390;
  _BatchNsteps = 10;
  _EXPSIGr = 0;
  _AGRAZr = 2;
  _ELOSS =0;
  _LMINN = 0;
  _JCMAXr = 0;

  _IZP = 20;
  _IAP = 48;
  _IZT = 50;
  _IAT = 124;
  _SPr = 0;
  _ST = 0;
  _QCNr = 0;

  //if(_INPUT != 1){
  _IZC = 85;
  _IAC = 215;
  _LMINN= 0;
  _EEXCN= 50;
  _EEXCN_MAX=150;
  _BatchNsteps = 10;
  _EREC= 0;
  _AJNUC= 40;

  for(int page=1; page<=4; page++)  {
      p->_ITAKE[page]= 1;
      _V0D[page]  = 47;
      _V01D[page] =  -0.267;
      _V02D[page] =  -0.002;
      _R0RD[page] =  1.264;
      _ARD[page]  =  0.660;
      _R0CD[page] =  0;
      _W0D[page]  =  9.520;
      _W01D[page] =  -0.053;
      _W02D[page] =   0;
      _R0ID[page] =  1.332;
      _AID[page]  =  0.480;
      _RMCHD[page]=  0;
      _NPD[page]  = 250;
      _IMAG[page] = 0;
      _IRAD[page] = 1;
    }


  QString labels[10] = {"NCASC","INPUT","FYRST","BARFAC","ARATIO","FACLA","IDIST","MDIR","ITRAC","NOSHL"};
  QString descriptions[10] = {"number of cascades (events in Monte Carlo calculation &lt; 1e6).", "= 1  projectile + target input. AGRAZ parameter determines diffuseness of partial wave "
                              "distribution <br>"
                              " = 2  compound nucleus input for single spin.<br>"
                              " = 3  compound nucleus input. Spin distribution read in.<br>"
                              " = 4  compound nucleus input. Spin distribution calculated taking spin-cutoff paramater at given Ex. <br>"
                              " = 5  triangular  (sigma = 2l+1 ) cross section between LMINN and maximum spin parameter determining<br> yrast line to be used.",
                              " &lt; 0 provides the G-C yrast line.<br> "
                              " &lt; 0  Gilbert-Cameron spin cutoff parameter. EROT = (SPIN)**2/(2.*SIGSQ)<br> "
                              " != 0  EROT = rotating liquid drop rotational energy,multiplied by factor of FYRST.<br> "
                              " == 0  value changed to FYRST = 1.  In both cases level density calculated aT E = EX-EROT.",
                              "The program assumes the A.J.Sierk modified rotating liquid drop barrier if this is equal to 0. <br>If "
                              "you provide a fission barrier of your own, the Sierk barrier will be renormalized accordingly. <br>"
                              "If BarFac is positive it will be taken as the desired zero spin fission barrier. <br>"
                              "If BarFac is negative, its absolute value will be taken as a factor to multiply the Sierk barrier.",
                              "Ratio of the Fermi gas level density parameter 'LITTLE-A' at the saddle point to the ground <br>"
                              "state value. The saddle point level density is determined by g.s. 'LITTLE-A' * ARATIO.",
                              "level density parameter = MASS/FACLA if not zero. <br>"
                              "if == 0  Gilbert and Cameron value used.","","",""};

  stackedw = new QStackedWidget;

  QWidget *page_1 = new QWidget;
  page_2 = new tfusion_2;
  page_2a = new tfusion_2a;
  page_2a2 = new tfusion_2a2;
  page3 = new tfusion3(p);
  QGridLayout *gridLayout1 = new QGridLayout;

  for(int i=0;i<6;i++){
      QLabel *l = new QLabel(labels[i]);
      QLabel *d = new QLabel(descriptions[i]);

      //if(i !=1){
      gridLayout1->addWidget(l,i,0);
      gridLayout1->addWidget(d,i,2,1,4);
      // }
    }

  lineEdit[0] = new QLineEdit;
  gridLayout1->addWidget(lineEdit[0],0,1,1,1);

  inp_combo = new QComboBox;
  inp_combo->addItem("1");
  inp_combo->addItem("2");
  inp_combo->addItem("3");
  inp_combo->addItem("4");
  inp_combo->addItem("5");
  gridLayout1->addWidget(inp_combo,1,1);
  // gridLayout1->addWidget(inp_box,1,0,1,6);

  for(int i=1;i<5;i++){
      lineEdit[i] = new QLineEdit;
      gridLayout1->addWidget(lineEdit[i],i+1,1);
    }

  QGroupBox *idist_box = new QGroupBox(labels[6]);
  QVBoxLayout *idist_layout = new QVBoxLayout;
  //idist_layout->setMargin(5);
  QString idist_text[3] = {" = 0 brief, schematic results of particle spectra and list of evaporated (residual) nuclei",
                           " = 1 detailed angular and energy distribution of residual nuclei and evaporated particles",
                           " = 2 detailed(#1) + transmission coefficients for particle emission"};

  for(int i=0;i<3;i++){
      idist[i] = new QRadioButton(idist_text[i]);
      connect(idist[i],SIGNAL(clicked()),this,SLOT(idist_clicked()));
      idist_layout->addWidget(idist[i]);
    }

  idist[_IDIST]->setChecked(true);
  idist_box->setLayout(idist_layout);
  // idist_box->setFixedWidth(530);
  gridLayout1->addWidget(idist_box,6,0,1,5);
  QWidget *lim_box = new QWidget;
  QGridLayout *lim_layout = new QGridLayout;
  // lim_layout->setMargin(0);
  QLabel *lim_descr = new QLabel("<p>Limits of residual<br> yields (in %) to show <br> angular and energy <br>distributions </p>");
  lim_descr->setAlignment(Qt::AlignBottom);
  lim_layout->addWidget(lim_descr,0,0,1,2);
  QLabel *lim_low = new QLabel("Low limit");
  lim_layout->addWidget(lim_low,1,0);
  limits[0] = new QLineEdit;
  lim_layout->addWidget(limits[0],1,1);
  QLabel *lim_high = new QLabel("High limit");
  lim_layout->addWidget(lim_high,2,0);
  limits[1] = new QLineEdit;
  lim_layout->addWidget(limits[1],2,1);
  lim_box->setLayout(lim_layout);
  gridLayout1->addWidget(lim_box,6,5,2,1);
  QGroupBox *mdir_box = new QGroupBox(labels[7]+" -  (MDIR=1 is appropriate for deep inelastic fragment de-excitation)");
  QVBoxLayout *mdir_layout = new QVBoxLayout;
  //  mdir_layout->setMargin(0);

  QString mdir_text[2] = {" = 0 compound nucleus is initially M=0 states and the Z-axis is the recoil axis",
                          " = 1 compound nucleus is initially M=J states and the Z-axis is perp. to recoil direction"};
  for(int i=0;i<2;i++){
      mdir[i] = new QRadioButton(mdir_text[i]);
      mdir_layout->addWidget(mdir[i]);
    }
  mdir[_MDIR]->setChecked(true);
  mdir_box->setLayout(mdir_layout);
  //mdir_box->setFixedWidth(530);
  gridLayout1->addWidget(mdir_box,7,0,1,5);
  gridLayout1->setColumnMinimumWidth(4,130);

  QGroupBox *itrac_box = new QGroupBox(labels[8]+" -  (controls the degree of event traceback)");
  QVBoxLayout *itrac_layout = new QVBoxLayout;
  // itrac_layout->setMargin(0);
  QString itrac_text[2] = {" = 0 produces compact traceback, summed over all residues",
                           " = 1 detailed traceback leading to each individual isotope separately"};
  for(int i=0;i<2;i++){
      itrac[i] = new QRadioButton(itrac_text[i]);
      itrac_layout->addWidget(itrac[i]);
    }
  itrac[_ITRAC]->setChecked(true);
  itrac_box->setLayout(itrac_layout);
  gridLayout1->addWidget(itrac_box,8,0,1,3);

  //-----------------------------------------------
  QGroupBox *noshl_box = new QGroupBox(labels[9]);
  QVBoxLayout *noshl_layout = new QVBoxLayout;
  QString noshl_text[2] = {" = 0 uses AME2016 values Chinese Physics C41 (2017) 030003",
                           " = 1 uses Lysekil massese with shell correction"};
  for(int i=0;i<2;i++){
      noshl[i] = new QRadioButton(noshl_text[i]);
      noshl_layout->addWidget(noshl[i]);
    }

  noshl[_NOSHL]->setChecked(true);
  noshl_box->setLayout(noshl_layout);
  gridLayout1->addWidget(noshl_box,9,0,1,3);
  //-----------------------------------------------

  QGroupBox *partBox = new QGroupBox("Particle analysis");
  QGridLayout *partLayout = new QGridLayout;
  partOutput[0] = new QCheckBox("Create output file");
  q = new QPushButton(" ? ");
  q->setMaximumWidth(20);
  q->setStyleSheet("margin-top: 0px;");
  connect(q,SIGNAL(clicked()),this, SLOT(q_clicked()));
  connect(partOutput[0],SIGNAL(clicked()),this,SLOT(createFile_clicked()));
  partOutput[1] = new QCheckBox("neutron");
  partOutput[2] = new QCheckBox("proton");
  partOutput[3] = new QCheckBox("alpha");
  partOutput[4] = new QCheckBox("gamma");
  for(int i=0;i<5;i++){
      if(i==0){
          partLayout->addWidget(partOutput[i],i,0,1,2);
        } else {
          partLayout->addWidget(partOutput[i],i,0);
        }
      if(i>0){
          partOutput[i]->setStyleSheet("margin-left: 7px;");
          partOutput[i]->setDisabled(true);
        }
      if(i == 1){
          q->setFixedWidth(30);
          q->setStyleSheet("padding:3px;");
          partLayout->addWidget(q,i,1,1,1);
        }

    }

  QGroupBox *gateBox = new QGroupBox("Nucleus gate");
  QGridLayout *gateLayout = new QGridLayout;

  use_gate = new QCheckBox("Use");
  connect(use_gate,SIGNAL(clicked()),this,SLOT(use_gate_clicked()));
  gateLayout->addWidget(use_gate,0,0,1,2);

  Zl = new QLabel("Z = ");
  gateAZ[0] = new QLineEdit(QString::number(p->ParticleFlags.Zgate)); 
  gateLayout->addWidget(Zl,1,0);
  gateLayout->addWidget(gateAZ[0],1,1);

  Al = new QLabel("A = ");
  gateAZ[1] = new QLineEdit(QString::number(p->ParticleFlags.Ngate+p->ParticleFlags.Zgate));
  gateLayout->addWidget(Al,2,0);
  gateLayout->addWidget(gateAZ[1],2,1);

  gateBox->setLayout(gateLayout);
  use_gate->setEnabled(false);
  gateAZ[0]->setEnabled(false);
  gateAZ[1]->setEnabled(false);

  Al->setEnabled(false);
  Zl->setEnabled(false);
  partLayout->addWidget(gateBox,0,2,5,1);
  partBox->setLayout(partLayout);
  gridLayout1->addWidget(partBox,8,4,2,2);
  // gridLayout1->setSpacing(10);
  page_1->setLayout(gridLayout1);
  setVariables();

  //------------------------------------------------------
  gridLayout1->setColumnStretch(2,1);
  stackedw->addWidget(page_1);
  stackedw->addWidget(page_2);
  stackedw->addWidget(page_2a);
  stackedw->addWidget(page_2a2);
  stackedw->addWidget(page3);

  page_1->setStyleSheet("QRadioButton { spacing: 0px;} ");
  // page_2->setStyleSheet("QRadioButton { spacing: 0px;}");
  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget(stackedw);
  this->setCentralWidget(stackedw);

  ui->action_Previous->setEnabled(false);
  stackedw->setCurrentIndex(0);
  //------------------------------------------------------

  FFileName  = FileArg.size()>0? FileArg : localPATH + FileNameAbsent;
  setFileName(FFileName);
  if(FileArg.size()>0) readFile(FFileName);

  QSize size = this->minimumSizeHint();
  this->resize(size);

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
tfusion::~tfusion()
{
  if(db.isOpen() != false) db.close();
  qApp->quit();
  delete ui;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::getVariables()
{
  QString val;

  for(int i=0;i<6;i++)
    {
      if(i != 1 && i != 0)    val = lineEdit[i-1]->text();

      if(i== 0)
        {
          val = lineEdit[i]->text();
          _NCASC = val.toInt();
        }
      else if (i ==1) {            _INPUT = inp_combo->currentIndex()+1; }
      else if (i == 2){            _FYRST = val.toFloat();}
      else if (i == 3){            _BARFACr = val.toFloat();}
      else if (i == 4){            _ARATIO = val.toFloat();}
      else if (i == 5){            _FACLA = val.toFloat();}
    }

  for(int j=0;j<3;j++)
    if(idist[j]->isChecked()){            _IDIST = j; }

  for(int j=0;j<2;j++){
      if(mdir[j]->isChecked()){            _MDIR = j;      }
      if(itrac[j]->isChecked()){            _ITRAC = j;    }
      if(noshl[j]->isChecked()){            _NOSHL = j;    }
    }

  if(partOutput[0]->isChecked())
    {
      p->ParticleFlags.SetParticleStateValue(ep_permit,true);
      p->ParticleFlags.SetParticleStateValue(ep_n,partOutput[1]->isChecked());
      p->ParticleFlags.SetParticleStateValue(ep_p,partOutput[2]->isChecked());
      p->ParticleFlags.SetParticleStateValue(ep_a,partOutput[3]->isChecked());
      p->ParticleFlags.SetParticleStateValue(ep_g,partOutput[4]->isChecked());
    }

  if(use_gate->isChecked())
    {
      p->ParticleFlags.UseGate = true;
      p->ParticleFlags.Zgate = gateAZ[0]->text().toInt();
      p->ParticleFlags.Ngate = gateAZ[1]->text().toInt() - p->ParticleFlags.Zgate;
    }
  else { p->ParticleFlags.UseGate = false;}

  _LowLimit = limits[0]->text().toInt();
  _HighLimit = limits[1]->text().toInt();

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::RunPACE3()
{

  if(_BatchMode)
    {
      QString FileNameSave;//, FileNameTemp;
      double _ELAB_save = _ELAB;
      double _EEXCN_save = _EEXCN;

      FileNameSave = FFileName;

      QString FN = FFileName;
      QString FP = FN.split(".",Qt::SkipEmptyParts).at(0);

      QString Fext = "txt"; //QString::fromStdString(FFileName).split(".",Qt::SkipEmptyParts).at(1);

      double h    = (_ELAB_MAX - _ELAB) / (_BatchNsteps-1);
      double h_CN = (_EEXCN_MAX - _EEXCN) / (_BatchNsteps-1);

      for(_CurrentBatchStep=0; _CurrentBatchStep < _BatchNsteps; _CurrentBatchStep++)
        {
          QString FileNameTemp =  FP+ "_step" + QString::number(_CurrentBatchStep+1) + "." + Fext;

          setFileName(FileNameTemp);
          if(_INPUT==1)
            {
              _ELAB =  _ELAB_save + double(_CurrentBatchStep)*h;
              extern double _QcLoc2;
              double ECM=_ELAB*double(_IAT)/double(_IAC);
              _EEXCN=ECM+_QcLoc2;
            }
          else    _EEXCN=  _EEXCN_save + double(_CurrentBatchStep)*h_CN;


          RunPACE3local();
        }

      _ELAB = _ELAB_save;
      _EEXCN = _EEXCN_save;
      setFileName(FileNameSave);
    }
  else
    RunPACE3local();

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void tfusion::RunPACE3local()
{

  int res = p->PACE3(FFileNameOut.toStdString().c_str(),FFileNameEve.toStdString().c_str(),
                     FFileNameCS.toStdString().c_str(), FFileNameParticles.toStdString().c_str(),
                     FFileNameHtml.toStdString().c_str());

  if(p->ParticleFlags.IsParticleStateSet(ep_opened))
    {
      //      FormParticles->Show();
      p->ParticleFlags.ClearParticleState(ep_opened);
    }
  if(p->ParticleFlags.IsParticleStateSet(ep_permit))
    {
      part_page = new particles(FFileNameParticles);
      part_page->setWindowFlags(Qt::CustomizeWindowHint |
                                Qt::WindowTitleHint | Qt::WindowCloseButtonHint);
      part_page->setModal(false);
      part_page->show();
    }

  // res_page->exec();
  if(res != 10){
      res_page = new results(FFileNameHtml);
      res_page->show();
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::on_actionExecute_triggered()
{
  getVariables();
  if(_INPUT == 1){        page_2->getVariables();    }
  else           {        page_2a->getVariables();   }

  for(int i=0;i<4;i++)
    page3->readPage(i+1,p);

  RunPACE3();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::on_actionNext_triggered()
{
  int i = stackedw->currentIndex();
  getVariables();

  if(i == 0 && _INPUT == 1)
    {
      stackedw->setCurrentIndex(1);
      page_2->setVariables();    //  Qt-Oleg
      ui->action_Previous->setEnabled(true);
    }
  else if(i == 0) {
      stackedw->setCurrentIndex(2);
      ui->action_Previous->setEnabled(true);
    }
  else if(i == 1)
    {
      page_2->getVariables();
      stackedw->setCurrentIndex(4);
      page3->CurrentPage = 1;
      page3->setPage(page3->CurrentPage,p);
    }
  else if(i == 2 && _INPUT != 3)
    {
      page_2a->getVariables();
      stackedw->setCurrentIndex(4);
    }
  else if(i == 2)
    {
      page_2a->getVariables();
      stackedw->setCurrentIndex(3);
    }
  else if(i == 3)
    {
      stackedw->setCurrentIndex(4);
      page3->CurrentPage = 1;
      page3->setPage(page3->CurrentPage,p);
    }
  else if(i == 4)
    {
      int j = page3->CurrentPage;
      page3->readPage(j,p);
      if(j<4) {
          ui->actionNext->setEnabled(true);
          page3->CurrentPage++;
          page3->setPage(page3->CurrentPage,p);
        }
      if(j == 3){ ui->actionNext->setEnabled(false); }
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::on_action_Previous_triggered()
{
  int i = stackedw->currentIndex();

  if(i == 1){
      page_2->getVariables();
      stackedw->setCurrentIndex(0);
      ui->action_Previous->setEnabled(false);
    }
  else if(i == 2){
      page_2a->getVariables();
      stackedw->setCurrentIndex(0);
      ui->action_Previous->setEnabled(false);
    }
  else if(i == 3) {
      stackedw->setCurrentIndex(2);
    }
  else if(i == 4){
      int j = page3->CurrentPage;
      page3->readPage(j,p);
      if(j==1){
          if(_INPUT==1){       stackedw->setCurrentIndex(1);  }
          else if(_INPUT!=3){       stackedw->setCurrentIndex(2);  }
          else {                    stackedw->setCurrentIndex(3);  }
        }
      else if(j<=4 && j >1) {
          ui->actionNext->setEnabled(true);
          page3->CurrentPage--;
          page3->setPage(page3->CurrentPage,p);
        }
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::createFile_clicked()
{
  if(partOutput[0]->isChecked())
    {
      for(int i=1;i<5;i++)
        {partOutput[i]->setEnabled(true);}

      use_gate->setEnabled(true);
      p->ParticleFlags.SetParticleStateValue(ep_permit,true);    }
  else {
      for(int i=1;i<5;i++){
          partOutput[i]->setChecked(false);
          partOutput[i]->setEnabled(false);
        }
      use_gate->setEnabled(false);

      p->ParticleFlags.SetParticleStateValue(ep_permit,false);
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::idist_clicked(){
  for(int i=0;i<3;i++){
      if(idist[i]->isChecked()){
          if(i>0){
              partOutput[0]->setEnabled(true);
              q->setEnabled(true);
            } else {
              partOutput[0]->setChecked(false);
              createFile_clicked();
              partOutput[0]->setDisabled(true);
              q->setEnabled(false);
            }
        }
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::use_gate_clicked()
{
  if(use_gate->isChecked())
    {
      p->ParticleFlags.UseGate = true;
      gateAZ[0]->setEnabled(true);
      Zl->setEnabled(true);
      Al->setEnabled(true);
      gateAZ[1]->setEnabled(true);
    }
  else {
      gateAZ[0]->setEnabled(false);
      gateAZ[1]->setEnabled(false);
      Zl->setEnabled(false);
      Al->setEnabled(false);
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::on_actionWeb_Documentation_triggered()
{
  QDesktopServices::openUrl(QUrl("http://lise.nscl.msu.edu/pace4.html"));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::on_actionOpen_triggered()
{

  QString file = QFileDialog::getOpenFileName(this,tr("Open File"),
                                              FFileName,tr("PACE files (*.pace);;All files (*.*)"));

  if(file.size()<=0) return;

  setFileName(file);
  readFile(file);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::readFile(const QString& file)
{
  QFile f(file);

  if(!f.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      QMessageBox MB;
      MB.setWindowTitle("Warning!");
      MB.setText(tr("Could not open file.\n%1").arg(file));
      MB.setIcon(QMessageBox::Warning);
      MB.exec();
      return;
    }

  QTextStream in(&f);

  QString line;
  line=in.readLine();
  QStringList fields = line.split(' ',Qt::SkipEmptyParts);
  if(fields.size() != 10 && fields.size() != 11) {format_error();}
  else {
      _NCASC = fields[0].toInt();
      _IDIST = fields[1].toInt();
      _INPUT = fields[2].toInt();
      _FYRST = fields[3].toDouble();
      _BARFACr = fields[4].toDouble();
      _ARATIO = fields[5].toDouble();
      _MDIR = fields[6].toInt();
      _FACLA = fields[7].toDouble();
      _ITRAC = fields[8].toInt();
      _NOSHL = fields[9].toInt();
      if(fields.size() == 11){
          int k = fields[10].toInt();
          p->ParticleFlags.SetParticleFlags(k);
        }
    }

  setVariables();
  line=in.readLine();
  fields = line.split(' ',Qt::SkipEmptyParts);

  if(_INPUT == 1){
      if(fields.size() != 7 && fields.size() != 8) {
          format_error();
          return;
        }

      _IZP = fields[0].toInt();
      _IAP = fields[1].toInt();
      _IZT = fields[2].toInt();
      _IAT = fields[3].toInt();
      _SPr = fields[4].toDouble();
      _ST = fields[5].toDouble();
      _QCNr = fields[6].toDouble();

      line=in.readLine();
      fields = line.split(' ',Qt::SkipEmptyParts);
      if(fields.size() != 6 && fields.size() != 8) {format_error(); }
      else {
          _ELAB = fields[0].toDouble();
          _EXPSIGr = fields[1].toDouble();
          _JCMAXr = fields[2].toInt();
          _AGRAZr = fields[3].toDouble();
          _ELOSS = fields[4].toDouble();
          _LMINN = fields[5].toInt();
          if(fields.size() == 8){
              _BarrierMode = fields[6].toInt();
              _h_omega = fields[7].toDouble();
            }
        }
    }else {
      if(fields.size() != 6) {
          format_error();
          return;
        }
      _IZC = fields[0].toInt();
      _IAC = fields[1].toInt();
      _EEXCN = fields[2].toDouble();
      _EREC = fields[3].toDouble();
      _AJNUC = fields[4].toDouble();
      _LMINN = fields[4].toInt();
    }

  if(_INPUT == 3){
      line=in.readLine();
      fields = line.split(' ',Qt::SkipEmptyParts);
      for(int i=1;i<fields.size();i++){
          _ALEVB[i] = fields[i].toDouble();
        }

    }

  if(_INPUT == 1){
      line=in.readLine();
      fields = line.split(' ',Qt::SkipEmptyParts);
      if(fields.size() != 15 && fields.size() != 0) {
          format_error();
          return;
        } else if(fields.size() == 0) {
          p->_ITAKE[4] = 0;
        } else {
          int I = 4;
          p->_ITAKE[I] = 1;
          _V0D[I] = fields[0].toDouble();
          _V01D[I] = fields[1].toDouble();
          _V02D[I] = fields[2].toDouble();
          _R0RD[I] = fields[3].toDouble();
          _ARD[I] = fields[4].toDouble();
          _R0CD[I] = fields[5].toDouble();
          _W0D[I] = fields[6].toDouble();
          _W01D[I] = fields[7].toDouble();
          _W02D[I] = fields[8].toDouble();
          _R0ID[I] = fields[9].toDouble();
          _AID[I] = fields[10].toDouble();
          _RMCHD[I] = fields[11].toDouble();
          _NPD[I] = fields[12].toInt();
          _IMAG[I] = fields[13].toInt();
          _IRAD[I] = fields[14].toInt();

        }
    }
  for(int I=1;I<=3;I++){
      line=in.readLine();
      fields = line.split(' ',Qt::SkipEmptyParts);
      if(fields.size() != 15 && fields.size() != 0) {
          format_error();
          return;
        } else if(fields.size() == 0) {
          p->_ITAKE[I] = 0;
        }else {
          p->_ITAKE[I] = 1;
          _V0D[I] = fields[0].toDouble();
          _V01D[I] = fields[1].toDouble();
          _V02D[I] = fields[2].toDouble();
          _R0RD[I] = fields[3].toDouble();
          _ARD[I] = fields[4].toDouble();
          _R0CD[I] = fields[5].toDouble();
          _W0D[I] = fields[6].toDouble();
          _W01D[I] = fields[7].toDouble();
          _W02D[I] = fields[8].toDouble();
          _R0ID[I] = fields[9].toDouble();
          _AID[I] = fields[10].toDouble();
          _RMCHD[I] = fields[11].toDouble();
          _NPD[I] = fields[12].toInt();
          _IMAG[I] = fields[13].toInt();
          _IRAD[I] = fields[14].toInt();
        }
    }

  stackedw->setCurrentIndex(0);           // v.4.34.7

  if(_INPUT==1) page_2->setVariables();   // v.4.34.6
  else          page_2a->setVariables();

  page_2a2->setPage();
  page3->setVariables(p);

  setFileName(file);
  f.close();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::q_clicked()
{
  QDesktopServices::openUrl(QUrl("http://lise.nscl.msu.edu/9_2/9_2_75/9_2_75.pdf"));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::format_error()
{
  QMessageBox MB1;
  MB1.setWindowTitle("Error!");
  MB1.setText("File not in correct format.");
  MB1.setIcon(QMessageBox::Warning);
  MB1.exec();
  return;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::setVariables()
{
  lineEdit[0]->setText(QString::number(_NCASC));
  //inp[_INPUT-1]->setChecked(true);
  inp_combo->setCurrentIndex(_INPUT-1);
  lineEdit[1]->setText(QString::number(_FYRST));
  lineEdit[2]->setText(QString::number(_BARFACr));
  lineEdit[3]->setText(QString::number(_ARATIO));
  lineEdit[4]->setText(QString::number(_FACLA));
  idist[_IDIST]->setChecked(true);
  mdir[_MDIR]->setChecked(true);
  itrac[_ITRAC]->setChecked(true);
  noshl[_NOSHL]->setChecked(true);
  limits[0]->setText(QString::number(_LowLimit));
  limits[1]->setText(QString::number(_HighLimit));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::keyPressEvent(QKeyEvent *e)
{
  // qDebug() << "keyPressEvent";
  if(e->key()==Qt::Key_Enter || e->key()==Qt::Key_Return){
      //  qDebug() << "Enter pressed";
      on_actionExecute_triggered();
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::on_actionSave_triggered()
{
  getVariables();
  QFile qf(FFileName);
  QFileInfo qfinfo(qf.fileName());
  if(qfinfo.fileName().contains(&FileNameAbsent[1]))
    { emit on_actionSave_As_triggered(); return;}

  f10=mfopen(FFileName,"wt");
  if(!f10) {
      qDebug() << "NOPE!";
      return;
    }

  Modified = false;

  fprintf(f10,"%5d%5d%5d%5.0f%5.0f%5.0f%5d%5.0f%5d%5d%5d\n",
          _NCASC,_IDIST,_INPUT,_FYRST, _BARFACr, _ARATIO, _MDIR, _FACLA, _ITRAC, _NOSHL, p->ParticleFlags.GetParticleFlags());

  // level2

  if( abs(_INPUT)==1) {
      /*
            1400  READ (10,1100)IZP,IAP,IZT,IAT,SP,ST,QCN
            1100 FORMAT(4I5,2F5.1,2F10.5)
            READ(10,1101)ELAB,EXPSIG,JCMAX,AGRAZ,ELOSS,
            1101 FORMAT(2F10.5,I5,2F5.0,I5)
           */
      page_2->getVariables();
      fprintf(f10,"%5d%5d%5d%5d%5.1f%5.1f%10.5f\n",
              _IZP, _IAP, _IZT, _IAT, _SPr, _ST, _QCNr);

      fprintf(f10, "%10.5f%10.5f%5d%5.0f%5.0f%5d %1d %.2f\n",
              _ELAB, _EXPSIGr, _JCMAXr, _AGRAZr, _ELOSS, _LMINN, _BarrierMode, _h_omega);
    }
  else {
      page_2a->getVariables();
      // READ(10,1402)IZ,IA,EEXCN,EREC,AJNUC ,LMINN

      fprintf(f10,"%5d%5d%10.4f%10.4f%10.4f%5d\n",
              _IZC, _IAC, _EEXCN, _EREC, _AJNUC, _LMINN);
    }

  //  WWWWWWWWWWWWWWWWWWWWWWWW-- level 2a INPUT=3
  if(_INPUT==3){
      // (I5,5X,7F10.5/(8F10.5));
      fprintf(f10,"%5d     ",int(_AJNUC));

      for(int i=1; i<=_AJNUC; i++)
        fprintf(f10,"%10.5f",_ALEVB[i]);

      fprintf(f10,"\n");
    }
  //  WWWWWWWWWWWWWWWWWWWWWWWW-- level 3   opt potentials
  char optic_f[70]="%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%5d%3d%3d";
  int I;

  if(_INPUT==1){
      I=4;
      if(p->_ITAKE[I])
        fprintf(f10,optic_f,_V0D[I],_V01D[I],_V02D[I],_R0RD[I],_ARD[I],_R0CD[I],
                _W0D[I],_W01D[I],_W02D[I],_R0ID[I],_AID[I],_RMCHD[I],
                _NPD[I],_IMAG[I],_IRAD[I]);
      fprintf(f10,"\n");
    }

  for(I=1; I<=3; I++) {
      if(p->_ITAKE[I])
        fprintf(f10,optic_f,_V0D[I],_V01D[I],_V02D[I],_R0RD[I],_ARD[I],_R0CD[I],
                _W0D[I],_W01D[I],_W02D[I],_R0ID[I],_AID[I],_RMCHD[I],
                _NPD[I],_IMAG[I],_IRAD[I]);
      fprintf(f10,"\n");
    }

  fprintf(f10,"\n");
  fprintf(f10,"\n");
  fclose(f10);

  // ameame(file);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::on_actionExit_triggered(){    qApp->quit();}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::on_actionSave_As_triggered()
{
  QString file = QFileDialog::getSaveFileName(this, tr("Save File"),
                                              FFileName,
                                              tr("PACE4 files (*.pace);;All files (*.*)"));
  if(file.size()>0)
    {
      FFileName = file;
      on_actionSave_triggered();
      setFileName(FFileName);
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::on_actionAbout_triggered()
{
  About *about_page = new About;
  about_page->setWindowFlags(Qt::CustomizeWindowHint |
                             Qt::WindowTitleHint | Qt::WindowCloseButtonHint);
  about_page->show();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
