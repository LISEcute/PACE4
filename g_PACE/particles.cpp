#include "particles.h"
#include "tfusion.h"
#include <QBoxLayout>
#include <QFileDialog>
#include <QTextStream>
#include <QPlainTextEdit>
#include <QDir>
#include <QMenuBar>
#include <QPrinter>
#include <QPrintDialog>
#include <QToolBar>
#include <QIcon>
//particles::particles(QWidget *parent) :
//    QWidget(parent)
//{
//    QVBoxLayout *layout = new QVBoxLayout;
//    QFile partFile(tfusion::FFileNameParticles);
//    partFile.open(QIODevice::ReadOnly);
//    QTextStream stream(&partFile);
//    QString content = stream.readAll();
//    partFile.close();
//}

extern QString localPATH;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
particles::particles(QString part_filename, QWidget *parent) :
    QDialog(parent)
{
    QToolBar *toolbar = new QToolBar;
    //QMenuBar *menubar = new QMenuBar;
    //QMenu *fileMenu = menubar->addMenu(tr("&File"));
    toolbar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    QAction *saveAction = new QAction(QIcon(":/save.png"),"&Save",this);
    //fileMenu->addAction(saveAction);
   // toolbar->setToolButtonStyle(QToolBar::toolButtonStyle().ToolButtonTextBesideIcon);
    this->setWindowTitle("PACE4: " + part_filename);
    this->setWindowIcon(QIcon(":/pace4.png"));
    this->setModal(false);
    toolbar->addAction(saveAction);
    connect(saveAction, SIGNAL(triggered()),this, SLOT(save_clicked()));
    saveAction->setIcon(QIcon(":/save.png"));
    saveAction->setIconText("Save");
    QAction *printAction = new QAction(tr("&Print"), this);
    printAction->setIcon(QIcon(":/print.gif"));
    printAction->setIconText("Print");
    toolbar->addAction(printAction);
    connect(printAction, SIGNAL(triggered()),this, SLOT(print_clicked()));

    QString filename(part_filename);
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(toolbar);
    partFile = new QFile(filename);
    partFile->open(QIODevice::ReadOnly);
    QTextStream stream(partFile);
    content = stream.readAll();
    partFile->close();

    text = new QPlainTextEdit(content);
    text->setReadOnly(true);
    layout->addWidget(text);
    this->setMinimumSize(1300,600);
    layout->setContentsMargins(0,0,0,0); //Qt6
    layout->setSpacing(0);
    setLayout(layout);
    //this->layout()->setMenuBar(menubar);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void particles::save_clicked()
{
 QString file = QFileDialog::getSaveFileName(this, tr("Save File"), localPATH, tr("PACE4 particle files (*.particles)"));
    if(file.size()>0){
        partFile->copy(file);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void particles::print_clicked()
{

    QPrinter printer;
    QPrintDialog *printDialog = new QPrintDialog(&printer, this);
    printDialog->setWindowTitle(tr("Print Results File"));

    if (printDialog->exec() != QDialog::Accepted)
        return;

    printer.setFullPage(true);
    QTextDocument textDoc;
    textDoc.setPlainText(content);
    textDoc.print(&printer);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
