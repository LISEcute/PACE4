#include "results.h"
#include <QBoxLayout>
#include <QFile>
#include <QTextStream>
#include <QMenuBar>
#include <QFileDialog>
#include <QStandardPaths>
#include <QPrinter>
#include <QPrintDialog>
#include <QToolBar>
#include <QIcon>
#include <QMainWindow>


extern QString localPATH;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
results::results(const QString &filename, QWidget *parent) :
    QDialog(parent)
{
   // QIcon ico(":icons/pace.png");
    this->setWindowIcon(QIcon(":/pace4.png"));
    //Qt::WindowFlags flags = windowFlags();
    this->setWindowFlags(Qt::Window
                         | Qt::WindowMinimizeButtonHint
                         | Qt::WindowMaximizeButtonHint
                         | Qt::WindowCloseButtonHint);
    //this->setModal(false);

    QToolBar *toolbar = new QToolBar;

    toolbar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    QAction *saveAction = new QAction(QIcon(":/save.png"),"&Save",this);
    toolbar->addAction(saveAction);
    connect(saveAction, SIGNAL(triggered()),this, SLOT(save_clicked()));
    saveAction->setIcon(QIcon(":/save.png"));
    saveAction->setIconText("Save");
    QAction *printAction = new QAction(tr("&Print"), this);
    printAction->setIcon(QIcon(":/print.gif"));
    printAction->setIconText("Print");
    toolbar->addAction(printAction);
    connect(printAction, SIGNAL(triggered()),this, SLOT(print_clicked()));

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(toolbar);

    resFile= new QFile(filename);
    resFile->open(QIODevice::ReadOnly);
    QTextStream stream(resFile);
    QString content = stream.readAll();
    resFile->close();

    text = new QTextEdit(content);
    text->setReadOnly(true);
    layout->addWidget(text);
    this->setMinimumSize(1000,800);
    layout->setContentsMargins(0,0,0,0); //Qt6
    layout->setSpacing(0);
    setLayout(layout);

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void results::save_clicked()
{
QString file = QFileDialog::getSaveFileName(this, tr("Save File"), localPATH, tr("PACE result files (*.html)"));

if(file.size()>0)
        resFile->copy(file);

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void results::print_clicked()
{
     QPrinter printer;
     QPrintDialog *printDialog = new QPrintDialog(&printer, this);
     printDialog->setWindowTitle(tr("Print Results File"));

     if (printDialog->exec() != QDialog::Accepted)
         return;

     QString htmlToPrint(text->toHtml());
     printer.setFullPage(true);
     QTextDocument textDoc;
     textDoc.setHtml(htmlToPrint);
     textDoc.print(&printer);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
