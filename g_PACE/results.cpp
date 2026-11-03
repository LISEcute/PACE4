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
#include <QDesktopServices>


extern QString localPATH;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
results::results(const QString &filename, QWidget *parent) :
    QDialog(parent)
{
    plotData = pacePlotDataSnapshot();
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

    text = new QTextBrowser;
    text->setReadOnly(true);
    text->setOpenLinks(false);
    text->setOpenExternalLinks(false);
    text->setHtml(content);
    connect(text, &QTextBrowser::anchorClicked, this, &results::link_clicked);
    layout->addWidget(text);
    this->setMinimumSize(1000,800);
    layout->setContentsMargins(0,0,0,0); //Qt6
    layout->setSpacing(0);
    setLayout(layout);

}

void results::link_clicked(const QUrl &url)
{
    if (url.scheme() != "pace")
    {
        QDesktopServices::openUrl(url);
        return;
    }

    if (url.host() == "plot_cm_spectra")
    {
        openPaceCMSpectraPlotWindow(this, plotData);
        return;
    }

    PaceAngularPlotKind kind = PacePlotEnergyVsTheta;
    if (url.host().endsWith("ntheta")) kind = PacePlotCountsVsTheta;
    else if (url.host().endsWith("ne")) kind = PacePlotCountsVsEnergy;
    else if (url.host().endsWith("dcsdtheta")) kind = PacePlotCrossSectionVsTheta;

    if (url.host().startsWith("plot_all"))
    {
        openPaceAngularPlotWindow(this, plotData, true, -1, kind);
        return;
    }

    if (url.host() == "plot_table" || url.host() == "plot_ntheta" ||
        url.host() == "plot_ne" || url.host() == "plot_dcsdtheta")
    {
        bool ok = false;
        const int index = url.path().mid(1).toInt(&ok);
        if (ok) openPaceAngularPlotWindow(this, plotData, false, index, kind);
    }
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
