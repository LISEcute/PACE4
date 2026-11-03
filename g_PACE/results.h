#ifndef RESULTS_H
#define RESULTS_H

#include <QDialog>
#include <QFile>
#include <QTextBrowser>
#include <QUrl>

#include "pace_angular_plots.h"

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class results : public QDialog
{
    Q_OBJECT
public:
    explicit results(const QString &filename, QWidget *parent = 0);

signals:
    //void results_closed();

public slots:
    void save_clicked();
    void print_clicked();
    void link_clicked(const QUrl &url);

private:
    QFile *resFile;
    QTextBrowser *text;
    PacePlotData plotData;
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#endif // RESULTS_H
