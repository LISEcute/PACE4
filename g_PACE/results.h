#ifndef RESULTS_H
#define RESULTS_H

#include <QDialog>
#include <QFile>
#include <QTextEdit>

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

private:
    QFile *resFile;
    QTextEdit *text;
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#endif // RESULTS_H
