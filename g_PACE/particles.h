#ifndef PARTICLES_H
#define PARTICLES_H

#include <QDialog>
#include <QString>
#include <QFile>
#include <QPlainTextEdit>
class particles : public QDialog
{
    Q_OBJECT
public:
    explicit particles(QString part_filename, QWidget *parent = 0);

public slots:
    void save_clicked();
    void print_clicked();

private:
     QPlainTextEdit *text;
     QFile *partFile;
     QString content;

};

#endif // PARTICLES_H
