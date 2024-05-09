#ifndef TFUSION_2A2_H
#define TFUSION_2A2_H

#include <QWidget>
#include <QLineEdit>
#include <QPushButton>
#include <QSignalMapper>

class tfusion_2a2 : public QWidget
{
    Q_OBJECT

public:
    explicit tfusion_2a2(QWidget *parent = 0);
    void setPage();

signals:

public slots:
    void up_clicked();
    void down_clicked();
    void edited(int);

private:
    QLineEdit *edit[4][2];
    QPushButton *up;
    QPushButton *down;
    int Current;
    QSignalMapper *signalMapper;
};

#endif // TFUSION_2A2_H
