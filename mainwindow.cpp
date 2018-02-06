#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->addVertexButton, SIGNAL(released()), ui->canvasWidget, SLOT(onAddVertexButton()));
    connect(ui->optimizeButton, SIGNAL(released()), ui->canvasWidget, SLOT(onOptimizeButton()));
    connect(ui->resetButton, SIGNAL(released()), ui->canvasWidget, SLOT(onResetButton()));
}

MainWindow::~MainWindow()
{
    delete ui;
}
