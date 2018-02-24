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
    connect(ui->delaunayButton, SIGNAL(released()), ui->canvasWidget, SLOT(onDelaunayButton()));

    connect(ui->voronoiButton, SIGNAL(released()), ui->canvasWidget, SLOT(onShowVoronoi()));
    connect(ui->hideVoronoiButton, SIGNAL(released()), ui->canvasWidget, SLOT(onHideVoronoi()));
}

MainWindow::~MainWindow()
{
    delete ui;
}
