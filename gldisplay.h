#ifndef GLDISPLAY_H
#define GLDISPLAY_H

#include <QGLWidget>
#include <mesh.h>
#include <QMouseEvent>

class GLDisplay : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLDisplay(QWidget *parent = 0);

    virtual void initializeGL();

    virtual void paintGL();

    virtual void resizeGL(int w, int h);

protected:
    virtual void mouseMoveEvent ( QMouseEvent * event );
    virtual void mousePressEvent ( QMouseEvent * event );

private:
    void drawSierpinski();
    Mesh mesh;

    float _angle;
    QPoint _position;
    bool voronoi = false;

signals:

public slots:
    void onAddVertexButton();
    void onOptimizeButton();
    void onResetButton();
    void onDelaunayButton();
    void onShowVoronoi();
    void onHideVoronoi();
};

#endif // GLDISPLAY_H
