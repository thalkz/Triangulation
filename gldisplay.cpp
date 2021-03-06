#include "gldisplay.h"

#define FRUSTUM_SIZE 1.0f

GLDisplay::GLDisplay(QWidget *parent) :
    QGLWidget(parent),
    _angle(0.0f)
{
}

void GLDisplay::initializeGL()
{
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_CULL_FACE);

    glFrontFace(GL_CCW);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    glColor3f(1.0, 1.0, 0.0);
}

void GLDisplay::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    glRotatef(-45.0f, 1.0f, 0.0f, 0.0f);
    glRotatef(_angle, 0.0f, 0.0f, 1.0f);

    mesh.draw();

    if (voronoi)
    {
        mesh.draw_voronoi();
    }
}

void GLDisplay::resizeGL(int w, int h)
{
    glMatrixMode(GL_PROJECTION);

    glViewport(0, 0, w, h);

    glOrtho(-FRUSTUM_SIZE, FRUSTUM_SIZE,
            -FRUSTUM_SIZE, FRUSTUM_SIZE,
            -FRUSTUM_SIZE, FRUSTUM_SIZE);

    glMatrixMode(GL_MODELVIEW);
}

void GLDisplay::mouseMoveEvent(QMouseEvent *event)
{
    if( event != NULL ) {
        QPoint position = event->pos();

        _angle += (position.x() - _position.x());

        _position = position;

        updateGL();
    }
}

void GLDisplay::mousePressEvent(QMouseEvent *event)
{
    if( event != NULL )
        _position = event->pos();
}

void GLDisplay::onAddVertexButton()
{
    mesh.add_random_vertex();
    updateGL();
}

void GLDisplay::onOptimizeButton()
{
    //mesh.flip_algorithm();
    mesh.lawson_algorithm();
    updateGL();
}

void GLDisplay::onResetButton()
{
    mesh.build_mesh();
    updateGL();
}

void GLDisplay::onDelaunayButton()
{
    mesh.add_delaunay_vertex();
    updateGL();
}

void GLDisplay::onShowVoronoi()
{
    voronoi = true;
    updateGL();
}

void GLDisplay::onHideVoronoi()
{
    voronoi = false;
    updateGL();
}
