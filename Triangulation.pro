#-------------------------------------------------
#
# Project created by QtCreator 2013-09-12T03:27:41
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Triangulation
TEMPLATE = app


SOURCES += main.cpp\
		mainwindow.cpp \
	gldisplay.cpp \
	mesh.cpp \
	triangle.cpp \
	vertex.cpp \
	colors.cpp

HEADERS  += mainwindow.h \
	gldisplay.h \
	mesh.h \
	triangle.h \
	vertex.h \
	colors.h

FORMS    += mainwindow.ui
