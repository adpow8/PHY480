TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib

SOURCES += \
    Project_1.cpp \
    ../../../Downloads/project1.cpp \
    project_2.cpp

LIBS += -larmadillo -llapack -lblas
