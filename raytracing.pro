TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -O3 -fopenmp
LIBS += -lgomp

SOURCES += main.cpp

HEADERS += \
    scene.hpp \
    object.hpp \
    segment.hpp \
    operation.hpp \
    light.hpp \
    color.hpp \
    constants.hpp
