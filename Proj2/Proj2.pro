TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CFLAGS += -march=native
QMAKE_CXXFLAGS += -march=native
QMAKE_CXXFLAGS_RELEASE += -march=native

SOURCES += \
    main.cpp \
    maxoffdiag.cpp \
    jacobi_rotate.cpp \
    analytic_eigvals.cpp \
    jacobi_method.cpp \
    armadillo_eigpair.cpp \
    test.cpp

HEADERS += catch.hpp \
    maxoffdiag.h \
    jacobi_rotate.h \
    analytic_eigvals.h \
    jacobi_method.h \
    armadillo_eigpair.h \
    test.h


INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -llapack -lblas -larmadillo
