QT = core
CONFIG += console c++17
CONFIG -= app_bundle
TEMPLATE = app
TARGET = random_regression
INCLUDEPATH += ../g_PACE
SOURCES += random_regression.cpp ../g_PACE/util.cpp
HEADERS += ../g_PACE/pace_random.h
