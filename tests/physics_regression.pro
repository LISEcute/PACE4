include(../PACE.pro)
TARGET = pace_physics_regression
CONFIG += console
CONFIG -= app_bundle
DESTDIR = $$OUT_PWD/bin
SOURCES -= g_PACE/p_main.cpp
PACE_ROOT = $$clean_path($$PWD/..)
for(source, SOURCES): ABS_SOURCES += $$PACE_ROOT/$$source
SOURCES = $$ABS_SOURCES $$PWD/physics_regression.cpp
for(header, HEADERS): ABS_HEADERS += $$PACE_ROOT/$$header
HEADERS = $$ABS_HEADERS
for(form, FORMS): ABS_FORMS += $$PACE_ROOT/$$form
FORMS = $$ABS_FORMS
for(resource, RESOURCES): ABS_RESOURCES += $$PACE_ROOT/$$resource
RESOURCES = $$ABS_RESOURCES
RC_ICONS = $$PACE_ROOT/g_PACE/icons/pace4.ico
INCLUDEPATH += $$PACE_ROOT $$PACE_ROOT/g_PACE
