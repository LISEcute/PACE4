TARGET = PACE4
TEMPLATE = app
CONFIG += c++17
QT       += widgets core gui printsupport sql


UI_DIR = ui
MOC_DIR = moc
RCC_DIR = rcc
OBJECTS_DIR = obj


win32-g++ {
DESTDIR = C:/PACE4/_install
#DESTDIR = ../../_install
}
win32-msvc {
DESTDIR = ../../_install_MSVC
}

macx {
DESTDIR = /Users/arjunray/Desktop/FRIB/PACE4/_install_MSVC/lisecfg/AME_DB.sqlite
}

win32:VERSION = 4.34.15.1 # major.minor.patch.build
else:VERSION = 4.34.15    # major.minor.patch

win32 {
       QMAKE_TARGET_COPYRIGHT = "LISE group at FRIB/MSU"
	QMAKE_TARGET_COMPANY   = "LISE group at FRIB/MSU"
	}

#==================================================
# The following define makes your compiler warn you if you use any
# feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Input

SOURCES += \
    g_PACE/lise_mass.cpp \
    g_PACE/p_fileReadWrite.cpp \
    g_PACE/p_main.cpp\
    g_PACE/tfusion.cpp \
    g_PACE/PACE3main.cpp \
    g_PACE/util.cpp \
    g_PACE/PACE3util.cpp \
    g_PACE/PACE3compound.cpp \
    g_PACE/PACE3Levden.cpp \
    g_PACE/PACE3product.cpp \
    g_PACE/PACE3statis.cpp \
    g_PACE/tfusion_2.cpp \
    g_PACE/tfusion_2a.cpp \
    g_PACE/tfusion_2a2.cpp \
    g_PACE/tfusion3.cpp \
    g_PACE/results.cpp \
    g_PACE/particles.cpp \
    g_PACE/p_about.cpp \
    w_Stuff/w_Label_clickable.cpp \
    w_Stuff/win_utilString.cpp

HEADERS  += \
    g_PACE/pace.h \
    g_PACE/tfusion.h \
    g_PACE/ftype.h \
    g_PACE/Constant.h \
    g_PACE/ParticleStates.h \
    g_PACE/WinLise_Constant.h \
    g_PACE/fisrot.h \
    g_PACE/yrast.h \
    g_PACE/barfit.h \
    g_PACE/tfusion_2.h \
    g_PACE/tfusion_2a.h \
    g_PACE/tfusion_2a2.h \
    g_PACE/tfusion3.h \
    g_PACE/results.h \
    g_PACE/particles.h \
    g_PACE/p_about.h \
    w_Stuff/w_Label_clickable.h


FORMS    += \
	g_PACE/tfusion.ui \
    	g_PACE/p_about.ui

RESOURCES += \
       lise.qrc \
	g_PACE/icons/pace.qrc


RC_ICONS += g_PACE/icons/pace4.ico

# Check if the platform is macOS
macx {
   DEFINES += __APPLE_
   ICON = ./Icons_macos/pace4.icns
   QMAKE_INFO_PLIST = ./Info.plist
}
