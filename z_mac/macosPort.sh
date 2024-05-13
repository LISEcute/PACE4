#!/bin/zsh
sed -i'' -e '$a\
DEFINES += __APPLE_
' $1/Charge.pro
sed -i'' -e '$a\
ICON = ./macosIcons/charge.icns
' $1/Charge.pro
sed -i'' -e '$a\
QMAKE_INFO_PLIST = ./Info.plist
' $1/Charge.pro

sed -i'' -e '$a\
DEFINES += __APPLE_
' $1/ETACHA4.pro
sed -i'' -e '$a\
ICON = ./macosIcons/etacha.icns
' $1/ETACHA4.pro
sed -i'' -e '$a\
QMAKE_INFO_PLIST = ./Info.plist
' $1/ETACHA4.pro

sed -i'' -e '$a\
DEFINES += __APPLE_
' $1/Gemini.pro
sed -i'' -e '$a\
ICON = ./macosIcons/gemini.icns
' $1/Gemini.pro
sed -i'' -e '$a\
QMAKE_INFO_PLIST = ./Info.plist
' $1/Gemini.pro

sed -i'' -e '$a\
DEFINES += __APPLE_
' $1/Global.pro
sed -i'' -e '$a\
ICON = ./macosIcons/global.icns
' $1/Global.pro
sed -i'' -e '$a\
QMAKE_INFO_PLIST = ./Info.plist
' $1/Global.pro

sed -i'' -e '$a\
DEFINES += __APPLE_
' $1/Kantele.pro
sed -i'' -e '$a\
ICON = ./macosIcons/kantele.icns
' $1/Kantele.pro
sed -i'' -e '$a\
QMAKE_INFO_PLIST = ./Info.plist
' $1/Kantele.pro

sed -i'' -e '$a\
DEFINES += __APPLE_
' $1/PACE.pro
sed -i'' -e '$a\
ICON = ./macosIcons/pace4.icns
' $1/PACE.pro
sed -i'' -e '$a\
QMAKE_INFO_PLIST = ./Info.plist
' $1/PACE.pro

sed -i'' -e '$a\
DEFINES += __APPLE_
' $1/LISEcute.pro
sed -i'' -e '$a\
ICON = ./macosIcons/lise.icns
' $1/LISEcute.pro
sed -i'' -e '$a\
QMAKE_INFO_PLIST = ./Info.plist
' $1/LISEcute.pro
sed -i'' -e '$a\
installFolder.files = _install
' $1/LISEcute.pro
sed -i'' -e '$a\
installFolder.path = Contents
' $1/LISEcute.pro
sed -i'' -e '$a\
QMAKE_BUNDLE_DATA += installFolder
' $1/LISEcute.pro