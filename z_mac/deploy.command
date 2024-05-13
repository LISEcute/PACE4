#!/bin/zsh
cd /Users/bazin/Documents/Github/LISEcute
/Users/bazin/Qt/6.5.0/macos/bin/macdeployqt ../build-LISE_Package-Qt_6_5_0_for_macOS-Release/LISEcute.app -sign-for-notarization=3B2C1D6D55CF593C716445BC4AD23C4C37DC68A3
#codesign --deep --timestamp --options runtime -s "3B2C1D6D55CF593C716445BC4AD23C4C37DC68A3" LISE++.app
/Users/bazin/Qt/6.5.0/macos/bin/macdeployqt ../build-LISE_Package-Qt_6_5_0_for_macOS-Release/Charge.app -sign-for-notarization=3B2C1D6D55CF593C716445BC4AD23C4C37DC68A3

/Users/bazin/Qt/6.5.0/macos/bin/macdeployqt ../build-LISE_Package-Qt_6_5_0_for_macOS-Release/Gemini.app -sign-for-notarization=3B2C1D6D55CF593C716445BC4AD23C4C37DC68A3

/Users/bazin/Qt/6.5.0/macos/bin/macdeployqt ../build-LISE_Package-Qt_6_5_0_for_macOS-Release/Global.app -sign-for-notarization=3B2C1D6D55CF593C716445BC4AD23C4C37DC68A3

/Users/bazin/Qt/6.5.0/macos/bin/macdeployqt ../build-LISE_Package-Qt_6_5_0_for_macOS-Release/KanteleHandbook.app -sign-for-notarization=3B2C1D6D55CF593C716445BC4AD23C4C37DC68A3

/Users/bazin/Qt/6.5.0/macos/bin/macdeployqt ../build-LISE_Package-Qt_6_5_0_for_macOS-Release/PACE4.app -sign-for-notarization=3B2C1D6D55CF593C716445BC4AD23C4C37DC68A3

/Users/bazin/Qt/6.5.0/macos/bin/macdeployqt ../build-LISE_Package-Qt_6_5_0_for_macOS-Release/ETACHA4.app -sign-for-notarization=3B2C1D6D55CF593C716445BC4AD23C4C37DC68A3
