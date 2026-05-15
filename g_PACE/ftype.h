#define PACE_version "Version 4.34.14"
#define PACE_date    "Last revision 16-AUG-2023"

#define Max_NCASC 1000000 //32767
#define Max_EBIN  4000
#define Max_MOM   200
#define DIM_EBIN  (Max_EBIN+1)
#define DIM_RLEV  99999
//#define N_fusion_event 6000 




// 4.17       10-08-2010   Jenna's corrections
// 4.18       06-03-2011   Detail output for particles
// 4.18.2     13-03-2011   New Random generator
// 4.18.3     14-03-2011   Output or gamma
// 4.18.4     17-03-2011   Gate for detail output
// 4.19.1     13-11-2011   Batch mode for energy
// 4.19.2     13-11-2011   FACLA, ALIT - avoid crashes for light guys at FACLA=0
// 4.19.3     20-11-2011   Batch mode for energy in compound case

// 4.20.      23-05-2013   QM -- passed also Classic L-value procedure.
//                           AGRAZ - default 2
// 4.21 bug with A = N+Z  where Z can be negative . new one is  A = N + abs(Z)
// 4.22 Michelle Qt
// 4.23 new about
// 4.24 Oleg Qt
// 4.31 -- jump to LISE-Qt-utilities version (4.30)
// 4.31.2 Correction in tfusion2 for value input
// 4.31.3 03/19/2021  labels modification in tfusion2
// 4.31.4 03/19/2021  modifications with lise2016 path
// 4.31.5 04/18/2021  fixed issues with edit cells for Elab and Qg, fixed reading file
// 4.31.6 04/18/2021  windows state change, layout of tfusion2, and other cosmetic issues
// 4.31.7 04/18/2021  modifications with lise2016 path : backup version with bathPATH2 from root

// 4.31.8  04/25/2021
// revision of paths, creation of arg-file reading utiltiy

// 4.31.9  04/25/2021
// revision of open and save dialogs

// 4.31.10  04/26/2021
// LISErootPATH = QCoreApplication::applicationDirPath();

//============================================================================
// 4.31.11  05/21/21
// corrections in HighDpiScaling attribute at start

//============================================================================
// 4.32.1  06/08/21
// Bug with Data filename was corrected

//============================================================================
// 4.32.2  06/09/21
// filename corrections in the case of net names //****.***.**/**
//============================================================================
// 4.32.3  06/10/21  adaptation of global for Qt 6 -- Oleg
//============================================================================
// 4.32.4  07/12/21  corrrection for lise.ini path
//============================================================================
// 4.32.5  12/03/21  modified for non-latin ame2016 path
//============================================================================
// 4.32.6  12/04/21  modified for const QString &
//============================================================================
// 4.32.7  12/29/21  compatibility with linux
//============================================================================
// 4.32.8  07/13/22  migrating to Qt63
//============================================================================
// 4.32.9  09/07/22  Fixed: bug with mapping in Qt63  (mapped(int) --> mappedInt(int))
//============================================================================
// 4.33.1  10/19/22  migrating to Qt64
// 4.33.2  10/19/22  QString::number to qbuf.number

//============================================================================
// 4.34.1    04/28/23
// Correction for ARATIO (int -> double)
// cleaning declaration of functions which moved to class PACE

//============================================================================
// 4.34.2    04/28/23
// modification for date/time  in output

//============================================================================
// 4.34.5    05/15/23
// corrections in tfusion_2a2 layout

//============================================================================
// 4.34.6    05/15/23
// corrections for file read in INPUT > 1

//============================================================================
// 4.34.7    05/15/23
// set current page==0 after file reading
//============================================================================
// 4.34.8    06/18/23 DK
// Changed the AME database from a text file to an ACCESS .accdb file
// Changed lise_mass.cpp to use sql to real the desired masses
//============================================================================
// 4.34.9    07/01/23
// compiled with MSVC
//============================================================================
// 4.34.10    07/04/23
// modifications in LISE mass
//============================================================================
// 4.34.11    07/04/23
// no dbf c-file in the project
//============================================================================
// 4.34.12    07/25/23
// moving to sqlite format
//============================================================================
// 4.34.13    08/16/23
// Bug: Corrections with database path
//============================================================================
// 4.34.14    08/16/23
// Bug: Corrections with AZ-gate (DB request)


//============================================================================

//
//
//     INPUT = 1  - PROJECTILE + TARGET
//           = 2  - COMP. NUC. , EX AND J FIXED
//           = 3  - COMP. NUC. , EX FIXED, J DISTRIBUTED
//           = 4  - COMP. NUC. ,J DISTIBUTION CALCULATED INTERNALLY
//   THE SPIN DISTRIBUTION OF THE STARTING EXCITED NUCLIDE HAS BEEN
//          DETERMINED.  THE DETERMINATION OF THE PROBABILITIES OF
//          THE PATHS OF DEEXCITATION WILL BE DETERMINED IN THE FOLLOWING
//          WAY
//         1) USING RANDOM NUMBERS, NCASC CHANNELS ARE ESTABLISHED FOR THE COMPOUND NUCLEUS.
//         2) THE HIGHEST A IN THE ARRAY IS DETERMINED.
//         3) THE FIRST NUCLEUS OF THIS A IN THE ARRAY IS DETERMINED.
//         4) THE HIGHEST ENERGY FOR THIS NUCLEUS IS DETERMINED.
//         5) THE FIRST SPIN OF THIS NUCLEUS AND ENERGY IS FOUND.
//         6) THE PROBABILITY MATRIX FOR THIS NUCLEUS ENERGY AND SPIN IS DETERMINED.
//         7) BY MEANS OF RANDOM NUMBERS, THE DEEXCITATION PATH IS CHOSEN
//            FOR ALL MEMBERS OF THE ARRAY OF THAT NUCLIDE, ENERGY AND
//            SPIN.  THE NEW Z, N, ENERGY, AND SPIN ARE SUBSTITUTED FOR
//            THE OLD VALUES IN THE ARRAY.
//         8) THE NEXT SPIN FOR THE (OLD) NUCLIDE AND ENERGY IS LOCATED.
//            IF IT EXISTS, GO TO STEP 7.
//         9) THE NEXT HIGHER ENERGY FOR THAT (OLD) NUCLIDE IS LOCATED.
//            IF IT EXISTS, GOT TO STEP 5.
//        10) THE NEXT NUCLIDE FOR THAT (OLD) A IS LOCATED.
//            IF IT EXISTS, GO TO STEP 4.
//        11) GO TO STEP 2.
//
//
//   SPROB   IS UNNORMALIZED SUM OF PROBABILITIES
//   PROB    (EQUIVALENT - PROB) IS UNNORMALIZED PROBABILITIES
//   MODES   IS NUMBER OF DIFFERENT PARTICLES (AND GAMMAS)
//           CONSIDERED FOR EVAPORATION
//   IZPART  AND  INPART  ARE Z AND N OF EMITTED PARTICLES
//     THE VECTOR ASSOCIATED WITH EACH EVENT IS AS FOLLOWS
//
//     ERGCS(J)- THE EXCITATION ENERGY
//     IZCS(J) - THE ATOMIC NUMBER ( Z )
//     INCS(J) - THE NUMBER OF NEUTRONS ( N )
//     JCS(J)  - THE SPIN INDEX  I.E.  1,2,3,   ETC.
//               SPIN = SPIN INDEX - 1  FOR EVEN A
//               SPIN = SPIN INDEX - 1/2 FOR ODD A
//     MJCS(J) - THE SPIN PROJECTION  =  0,1,2,   SPIN.
//              WE NEGLEGT 1/2 INTEGER VALUES FOR M *******
//  *** NOTE ***  THESE INTEGERS ARE   INTEGER*2  TO SAVE SPACE
//
//
//   SPIN DISTRIBUTION OF INITIAL NUCLIDE NORMALIZED TO SUM=1.
//     FOR J.GT.100 THERE IS NO COMPOUND NUCLEUS IN PRACTICE.
//     ALEV  IS PARTIAL XSECTION FOR SPIN INDEX I
//           (FROM  COMPOSE  THRU COMMON)


