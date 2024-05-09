#if !defined(constant_h)
	#define  constant_h


# define Number_CrossSections 10000      //TO
# define Number_tab 		20000      //TO
# define Number_tab_MC 		200        //TO
#define  Number_tab_SR	      10000
#define  Number_tab_iso		4000       // isomer
# define Number_EPAX		4
# define Number_Fusion		1
# define NumberTabELOSS		65
# define NumberTabCHARGE      65
# define Number_ELOSS		4
# define Number_ELOSS_PLOT	(Number_ELOSS+2)
# define Number_RangePoints  217
# define Number_StragEne	3
# define Number_StragAng	4
# define Number_LISE_web	1
# define Number_LISE_ftp	1
# define Number_time_new_versions	8     // Ntnv
# define Number_reactions	   9    // 3 - Abrasion-Fissuion
# define Number_SecondaryTargets 3    // 3 - Secondary targets
# define Number_TI 17    // Time interval
# define Number_T12_alpha	   5
# define Number_CatcherMaterials    6 

//# define Number_Button_reactions  Number_reactions
# define Number_AxisElementModes 9       // 5 regular, 3sum+NdZ
# define Number_ProfilesPhys	3
# define Number_Profiles	4
# define Number_MassFormula   3
# define Number_MassDataBase  2
# define Number_FisBar		6
# define Number_T12		4
# define Number_T12_compile	3
# define Number_SecR_filtres  4
# define Number_RayFields	10

#define FlagLogTrue    true
#define FlagLogFalse   false
#define LengthNameFile _MAX_PATH
//#define NumberBLOCKs 		13    // v76
//#define NumberBLOCKs 		14    // v83
#define NumberBLOCKs 		15    // v83.141
#define NumberDIfields		11
#define NumberGanilIsomers    500
#define Number_MCgates		4

# define Number_Charges		6      ///  charge state methods

# define Number_Proton		200
# define Number_Neutron		500
# define Number_Nuclon		700
# define Number_NP		7		//16,32,64,128,256,512,1024
# define MaxZtarget		100
#define  NumberMethodSeparation 11   // 8 method + fission + breakup + gamma&residual
//#define  MC_RB 13                    // MC transmission parameters
//#define  MC_RB 17                    // MC transmission parameters     08/14/07
#define  MC_RB 18                    // MC transmission parameters     03/10/10  + Envelope
//#define  MC_RB_files 17              // MC transmission parameters for ray files    04/07/08
#define  MC_RB_files 19              // MC transmission parameters for ray files     02/28/2010
#define  MC_plot_from_files 12       // MC transmission parameters to plot in File dialog


#define SR_NpickupShift 2		// for fragmentation -- neutron pickup secondary reactions
#define SR_PpickupShift 1           // for fragmentation -- proton pickup secondary reactions


#define  max_Matter_into_Compound 5
#define  DecayString           12
#define  LengthNameCompound (max_Matter_into_Compound*4+1)
#define  LengthNameCompoundPlusMass (LengthNameCompound+17)

#define  N_sigma_corr		   2
#define  N_sigma_distr		   4
#define  N_velocity_distr	   6
#define  N_mode_separation	   3
#define  NumberDriftModes	   4

// secondary reactions
#define  n_draw_SecReact	   5     // and +1 for thickness target
#define  coef_decrease_SecReact  4
#define  Threshold_Alpha 	   6.0



//# define StrLenFileName		50


//   PHYS


# define PI				3.14159265358979
# define PI2			6.28318530717959
# define PI_L			3.1415926535897932384626433832795L
# define aem_MeV			931.49432  //TO  MeV
# define vc_mm 			299.792458    //TO  mm/ns
# define vc_cm 			29.9792458    //TO  cm/ns
# define vc_m 			0.299792458   //TO  m/ns
# define vc_m_s 			2.99792458e+8 //TO  m/s
#define RAD 	  		0.017453292519943		// (PI/180.)
#define RAD_L	 	  	0.0174532925199432957692369076848861L
#define GRAD		  	57.2957795130823	// (180./PI)
#define GRAD_L			57.2957795130823208767981548141052L

#define RADM 	  		17.453292519943		// (PI/180.)
#define RAD_LM	 	  	17.4532925199432957692369076848861L
#define GRADM		  	0.0572957795130823	// (180./PI)
#define GRAD_LM			0.0572957795130823208767981548141052L

#define aem_mg           (1.6605402e-21)   // mg
#define _joule           1.             //system value of Energy joule=1
#define _eV              (1.60217733e-19*_joule)
#define _MeV             (1.e6*_eV)
#define _cm              1.             //system value of length cm=1
#define _fermi           (1.e-13*_cm)
#define _h_bar           6.58211889e-22   // MeV/s
#define _hc_MeV_fermi    197.327053
#define _hc              (_hc_MeV_fermi*_MeV*_fermi)
#define _e2              _hc                  // costant of fine structure
#define FineStructConst  7.297352521e-03       // fine-structure constant
#define e2_MeV_fermi    (_hc_MeV_fermi*FineStructConst)   // costant of fine structure
#define barn            (1.e-24*_cm*_cm)
#define mbarn           (1.e-3*barn)
//#define U               aem_MeV
#define _vc              2.9997924582e10       // light velocity in cm/sec
#define _sec             1.             //system value of time sec=1
#define _ns              (1.e-9*_sec)
#define _vcn             (_vc*_ns)           // light velocity in cm/ns
#define Brho_const	3.107129578
#define Erho_const      9.32764e8   // Joule/coul

//from the code "charge"
#define BohrRadius_cm      0.529177249e-08       // Bohr radius in cm;
#define Rydberg_eV         13.605698100000D00    // Rydberg in eV;
// for Schiwietz
#define Bohr_Velocity    2.19e6  // m/s


#define Mass_p          938.271998          // MeV
#define Mass_n          939.565330722       // MeV
#define Mass_e          0.510998902          // MeV
//#define ME_p            6.777985             // MeV   it is real value
#define ME_p            7.28897              // MeV     it is from database
#define ME_n            8.07132	         // MeV
#define ME_a            2.4249	        // MeV
#define ME_t            14.9498	        // MeV
#define ME_3he          14.9312	        // MeV
#define ME_d            13.1357	        // MeV

#define DontCalculateCSdegrader   -7711
//-----------------------------------------------------------
#define  o1_div_32767	(3.05185e-5)


#define calc_index(Z,N)          (Z*1000+N)
#define calc_indexG(Z,N,G)       ((Z*1000+N)*1000 + int(G))
#define Z_from_index(index)   (index/1000)
#define N_from_index(index)   (index%1000)


#define exp_ln_limit  1e-43 //   exp(-99) =  1.011e-43
#define min_gauss  1e-6

#endif
