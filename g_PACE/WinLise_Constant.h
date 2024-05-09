#if !defined(constant_h)
	#define  constant_h


# define Number_Calc          25        //TO    
# define Number_CrossSections 1000       //TO
# define Number_tab 		8000      //TO
# define NRD	 		25        //TO     number of result's distribution  12-d6       +4
//# define NP				128
# define NumberMaterial		7
# define Number_EPAX		3
# define Number_Fusion		1
# define NumberTabELOSS		10
# define Number_ELOSS		4
# define Number_ELOSS_PLOT	(Number_ELOSS+2)
# define Number_RangePoints  217
# define Number_StragEne	2
# define Number_StragAng	2
# define Number_LISE_web	2
# define Number_LISE_ftp	2
# define Number_time_new_versions	8     // Ntnv
# define Number_reactions	6
# define Number_AxisElementModes 4


#define FlagLogTrue    true
#define FlagLogFalse   false
#define LengthNameFile _MAX_PATH


# define Number_Charges		3
# define Number_Proton		130
# define Number_Neutron		200
# define Number_NP		7		//16,32,64,128,256,512,1024
# define NumberPlaceEnergyCalculation 12
# define NumberPlaceTOFCalculation (2+NumberMaterial)
#define  NumberMethodSeparation 8

#define  max_Matter_into_Compound 5
#define  DecayString           11
#define  LengthNameCompound (max_Matter_into_Compound*4+1)
#define  LengthNameCompoundPlusMass (LengthNameCompound+17)

#define  N_sigma_corr		   2
#define  N_sigma_distr		   4
#define  N_velocity_distr	   6
#define  N_mode_separation	   3

// secondary reactions
#define  n_draw_SecReact	   5     // and +1 for thickness target
#define  coef_decrease_SecReact  4


//# define StrLenFileName		50


//   PHYS


# define PI				3.14159265358979
# define PI_L			3.1415926535897932384626433832795L
# define aem_MeV			931.494013  //TO  MeV
# define vc_cm 			29.99792458    //TO  cm/ns
# define vc_m 			0.2999792458   //TO  m/ns
# define vc_m_s 			2.999792458e+8 //TO  m/s
#define RAD 	  		0.017453292519943		// (PI/180.)
#define RAD_L	 	  	0.0174532925199432957692369076848861L
#define GRAD		  	57.2957795130823	// (180./PI)
#define GRAD_L			57.2957795130823208767981548141052L
#define aem_mg          (1.66053873e-21)   // mg
#define joule           1.             //system value of Energy joule=1
#define eV              (1.602176462e-19*joule)
#define MeV             (1.e6*eV)
#define cm              1.             //system value of length cm=1
#define fermi           (1.e-13*cm)
#define hc_MeV_fermi    197.327053
#define hc              (hc_MeV_fermi*MeV*fermi)
#define e2              hc                  // costant of fine structure
#define e2_MeV_fermi    (hc_MeV_fermi/137.)   // costant of fine structure
#define barn            (1.e-24*cm*cm)
#define mbarn           (1.e-3*barn)
#define U               aem_MeV
#define vc              2.9997924582e10       // light velocity in cm/sec
#define sec             1.             //system value of time sec=1
#define ns              (1.e-9*sec)
#define vcn             (vc*ns)           // light velocity in cm/ns
#define Brho_const	3.107129578


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


//-----------------------------------------------------------
#define  o1_div_32767	(3.05185e-5)
#define Z_Al                     13.
#define Z_U                      92.
#define Z_He                      2.
#define A_Al                     27.
#define A_U                      238.
#define A_6He                       6.


#define calc_index(Z,N)   (Z*1000+N)


#endif
