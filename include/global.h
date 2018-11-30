#ifndef INCLUDED_GLOBAL_SUNA
#define INCLUDED_GLOBAL_SUNA

const int g_rng_cg            =20;    //Search range for center of gracity[ns]
const int g_num_slice         =1024;  //Total number of drs4 cell
const int g_edg_cut           =10;    //Deviation of the center of gravity[cell]
const float g_chrg_rng        =10;    //Search range for chrage 
const float g_chrg_per_up     =0.9;   //
const float g_chrg_per_dwn    =0.1;   //  
const float g_chrg_per =g_chrg_per_up - g_chrg_per_dwn; //
const float g_echrg = 1.602176e-19;   //Elementary Charge
const float g_impd  = 1.2e3;          //Impedance of PACTA
const TString Referance = "AA3300";   //Reference PMT Serial
float af = 1.0;  //Attennuation factor


#endif
