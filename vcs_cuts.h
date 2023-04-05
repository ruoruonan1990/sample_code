#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TF1.h"

double m_e = 0.0005109989461;
double m_p = 0.9382720813;
double m_pi0 = 0.1349766;
double type_cut = 4;
double maxtime = -5.;
double mintime = -21.;
//double Ngen = 3000000.;
//double Ngen = 2000000.;
double Ngen = 500000.;
double Ntried = 139452081.;

//----- pion

double hms_dp_cut=8.;//7
double shms_dp_min=-8.;//6
double shms_dp_max=4.;//6
double shms_dp_cut=10.;//6
double hms_theta_cut=0.1;//.07
double hms_phi_cut=0.04;
double shms_theta_cut=0.06;//0.04
double shms_phi_cut=0.04;
double entrth = 0.232126;
double exitth = 0.295711;
double thickrate = 0.263911;
double stheta_off = 0.;
double htheta_off = 0.;
double sphi_off = 0.;
double hphi_off = 0.;
double hchisq_eff = 0.;
double pchisq_eff = 0.;
double dumhchisq_eff = 0.;
double dumpchisq_eff = 0.;
double acchchisq_eff = 0.;
double accpchisq_eff = 0.;
double accdumhchisq_eff = 0.;
double accdumpchisq_eff = 0.;
double hchisq_fac = 0.;
double pchisq_fac = 0.;
double dumhchisq_fac = 0.;
double dumpchisq_fac = 0.;
double acchchisq_fac = 0.;
double accpchisq_fac = 0.;
double accdumhchisq_fac = 0.;
double accdumpchisq_fac = 0.;
double datachi_fac = 0.;
double dumchi_fac = 0.;
double dataaccchi_fac = 0.;
double dumaccchi_fac = 0.;


/*
double hms_dp_cut=5.;//7
double shms_dp_min=-8.;//6
double shms_dp_max=4.;//6
double shms_dp_cut=2.;//6
double hms_theta_cut=0.04;//.07
double hms_phi_cut=0.015;
double shms_theta_cut=0.025;//0.04
double shms_phi_cut=0.015;*/

double cointime_cut = 3.;
double mm_min = -0.01;
double mm_max = 0.01;


//----- vcs
/*double hms_dp_cut=8.;//7
double shms_dp_cut=8.;//6
double hms_theta_cut=0.1;//.07
double shms_theta_cut=0.06;//0.04
double hms_phi_cut=0.04;
double shms_phi_cut=0.04;
double cointime_cut = 2.5;
double mm_min = -0.005;
double mm_max = 0.005;
*/
int t = 0;
int f = 0;
int k = 0;
int b = 0;
double maxcut=0.;
double mincut=0.;
double factor=0.;
double dumfactor=0.;
double fac=0.;
double Data_eff = 0.;
double Dum_eff = 0.;
double beam_E = 0.;
double hms_p = 0.;
double hms_pp = 0.;
double hms_deg = 0.;
double shms_p = 0.;
double shms_pp = 0.;
double shms_deg = 0.;
double normfac = 0.;
double Exp_charge = 0.;
double Exp_dumcharge = 0.;
double dataps6 = 0.;
double dumps6 = 0.;
double timeshift = 0.;
double chargerate = 0.;
double theta_gg_rad = 0.;
double theta_gg = 0.;
double phi_gg_rad = 0.;
double phi_gg = 0.;
double dumtheta_gg_rad = 0.;
double dumtheta_gg = 0.;
double dumphi_gg_rad = 0.;
double dumphi_gg = 0.;
double simtheta_gg_rad = 0.;
double simtheta_gg = 0.;
double simphi_gg_rad = 0.;
double simphi_gg = 0.;
double theggmin = 0.;
double theggmax = 0.;
double h2_entr = 0.15; // mm
double h2_exit = 0.191; // mm
double dum_entr = 0.1815; // g/cm^2
double dum_exit = 0.1816; // g/cm^2
double dum_density = 2.699; // g/cm^3
double thick_entr = (h2_entr + h2_exit) / ((dum_entr + dum_entr)/dum_density * 10.);
double thick_exit = (h2_entr + h2_exit) / ((dum_exit + dum_exit)/dum_density * 10.);


// define LorentzVecotors
TLorentzVector k_in;
TLorentzVector k_out;
TLorentzVector p_in;
TLorentzVector q_in;
TLorentzVector p_out;
TLorentzVector q_out;
TLorentzVector delta;
TLorentzVector missing_mass;

TLorentzVector dumk_in;
TLorentzVector dumk_out;
TLorentzVector dump_in;
TLorentzVector dumq_in;
TLorentzVector dump_out;
TLorentzVector dumq_out;
TLorentzVector dumdelta;
TLorentzVector dummissing_mass;

TLorentzVector simk_in;
TLorentzVector simk_out;
TLorentzVector simp_in;
TLorentzVector simq_in;
TLorentzVector simp_out;
TLorentzVector simq_out;
TLorentzVector simdelta;
TLorentzVector simmissing_mass;

TLorentzVector gimk_in;
TLorentzVector gimk_out;
TLorentzVector gimp_in;
TLorentzVector gimq_in;
TLorentzVector gimp_out;
TLorentzVector gimq_out;
TLorentzVector gimdelta;
TLorentzVector gimmissing_mass;



double cotime    ;  
double hms_dp    ;
double hms_phi   ;
double hms_theta ;
double hytar     ;
double hztar     ;
double shms_dp   ;
double shms_phi  ;
double shms_theta;
double sytar     ;
double sztar     ;
double W         ;
double Q2        ;
double datafac   ;
double dataaccfac;
double missmass2 ;
double phigg     ;
double thegg     ;
double hxfp ;
double hxpfp;
double hyfp ;
double hypfp;
double sxfp ;
double sxpfp;
double syfp ;
double sypfp;
double hxver;
double hyver;
double pxver;
double pyver;
double hchisq;
double pchisq;;
double pntrack	 ;
double hntrack	 ;
double phodstart ;
double hhodstart ;

double dumcotime  ;
double dumh_dp    ;
double dumh_phi   ;
double dumh_theta ;
double dumhytar   ;
double dumhztar   ;
double dums_dp    ;
double dums_phi   ;
double dums_theta ;
double dumsytar   ;
double dumsztar   ;
double dumW       ;
double dumQ2      ;
double dumfac     ;
double dumaccfac  ;
double dummiss2   ;
double dumphigg   ;
double dumthegg   ;
double dumhxfp	;	
double dumhxpfp ;
double dumhyfp  ;
double dumhypfp ;
double dumsxfp  ;
double dumsxpfp ;
double dumsyfp  ;
double dumsypfp ;
double dumhxver;
double dumhyver;
double dumpxver;
double dumpyver;
double dumhchisq;
double dumpchisq;
double dumhntrack ;
double dumpntrack ;
double dumhhodstart;
double dumphodstart;

double shdp	   ;
double shphi   ;
double shtheta ;
double shytar  ;
double shztar  ;
double ssdp	   ;
double ssphi   ;
double sstheta ;
double ssytar  ;
double ssztar  ;
double w	   ;
double simm2   ;
double simW	   ;
double simQ2   ;
double simcfac ;
double simphigg;
double simthegg;
double shxfp ; 
double shxpfp;
double shyfp ;
double shypfp;
double ssxfp ;
double ssxpfp;
double ssyfp ;
double ssypfp;

double sphdp    ;
double sphphi   ;
double sphtheta ;
double sphytar  ;
double sphztar  ;
double spsdp    ;
double spsphi   ;
double spstheta ;
double spsytar  ;
double spsztar  ;
double wp	    ;
double spmm2    ;
double spmW	    ;
double spmQ2    ;
double spmfac   ;
double spmphigg ;
double spmthegg ;
double sphxfp  ;              
double sphxpfp ;
double sphyfp  ;
double sphypfp ;
double spsxfp  ;
double spsxpfp ;
double spsyfp  ;
double spsypfp ;
