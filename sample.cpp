#include <iostream>
#include <string.h>
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "vcs_cuts.h"
#include "findPara.h"
#include "vcs_runlist.h"

using namespace std;


void sample(Int_t RunNumber=8916, Int_t dumNumber=9034, string kinType = "kin1a"){

	//	cout << "Enter kin type: kin1a/kin1b..." << endl;
	//	cin >> kinType;
	kinType = "kin1a";
	vector<int> kinvec;
	vector<int> dumvec;
	kinvec = list1a;
	dumvec = list1adummy;
	gStyle->SetTitleOffset(.7,"Y");
	gStyle->SetTitleOffset(.7,"X");
	gStyle->SetLabelSize(0.03,"XY");
	gStyle->SetTitleSize(0.05,"XY");
	TH1::SetDefaultSumw2(kTRUE);

	////////////////////////////////////////////////////find ps6 and kinematics//////////////////////////////////////////

	TObjArray HList(0);
	Functions data;
	Functions dummy;
	Functions simc;
	dataps6 = data.findps6(RunNumber);
	dumps6 = dummy.findps6(dumNumber);
	double* datakin = data.readDbase(RunNumber);
	hms_pp = datakin[4] *0.9987;
	hms_deg = TMath::Abs(datakin[2]);
	shms_pp = datakin[5] *0.9967;
	shms_deg = TMath::Abs(datakin[3]) - 0.017;
	beam_E = datakin[0];
	double* datapara = data.findeff(RunNumber);
	double* dumpara = dummy.findeff(dumNumber);
	Data_eff = datapara[0];
	Dum_eff = dumpara[0];
	normfac = datapara[2];
	timeshift = datapara[1];
	factor = 1./datapara[3];
	dumfactor = 1./dumpara[3];
	double scal_e = shms_pp/shms_p;
	double scal_p = hms_pp/hms_p;

	int theggnum = theggmax-theggmin;
	double theggcm[theggnum];
	for(int i = 0; i < theggnum; i++) {
		theggcm[i] = theggmin + 0.5 + i;

	}

	double theggcut = 0.5;
	double phicut = 70.;
	mm_max = 0.06;	

	cout << "factor is " << factor << endl;
	cout << "normfac is " << normfac << endl;
	cout << "mintime is " << mintime << endl;
	cout << "shms angle is " << shms_deg << endl;
	cout << "hms angle is " << hms_deg << endl;

	////////////////////////////////////////////////////read rootfiles and define variables &  TH1F////////////////////////////////////////

	TString coinFile, dummyFile, simcrootFile, labrootFile, outFile;

	simcrootFile = "/work/hallc/alphaE/ruonanli/simc_michael/worksim/coin_vcs_"+kinType+"_newtarget_muscext.root";
	outFile = Form("%s_vcs_wq2_multi_thephiggcut.root",kinType.c_str());

	TFile *simcFile = new TFile(simcrootFile);

	TChain *ch = new TChain("T");
	TChain *dumch = new TChain("T");

	TTreeReader myReader;
	TTreeReader dumReader;
	myReader.SetTree(ch);
	dumReader.SetTree(dumch);

	for (int i=0; i<kinvec.size(); i++) {  // kin1a 8909

		coinFile = Form("$rootfile/coin_replay_production_%d_-1.root",kinvec[i]);
		TFile *myFile = new TFile(coinFile);
		if (myFile->IsOpen()) {
			cout << "Run " << kinvec[i] << " exists" << endl;
		}else{
			cout << "Run " << kinvec[i] << " does not exist" << endl;
			continue;
		}

		ch->Add(coinFile);
		TTreeReader eReader;
		TChain *chp = new TChain("TSP");
		chp->Add(Form("$rootfile/coin_replay_production_%d_-1.root",kinvec[i]));
		eReader.SetTree(chp);
		TTreeReaderValue<Double_t> charge(eReader, "P.BCM1.scalerChargeCut");
		TTreeReaderValue<Double_t> current(eReader, "P.BCM1.scalerCurrent");
		cout << eReader.GetEntries(1) << endl;
		while(eReader.Next()) {
			if(*current>=5.) {
				charge_cut += (*charge-prev_charge);
			}
			prev_charge = *charge;
		}
		totcharge += charge_cut;
		charge_cut = 0.;
		prev_charge = 0.;
		cout << kinvec[i] << "	charge is	" << totcharge << endl;

	}

	for (int j=0;j<dumvec.size(); j++) {
		dummyFile = Form("$rootfile/coin_replay_production_%d_-1.root",dumvec[j]);
		dumch->Add(dummyFile);
		TFile *dumFile = new TFile(dummyFile);
		if (dumFile->IsOpen()) {
			cout << "Run " << dumvec[j] << " exists" << endl;
		}else{
			cout << "Run " << dumvec[j] << " does not exist" << endl;
			continue;
		}

		TTreeReader dumeReader;
		TChain *dumchp = new TChain("TSP");
		dumchp->Add(Form("$rootfile/coin_replay_production_%d_-1.root",dumvec[j]));
		dumeReader.SetTree(dumchp);
		TTreeReaderValue<Double_t> dumcharge(dumeReader, "P.BCM1.scalerChargeCut");
		TTreeReaderValue<Double_t> dumcurrent(dumeReader, "P.BCM1.scalerCurrent");
		TTreeReaderValue<Double_t> dumScal_evNumber(dumeReader, "evNumber");
		cout << dumeReader.GetEntries(1) << endl;
		while (dumeReader.Next()) {

			if(*dumcurrent>=20.) {
				dumcharge_cut += (*dumcharge - dumprev_charge);

			}
			dumprev_charge = *dumcharge;
		}
		dumtotcharge += dumcharge_cut;
		dumcharge_cut = 0.;
		dumprev_charge = 0.;
		cout << "dummy	" << dumvec[j] << "	  charge is  " << dumtotcharge << endl;
	}

	TTreeReader simcReader("h666",simcFile);

	//  get data values
	TTreeReaderValue<Double_t> hms_dp(myReader, "H.gtr.dp");
	TTreeReaderValue<Double_t> hms_phi(myReader, "H.gtr.yptar");
	TTreeReaderValue<Double_t> hms_theta(myReader, "H.gtr.xptar");
	TTreeReaderValue<Double_t> hytar(myReader, "H.gtr.y");
	TTreeReaderValue<Double_t> shms_dp(myReader, "P.gtr.dp");
	TTreeReaderValue<Double_t> shms_phi(myReader, "P.gtr.yptar");
	TTreeReaderValue<Double_t> shms_theta(myReader, "P.gtr.xptar");
	TTreeReaderValue<Double_t> sytar(myReader, "P.gtr.y");

	TTreeReaderValue<Double_t> hms_px(myReader, "H.gtr.px");
	TTreeReaderValue<Double_t> hms_py(myReader, "H.gtr.py");
	TTreeReaderValue<Double_t> hms_pz(myReader, "H.gtr.pz");
	TTreeReaderValue<Double_t> shms_px(myReader, "P.gtr.px");
	TTreeReaderValue<Double_t> shms_py(myReader, "P.gtr.py");
	TTreeReaderValue<Double_t> shms_pz(myReader, "P.gtr.pz");

	TTreeReaderValue<Double_t> shms_tracknorm(myReader, "P.cal.etottracknorm");
	TTreeReaderValue<Double_t> hms_beta(myReader, "H.gtr.beta");
	TTreeReaderValue<Double_t> cointime(myReader,"CTime.epCoinTime_ROC2");
	TTreeReaderValue<Double_t> t1r2(myReader,"T.coin.pTRIG1_ROC2_tdcTime");
	TTreeReaderValue<Double_t> t4r2(myReader,"T.coin.pTRIG4_ROC2_tdcTime");
	TTreeReaderValue<Double_t> type(myReader,"g.evtyp");
	TTreeReaderValue<Double_t> gevnum(myReader,"g.evnum");

	TTreeReaderValue<Double_t> h1ytdc(myReader, "H.hod.1y.NegTdcRefTime");
	TTreeReaderValue<Double_t> h1ytdcdiff(myReader, "H.hod.1y.NegTdcRefDiffTime");
	TTreeReaderValue<Double_t> h1yadc(myReader, "H.hod.1y.NegAdcRefTime");
	TTreeReaderValue<Double_t> h1yadcdiff(myReader, "H.hod.1y.NegAdcRefDiffTime");
	TTreeReaderValue<Double_t> p1ytdc(myReader, "P.hod.1y.NegTdcRefTime");
	TTreeReaderValue<Double_t> p1ytdcdiff(myReader, "P.hod.1y.NegTdcRefDiffTime");
	TTreeReaderValue<Double_t> p1yadc(myReader, "P.hod.1y.NegAdcRefTime");
	TTreeReaderValue<Double_t> p1yadcdiff(myReader, "P.hod.1y.NegAdcRefDiffTime");
	TTreeReaderValue<Double_t> h1x1dcref(myReader,"H.dc.1x1.RefTime");
	TTreeReaderValue<Double_t> p1x1dcref(myReader,"P.dc.1x1.RefTime");
	TTreeReaderValue<Double_t> hztar(myReader, "H.react.z");
	TTreeReaderValue<Double_t> pztar(myReader,"P.react.z");
	TTreeReaderValue<Double_t> htime(myReader, "H.hod.starttime");
	TTreeReaderValue<Double_t> ptime(myReader, "P.hod.starttime");

	// get dummy values
	TTreeReaderValue<Double_t> dumh_dp(dumReader, "H.gtr.dp");
	TTreeReaderValue<Double_t> dumh_phi(dumReader, "H.gtr.yptar");
	TTreeReaderValue<Double_t> dumh_theta(dumReader, "H.gtr.xptar");
	TTreeReaderValue<Double_t> dumhytar(dumReader, "H.gtr.y");
	TTreeReaderValue<Double_t> dums_dp(dumReader, "P.gtr.dp");
	TTreeReaderValue<Double_t> dums_phi(dumReader, "P.gtr.yptar");
	TTreeReaderValue<Double_t> dums_theta(dumReader, "P.gtr.xptar");
	TTreeReaderValue<Double_t> dumsytar(dumReader, "P.gtr.y");

	TTreeReaderValue<Double_t> dumh_px(dumReader, "H.gtr.px");
	TTreeReaderValue<Double_t> dumh_py(dumReader, "H.gtr.py");
	TTreeReaderValue<Double_t> dumh_pz(dumReader, "H.gtr.pz");
	TTreeReaderValue<Double_t> dums_px(dumReader, "P.gtr.px");
	TTreeReaderValue<Double_t> dums_py(dumReader, "P.gtr.py");
	TTreeReaderValue<Double_t> dums_pz(dumReader, "P.gtr.pz");

	TTreeReaderValue<Double_t> dumcointime(dumReader,"CTime.epCoinTime_ROC2");
	TTreeReaderValue<Double_t> dumt1r2(dumReader,"T.coin.pTRIG1_ROC2_tdcTimeRaw");
	TTreeReaderValue<Double_t> dumt4r2(dumReader,"T.coin.pTRIG4_ROC2_tdcTimeRaw");
	TTreeReaderValue<Double_t> dumtype(dumReader,"g.evtyp");
	TTreeReaderValue<Double_t> dumgevnum(dumReader,"g.evnum");
	TTreeReaderValue<Double_t> dum_tracknorm(dumReader, "P.cal.etottracknorm");
	TTreeReaderValue<Double_t> dum_beta(dumReader, "H.gtr.beta");

	TTreeReaderValue<Double_t> dumh1ytdc(dumReader, "H.hod.1y.NegTdcRefTime");
	TTreeReaderValue<Double_t> dumh1ytdcdiff(dumReader, "H.hod.1y.NegTdcRefDiffTime");
	TTreeReaderValue<Double_t> dumh1yadc(dumReader, "H.hod.1y.NegAdcRefTime");
	TTreeReaderValue<Double_t> dumh1yadcdiff(dumReader, "H.hod.1y.NegAdcRefDiffTime");
	TTreeReaderValue<Double_t> dump1ytdc(dumReader, "P.hod.1y.NegTdcRefTime");
	TTreeReaderValue<Double_t> dump1ytdcdiff(dumReader, "P.hod.1y.NegTdcRefDiffTime");
	TTreeReaderValue<Double_t> dump1yadc(dumReader, "P.hod.1y.NegAdcRefTime");
	TTreeReaderValue<Double_t> dump1yadcdiff(dumReader, "P.hod.1y.NegAdcRefDiffTime");
	TTreeReaderValue<Double_t> dumh1x1dcref(dumReader,"H.dc.1x1.RefTime");
	TTreeReaderValue<Double_t> dump1x1dcref(dumReader,"P.dc.1x1.RefTime");
	TTreeReaderValue<Double_t> dumhztar(dumReader, "H.react.z");
	TTreeReaderValue<Double_t> dumpztar(dumReader,"P.react.z");
	TTreeReaderValue<Double_t> dumhtime(dumReader, "H.hod.starttime");
	TTreeReaderValue<Double_t> dumptime(dumReader, "P.hod.starttime");

	// get simc values
	TTreeReaderValue<Float_t> shdp(simcReader, "hsdelta");
	TTreeReaderValue<Float_t> shphi(simcReader, "hsyptar");
	TTreeReaderValue<Float_t> shtheta(simcReader, "hsxptar");
	TTreeReaderValue<Float_t> ssdp(simcReader, "ssdelta");
	TTreeReaderValue<Float_t> ssphi(simcReader, "ssyptar");
	TTreeReaderValue<Float_t> sstheta(simcReader, "ssxptar");
	TTreeReaderValue<Float_t> w(simcReader, "Weight");
	TTreeReaderValue<Float_t> shxfp(simcReader, "hsxfp");
	TTreeReaderValue<Float_t> shxpfp(simcReader, "hsxpfp");
	TTreeReaderValue<Float_t> shyfp(simcReader, "hsyfp");
	TTreeReaderValue<Float_t> shypfp(simcReader, "hsypfp");
	TTreeReaderValue<Float_t> shytar(simcReader, "hsytar");
	TTreeReaderValue<Float_t> ssytar(simcReader, "ssytar");
	TTreeReaderValue<Float_t> ssztar(simcReader, "ssztar");
	TTreeReaderValue<Float_t> ssxfp(simcReader, "ssxfp");
	TTreeReaderValue<Float_t> ssxpfp(simcReader, "ssxpfp");
	TTreeReaderValue<Float_t> ssyfp(simcReader, "ssyfp");
	TTreeReaderValue<Float_t> ssypfp(simcReader, "ssypfp");
	TTreeReaderValue<Float_t> smiss_E(simcReader, "Em");
	TTreeReaderValue<Float_t> smissmass(simcReader, "Mm2");
	TTreeReaderValue<Float_t> simW(simcReader, "W");
	TTreeReaderValue<Float_t> simQ2(simcReader, "Q2");
	TTreeReaderValue<Float_t> siglab(simcReader, "siglab");
	TTreeReaderValue<Float_t> simhztar(simcReader, "hsztar");


	//	define array of histograms
	TH1F *h_sW[theggnum];
	TH1F *h_W[theggnum];
	TH1F *h_W_dum[theggnum];
	TH1F *h_W_acc[theggnum];
	TH1F *h_W_dumacc[theggnum];
	TH1F *h_sQ2[theggnum];
	TH1F *h_Q2[theggnum];
	TH1F *h_Q2_dum[theggnum];
	TH1F *h_Q2_acc[theggnum];
	TH1F *h_Q2_dumacc[theggnum];
	TH1F *h_smm2[theggnum];
	TH1F *h_mm2[theggnum];
	TH1F *h_mm2_dum[theggnum];
	TH1F *h_mm2_acc[theggnum];
	TH1F *h_mm2_dumacc[theggnum];
	TH1F* h_W_dum_new[theggnum];
	TH1F* h_Q2_dum_new[theggnum];
	TH1F* h_mm2_dum_new[theggnum];

	TString temp;

	temp=Form("%s ; #phi_{#gamma#gamma^{*}}; Counts", kinType.c_str());
	h_phigg = new TH1F("h_phigg", temp,nbins, 0., 360);
	h_phigg_acc = new TH1F("h_phigg_acc", temp,nbins, 0., 360);
	h_phigg_dum = new TH1F("h_phigg_dum", temp,nbins, 0., 360);
	h_phigg_dumacc = new TH1F("h_phigg_dumacc", temp,nbins, 0., 360);
	h_sphigg = new TH1F("h_sphigg", temp,nbins, 0., 360);

	temp=Form("%s ; #theta_{#gamma#gamma^{*}}; Counts", kinType.c_str());
	TH1F *h_sthegg = new TH1F("h_sthegg",temp,nbins, 80, 180);
	TH1F *h_thegg = new TH1F("h_thegg",temp,nbins, 80, 180);
	TH1F *h_thegg_dum = new TH1F("h_thegg_dum",temp,nbins, 80, 180);
	TH1F *h_thegg_acc = new TH1F("h_thegg_acc",temp,nbins, 80, 180);
	TH1F *h_thegg_dumacc = new TH1F("h_thegg_dumacc",temp,nbins, 80, 180);

	temp=Form("%s ; W(GeV); Counts", kinType.c_str());

	for(int i=0;i<theggnum;i++){
		h_sW[i] = new TH1F(Form("h_sW%d",i),temp,100,1.0,1.4);
		h_W[i] = new TH1F(Form("h_W%d",i),temp,100,1.0,1.4);
		h_W_dum[i] = new TH1F(Form("h_W_dum%d",i),temp,100,1.0,1.4);
		h_W_acc[i] = new TH1F(Form("h_W_acc%d",i),temp,100,1.0,1.4);
		h_W_dumacc[i] = new TH1F(Form("h_W_dumacc%d",i),temp,100,1.0,1.4);

	}
	temp=Form("%s ; Q2(GeV^{2}); Counts", kinType.c_str());

	for(int i=0;i<theggnum;i++){
		h_sQ2[i] = new TH1F(Form("h_sQ2%d",i),temp,100,0.15,0.55);
		h_Q2[i] = new TH1F(Form("h_Q2%d",i),temp,100,0.15,0.55);
		h_Q2_dum[i] = new TH1F(Form("h_Q2_dum%d",i),temp,100,0.15,0.55);
		h_Q2_acc[i] = new TH1F(Form("h_Q2_acc%d",i),temp,100,0.15,0.55);
		h_Q2_dumacc[i] = new TH1F(Form("h_Q2_dumacc%d",i),temp,100,0.15,0.55);

	}
	temp=Form("%s ; missing_mass^{2}; Counts", kinType.c_str());
	for(int i=0;i<theggnum;i++){

		h_smm2[i] = new TH1F(Form("h_smm2%d",i),temp,nbins,-0.02, 0.03);
		h_mm2[i] = new TH1F(Form("h_mm2%d",i),temp,nbins,-0.02, 0.03);
		h_mm2_dum[i] = new TH1F(Form("h_mm2_dum%d",i),temp,nbins,-0.02, 0.03);
		h_mm2_dum_new[i] = new TH1F(Form("h_mm2_dum_new%d",i),temp,nbins,-0.02, 0.03);
		h_mm2_acc[i] = new TH1F(Form("h_mm2_acc%d",i),temp,nbins,-0.02, 0.03);
		h_mm2_dumacc[i] = new TH1F(Form("h_mm2_dumacc%d",i),temp,nbins,-0.02, 0.03);

	};



	///////////////////////////////////////////////// calcuatlate charge rate///////////////////////////////////////////////////////


	Exp_charge = totcharge/1000.;
	Exp_dumcharge = dumtotcharge/1000.;
	cout << Exp_charge << "\n";
	cout << Exp_dumcharge << "\n";
	chargerate = Exp_charge/Exp_dumcharge;

	double timediff = maxtime - mintime;
	double timefac = 1./timediff*2.*cointime_cut;
	double simc_fac = normfac * Exp_charge * Data_eff /dataps6 / Ngen;
	double dum_fac = thickrate * chargerate * dumps6/dataps6 * Data_eff/Dum_eff;
	double dum_fac_nothk = chargerate * dumps6/dataps6 * Data_eff/Dum_eff;
	double datafac = factor;
	double dataaccfac = timefac*datafac;
	double dumfac, dumaccfac;
	double dumfac_nothk = dum_fac_nothk * dumfactor;
	double dumaccfac_nothk = timefac * dumfac_nothk;
	cout << "charge rate " << chargerate << endl;
	cout << "simc_fac " << simc_fac << endl;
	double* param = data.elosspara(kinType);
	double eeap0 = param[0];
	double eeap1 = param[1];
	double eebp0 = param[2];
	double eebp1 = param[3];
	double eecp0 = param[4];
	double eecp1 = param[5];
	double ebap0 = param[6];
	double ebap1 = param[7];
	double ebbp0 = param[8];
	double ebbp1 = param[9];
	double ebcp0 = param[10];
	double ebcp1 = param[11];
	double maxpX = param[12];
	double beloss,eeloss,dumbeloss,dumeeloss,simbeloss,simeeloss;
	double stheta_off = 0.0005;
	double htheta_off = 0.0;

	//////////////////////////////////////////////////////////loop over data rootfile/////////////////////////////////////////////

	if (beam_E < 0) {
		beam_E = -beam_E;
	}

	while(myReader.Next()){

		if (t % 100000000 == 0) cout << "experiment  " << t << endl;
		//	if (t == 1000) return;
		double cotime = *cointime - timeshift;
		double hadcns= *h1yadc*0.0625;
		double padcns= *p1yadc*0.0625;
		double hadcdiffns= *h1yadcdiff*0.0625;
		double padcdiffns= *p1yadcdiff*0.0625;
		double htdcns= *h1ytdc*0.09776;
		double ptdcns= *p1ytdc*0.09776;
		double htdcdiffns= *h1ytdcdiff*0.09776;
		double ptdcdiffns= *p1ytdcdiff*0.09776;
		bool diffcut = (hadcdiffns<130. || hadcdiffns>170.) && (htdcdiffns<130. || htdcdiffns>170.) && (padcdiffns<160. || padcdiffns>200.) && (ptdcdiffns<160. || ptdcdiffns>200.);
		bool newvarcut = (*h1ytdc>4200. && *h1ytdc<4800.) && (*h1yadc>6630. && *h1yadc<6950.) && (*h1x1dcref>17500. && *h1x1dcref<17800.) && (*p1ytdc>4750. && *p1ytdc<5500.) && (*p1yadc>6040. && *p1yadc<7100.) && (*p1x1dcref>15680. && *p1x1dcref<16400.);
		bool data_cut =  *type == type_cut && diffcut && newvarcut && abs(*hms_dp)<8. && *shms_dp>-10. && *shms_dp < 22.;
		if (data_cut) {

			if(*hztar<= -4.5){
				beloss = ebap0 + *hztar * ebap1;
				eeloss = eeap0 + *hztar * eeap1;
			} else if (*hztar > -4.5 && *hztar < 4.5) {
				beloss = ebbp0 + *hztar * ebbp1;
				eeloss = eebp0 + *hztar * eebp1;

			} else {
				beloss = ebcp0 + *hztar * ebcp1;
				eeloss = eecp0 + *hztar * eecp1;

			}

			// energy loss
			double Eescat1 = sqrt(shms_p*shms_p + m_e*m_e) + eeloss/1000.;  ////////////
			double Pescat1 = sqrt(Eescat1*Eescat1 - m_e*m_e);                 // 0.997  //
			double Eprecoil1 = sqrt(hms_p*hms_p + m_p*m_p) + maxpX/1000.;   // 0.997  //
			double Pprecoil1 = sqrt(Eprecoil1*Eprecoil1 - m_p*m_p);           ////////////

			double Eescat2 = sqrt(shms_pp*shms_pp + m_e*m_e) + eeloss/1000.;////////////
			double Pescat2 = sqrt(Eescat2*Eescat2 - m_e*m_e);                 // 0.9967 //
			double Eprecoil2 = sqrt(hms_pp*hms_pp + m_p*m_p) + maxpX/1000.; // 0.9987 //
			double Pprecoil2 = sqrt(Eprecoil2*Eprecoil2 - m_p*m_p);           ////////////
			// calculate 4-vectors
			double Ebeam = beam_E - beloss/1000.;
			k_in.SetXYZM(0.,0.,Ebeam, m_e);  //k_in with energy loss correction
			//	k_in.SetXYZM(0.,0.,beam_E, m_e); //k_in without energy loss correction
			p_in.SetXYZM(0.,0.,0.,m_p);

			p_out = data.converter(Pprecoil2, m_p, -hms_deg, *hms_dp, *hms_theta+htheta_off, *hms_phi); // with energy loss correction
			k_out = data.converter(Pescat2, m_e, shms_deg, *shms_dp, *shms_theta +8.681269905E-4+stheta_off, *shms_phi); // with energy loss correction 
			q_out = k_in + p_in - k_out - p_out;
			q_in = k_in - k_out;
			delta = q_in + p_in;

			// calculate w q2 mm2
			double W = delta.Mag();
			double Q2 = -q_in.M2();
			missmass2 = q_out.M2();
			double mm1 = sqrt(abs(q_out.M2()) * 1e6);
			double emiss = q_out.E();
			TRotation rot_to_q;
			rot_to_q.SetZAxis(q_in.Vect(),k_out.Vect()).Invert();
			TVector3 bq = q_out.Vect();
			bq *= rot_to_q;
			TVector3 p_miss = -bq;
			double pmiss = p_miss.Mag();

			// calculate thetagg			
			double pin = m_p * q_in.P()/W;
			double t = - (p_in - p_out).M2();
			double egamma = (W*W - m_p*m_p)/(2.*W);
			double pgamma = sqrt(egamma*egamma);
			double costhgamma = (-t + Q2 + 2.* egamma * sqrt(pin*pin-Q2))/(2.*pin*pgamma);
			double thegg_rad, thegg;
			thegg_rad = acos(costhgamma);
			thegg = thegg_rad * TMath::RadToDeg();
			//	calculate phigg
			TVector3 normal_mu, normal_out;
			double phigg_rad, phigg;
			normal_mu = (k_in.Vect().Cross(k_out.Vect())).Unit();
			normal_out = (q_in.Vect().Cross(p_out.Vect())).Unit();
			phigg_rad = TMath::ACos(normal_mu*normal_out);
			phigg = 360. - phigg_rad * TMath::RadToDeg();
			if((normal_mu * q_out.Vect()) < 0.) {
				phigg = 360. - phigg;
			}	
			bool MM2_cut =  missmass2>mm_min && missmass2<mm_max && W>wmin && W<wmax;
			if(MM2_cut && abs(Q2-q2_cent)<q2cut && (phigg < phicut || phigg > (360.- phicut))) {

				if (abs(cotime) < cointime_cut) {

					if(abs(W-w_cent)<wcut){
						h_thegg->Fill(thegg, datafac);
						h_phigg->Fill(phigg, datafac);
						for(int ct =0;ct<theggnum;ct++){
							if(abs(thegg-theggcm[ct])<=theggcut) {

								h_mm2[ct]->Fill(missmass2, datafac);
								h_W[ct]->Fill(W, datafac);
								h_Q2[ct]->Fill(Q2, datafac);
							}
						}
					}
				}

				if (cotime> mintime && cotime < maxtime) {
					if(abs(W-w_cent)<wcut){	
						h_thegg_acc->Fill(thegg, dataaccfac);
						h_phigg_acc->Fill(phigg, dataaccfac);
						for(int ct =0;ct<theggnum;ct++){
							if(abs(thegg-theggcm[ct])<=theggcut) {

								h_mm2_acc[ct]->Fill(missmass2, dataaccfac);
								h_W_acc[ct]->Fill(W, dataaccfac);
								h_Q2_acc[ct]->Fill(Q2, dataaccfac);
							}
						}					
					}
				}

			}
		}

		t++;
	}	

	int k = 0;
	cout << " start to loop dummy " << endl;
	//////////////////////////////////////////////////////////loop over dummy rootfile/////////////////////////////////////////////
	while(dumReader.Next()){

		if (k % 10000000 == 0) cout << "dummy  " << k << endl;

		double dumcotime = *dumcointime - timeshift;
		double dumhadcns= *dumh1yadc*0.0625;
		double dumpadcns= *dump1yadc*0.0625;
		double dumhadcdiffns= *dumh1yadcdiff*0.0625;
		double dumpadcdiffns= *dump1yadcdiff*0.0625;
		double dumhtdcns= *dumh1ytdc*0.09776;
		double dumptdcns= *dump1ytdc*0.09776;
		double dumhtdcdiffns= *dumh1ytdcdiff*0.09776;
		double dumptdcdiffns= *dump1ytdcdiff*0.09776;
		bool dumdiffcut = (dumhadcdiffns<130. || dumhadcdiffns>170.) && (dumhtdcdiffns<130. || dumhtdcdiffns>170.) && (dumpadcdiffns<160. || dumpadcdiffns>200.) && (dumptdcdiffns<160. || dumptdcdiffns>200.);
		bool dumnewvarcut = (*dumh1ytdc>4200. && *dumh1ytdc<4800.) && (*dumh1yadc>6630. && *dumh1yadc<6950.) && (*dumh1x1dcref>17500. && *dumh1x1dcref<17800.) && (*dump1ytdc>4750. && *dump1ytdc<5500.) && (*dump1yadc>6040. && *dump1yadc<7100.) && (*dump1x1dcref>15680. && *dump1x1dcref<16400.);
		bool dummy_cut = *dumtype == type_cut && dumdiffcut && dumnewvarcut && abs(*dumh_dp)<8. && *dums_dp>-10. && *dums_dp < 22. ;
		if (dummy_cut) {

			if(*dumhztar<= -4.5){
				dumbeloss = ebap0 + *dumhztar * ebap1;
				dumeeloss = eeap0 + *dumhztar * eeap1;
			} else if (*dumhztar > -4.5 && *dumhztar < 4.5) {
				dumbeloss = ebbp0 + *dumhztar * ebbp1;
				dumeeloss = eebp0 + *dumhztar * eebp1;

			} else {
				dumbeloss = ebcp0 + *dumhztar * ebcp1;
				dumeeloss = eecp0 + *dumhztar * eecp1;

			}

			double dumEescat1 = sqrt(shms_p*shms_p + m_e*m_e) + dumeeloss/1000.;		////////////
			double dumPescat1 = sqrt(dumEescat1*dumEescat1 - m_e*m_e);                   // 0.997  //
			double dumEprecoil1 = sqrt(hms_p*hms_p + m_p*m_p) + maxpX/1000.;         // 0.997  //
			double dumPprecoil1 = sqrt(dumEprecoil1*dumEprecoil1 - m_p*m_p);             ////////////

			double dumEescat2 = sqrt(shms_pp*shms_pp + m_e*m_e) + dumeeloss/1000.;   ////////////
			double dumPescat2 = sqrt(dumEescat2*dumEescat2 - m_e*m_e);                   // 0.9967 //
			double dumEprecoil2 = sqrt(hms_pp*hms_pp + m_p*m_p) + maxpX/1000.;       // 0.9987 //
			double dumPprecoil2 = sqrt(dumEprecoil2*dumEprecoil2 - m_p*m_p);             ////////////
			// calculate 4-vectors
			double dumEbeam = beam_E - dumbeloss/1000.;
			dumk_in.SetXYZM(0.,0.,dumEbeam, m_e); // with energy loss correction
			//	dumk_in.SetXYZM(0.,0.,beam_E, m_e); // w/o energy loss correction
			dump_in.SetXYZM(0.,0.,0.,m_p);

			dump_out = dummy.converter(dumPprecoil2, m_p, -hms_deg, *dumh_dp, *dumh_theta+htheta_off, *dumh_phi); // with energy loss correction
			dumk_out = dummy.converter(dumPescat2, m_e, shms_deg, *dums_dp, *dums_theta +8.681269905E-4+stheta_off, *dums_phi); // with energy loss correction
			//	dump_out = dummy.converter(hms_p, m_p, -hms_deg, *dumh_dp, *dumh_theta, *dumh_phi); // w/o energy loss
			//	dumk_out = dummy.converter(shms_p, m_e, shms_deg, *dums_dp, *dums_theta +8.681269905E-4, *dums_phi); // w/o energy loss
			dumq_out = dumk_in + dump_in - dumk_out - dump_out;
			dumq_in = dumk_in - dumk_out;
			dumdelta = dumq_in + dump_in;
			// calculate w q2 mm2
			double dumW = dumdelta.Mag();
			double dumQ2 = -dumq_in.M2();
			dummiss2 = dumq_out.M2();
			double dumm1 = sqrt(abs(dumq_out.M2()) * 1e6);
			double dumemiss = dumq_out.E();
			TRotation dumrot_to_q;
			dumrot_to_q.SetZAxis(dumq_in.Vect(),dumk_out.Vect()).Invert();
			TVector3 dumbq = dumq_out.Vect();
			dumbq *= dumrot_to_q;
			TVector3 dump_miss = - dumbq;
			double dumpmiss = dump_miss.Mag();

			//	calculate thetagg
			double dumpin = m_p * dumq_in.P()/dumW;
			double dumt = - (dump_in - dump_out).M2();
			double dumegamma = (dumW*dumW - m_p*m_p)/(2.*dumW);
			double dumpgamma = sqrt(dumegamma*dumegamma);
			double dumcosthgamma = (-dumt + dumQ2 + 2.* dumpgamma * sqrt(dumpin*dumpin-dumQ2))/(2.*dumpin*dumpgamma);
			double dumthegg_rad, dumthegg, dumphigg_rad, dumphigg;
			dumthegg_rad = acos(dumcosthgamma);
			dumthegg = dumthegg_rad * TMath::RadToDeg();
			//	calculate phigg
			TVector3 dumnormal_mu, dumnormal_out;
			dumnormal_mu = (dumk_in.Vect().Cross(dumk_out.Vect())).Unit();
			dumnormal_out = (dumq_in.Vect().Cross(dump_out.Vect())).Unit();
			dumphigg_rad = TMath::ACos(dumnormal_mu*dumnormal_out);
			dumphigg = 360. - dumphigg_rad * TMath::RadToDeg();
			if((dumnormal_mu * dumq_out.Vect()) < 0.) {
				dumphigg = 360. - dumphigg;
			}

			bool dumMM2_cut = dummiss2>mm_min && dummiss2<mm_max && dumW>wmin && dumW<wmax;
			if (dumMM2_cut && abs(dumQ2-q2_cent)<q2cut && (dumphigg < phicut || dumphigg > (360.- phicut))) {
				if (*dumpztar >= 0) {
					dumfac = dumfac_nothk * exitth;
					dumaccfac = dumaccfac_nothk * exitth;
					if (abs(dumcotime) < cointime_cut) {

						if(abs(dumW-w_cent)<wcut){
							h_thegg_dum->Fill(dumthegg, dumfac);
							h_phigg_dum->Fill(dumphigg, dumfac);
							for(int ct =0;ct<theggnum;ct++){
								if(abs(dumthegg-theggcm[ct])<=theggcut) {

									h_mm2_dum[ct]->Fill(dummiss2, dumfac);
									h_W_dum[ct]->Fill(dumW, dumfac);
									h_Q2_dum[ct]->Fill(dumQ2, dumfac);
								}
							}					
						}
					}
					if (dumcotime > mintime && dumcotime < maxtime) {

						if(abs(dumW-w_cent)<wcut){
							h_thegg_dumacc->Fill(dumthegg, dumaccfac);
							h_phigg_dumacc->Fill(dumphigg, dumaccfac);
							for(int ct =0;ct<theggnum;ct++){
								if(abs(dumthegg-theggcm[ct])<=theggcut) {

									h_mm2_dumacc[ct]->Fill(dummiss2, dumaccfac);
									h_W_dumacc[ct]->Fill(dumW, dumaccfac);
									h_Q2_dumacc[ct]->Fill(dumQ2, dumaccfac);
								}
							}					

						}

					}
				}	
				if (*dumpztar < 0) {
					dumfac = dumfac_nothk * entrth;
					dumaccfac = dumaccfac_nothk * entrth;

					if (abs(dumcotime) < cointime_cut) {

						if(abs(dumW-w_cent)<wcut){
							h_thegg_dum->Fill(dumthegg, dumfac);
							h_phigg_dum->Fill(dumphigg, dumfac);
							for(int ct =0;ct<theggnum;ct++){
								if(abs(dumthegg-theggcm[ct])<=theggcut) {

									h_mm2_dum[ct]->Fill(dummiss2, dumfac);
									h_W_dum[ct]->Fill(dumW, dumfac);
									h_Q2_dum[ct]->Fill(dumQ2, dumfac);
								}
							}					
						}
					}
					if (dumcotime > mintime && dumcotime < maxtime) {

						if(abs(dumW-w_cent)<wcut){
							h_thegg_dumacc->Fill(dumthegg, dumaccfac);
							h_phigg_dumacc->Fill(dumphigg, dumaccfac);
							for(int ct =0;ct<theggnum;ct++){
								if(abs(dumthegg-theggcm[ct])<=theggcut) {

									h_mm2_dumacc[ct]->Fill(dummiss2, dumaccfac);
									h_W_dumacc[ct]->Fill(dumW, dumaccfac);
									h_Q2_dumacc[ct]->Fill(dumQ2, dumaccfac);
								}
							}					

						}

					}
				}		

			}
		}
		k++;
	}
	h_phigg_dum->Add(h_phigg_dumacc,-1);
	h_thegg_dum->Add(h_thegg_dumacc,-1);
	h_phigg->Add(h_phigg_acc,-1);
	h_thegg->Add(h_thegg_acc,-1);
	h_phigg->Add(h_phigg_dum,-1);
	h_thegg->Add(h_thegg_dum, -1);

	Functions data_new[theggnum];
	Functions data_new_low[theggnum];

	cout << "end definition" << endl;

	for(int t =0;t<theggnum;t++){

		h_mm2_dum[t]->Add(h_mm2_dumacc[t],-1);

		double wbin_min = h_mm2_dum[t]->GetXaxis()->FindBin(-0.01);
		double wbin_max = h_mm2_dum[t]->GetXaxis()->FindBin(0.01);
		int wbinnum = wbin_max - wbin_min;
		double sumw = 0.;
		for(int i = 0; i<wbinnum; i++) {

			sumw += h_mm2_dum[t]->GetBinContent(wbin_min+i);
		}
		for(int j=1; j<h_mm2_dum[t]->GetNbinsX(); j++) {
			if(j>=wbin_min && j<=wbin_max) {
				h_mm2_dum_new[t]->SetBinContent(j,sumw/wbinnum);
			}else{
				h_mm2_dum_new[t]->SetBinContent(j,0);
			}
		}


		h_mm2[t]->Add(h_mm2_acc[t],-1);
		h_mm2[t]->Add(h_mm2_dum_new[t], -1);

		h_W_dum[t]->Add(h_W_dumacc[t],-1);
		h_W[t]->Add(h_W_acc[t], -1);
		h_W[t]->Add(h_W_dum[t], -1);

		h_Q2_dum[t]->Add(h_Q2_dumacc[t],-1);
		h_Q2[t]->Add(h_Q2_acc[t], -1);
		h_Q2[t]->Add(h_Q2_dum[t], -1);
	}



	HList.Add(h_phigg);
	HList.Add(h_thegg);
	for(int t =0;t<theggnum;t++){

		HList.Add(h_mm2[t]);
		HList.Add(h_W[t]);
		HList.Add(h_Q2[t]);
	}

	//////////////////////////////////////////////////////////loop over simulation rootfile/////////////////////////////////////////////////
	//resolution test
	TRandom3 *rand3 = new TRandom3();
	while(simcReader.Next()) {
		if (f % 100000 == 0) cout << "simc  " << f << endl;
		//	bool simc_cut = abs(*shdp-0) < hms_dp_cut && abs(*ssdp-0) < shms_dp_cut && abs(*shtheta - 0) < hms_theta_cut && abs(*shphi - 0) < hms_phi_cut && abs(*sstheta - 0) < shms_theta_cut && abs(*ssphi - 0) < shms_phi_cut;
		bool simc_cut = abs(*shdp-0) < 8. && *ssdp >-10.;
		if (simc_cut) {
			if(*simhztar<= -4.5){
				simbeloss = ebap0 + *simhztar * ebap1;
				simeeloss = eeap0 + *simhztar * eeap1;
			} else if (*simhztar > -4.5 && *simhztar < 4.5) {
				simbeloss = ebbp0 + *simhztar * ebbp1;
				simeeloss = eebp0 + *simhztar * eebp1;

			} else {
				simbeloss = ebcp0 + *simhztar * ebcp1;
				simeeloss = eecp0 + *simhztar * eecp1;

			}
			double simEescat1 = sqrt(shms_p*shms_p + m_e*m_e) + simeeloss/1000.;			////////////
			double simPescat1 = sqrt(simEescat1*simEescat1 - m_e*m_e);                          // 0.997  //
			double simEprecoil1 = sqrt(hms_p*hms_p + m_p*m_p) + maxpX/1000.;            // 0.997  //
			double simPprecoil1 = sqrt(simEprecoil1*simEprecoil1 - m_p*m_p);                    ////////////

			double simEescat2 = sqrt(shms_pp*shms_pp + m_e*m_e) + simeeloss/1000.;        ////////////
			double simPescat2 = sqrt(simEescat2*simEescat2 - m_e*m_e);                          // 0.9967 //
			double simEprecoil2 = sqrt(hms_pp*hms_pp + m_p*m_p) + maxpX/1000.;	         // 0.9987 //
			double simPprecoil2 = sqrt(simEprecoil2*simEprecoil2 - m_p*m_p);                    ////////////
			// calculate 4-vectors
			double simEbeam = beam_E - simbeloss/1000.;
			simk_in.SetXYZM(0.,0.,simEbeam, m_e); // with energy loss correction
			//	simk_in.SetXYZM(0.,0.,beam_E, m_e); // w/o energy loss correction
			simp_in.SetXYZM(0.,0.,0.,m_p);

			double shms_xptar = *sstheta + rand3->Gaus(0.0, 0.001);
			double shms_yptar = *ssphi + rand3->Gaus(0.0, 0.001);
			double hms_xptar = *shtheta + rand3->Gaus(0.0, 0.004);
			double hms_yptar = *shphi + rand3->Gaus(0.0, 0.003);
			//	simk_out = simc.converter(simPescat2,m_e,shms_deg,*ssdp,shms_xptar,shms_yptar); // with energy loss correction
			//	simp_out = simc.converter(simPprecoil2,m_p,-hms_deg,*shdp,hms_xptar,hms_yptar); // with energy loss correction
			simk_out = simc.converter(simPescat2,m_e,shms_deg,*ssdp,*sstheta,*ssphi); // w/o energy loss
			simp_out = simc.converter(simPprecoil2,m_p,-hms_deg,*shdp,*shtheta,*shphi); //w/o energy loss
			simq_out = simk_in+simp_in-simk_out-simp_out;
			simq_in = simk_in - simk_out;
			simdelta = simq_in + simp_in;
			// calculate w q2 mm2
			double simm2 = simq_out.M2();
			double simm1 = sqrt(abs(simq_out.M2())*1e6);
			double simQ2, simW;
			simQ2 = - simq_in.Mag2();
			simW = simdelta.Mag();
			double semiss = simq_out.E();
			TRotation simrot_to_q;
			simrot_to_q.SetZAxis(simq_in.Vect(),simk_out.Vect()).Invert();
			TVector3 simbq = simq_out.Vect();
			simbq *= simrot_to_q;
			TVector3 simp_miss = -simbq;
			double spmiss = simp_miss.Mag();

			//	calculate thetagg			
			double simpin = m_p * simq_in.P()/ simW;
			double simt = - (simp_in - simp_out).M2();
			double simegamma = (pow(simW,2.) - m_p*m_p)/(simW * 2.);
			double simpgamma = sqrt(simegamma*simegamma);
			double simcosthgamma = (-simt + simQ2 + 2.* simegamma * sqrt(simpin*simpin- simQ2))/(2.*simpin*simpgamma);
			double simthegg_rad, simthegg, simphigg_rad, simphigg;
			simthegg_rad = acos(simcosthgamma);
			simthegg = simthegg_rad * TMath::RadToDeg();
			//	calculate phigg
			TVector3 simnormal_mu, simnormal_out;
			simnormal_mu = (simk_in.Vect().Cross(simk_out.Vect())).Unit();
			simnormal_out = (simq_in.Vect().Cross(simp_out.Vect())).Unit();
			simphigg_rad = TMath::ACos(simnormal_mu*simnormal_out);
			simphigg = 360. - simphigg_rad * TMath::RadToDeg();
			if((simnormal_mu * simq_out.Vect()) < 0.) {
				simphigg = 360. - simphigg;
			}

			bool simMM2_cut = simm2>mm_min && simm2<mm_max && simW>wmin && simW<wmax;
			if(simMM2_cut && abs(simQ2-q2_cent)<q2cut && (simphigg < phicut || simphigg > (360.- phicut))) {
				if(abs(simW-w_cent)<wcut){
					h_sthegg->Fill(simthegg,*w * simc_fac);
					h_sphigg->Fill(simphigg, *w * simc_fac);
					for(int ct =0;ct<theggnum;ct++){
						if(abs(simthegg-theggcm[ct])<=theggcut) {

							h_smm2[ct]->Fill(simm2,*w * simc_fac);
							h_sW[ct]->Fill(simW, *w * simc_fac );
							h_sQ2[ct]->Fill(simQ2, *w * simc_fac);
						}
					}					
				}
			}
		}
		f++;
	}

	HList.Add(h_sphigg);
	HList.Add(h_sthegg);
	for(int t =0;t<theggnum;t++){

		HList.Add(h_smm2[t]);
		HList.Add(h_sW[t]);
		HList.Add(h_sQ2[t]);
	}

	//////////////////////////////////////////////find integrals and uncertainties//////////////////////////////////////////////////////////
	//	double* datastat = data.find_uncertain(h_W, h_sW, 1.0, 1.3, 1);
	//	double ratio = datastat[0];
	//	double uncertain = datastat[1];
	//	double numperc= datastat[2]/Data_eff/Exp_charge; 
	//	double uncer = sqrt(datastat[2]) * 1/Data_eff/Exp_charge;
	//	cout << "counts per charge is: " << numperc << " +/- " << uncer << endl;
	TFile hhist(outFile,"recreate");
	HList.Write();
	/////////////////////////////////////////////////////////draw plots////////////////////////////////////////////////////////////////////
	//

	TLine* line[theggnum];
	TPaveText* text[theggnum];
	TCanvas *c1 = new TCanvas("c1","c1");
	c1->Divide(4,2);
	for(int i = 0; i < 8; i++) {

		c1->cd(i+1);	
		h_mm2[i]->SetLineColor(2); h_mm2[i]->Draw("ep");
		h_smm2[i]->SetLineColor(4);h_smm2[i]->Draw("HIST,same");
		h_mm2_dum_new[i]->Draw("epsame");h_mm2_dum_new[i]->SetLineColor(8);
		h_mm2_acc[i]->Draw("epsame");h_mm2_acc[i]->SetLineColor(1);

		text[i]= new TPaveText(0.7, 0.5, 0.9,0.7,"b1NDC");
		text[i]->SetBorderSize(1);
		text[i]->SetFillStyle(0);
		text[i]->SetLineColor(0);
		text[i]->AddText(Form("#frac{dummy}{data}=%.3f",h_mm2_dum_new[i]->Integral()/h_mm2[i]->Integral()));
		text[i]->AddText(Form("#frac{data}{SIMC}=%.3f",h_mm2[i]->Integral()/h_smm2[i]->Integral()));
		text[i]->Draw();

		line[i] = new TLine(0.0,0.0,0.0,h_smm2[i]->GetMaximum());
		line[i]->SetLineStyle(10);
		line[i]->Draw();
	}
	c1->SaveAs(Form("vcs_xs_%s_multibins.pdf",kinType.c_str()));


}

