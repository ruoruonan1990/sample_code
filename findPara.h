#include <iostream>
#include <string>
#include <sstream>
#include "TProfile.h"

using namespace std;


class Functions{

	public:

		//------------------------------------------------------------------find ps6--------------------------------------------------------------//
		double findps6(int nRun){

			ifstream infile1(Form("/work/hallc/alphaE/ruonanli/analyzer/hallc_replay_vcs/REPORT_OUTPUT/COIN/PRODUCTION/replay_coin_production_%d_-1.report",nRun));
			stringstream ss1;
			string line1;
			char const* key1 = "Ps6_factor";
			string junk1;

			if(infile1.is_open()) {

				while(getline(infile1,line1)) {
					ss1.clear();
					ss1.str(line1);

					if(line1.find(key1, 0) != string::npos) {
						ss1 >> junk1 >> junk1 >> ps6;
					}
				}
			}else{
				cout << "report File doesn't exist! " << "\n";
				return 0;
			}
			return ps6;

		}
	//-----------------------------------------------------------find efficiencies-------------------------------------------------------------//

		double* findeff(int Number) {


			ifstream infile3("/work/hallc/alphaE/ruonanli/vcs_ruonan/scripts/effciencies/efficiency.txt");
			stringstream ss3;
			string line3;
			string junk3;
			int key3;
			double CLT, ELT, PTRK, PTRI, HTRK, HTRI, PROTON,EDTM,CUR,CTIME,NORM,FAC;

			if(infile3.is_open()) {

				while(getline(infile3,line3)) {
				key3 = 0;
				ss3.clear();
				ss3.str(line3);
				ss3 >> key3;

					if(key3 == Number) {
						ss3 >> CLT >> ELT >> PTRK >> PTRI >> HTRK >> HTRI>> PROTON >> EDTM >> CUR >> CTIME >> NORM >> FAC;
						eff[0] = CLT*ELT*PTRK*PTRI*HTRK*HTRI*PROTON;
						eff[1] = CTIME;
						eff[2] = NORM;
						eff[3] = FAC;
						eff[4] = CLT*ELT*PTRK*HTRK;
					}
				}
			}
			else{
				cout << "EFF File doesn't exist! " << "\n";
				return 0;
			}
			cout << CLT <<"\t" << ELT <<"\t" << PTRK <<"\t" << PTRI <<"\t" << HTRK <<"\t" << HTRI <<"\t" << PROTON << endl;
			return eff;
		}
		//---------------------------------------------------------find kinematics-------------------------------------------------------------//
		double* readDbase(int Number) {

			ifstream file("/work/hallc/alphaE/ruonanli/analyzer/hallc_replay_vcs/DBASE/COIN/standard.kinematics");

			stringstream ss;
			string line;

			int record_values = 0;  //this is switch, it starts in the off position.
			int key;
			string junk;  //junk don't want to save

			if(file.is_open()){

				while(getline(file,line)){

					key = 0;
					ss.clear();
					ss.str(line);

					if(record_values == 6){
						ss >> junk >> junk >> arr[5];  //SHMS p
						record_values = 0;
					}
					if(record_values == 5){
						ss >> junk >> junk >> arr[4]; // HMS p
						record_values++;
					}
					if(record_values == 4){
						ss >> junk >> junk >> arr[3]; //SHMS deg 
						record_values++;
					}
					if(record_values == 3){
						ss >> junk >> junk >> arr[2];  //HMS deg
						record_values++;
					}
					if(record_values == 2){
						ss >> junk >> junk >> arr[1]; //target mass
						record_values++;
					}

					if(record_values == 1){
						ss >> junk >> junk >> arr[0];  //beam energy
						record_values++;
					}
					if(!record_values) {

						ss >> key;

						if(key == Number) 	record_values = 1;  //next line iteration, will record something
					}
				}


			}
			else{
				cout << "dbase File doesn't exist" << "\n";
				return 0;
			}

			return arr;
		}

	//--------------------------------------------------------------converter-from spectrometer frame to C.M. frame---------------------------------//
		TLorentzVector converter(double central_p, double mass, double angle, double dp, double theta, double phi) {

			double total_p = central_p * (1.+ dp/100.);
			double pz = sqrt(total_p * total_p / (theta * theta + phi * phi + 1.));
			double px = pz * theta;
			double py = pz * phi;

			ret.SetXYZM(px,py,pz,mass);

			ret.RotateZ(TMath::Pi()/2.);
			ret.RotateY(-angle * TMath::DegToRad());

			return ret;
		}

		void newconverter(double central_p, double mass, double angle, double dp, double theta, double phi,double eloss,TLorentzVector &ret) {

			double total_p = central_p * (1.+ dp/100.);
			double energy = sqrt(total_p*total_p + mass*mass) + eloss/1000;
			double new_total_p = sqrt(energy*energy - mass*mass);
			double pz = sqrt(new_total_p * new_total_p / (theta * theta + phi * phi + 1.));
			double px = pz * theta;
			double py = pz * phi;

			ret.SetXYZM(px,py,pz,mass);

			ret.RotateZ(TMath::Pi()/2.);
			ret.RotateY(-angle * TMath::DegToRad());

		}

		//-------------------------------------------------------find averaged dummy backgroud-----------------------------------------------------//
		TH1F* ave_dum(TH1F* h_dum) {
		
			double wbin_min = h_dum->GetXaxis()->FindBin(-0.01);
			double wbin_max = h_dum->GetXaxis()->FindBin(0.01);
			int wbinnum = wbin_max - wbin_min;
			double sumw = 0.;
			for(int i = 0; i<wbinnum; i++) {

				sumw += h_dum->GetBinContent(wbin_min+i);
			}
			
			for(int j=1; j<h_dum->GetNbinsX(); j++) {
				if(j>=wbin_min && j<=wbin_max) {
					h_ave->SetBinContent(j,sumw/wbinnum);
				}else{
					h_ave->SetBinContent(j,0);
				}
			}
			return h_ave;
		}
		
		//-----------------------------------------------------find uncertainties-----------------------------------------------------------------//
		double* find_uncertain(TH1F* h_data, TH1F* h_simc, double min, double max) {

			double sum_err2= 0.;
			double datamin = h_data->FindBin(min);
			double datamax = h_data->FindBin(max);
			double simcmin = h_simc->FindBin(min);
			double simcmax = h_simc->FindBin(max);
			
			// err calculated with IntegralAndError
			double dataerr, simcerr;
			double integral_data, integral_simc;

			integral_data = h_data->IntegralAndError(datamin, datamax, dataerr,"");
			integral_simc = h_simc->IntegralAndError(simcmin, simcmax, simcerr,"");
			stat[0] = integral_data/integral_simc; 
			stat[2] = integral_data;
			stat[3] = dataerr;
			//uncertainty for A/B = A/B * sqrt(pow(sigmaA/A,2)+pow(sigmaB/B,2))
			stat[1] = stat[0] * sqrt(pow(dataerr/integral_data,2) + pow(simcerr/integral_simc,2));
			cout << "Data IntegralAndError: " << integral_data << " +/- "<< dataerr << endl;
			cout << "Simc IntegralAndError: " << integral_simc << " +/- " << simcerr << endl;
			cout << "IntegralAndError:data/simc = " << stat[0]  << " +/- " << stat[1] << endl;
		
			// err calculated with error propogration from bin errors(simc) - cross check
	/*		double uncertainty_data, uncertainty_simc;

			for (int d = simcmin; d<= simcmax; d++) {
				sum_err2 += pow(h_simc->GetBinError(d), 2);
			}
			uncertainty_simc = sqrt(sum_err2);
			uncertainty_data = dataerr;
			double uncertainty_ratio = stat[0] * sqrt(pow(uncertainty_data/integral_data,2) + pow(uncertainty_simc/integral_simc,2));*/
		//	cout << "Propogation of bin errors:simc error " << uncertainty_simc << endl;
		//	cout << "Propogation of bin errors:data/simc = " << stat[0]  << "	+/-		" << uncertainty_ratio << endl;

			return stat;
		}

		//-----------------------------------------------------Gaussian fitting-----------------------------------------------------------------//
		double* findpeak(TH1F* h1, double width, TF1* fh,int linecol) {

			double a = h1->GetBinCenter(h1->GetMaximumBin());
			fh = new TF1("fh","gaus",a-width,a+width);
			h1->Fit("fh","QR");
			peak[0] = fh->GetParameter(1);
			peak[1] = fh->GetParError(1);
		//	cout << "peak is " << peak[0] << "  +/-  " << peak[1] << endl;
			fh->SetLineColor(linecol);
			fh->Draw("same");
			return peak;
		}
		//----------------------------------------------------draw TProfile----------------------------------------------------------------------//

		TProfile* drawprofile(TH2F* h1, double xmin, double xmax, double ymin, double ymax) {

			TProfile* h2;
			h2 = h1->ProfileX("h2");
			h2->GetXaxis()->SetRange(h1->GetXaxis()->FindBin(xmin), h1->GetXaxis()->FindBin(xmax));
			h2->SetMinimum(ymin);
			h2->SetMaximum(ymax);

			return h2;
		}


		//----------------------------------------------------read energy loss correction vals-----------------------------------------------------//
		double* elosspara(string kinType) {

			if(kinType=="kin1a") {
			//	param[0] = 3.21855;  // newtarget_muscext
			//	param[1] = 0.0504515;
			//	param[2] = 2.07455;
			//	param[3] = -0.231156;
			//	param[4] = 0.959132;
			//	param[5] = 0.0392848;
			//	param[6] = 0.417479;
			//	param[7] = -0.0025103;
			//	param[8] = 1.60672;
			//	param[9] = 0.277392;
			//	param[10] = 2.64215; 
			//	param[11] = 0.0361774; 
			//	param[12] = 3.09;  
				param[0] = 3.1769; // new vovilt 3 runs
				param[1] = 0.0427051;
				param[2] = 2.07606;
				param[3] = -0.232071;
				param[4] = 0.913448;
				param[5] = 0.047339;
				param[6] = 0.317981;
				param[7] = -0.0230869;
				param[8] = 1.60531;
				param[9] = 0.277375;
				param[10] = 2.67739; 
				param[11] = 0.029935; 
				param[12] = 3.07; 
			}
			if(kinType=="kin1b") {
			//	param[0] = 3.24916;  // newtarget_muscext
			//	param[1] = 0.0282571;
			//	param[2] = 2.0741;
			//	param[3] = -0.258393;
			//	param[4] = 0.814932;
			//	param[5] = 0.0404905;
			//	param[6] = 0.281178;
			//	param[7] = -0.010876;
			//	param[8] = 1.59693;
			//	param[9] = 0.290196;
			//	param[10] =2.77276;
			//	param[11] =0.0268079;
			//	param[12] = 2.65;
				param[0] = 3.34923;  // new vovilt 3 runs
				param[1] = 0.0524643;
				param[2] = 2.07899;
				param[3] = -0.256065;
				param[4] = 0.835706;
				param[5] = 0.0394785;
				param[6] = 0.255668;
				param[7] = -0.0200462;
				param[8] = 1.59625;
				param[9] = 0.288017;
				param[10] = 2.73451;
				param[11] = 0.0309081;
				param[12] = 2.65;

			}
			if(kinType=="kin2a") {
			//	param[0] =  3.1797;  // newtarget_muscext
			//	param[1] =  0.0505098;
			//	param[2] =  2.07386;
			//	param[3] =  -0.221387;
			//	param[4] =  1.00431;
			//	param[5] =  0.0365711;
			//	param[6] =  0.349406;
			//	param[7] =  -0.0215562;
			//	param[8] =  1.60841;
			//	param[9] =  0.271537;
			//	param[10] = 2.62987;
			//	param[11] = 0.033551;
			//	param[12] = 3.41;
				param[0] =  3.14062;  // new vovilt 3 runs
				param[1] =  0.0454896;
				param[2] =  2.0738;
				param[3] =  -0.219954;
				param[4] =  1.0502;
				param[5] =  0.0296818;
				param[6] =  0.387476;
				param[7] =  -0.0195034;
				param[8] =  1.60972;
				param[9] =  0.2686;
				param[10] = 2.58439;
				param[11] = 0.0395248;
				param[12] = 3.33;

			}
			if(kinType=="kin2b") {
			//	param[0] =  3.31985;  // newtarget_muscext
			//	param[1] =  0.0402951;
			//	param[2] =  2.08511;
			//	param[3] =  -0.260999;
			//	param[4] =  0.723082;
			//	param[5] =  0.0565357;
			//	param[6] =  0.261602;
			//	param[7] =  -0.0160132;
			//	param[8] =  1.59057;
			//	param[9] =  0.28935;
			//	param[10] = 2.773;
			//	param[11] = 0.0238289;
			//	param[12] = 2.65;
				param[0] =  3.36696;  // new vovilt 3 runs
				param[1] =  0.0503432;
				param[2] =  2.08082;
				param[3] =  -0.257484;
				param[4] =  0.875552;
				param[5] =  0.030451;
				param[6] =  0.218867;
				param[7] =  -0.0258478;
				param[8] =  1.59209;
				param[9] =  0.286585;
				param[10] = 2.73218;
				param[11] = 0.0282284;
				param[12] = 2.61;
			
			}
			if(kinType=="kin3b") {
			//	param[0] =  3.52694;  // newtarget_muscext
			//	param[1] =  0.0786856;
			//	param[2] =  2.10858;
			//	param[3] =  -0.255974;
			//	param[4] =  1.02366;
			//	param[5] =  0.00821375;
			//	param[6] =  0.174108;
			//	param[7] =  -0.0378911;
			//	param[8] =  1.56547;
			//	param[9] =  0.277456;
			//	param[10] = 2.64894;
			//	param[11] = 0.027786;
			//	param[12] = 2.83;
				param[0] =  3.45638;  // new vovilt 3 runs
				param[1] =  0.0663388;
				param[2] =  2.10954;
				param[3] =  -0.252933;
				param[4] =  0.868459;
				param[5] =  0.0429854;
				param[6] =  0.175194;
				param[7] =  -0.0409396;
				param[8] =  1.56545;
				param[9] =  0.274313;
				param[10] = 2.79369;
				param[11] = -0.00638383;
				param[12] = 2.79;

			}
			return param;
		}

		//----------------------------------------------------calculate CM angles-----------------------------------------------------//

	//	void get_angles(TLorentzVector k1, TLorentzVector k2,TLorentzVector p1, TLorentzVector p2,TLorentzVector q1, TLorentzVector q2, double &a, double &b){
		void get_angles(TLorentzVector k1, TLorentzVector k2,TLorentzVector p1, TLorentzVector p2,TLorentzVector q1, TLorentzVector q2, vector<double>* angles){

			//		boost to CM frame	
			TLorentzVector tot;
			tot = q1 + p1;

			TVector3 boost;
			boost = tot.BoostVector();

			k1.Boost(-boost);
			q1.Boost(-boost);
			p1.Boost(-boost);
			k2.Boost(-boost);
			p2.Boost(-boost);


			//		rotate q_initial to make it along positive z direction

			TVector3 unit_z;
			unit_z.SetXYZ(0.,0.,1.);

			double angle1 = q1.Angle(unit_z);

			TVector3 rotate_axis;
			rotate_axis = (q1.Vect()).Cross(unit_z);

			k1.Rotate(angle1, rotate_axis);
			q1.Rotate(angle1, rotate_axis);
			p1.Rotate(angle1, rotate_axis);
			k2.Rotate(angle1, rotate_axis);
			p2.Rotate(angle1, rotate_axis);

			//		rotate scattered electron to x-z plane
			double angle2 = k2.Phi();
			k1.RotateZ(-angle2);
			q1.RotateZ(-angle2);
			p1.RotateZ(-angle2);
			k2.RotateZ(-angle2);
			p2.RotateZ(-angle2);

			//		set values to q_final
			//q2.SetXYZM(-(p2.X()), -(p2.Y()),-(p2.Z()),0.);

			TLorentzVector mm = k1 + p1 - k2 - p2;
			double mm2 = mm.M2();
			q2=mm;
			//		find theta_gg and phi_gg
			double thetagg_rad = q1.Angle(q2.Vect());
			double thetagg = thetagg_rad * TMath::RadToDeg();

			double phigg_rad = q2.Phi() + TMath::Pi();
			double phigg = phigg_rad * TMath::RadToDeg();
			
			angles->push_back(thetagg);
			angles->push_back(phigg);
			angles->push_back(mm2);
		//	a = thetagg;
		//	b = phigg;

		}
		vector<double> get_cmangles(TLorentzVector k1, TLorentzVector k2,TLorentzVector p1, TLorentzVector p2,TLorentzVector q1, TLorentzVector q2){
				
			vector<double> angles;
			//		boost to CM frame	
			TLorentzVector tot;
			tot = q1 + p1;

			TVector3 boost;
			boost = tot.BoostVector();

			k1.Boost(-boost);
			q1.Boost(-boost);
			p1.Boost(-boost);
			k2.Boost(-boost);
			p2.Boost(-boost);


			//		rotate q_initial to make it along positive z direction

			TVector3 unit_z;
			unit_z.SetXYZ(0.,0.,1.);

			double angle1 = q1.Angle(unit_z);

			TVector3 rotate_axis;
			rotate_axis = (q1.Vect()).Cross(unit_z);

			k1.Rotate(angle1, rotate_axis);
			q1.Rotate(angle1, rotate_axis);
			p1.Rotate(angle1, rotate_axis);
			k2.Rotate(angle1, rotate_axis);
			p2.Rotate(angle1, rotate_axis);

			//		rotate scattered electron to x-z plane
			double angle2 = k2.Phi();
			k1.RotateZ(-angle2);
			q1.RotateZ(-angle2);
			p1.RotateZ(-angle2);
			k2.RotateZ(-angle2);
			p2.RotateZ(-angle2);

			//		set values to q_final
			//q2.SetXYZM(-(p2.X()), -(p2.Y()),-(p2.Z()),0.);

			TLorentzVector mm = k1 + p1 - k2 - p2;
			double mm2 = mm.M2();
			q2=mm;
			//		find theta_gg and phi_gg
			double thetagg_rad = q1.Angle(q2.Vect());
			double thetagg = thetagg_rad * TMath::RadToDeg();

			double phigg_rad = q2.Phi() + TMath::Pi();
			double phigg = phigg_rad * TMath::RadToDeg();
			
			angles.push_back(thetagg);
			angles.push_back(phigg);
		//	a = thetagg;
		//	b = phigg;
			return angles;
		}
	private:
		double ps6;
		double arr[6];
		double eff[5];
		double stat[4];
		double peak[2];
		double param[13];
		double uncertain;
		double ratio;
		double counts;
		TLorentzVector ret;
		TH1F* h_ave;

};

