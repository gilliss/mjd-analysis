// Script to calculate exposure and to fit the added histograms produced by AddHists_KrisRootFile.cc
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TChain.h>
#include <TEventList.h>
#include <TEntryList.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>

#include <TApplication.h>  //This class creates the ROOT Application Environment that interfaces to the windowing system eventloop and eventhandlers

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <string>
//#include <stdlib.h> 	// atof, atoi
//#include <iomanip>      // std::setprecision

#include <GATBaseClassesDICT.h>   //#include <GATDataSet.hh>   //#include <GATPeakShape.hh>
#include <MJTChannelSettings.hh>
#include <MJTChannelMap.hh>
//#include <MGDOUtils.hh>

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Body
//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	TApplication *App = new TApplication("App", 0, NULL);

////////////////////////////////////
//////// CALCULATE EXPOSURES FIRST
////////////////////////////////////

	TFile * originalFile2 = TFile::Open("/global/u2/g/gilliss/dev/Analysis/Cryo1/SpectralAnalysis68Ge/AnalysisWithKrisRootFile/m1Diagnostic19Oct2015.root");

	char detname[200], histname[200]; 	
	int chan;
	double exposureTimeS, exposureTimeDy, massg, masskg; // http://mjwiki.npl.washington.edu/pub/Majorana/AnalysisReports/20151020_JGruszko_Alpha_Part3.pdf
	double exposureKgDyNat = 0;
	double exposureKgDyEnr = 0;

for(int i=0; i<29; i++)
{
	if(i==0){sprintf(histname,"EProjectionP1D2"); sprintf(detname,"P1D2"); massg=1055.4; chan =690;};
	if(i==1){sprintf(histname,"EProjectionP5D2"); sprintf(detname,"P5D2"); massg=796.0; chan =662;};
	if(i==2){sprintf(histname,"EProjectionP7D4"); sprintf(detname,"P7D4"); massg=1043.6; chan =1000;};
	if(i==3){sprintf(histname,"EProjectionP2D2"); sprintf(detname,"P2D2"); massg=809.3; chan =674;};
	if(i==4){sprintf(histname,"EProjectionP1D1"); sprintf(detname,"P1D1"); massg=562.0; chan =692;};
	if(i==5){sprintf(histname,"EProjectionP4D3"); sprintf(detname,"P4D3"); massg=633.0; chan =600;};
	if(i==6){sprintf(histname,"EProjectionP6D2"); sprintf(detname,"P6D2"); massg=749.8; chan =626;};
	if(i==7){sprintf(histname,"EProjectionP5D1"); sprintf(detname,"P5D1"); massg=617.0; chan =664;};
	if(i==8){sprintf(histname,"EProjectionP6D4"); sprintf(detname,"P6D4"); massg=633.7; chan =1000;};
	if(i==9){sprintf(histname,"EProjectionP2D3"); sprintf(detname,"P2D3"); massg=717.5; chan =576;};
	if(i==10){sprintf(histname,"EProjectionP6D1"); sprintf(detname,"P6D1"); massg=802.8; chan =628;};
	if(i==11){sprintf(histname,"EProjectionP5D3"); sprintf(detname,"P5D3"); massg=694.0; chan =656;};
	if(i==12){sprintf(histname,"EProjectionP4D5"); sprintf(detname,"P4D5"); massg=621.0; chan =592;};
	if(i==13){sprintf(histname,"EProjectionP4D1"); sprintf(detname,"P4D1"); massg=622.0; chan =608;};
	if(i==14){sprintf(histname,"EProjectionP2D4"); sprintf(detname,"P2D4"); massg=749.2; chan =680;};
	if(i==15){sprintf(histname,"EProjectionP3D3"); sprintf(detname,"P3D3"); massg=1025.2; chan =614;};
	if(i==16){sprintf(histname,"EProjectionP1D4"); sprintf(detname,"P1D4"); massg=1037.9; chan =640;};
	if(i==17){sprintf(histname,"EProjectionP3D1"); sprintf(detname,"P3D1"); massg=615.0; chan =1000;};
	if(i==18){sprintf(histname,"EProjectionP3D2"); sprintf(detname,"P3D2"); massg=962.0; chan =1000;};
	if(i==19){sprintf(histname,"EProjectionP7D2"); sprintf(detname,"P7D2"); massg=778.9; chan =644;};
	if(i==20){sprintf(histname,"EProjectionP1D3"); sprintf(detname,"P1D3"); massg=897.7; chan =688;};
	if(i==21){sprintf(histname,"EProjectionP3D4"); sprintf(detname,"P3D4"); massg=1051.7; chan =610;};
	if(i==22){sprintf(histname,"EProjectionP7D3"); sprintf(detname,"P7D3"); massg=625.7; chan =642;};
	if(i==23){sprintf(histname,"EProjectionP4D2"); sprintf(detname,"P4D2"); massg=629.0; chan =598;};
	if(i==24){sprintf(histname,"EProjectionP7D1"); sprintf(detname,"P7D1"); massg=626.0; chan =646;};
	if(i==25){sprintf(histname,"EProjectionP5D4"); sprintf(detname,"P5D4"); massg=1068.0; chan =696;};
	if(i==26){sprintf(histname,"EProjectionP2D1"); sprintf(detname,"P2D1"); massg=625.0; chan =1000;};
	if(i==27){sprintf(histname,"EProjectionP6D3"); sprintf(detname,"P6D3"); massg=764.0; chan =624;};
	if(i==28){sprintf(histname,"EProjectionP4D4"); sprintf(detname,"P4D4"); massg=608.0; chan =594;};
	//cout << histname << " " << detname << endl;

	TH1D *h = (TH1D*)originalFile2->Get("exposure"); // pull out the original TH1D

	////
	//// Natural
	////
	if(i==12 || i==28 || i==23 || i==5 || i==13 || i==24 || i==7)
	{
		TAxis *xaxis = h->GetXaxis();
		Int_t binChan = xaxis->FindBin(chan);
		exposureTimeS=h->Integral(binChan,binChan+1);
		exposureTimeDy = exposureTimeS/86400.0; // in days
		masskg=massg/1000;
		exposureKgDyNat += masskg*exposureTimeDy;	
		cout << detname << " " << exposureTimeS/3600.0 << " hrs " << exposureKgDyNat << " kg*dy" << " Natural" << endl;		
	}


	////
	//// Enriched
	////
	if(i==3 || i==20 || i==0 || i==4 || i==25 || i==9 || i==21 || i==27 || i==6 || i==16 || i==22 || i==19 || i==11 || i==1)
	{
		TAxis *xaxis = h->GetXaxis();
		Int_t binChan = xaxis->FindBin(chan);
		exposureTimeS=h->Integral(binChan,binChan+1);
		exposureTimeDy = exposureTimeS/86400.0; // in days
		masskg=massg/1000;
		exposureKgDyEnr += masskg*exposureTimeDy;
		cout << detname << " " << exposureTimeS/3600.0 << " hrs " << exposureKgDyEnr << " kg*dy" << " Enriched" << endl;				
	}
	
	delete h;		


} // End of loop over energy hists
originalFile2->Close();

////////////////////////////////////
//////// PLOT ADDED HISTS
////////////////////////////////////

	//TFile * originalFile3 = TFile::Open("/global/u2/g/gilliss/dev/Analysis/Cryo1/SpectralAnalysis68Ge/AnalysisWithKrisRootFile/AddedHists_4sigmaToE_halfkeVBins.root");
	//originalFile3->GetListOfKeys()->Print();
	TFile fff("/global/u2/g/gilliss/dev/Analysis/Cryo1/SpectralAnalysis68Ge/AnalysisWithKrisRootFile/AddedHists_4sigmaToE_halfkeVBins.root");
	
	////
	//// Natural
	////
	TCanvas *cNat = new TCanvas;
	cNat->cd();
	TH1D *hNat = (TH1D*)fff.Get("NaturalDetsAdded"); // //TH1D *hNat = (TH1D*)originalFile3->Get("NaturalDetsAdded"); // pull out the original TH1D
		char NaturalHistTitle[200]; sprintf(NaturalHistTitle,"Natural Dets, Hit Spectrum, %.1f kg*dy",exposureKgDyNat);
		hNat->SetTitle(NaturalHistTitle);			
		hNat->GetYaxis()->SetTitle("Cnts/kg/day/0.5keV");
		hNat->Scale(1/exposureKgDyNat); // puts histo in counts/kg/dy (work with binning for "per keV")
        	hNat->SetStats(0);
	hNat->Draw();
	cNat->SetLogy();

		double peakpos = 10.37; // guess peak

	        TF1 * fGausFit = new TF1("GausFit", "[0]*Gaus(x,[1],[2]) + [3] + [4]*x + [5]*Gaus(x,[6],[7])",0,1000);  // TMath::Gaus(x,[mean],[sigma]) , "0,1000" are
											   // https://root.cern.ch/root/html534/TMath.html#TMath:Gaus
                fGausFit->SetRange(0.53*peakpos,1.4*peakpos);//first slides: (0.6*peakpos,1.4*peakpos); // 0.8, 1.2
                fGausFit->SetLineColor(kRed);
                fGausFit->SetLineWidth(1);
                fGausFit->SetLineStyle(2);
                fGausFit->SetParameter(0,0.0); // gaus amplitude
                fGausFit->SetParameter(1,peakpos); // gaus mean
                fGausFit->SetParameter(2,peakpos/20.0); // gaus sigma // maybe use a FWHM=2.3548*sigma -> sigma=FWHM/2.3548 conversion for a guess?
                fGausFit->SetParameter(3,0); // flat BG
                fGausFit->SetParameter(4,0); // linear BG
                fGausFit->SetParameter(5,0.0); // gaus2 amplitude
                fGausFit->SetParameter(6,6.6); // gaus2 mean
                fGausFit->SetParameter(7,6.6/20.0); // gaus2 sigma


                
		hNat->Fit("GausFit","qr+"); // for fit options https://root.cern.ch/doc/master/classTH1.html#a63eb028df86bc86c8e20c989eb23fb2a
		
		cout << "--------Natural FIT RESULTS--------" << endl;                      
		cout << "peakpos input " << peakpos << " | fit mean " << fGausFit->GetParameter(1) << " +/- " << fGausFit->GetParError(1) << endl;
		cout << "peakheight input " << "0"<< " " << " | fit height " << fGausFit->GetParameter(0) << " +/- " << fGausFit->GetParError(0) << endl;
		cout << "peaksigma input " << peakpos/20.0 << " | fit sigma " << fGausFit->GetParameter(2) << " +/- " << fGausFit->GetParError(2) << endl;

		cout << "BACKGROUND AND INTEGRAL RESULTS" << endl;
		double sigmaNat = fGausFit->GetParameter(2);
		double meanNat = fGausFit->GetParameter(1);		
		TAxis *xaxis = hNat->GetXaxis();
		Int_t binx1 = xaxis->FindBin(meanNat-3*sigmaNat);
		Int_t binx2 = xaxis->FindBin(meanNat+3*sigmaNat);
		cout << "3sigma integral bounds " << meanNat-3*sigmaNat << " to " << meanNat+3*sigmaNat << endl;
		cout << "3 sigma Integral for gaussian portion = " << hNat->Integral(binx1,binx2) << endl;
		cout << "Linear BG: y="<< fGausFit->GetParameter(4) << " +/- " << fGausFit->GetParError(4) << "x+" <<  fGausFit->GetParameter(3)<< " +/- " << fGausFit->GetParError(3) << endl; 

		cout << "--------Natural SecondPeak-FIT RESULTS--------" << endl;                      
		cout << "peakpos input " << 6.6 << " | fit mean " << fGausFit->GetParameter(6) << " +/- " << fGausFit->GetParError(6) << endl;
		cout << "peakheight input " << "0"<< " " << " | fit height " << fGausFit->GetParameter(5) << " +/- " << fGausFit->GetParError(5) << endl;
		cout << "peaksigma input " << 6.6/20.0 << " | fit sigma " << fGausFit->GetParameter(7) << " +/- " << fGausFit->GetParError(7) << endl;

		cout << "SecondPeak-FIT BACKGROUND AND INTEGRAL RESULTS" << endl;
		sigmaNat = fGausFit->GetParameter(7);
		meanNat = fGausFit->GetParameter(6);		
		binx1 = xaxis->FindBin(meanNat-3*sigmaNat);
		binx2 = xaxis->FindBin(meanNat+3*sigmaNat);
		cout << "3sigma integral bounds " << meanNat-3*sigmaNat << " to " << meanNat+3*sigmaNat << endl;
		cout << "3 sigma Integral for gaussian portion = " << hNat->Integral(binx1,binx2) << endl;
		cout << "Linear BG: y="<< fGausFit->GetParameter(4) << " +/- " << fGausFit->GetParError(4) << "x+" <<  fGausFit->GetParameter(3)<< " +/- " << fGausFit->GetParError(3) << endl; 


	////
	//// Enriched
	////
	TCanvas *cEnr = new TCanvas;
	cEnr->cd();
	TH1D *hEnr = (TH1D*)fff.Get("EnrichedDetsAdded"); // //TH1D *hEnr = (TH1D*)originalFile3->Get("EnrichedDetsAdded"); // pull out the original TH1D
		char EnrichedHistTitle[200]; sprintf(EnrichedHistTitle,"Enriched Dets, Hit Spectrum, %.1f kg*dy",exposureKgDyEnr);
		hEnr->SetTitle(EnrichedHistTitle);		
		hEnr->GetYaxis()->SetTitle("Cnts/kg/day/0.5keV");
		hEnr->Scale(1/exposureKgDyEnr); // puts histo in counts/kg/dy (work with binning for "per keV")
        	hEnr->SetStats(0);
	hEnr->Draw();
	cEnr->SetLogy();
                
/* CONSIDER A PIECEWISE FIT HERE WITH A GAUSSIAN AND THEN ~LINEAR TO THE LEFT AND RIGHT. Have to think about peak height - BG...
                fGausFit->SetRange(9.0,12.0); // 0.8, 1.2
                fGausFit->SetParameter(0,.01); // gaus amplitude
                fGausFit->SetParameter(2,peakpos/20.0); // gaus sigma // maybe use a FWHM=2.3548*sigma -> sigma=FWHM/2.3548 conversion for a guess?

		hEnr->Fit("GausFit","qr+"); // for fit options https://root.cern.ch/doc/master/classTH1.html#a63eb028df86bc86c8e20c989eb23fb2a
*/
		//TF1 *g1 = new TF1("g1","gaus",9,12);
                //g1->SetLineColor(kRed);
                //g1->SetLineWidth(1);
                //g1->SetLineStyle(2);

	        TF1 * g1 = new TF1("g1Fit", "[0]*Gaus(x,[1],[2])",0,1000);  // TMath::Gaus(x,[mean],[sigma]) , "0,1000" are
									   // https://root.cern.ch/root/html534/TMath.html#TMath:Gaus
                	g1->SetRange(9,12); // 0.8, 1.2
                	g1->SetLineColor(kBlue);
                	g1->SetLineWidth(1);
                	g1->SetLineStyle(2);
                	//g1->SetParameter(0,0.05); // gaus amplitude
			g1->FixParameter(0,hEnr->GetBinContent(hEnr->FindBin(10.6)));
                	g1->SetParameter(1,peakpos); // gaus mean
                	g1->SetParameter(2,peakpos/20.0); // gaus sigma // maybe use a FWHM=2.3548*sigma -> sigma=FWHM/2.3548 conversion for a guess?
			hEnr->Fit(g1,"R");

	        TF1 * g2 = new TF1("g1Fit", "[0]",0,1000);  // TMath::Gaus(x,[mean],[sigma]) , "0,1000" are
									   // https://root.cern.ch/root/html534/TMath.html#TMath:Gaus
                	g2->SetRange(12,35); // 0.8, 1.2
                	g2->SetLineColor(kBlue);
                	g2->SetLineWidth(1);
                	g2->SetLineStyle(2);
                	g2->SetParameter(0,0.0); // gaus amplitude
			hEnr->Fit(g2,"R+");

		cout << "--------Enriched Gaus FIT RESULTS--------" << endl;                      
		cout << "peakpos input " << peakpos << " | fit mean " << g1->GetParameter(1) << " +/- " << g1->GetParError(1) << endl;
		cout << "peakheight input fixed" << hEnr->GetBinContent(hEnr->FindBin(10.6)) << " " << " | fit height " << g1->GetParameter(0) << " +/- " << g1->GetParError(0) << endl;
		cout << "peaksigma input " << peakpos/20.0 << " | fit sigma " << g1->GetParameter(2) << " +/- " << g1->GetParError(2) << endl;

		cout << "--------Enriched Flat FIT RESULTS--------" << endl;                      
		cout << "constant input " << "0" << " | fit constant " << g2->GetParameter(0) << " +/- " << g2->GetParError(0) << endl;

		double sigmaEnr = g1->GetParameter(2);
		double meanEnr = g1->GetParameter(1);
		double flatEnr = g2->GetParameter(0);		

		TAxis *xaxis_2 = hEnr->GetXaxis();
		Int_t binx1_2 = xaxis->FindBin(meanEnr-3*sigmaEnr);
		Int_t binx2_2 = xaxis->FindBin(meanEnr+3*sigmaEnr);
		cout << "3sigma integral bounds " << meanEnr-3*sigmaEnr << " to " << meanEnr+3*sigmaEnr << endl;
		cout << "3 sigma Integral for gaussian portion = " << hEnr->Integral(binx1_2,binx2_2) << endl;
		cout << "Linear BG: y="<< "0" << "x+" << flatEnr << "+/-"<< g2->GetParError(0) << endl; 
		cout << "3 sigma Integral for flat portion = " << flatEnr*(2*(3*sigmaEnr)) << endl;
		

/*			
		Double_t par[4];
		g1->GetParameters(&par[0]);
		g2->GetParameters(&par[3]);
*/
/*
		TF1 * TotalPiece = new TF1("TotalPieceFit","[0]*Gaus(x,[1],[2])+[3]",9,17); // https://root.cern.ch/root/html/tutorials/fit/multifit.C.html
                	//TotalPiece->SetRange(9,17); // 0.8, 1.2
			TotalPiece->SetParameters(par);
                	TotalPiece->SetLineColor(kRed);
                	TotalPiece->SetLineWidth(1);
                	TotalPiece->SetLineStyle(2);


		hEnr->Fit(TotalPiece,"R+");
*/
		/*
		cout << "--------Enriched FIT RESULTS--------" << endl;                      
		cout << "peakpos input " << peakpos << " | fit mean " << fGausFit->GetParameter(1) << " +/- " << fGausFit->GetParError(1) << endl;
		cout << "peakheight input " << "0"<< " " << " | fit height " << fGausFit->GetParameter(0) << " +/- " << fGausFit->GetParError(0) << endl;
		cout << "peaksigma input " << peakpos/20.0 << " | fit sigma " << fGausFit->GetParameter(2) << " +/- " << fGausFit->GetParError(2) << endl;
		*/
	


	//originalFile3->Close();






	cout << "F I N I S H E D" << endl;

	App->Run();
} // End of main()




