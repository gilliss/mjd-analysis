// A script to plot uniform 1D hists for all chans and them add them together. 
// Fit the ToE projection to get that cut first. 
// Save all hists in a root file. ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch03s11.html
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

	TFile * originalFile = TFile::Open("/global/u2/g/gilliss/dev/Analysis/Cryo1/SpectralAnalysis68Ge/AnalysisWithKrisRootFile/m1Diagnostic19Oct2015.root");
	//TH2F *h = (TH2F*)originalFile->Get(histname); // pull out the original TH2F // P7D2 is a good example
	
	TFile *f = new TFile("ToE_E_ProjectionHists_4sigmaToE_halfkeVBins_2.root","NEW"); // http://www-numi.fnal.gov/offline_software/srt_public_context/WebDocs/Companion/root_crib/tfile.html
	//TFile f("ToEProjectionHistsTest.root","new"); // create a new root file to hold the hists made below ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch03s11.html
	f->Close();
	
for(int i=0; i<29; i++)
{
	char detname[200], histname[200];
	if(i==0){sprintf(histname,"cutToEP1D2"); sprintf(detname,"P1D2");};
	if(i==1){sprintf(histname,"cutToEP5D2"); sprintf(detname,"P5D2");};
	if(i==2){sprintf(histname,"cutToEP7D4"); sprintf(detname,"P7D4");};
	if(i==3){sprintf(histname,"cutToEP2D2"); sprintf(detname,"P2D2");};
	if(i==4){sprintf(histname,"cutToEP1D1"); sprintf(detname,"P1D1");};
	if(i==5){sprintf(histname,"cutToEP4D3"); sprintf(detname,"P4D3");};
	if(i==6){sprintf(histname,"cutToEP6D2"); sprintf(detname,"P6D2");};
	if(i==7){sprintf(histname,"cutToEP5D1"); sprintf(detname,"P5D1");};
	if(i==8){sprintf(histname,"cutToEP6D4"); sprintf(detname,"P6D4");};
	if(i==9){sprintf(histname,"cutToEP2D3"); sprintf(detname,"P2D3");};
	if(i==10){sprintf(histname,"cutToEP6D1"); sprintf(detname,"P6D1");};
	if(i==11){sprintf(histname,"cutToEP5D3"); sprintf(detname,"P5D3");};
	if(i==12){sprintf(histname,"cutToEP4D5"); sprintf(detname,"P4D5");};
	if(i==13){sprintf(histname,"cutToEP4D1"); sprintf(detname,"P4D1");};
	if(i==14){sprintf(histname,"cutToEP2D4"); sprintf(detname,"P2D4");};
	if(i==15){sprintf(histname,"cutToEP3D3"); sprintf(detname,"P3D3");};
	if(i==16){sprintf(histname,"cutToEP1D4"); sprintf(detname,"P1D4");};
	if(i==17){sprintf(histname,"cutToEP3D1"); sprintf(detname,"P3D1");};
	if(i==18){sprintf(histname,"cutToEP3D2"); sprintf(detname,"P3D2");};
	if(i==19){sprintf(histname,"cutToEP7D2"); sprintf(detname,"P7D2");};
	if(i==20){sprintf(histname,"cutToEP1D3"); sprintf(detname,"P1D3");};
	if(i==21){sprintf(histname,"cutToEP3D4"); sprintf(detname,"P3D4");};
	if(i==22){sprintf(histname,"cutToEP7D3"); sprintf(detname,"P7D3");};
	if(i==23){sprintf(histname,"cutToEP4D2"); sprintf(detname,"P4D2");};
	if(i==24){sprintf(histname,"cutToEP7D1"); sprintf(detname,"P7D1");};
	if(i==25){sprintf(histname,"cutToEP5D4"); sprintf(detname,"P5D4");};
	if(i==26){sprintf(histname,"cutToEP2D1"); sprintf(detname,"P2D1");};
	if(i==27){sprintf(histname,"cutToEP6D3"); sprintf(detname,"P6D3");};
	if(i==28){sprintf(histname,"cutToEP4D4"); sprintf(detname,"P4D4");};
	//cout << histname << " " << detname << endl;

	char ToEHistTitle[200];
	sprintf(ToEHistTitle,"ToE Projection for %s, Hit Spectrum",detname);
	char ToEHistNameKey[200];
	sprintf(ToEHistNameKey,"ToEProjection%s",detname);
	char EHistTitle[200];
	sprintf(EHistTitle,"E Projection for %s, Hit Spectrum",detname);
	char EHistNameKey[200];
	sprintf(EHistNameKey,"EProjection%s",detname);


	TH2F *h = (TH2F*)originalFile->Get(histname); // pull out the original TH2F // P7D2 is a good example
	
	/////
	///// ToE HIST
	/////
	TAxis *xaxisToE = h->GetXaxis(); 
	//Int_t biny1ToE = xaxisToE->FindBin(0); 
	Int_t biny2ToE = xaxisToE->FindBin(100); 
	//Int_t bin1ToE = h->FindFirstBinAbove(0,1); // first x-axis(axis 1) bin with nonzero value
	//Int_t bin2ToE = h->FindLastBinAbove(0,1); // last x-axis(axis 1) bin with nonzero value 
	//h->ProjectionY("energyCal",biny1ToE,biny2ToE)->Draw();
	// rebin? 
	TH1D *hToE = h->ProjectionY(ToEHistNameKey,0,biny2ToE); // make a TH1F out of the projection
		hToE->SetTitle(ToEHistTitle);		
		hToE->GetXaxis()->SetTitle("ToE");
		hToE->GetYaxis()->SetTitle("Counts");
		hToE->GetXaxis()->SetRangeUser(0,5); // input axis values, not bin numbers like for TAxis::SetRange
	hToE->Draw();
        double peakpos = hToE->GetBinCenter(hToE->GetMaximumBin()); // guess peak by max bin in window
	double peakheight = hToE->GetBinContent(hToE->FindBin(peakpos)); // get the height of the max bin
	//cout << "peakpos " << peakpos << endl;
	//cout << "peakheight " << peakheight << endl;
	
        TF1 * fGausFit = new TF1("GausFit", "[0]*Gaus(x,[1],[2]) + [3] + [4]*x ",0,1000);  // TMath::Gaus(x,[mean],[sigma]) , "0,1000" are
											   // https://root.cern.ch/root/html534/TMath.html#TMath:Gaus
                fGausFit->SetRange(0.8*peakpos,1.2*peakpos);
                fGausFit->SetLineColor(kRed);
                fGausFit->SetLineWidth(1);
                fGausFit->SetLineStyle(2);
                fGausFit->FixParameter(0,peakheight*0.92); // gaus amplitude
                fGausFit->SetParameter(1,peakpos); // gaus mean
                fGausFit->SetParameter(2,peakpos/20.0); // gaus sigma // maybe use a FWHM=2.3548*sigma -> sigma=FWHM/2.3548 conversion for a guess?
                fGausFit->SetParameter(3,0); // flat BG
                fGausFit->SetParameter(4,0); // linear BG
                
	hToE->Fit("GausFit","qr+"); // for fit options https://root.cern.ch/doc/master/classTH1.html#a63eb028df86bc86c8e20c989eb23fb2a
/*
	cout << "--------"<<detname<<" FIT RESULTS--------" << endl;                      
	cout << "peakpos input " << peakpos << " | fit mean " << fGausFit->GetParameter(1) << " +/- " << fGausFit->GetParError(1) << endl;
	cout << "peakheight*.92 input " << peakheight<< "*0.92" << " | fit height " << fGausFit->GetParameter(0) << " +/- " << fGausFit->GetParError(0) << endl;
	cout << "peaksigma input " << peakpos/20.0 << " | fit sigma " << fGausFit->GetParameter(2) << " +/- " << fGausFit->GetParError(2) << endl;
*/
	double meanToE = fGausFit->GetParameter(1);
	double sigmaToE = fGausFit->GetParameter(2);

	if(sigmaToE < 0)
	{
		sigmaToE=-sigmaToE;
	}

	//f->cd();
	TFile ff("ToE_E_ProjectionHists_4sigmaToE_halfkeVBins_2.root","UPDATE"); // https://root.cern.ch/root/HowtoWrite.html
	hToE->Write(); // ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch03s11.html
	//f->Write();
	//f->Close();
	
		

	/////
	///// ENERGY HIST
	/////
	TAxis *yaxis = h->GetYaxis(); 
	Int_t biny1 = yaxis->FindBin(meanToE-4*sigmaToE);//(1.0);   // 1sigma=68.27%, 2sigma=95.45%, 3sigma=99.73% 
	Int_t biny2 = yaxis->FindBin(meanToE+4*sigmaToE);//(1.6);
	TH1D *hInter = h->ProjectionX("EProjectionIntermediate",biny1,biny2);
	TH1D *h1 = new TH1D(EHistNameKey, EHistTitle, 200, 0, 100);
	//TH1D *h1 = (TH1D*)hInter->Rebin(15,EHistNameKey); // hInter has 3000(3002) bins // https://root.cern.ch/root/html520/TH1.html#TH1:Rebin
	h1 = (TH1D*)hInter->Rebin(15,EHistNameKey); // hInter has 3000(3002) bins // https://root.cern.ch/root/html520/TH1.html#TH1:Rebin
		h1->SetTitle(EHistTitle);		
		h1->GetXaxis()->SetTitle("Energy (keV)");
		h1->GetYaxis()->SetTitle("Counts");
		h1->GetXaxis()->SetRangeUser(0,50); // input axis values, not bin numbers like for TAxis::SetRange
	h1->Draw();

	//TFile ff("ToEProjectionHistsTest.root","UPDATE"); // https://root.cern.ch/root/HowtoWrite.html
	h1->Write(); // ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch03s11.html


} // End of loop over Kris' TH2F hists






	cout << "F I N I S H E D" << endl;

	App->Run();
} // End of main()


 
/* the previous energy hist section
	/////
	///// ENERGY HIST
	/////
	TAxis *yaxis = h->GetYaxis(); 
	Int_t biny1 = yaxis->FindBin(meanToE-4*sigmaToE);//(1.0);   // 1sigma=68.27%, 2sigma=95.45%, 3sigma=99.73% 
	Int_t biny2 = yaxis->FindBin(meanToE+4*sigmaToE);//(1.6);
	h->ProjectionX("EProjectionIntermediate",biny1,biny2)->Draw(); // project the TH2F onto an axis
	TH2 *hnew = h->RebinX(15,"hnew"); // rebin the projected TH2F
	//hnew->Draw(); 
	TH1D *h1 = hnew->ProjectionX(EHistNameKey,biny1,biny2); // make a TH1F out of the projection
		h1->SetTitle(EHistTitle);		
		h1->GetXaxis()->SetTitle("Energy (keV)");
		h1->GetYaxis()->SetTitle("Counts");
		h1->GetXaxis()->SetRangeUser(0,50); // input axis values, not bin numbers like for TAxis::SetRange
	h1->Draw();

	//TFile ff("ToEProjectionHistsTest.root","UPDATE"); // https://root.cern.ch/root/HowtoWrite.html
	h1->Write(); // ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch03s11.html
*/
