// A script to add the standardized energy histograms produced by PlotFitSaveHists_KrisRootFile.cc
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

	TFile * originalFile = TFile::Open("/global/u2/g/gilliss/dev/Analysis/Cryo1/SpectralAnalysis68Ge/AnalysisWithKrisRootFile/ToE_E_ProjectionHists_4sigmaToE_halfkeVBins_2.root");
	////TFile originalFile("/global/u2/g/gilliss/dev/Analysis/Cryo1/SpectralAnalysis68Ge/AnalysisWithKrisRootFile/m1Diagnostic19Oct2015.root");	

	TFile *f = new TFile("AddedHists_4sigmaToE_halfkeVBins.root","NEW"); // http://www-numi.fnal.gov/offline_software/srt_public_context/WebDocs/Companion/root_crib/tfile.html
	f->Close();


	char NaturalHistTitle[200];
	sprintf(NaturalHistTitle,"Natural Dets, Hit Spectrum");
	char NaturalHistNameKey[200];
	sprintf(NaturalHistNameKey,"NaturalDetsAdded");
	char EnrichedHistTitle[200];
	sprintf(EnrichedHistTitle,"Enriched Dets, Hit Spectrum");
	char EnrichedHistNameKey[200];
	sprintf(EnrichedHistNameKey,"EnrichedDetsAdded");

	TH1D *hNatSum = new TH1D(NaturalHistNameKey,NaturalHistTitle,200,0,100);
		hNatSum->SetTitle(NaturalHistTitle);		
		hNatSum->GetXaxis()->SetTitle("Energy (keV)");
		hNatSum->GetYaxis()->SetTitle("Counts");
		hNatSum->GetXaxis()->SetRangeUser(0,50); // input axis values, not bin numbers like for TAxis::SetRange
	TH1D *hEnrSum = new TH1D(EnrichedHistNameKey,EnrichedHistTitle,200,0,100);
		hEnrSum->SetTitle(EnrichedHistTitle);		
		hEnrSum->GetXaxis()->SetTitle("Energy (keV)");
		hEnrSum->GetYaxis()->SetTitle("Counts");
		hEnrSum->GetXaxis()->SetRangeUser(0,50); // input axis values, not bin numbers like for TAxis::SetRange
	
for(int i=0; i<29; i++)
{
	char detname[200], histname[200];
	double massg; // http://mjwiki.npl.washington.edu/pub/Majorana/AnalysisReports/20151020_JGruszko_Alpha_Part3.pdf
	if(i==0){sprintf(histname,"EProjectionP1D2"); sprintf(detname,"P1D2"); massg=1055.4;};
	if(i==1){sprintf(histname,"EProjectionP5D2"); sprintf(detname,"P5D2"); massg=796.0;};
	if(i==2){sprintf(histname,"EProjectionP7D4"); sprintf(detname,"P7D4"); massg=1043.6;};
	if(i==3){sprintf(histname,"EProjectionP2D2"); sprintf(detname,"P2D2"); massg=809.3;};
	if(i==4){sprintf(histname,"EProjectionP1D1"); sprintf(detname,"P1D1"); massg=562.0;};
	if(i==5){sprintf(histname,"EProjectionP4D3"); sprintf(detname,"P4D3"); massg=633.0;};
	if(i==6){sprintf(histname,"EProjectionP6D2"); sprintf(detname,"P6D2"); massg=749.8;};
	if(i==7){sprintf(histname,"EProjectionP5D1"); sprintf(detname,"P5D1"); massg=617.0;};
	if(i==8){sprintf(histname,"EProjectionP6D4"); sprintf(detname,"P6D4"); massg=633.7;};
	if(i==9){sprintf(histname,"EProjectionP2D3"); sprintf(detname,"P2D3"); massg=717.5;};
	if(i==10){sprintf(histname,"EProjectionP6D1"); sprintf(detname,"P6D1"); massg=802.8;};
	if(i==11){sprintf(histname,"EProjectionP5D3"); sprintf(detname,"P5D3"); massg=694.0;};
	if(i==12){sprintf(histname,"EProjectionP4D5"); sprintf(detname,"P4D5"); massg=621.0;};
	if(i==13){sprintf(histname,"EProjectionP4D1"); sprintf(detname,"P4D1"); massg=622.0;};
	if(i==14){sprintf(histname,"EProjectionP2D4"); sprintf(detname,"P2D4"); massg=749.2;};
	if(i==15){sprintf(histname,"EProjectionP3D3"); sprintf(detname,"P3D3"); massg=1025.2;};
	if(i==16){sprintf(histname,"EProjectionP1D4"); sprintf(detname,"P1D4"); massg=1037.9;};
	if(i==17){sprintf(histname,"EProjectionP3D1"); sprintf(detname,"P3D1"); massg=615.0;};
	if(i==18){sprintf(histname,"EProjectionP3D2"); sprintf(detname,"P3D2"); massg=962.0;};
	if(i==19){sprintf(histname,"EProjectionP7D2"); sprintf(detname,"P7D2"); massg=778.9;};
	if(i==20){sprintf(histname,"EProjectionP1D3"); sprintf(detname,"P1D3"); massg=897.7;};
	if(i==21){sprintf(histname,"EProjectionP3D4"); sprintf(detname,"P3D4"); massg=1051.7;};
	if(i==22){sprintf(histname,"EProjectionP7D3"); sprintf(detname,"P7D3"); massg=625.7;};
	if(i==23){sprintf(histname,"EProjectionP4D2"); sprintf(detname,"P4D2"); massg=629.0;};
	if(i==24){sprintf(histname,"EProjectionP7D1"); sprintf(detname,"P7D1"); massg=626.0;};
	if(i==25){sprintf(histname,"EProjectionP5D4"); sprintf(detname,"P5D4"); massg=1068.0;};
	if(i==26){sprintf(histname,"EProjectionP2D1"); sprintf(detname,"P2D1"); massg=625.0;};
	if(i==27){sprintf(histname,"EProjectionP6D3"); sprintf(detname,"P6D3"); massg=764.0;};
	if(i==28){sprintf(histname,"EProjectionP4D4"); sprintf(detname,"P4D4"); massg=608.0;};
	//cout << histname << " " << detname << endl;

	TH1D *h = (TH1D*)originalFile->Get(histname); // pull out the original TH1D

	////
	//// Natural
	////
	if(i==12 || i==28 || i==23 || i==5 || i==13 || i==24 || i==7)
	{
		hNatSum->Add(h);		
	}


	////
	//// Enriched
	////
	if(i==3 || i==20 || i==0 || i==4 || i==25 || i==9 || i==21 || i==27 || i==6 || i==16 || i==22 || i==19 || i==11 || i==1)
	{
		hEnrSum->Add(h);		
	}
	
	delete h;		


} // End of loop over energy hists

	TFile ff("AddedHists_4sigmaToE_halfkeVBins.root","UPDATE"); // https://root.cern.ch/root/HowtoWrite.html
	hNatSum->Write(); // ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch03s11.html
	hEnrSum->Write(); // ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch03s11.html





	cout << "F I N I S H E D" << endl;

	App->Run();
} // End of main()


/*
	if(i==28) 
	{
		////TH1D *hexample = (TH1D*)originalFile.Get(histname);
		TH1D *hexample = (TH1D*)originalFile->Get(histname); // pull out the original TH1D
		hexample->Draw();
	}
	cout << histname << endl;
*/
