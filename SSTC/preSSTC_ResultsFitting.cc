/*************************************************
1. In pre-SSTC script, shift energy, apply trapETailMin … make a skimfile from this that would be run on by the SSTC script (skim: shifted trapEcal, trapETailMin good/bad, triFilter, startTime, stopTime, timestamp, channel, mH) …  plot ToEvsE … pull out ToE cut params
 
 Notes/Issues:
x-For proof of concept, start by plotting lowE, then ToE v E.
 -need to scan through values of mH to understand its useage: trapETailMin seems to work as Kris suggested, and mH is the number of HG channels hit in an event. Might need lots of data to overcome an mH==1 cut.

 Work Flow:
 

 *************************************************/

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TEventList.h>
#include <TEntryList.h>
#include <TCut.h>
#include <TMath.h>
#include <TApplication.h>  //This class creates the ROOT Application Environment that interfaces to the windowing system eventloop and eventhandlers

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h> 	// atof, atoi
#include <iomanip>      // std::setprecision
#include <utility>      // pair<Type1,Type2>
#include <string>
//#include <vector>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Mapping Functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// A FUNCTION TO CHECK FOR GOOD CHANNELS, RETURNS 1000 IF NOT A CHANNEL TO BE INCLUDED //
int CheckChannel(int ch)
{
    int cheasy = 1000;
    if(ch == 646) cheasy = 0; //HG channels, * * * BASED ON P3JDY * * *
    if(ch == 644) cheasy = 1;
    if(ch == 642) cheasy = 2;
    if(ch == 626) cheasy = 3;
    if(ch == 624) cheasy = 4;
    if(ch == 674) cheasy = 5;
    if(ch == 576) cheasy = 6;
    if(ch == 692) cheasy = 7;
    if(ch == 690) cheasy = 8;
    if(ch == 688) cheasy = 9;
    if(ch == 640) cheasy = 10;
    if(ch == 610) cheasy = 11;
    if(ch == 664) cheasy = 12;
    if(ch == 662) cheasy = 13;
    if(ch == 656) cheasy = 14;
    if(ch == 696) cheasy = 15;
    if(ch == 608) cheasy = 16;
    if(ch == 598) cheasy = 17;
    if(ch == 600) cheasy = 18;
    if(ch == 594) cheasy = 19;
    if(ch == 592) cheasy = 20;
    /*
     if(ch == 584) cheasy = 21; //Bad HV connection
     if(ch == 680) cheasy = 22; //High leakage current1
     if(ch == 676) cheasy = 23; //Bad signal connection
     if(ch == 616) cheasy = 24; //Bad signal connection
     if(ch == 614) cheasy = 25; //Gain oscillation
     if(ch == 628) cheasy = 26; //Sparking
     if(ch == 632) cheasy = 27; //High leakage current
     if(ch == 630) cheasy = 28; //Bad signal connection
     */
    return cheasy;
}
// A FUNCTION TO CHECK FOR GOOD CHANNELS, RETURNS 1000 IF NOT A CHANNEL TO BE INCLUDED //
int IndexToChannel(int ch)
{
    int cheasy = 1000;
    if(ch == 0) cheasy = 646; //HG channels, * * * BASED ON P3JDY * * *
    if(ch == 1) cheasy = 644;
    if(ch == 2) cheasy = 642;
    if(ch == 3) cheasy = 626;
    if(ch == 4) cheasy = 624;
    if(ch == 5) cheasy = 674;
    if(ch == 6) cheasy = 576;
    if(ch == 7) cheasy = 692;
    if(ch == 8) cheasy = 690;
    if(ch == 9) cheasy = 688;
    if(ch == 10) cheasy = 640;
    if(ch == 11) cheasy = 610;
    if(ch == 12) cheasy = 664;
    if(ch == 13) cheasy = 662;
    if(ch == 14) cheasy = 656;
    if(ch == 15) cheasy = 696;
    if(ch == 16) cheasy = 608;
    if(ch == 17) cheasy = 598;
    if(ch == 18) cheasy = 600;
    if(ch == 19) cheasy = 594;
    if(ch == 20) cheasy = 592;
     if(ch == 21) cheasy = 584; //Bad HV connection
     if(ch == 22) cheasy = 680; //High leakage current1
     if(ch == 23) cheasy = 676; //Bad signal connection
     if(ch == 24) cheasy = 616; //Bad signal connection
     if(ch == 25) cheasy = 614; //Gain oscillation
     if(ch == 26) cheasy = 628; //Sparking
     if(ch == 27) cheasy = 632; //High leakage current
     if(ch == 28) cheasy = 630; //Bad signal connection
    return cheasy;
}
// A FUNCTION TO SHIFT ENERGIES BASED ON ZERO-POINT CALIBRATION //
double EnergyShift(int ch)
{
    double shift = 0;
    if(ch == 646) shift = 0.022; // oldE - shift = correctedE, * * * BASED ON P3JDY * * *
    if(ch == 644) shift = 0.091;
    if(ch == 642) shift = 0.009;
    if(ch == 626) shift = -0.181;
    if(ch == 624) shift = -0.153;
    if(ch == 674) shift = -0.022;
    if(ch == 576) shift = 0.075;
    if(ch == 692) shift = 0.318;
    if(ch == 690) shift = 0.124;
    if(ch == 688) shift = 0.469;
    if(ch == 640) shift = 0.432;
    if(ch == 610) shift = -0.013;
    if(ch == 664) shift = 0.474;
    if(ch == 662) shift = 0.159;
    if(ch == 656) shift = 0.196;
    if(ch == 696) shift = 0.449;
    if(ch == 608) shift = 0.381;
    if(ch == 598) shift = 0.322;
    if(ch == 600) shift = -0.209;
    if(ch == 594) shift = 0.617;
    if(ch == 592) shift = 0.288;
     if(ch == 584) shift = 0.000; //Bad HV connection
     if(ch == 680) shift = 0.000; //High leakage current1
     if(ch == 676) shift = 0.000; //Bad signal connection
     if(ch == 616) shift = 0.000; //Bad signal connection
     if(ch == 614) shift = 0.448; //Gain oscillation
     if(ch == 628) shift = 0.587; //Sparking
     if(ch == 632) shift = 0.000; //High leakage current
     if(ch == 630) shift = 0.000; //Bad signal connection
    return shift;
}
// A FUNCTION TO CHECK FOR Enr VS Nat DETs //
int CheckEnriched(int ch)
{
    int cheasy = 1000;
    if(ch == 646) cheasy = 0; //Enr=1, Nat=0 * * * BASED ON P3JDY * * *
    if(ch == 644) cheasy = 1;
    if(ch == 642) cheasy = 1;
    if(ch == 626) cheasy = 1;
    if(ch == 624) cheasy = 1;
    if(ch == 674) cheasy = 1;
    if(ch == 576) cheasy = 1;
    if(ch == 692) cheasy = 1;
    if(ch == 690) cheasy = 1;
    if(ch == 688) cheasy = 1;
    if(ch == 640) cheasy = 1;
    if(ch == 610) cheasy = 1;
    if(ch == 664) cheasy = 0;
    if(ch == 662) cheasy = 1;
    if(ch == 656) cheasy = 1;
    if(ch == 696) cheasy = 1;
    if(ch == 608) cheasy = 0;
    if(ch == 598) cheasy = 0;
    if(ch == 600) cheasy = 0;
    if(ch == 594) cheasy = 0;
    if(ch == 592) cheasy = 0;
     if(ch == 584) cheasy = 0; //Bad HV connection
     if(ch == 680) cheasy = 1; //High leakage current1
     if(ch == 676) cheasy = 0; //Bad signal connection
     if(ch == 616) cheasy = 1; //Bad signal connection
     if(ch == 614) cheasy = 1; //Gain oscillation
     if(ch == 628) cheasy = 1; //Sparking
     if(ch == 632) cheasy = 1; //High leakage current
     if(ch == 630) cheasy = 1; //Bad signal connection
    return cheasy;
}

    
int main(int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Initialize Vars
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    // RUN LIST CONSIDERATIONS //
    ///int startrun = 4854;
    ///int endrun = 4855;
    int runNumber;
    
    // OtheR and FLOW CONTROL //
    char infile[200], infilename[200];
    int nentriest; // number of entries in the MJD tree
    int ch;
    int chIndex;
    int chanSize;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Initialize Plots and TFile
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    TFile * originalFile = TFile::Open("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/data/preSSTC_Plots_5126-5146.root"); // Local
    
    //TFile *newFile = new TFile("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/data/preSSTC_FittedPlots_5126-5146.root","NEW"); // http://www-numi.fnal.gov/offline_software/srt_public_context/WebDocs/Companion/root_crib/tfile.html
    //newFile->Close();
    
    int nGoodChannels = 21; // number of dets to be plotted (add 3 for "combined +2", "enr +1", and "nat +0" plots)
    
    TCanvas* cToE;
    cToE = new TCanvas();
    
    TCanvas* cLowE;
    cLowE = new TCanvas();
    
    TH1D *hToE[nGoodChannels+3];
    TH1D *hLowE[nGoodChannels+3];
    
    char hname_ToE[200];
    char hname_LowE[200];
    char specialName[200];

    for (int i=0; i<(nGoodChannels+3); i++)
    {
        if(i==nGoodChannels+0) sprintf(specialName,"Natural");
        if(i==nGoodChannels+1) sprintf(specialName,"Enriched");
        if(i==nGoodChannels+2) sprintf(specialName,"AllDets");
        
        if (i>=nGoodChannels)
        {
            sprintf(hname_ToE,"hist_ToE_%s",specialName);
            hToE[i] = (TH1D*)originalFile->Get(hname_ToE);
            
            sprintf(hname_LowE,"hist_LowE_%s",specialName);
            hLowE[i] = (TH1D*)originalFile->Get(hname_LowE);
        }
        if (i<nGoodChannels)
        {
            sprintf(hname_ToE,"hist_ToE_%d",IndexToChannel(i));
            hToE[i] = (TH1D*)originalFile->Get(hname_ToE);
            
            sprintf(hname_LowE,"hist_LowE_%d",IndexToChannel(i));
            hLowE[i] = (TH1D*)originalFile->Get(hname_LowE);
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Plot Plots, Fit and Write
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    double peakpos;
    double peakheight;
    TF1 * fit_ToE = new TF1("ToEFit", "[0]*Gaus(x,[1],[2]) + [3] + [4]*x ",0,1000);  // TMath::Gaus(x,[mean],[sigma]) , "0,1000" are
    TF1 * fit_LowE = new TF1("LowEFit", "[0]*Gaus(x,[1],[2]) + [3] + [4]*x ",0,1000);  // TMath::Gaus(x,[mean],[sigma]) , "0,1000" are


/*
    for (int i=0; i<(nGoodChannels+3); i++)
    {
        hToE[i]->Draw();
        peakpos = hToE[i]->GetBinCenter(hToE[i]->GetMaximumBin()); // guess peak by max bin in window
        peakheight = hToE[i]->GetBinContent(hToE[i]->FindBin(peakpos)); // get the height of the max bin
        fit_ToE->SetLineColor(kRed);
        fit_ToE->SetLineWidth(1);
        fit_ToE->SetLineStyle(2);
        fit_ToE->FixParameter(0,peakheight*0.92); // gaus amplitude
        fit_ToE->SetParameter(1,peakpos); // gaus mean
        fit_ToE->SetParameter(2,peakpos/20.0); // gaus sigma
        fit_ToE->SetParameter(3,0); // flat BG
        fit_ToE->SetParameter(4,0); // linear BG
        hToE[i]->Fit("ToEFit","qr+"); // for fit options https://root.cern.ch/doc/master/classTH1.html#a63eb028df86bc86c8e20c989eb23fb2a
        
        hLowE[i]->Draw();
        peakpos = 10.3; // guess peak by max bin in window
        peakheight = 0.0;; // get the height of the max bin
        LowEFit->SetRange(0.53*peakpos,1.4*peakpos);
        LowEFit->SetLineColor(kRed);
        LowEFit->SetLineWidth(1);
        LowEFit->SetLineStyle(2);
        LowEFit->SetParameter(0,0.0); // gaus amplitude
        LowEFit->SetParameter(1,peakpos); // gaus mean
        LowEFit->SetParameter(2,peakpos/20.0); // gaus sigma
        LowEFit->SetParameter(3,0); // flat BG
        LowEFit->SetParameter(4,0); // linear BG
        hLowE[i]->Fit("LowEFit","qr+"); // for fit options https://root.cern.ch/doc/master/classTH1.html#a63eb028df86bc86c8e20c989eb23fb2a
    }
*/
    
    /*
     cLowE->cd();
     hLowE[nGoodChannels+2]->Draw();
     cLowE->Update();
     */
    
    cLowE->cd();
    cLowE->SetLogy();
    hLowE[nGoodChannels+2]->Draw();
    peakpos = 10.3; // guess peak by max bin in window
    peakheight = 0.0;; // get the height of the max bin
    fit_LowE->SetRange(0.53*peakpos,1.4*peakpos);
    fit_LowE->SetLineColor(kRed);
    fit_LowE->SetLineWidth(1);
    fit_LowE->SetLineStyle(2);
    fit_LowE->SetParameter(0,0.0); // gaus amplitude
    fit_LowE->SetParameter(1,peakpos); // gaus mean
    fit_LowE->SetParameter(2,peakpos/20.0); // gaus sigma
    fit_LowE->SetParameter(3,0); // flat BG
    fit_LowE->SetParameter(4,0); // linear BG
    hLowE[nGoodChannels+2]->Fit("LowEFit","qr+"); // for fit options https://root.cern.ch/doc/master/classTH1.html#a63eb028df86bc86c8e20c989eb23fb2a
    hLowE[nGoodChannels+2]->GetXaxis()->SetRangeUser(0,20);
    cLowE->Update();

    double meanLowE = fit_LowE->GetParameter(1);
    double sigmaLowE = fit_LowE->GetParameter(2);
    
    cout<<"mean = " << meanLowE << " sigma = " << sigmaLowE << endl;

    originalFile->Close();
    
    
    cout << "F I N I S H E D" << endl;
    App->Run();
} // End of main()
