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

    // PROMPT-DELAYED DATA & SKIM TREE BRANCHES //
    vector<double>* globalTimestamp; // vector containing globalTimestampValues of each waveform of an entry
    vector<int>* promptTag; // a vector containing 0 for non-K-shell WFs, 1 for K-shell WFs
    vector<int>* delayedTag; // 0 for out of window, 1 for in window
    vector<double>* timeSincePrompt; // Time between a prompt (K-shell) and delayed event (in the SSTC window)
    vector<double>* shiftedE;
    vector<int>* delayedOutsideTag; // 0 for out of window, 1 for in window
    vector<int>* chan; // 0 for out of window, 1 for in window
    vector<double>* ToE;
    vector<double>* trapETM;
    vector<int>* mult;
    
    double E;
    int mH;
    double trapETailMin;
    double toe;
    
    // OtheR and FLOW CONTROL //
    char infile[200], infilename[200];
    int nentriesOriginal; // number of entries in the MJD tree
    int ch;
    int chIndex;
    int chanSize;
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initialize Plots and TFile
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    //TFile *f = new TFile("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/data/preSSTC_Plots_5610-5680.root","NEW"); // Local
        //http://www-numi.fnal.gov/offline_software/srt_public_context/WebDocs/Companion/root_crib/tfile.html
    
    int nGoodChannels = 21; // number of dets to be plotted (add 3 for "combined +2", "enr +1", and "nat +0" plots)
    
    int nbins_LowE = 100;
    double max_LowE = 50;
    int nbins_E = 3000;
    double max_E = 3000;
    
    TCanvas* cLowE_Cleaned = new TCanvas();
    TCanvas* cLowE_Cleaned_Nat = new TCanvas();
    TCanvas* cLowE_Cleaned_Enr = new TCanvas();
    TH1D *hLowE_Cleaned = new TH1D("hLowE_Cleaned", "hLowE_Cleaned", nbins_LowE, 0, max_LowE);
    TH1D *hLowE_Cleaned_Nat = new TH1D("hLowE_Cleaned_Nat", "hLowE_Cleaned_Nat", nbins_LowE, 0, max_LowE);
    TH1D *hLowE_Cleaned_Enr = new TH1D("hLowE_Cleaned_Enr", "hLowE_Cleaned_Enr", nbins_LowE, 0, max_LowE);
    TH1D *hLowE_Prompt_Cleaned  = new TH1D("hLowE_Prompt_Cleaned", "hLowE_Prompt_Cleaned", nbins_LowE, 0, max_LowE);
    TH1D *hLowE_Prompt_Cleaned_Nat  = new TH1D("hLowE_Prompt_Cleaned_Nat", "hLowE_Prompt_Cleaned_Nat", nbins_LowE, 0, max_LowE);
    TH1D *hLowE_Prompt_Cleaned_Enr  = new TH1D("hLowE_Prompt_Cleaned_Enr", "hLowE_Prompt_Cleaned_Enr", nbins_LowE, 0, max_LowE);
    
    TCanvas* cDelayed_Cleaned = new TCanvas();
    TCanvas* cDelayed_Cleaned_Nat = new TCanvas();
    TCanvas* cDelayed_Cleaned_Enr = new TCanvas();
    TH1D *hDelayed_Cleaned = new TH1D("hDelayed_Cleaned", "hDelayed_Cleaned", 300, 0, max_E);
    TH1D *hDelayed_Cleaned_Nat = new TH1D("hDelayed_Cleaned_Nat", "hDelayed_Cleaned_Nat", 300, 0, max_E);
    TH1D *hDelayed_Cleaned_Enr = new TH1D("hDelayed_Cleaned_Enr", "hDelayed_Cleaned_Enr", 300, 0, max_E);

    TCanvas* cFullE_Cleaned = new TCanvas();
    TCanvas* cFullE_Cleaned_Nat = new TCanvas();
    TCanvas* cFullE_Cleaned_Enr = new TCanvas();
    TH1D *hFullE_Cleaned = new TH1D("hFullE_Cleaned", "hFullE_Cleaned", nbins_E, 0, max_E);
    TH1D *hFullE_Cleaned_Nat = new TH1D("hFullE_Cleaned_Nat", "hFullE_Cleaned_Nat", nbins_E, 0, max_E);
    TH1D *hFullE_Cleaned_Enr = new TH1D("hFullE_Cleaned_Enr", "hFullE_Cleaned_Enr", nbins_E, 0, max_E);
    TH1D *hFullE_SansDelayed_Cleaned = new TH1D("hFullE_SansDelayed_Cleaned", "hFullE_SansDelayed_Cleaned", nbins_E, 0, max_E);
    TH1D *hFullE_SansDelayed_Cleaned_Nat = new TH1D("hFullE_SansDelayed_Cleaned_Nat", "hFullE_SansDelayed_Cleaned_Nat", nbins_E, 0, max_E);
    TH1D *hFullE_SansDelayed_Cleaned_Enr = new TH1D("hFullE_SansDelayed_Cleaned_Enr", "hFullE_SansDelayed_Cleaned_Enr", nbins_E, 0, max_E);


    
    /*
     TCanvas* cFullE_SansDelayed_Cleaned = new TCanvas();

    TCanvas* cDelayedOutside = new TCanvas();
    TCanvas* cDelayedOutsideTimeSince = new TCanvas();
    TCanvas* cFullE = new TCanvas();
    TCanvas* cFullE_SansDelayed = new TCanvas();
*/



    //TH1D *hDelayedOutside = new TH1D("hDelayedOutside", "hDelayedOutside", nbins_E, 0, max_E);
    //TH1D *hDelayedOutsideTimeSince = new TH1D("hDelayedOutsideTimeSince", "hDelayedOutsideTimeSince", 800, 0, 8000);
    
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Fill Plots
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    TChain* originalTree = new TChain("sstcTree", "sstcTree");
    originalTree->Add("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/SSTC_SkimFiles_P3JDYGolden_CleanedFixedBounds3Sigma.root"); // Local
    originalTree->SetBranchAddress("globalTimestamp",&globalTimestamp); // Declare branch for storing globalTimestamps
    originalTree->SetBranchAddress("promptTag",&promptTag); // Declare branch for storing K-shell tag
    originalTree->SetBranchAddress("delayedTag",&delayedTag); // Declare branch for storing SSTC tag
    originalTree->SetBranchAddress("timeSincePrompt",&timeSincePrompt); // Declare branch for storing timeSincePrompt
    originalTree->SetBranchAddress("shiftedE",&shiftedE);
    originalTree->SetBranchAddress("delayedOutsideTag",&delayedOutsideTag); // Declare branch for storing SSTC tag
    originalTree->SetBranchAddress("chan",&chan);
    originalTree->SetBranchAddress("trapETM",&trapETM);
    originalTree->SetBranchAddress("ToE",&ToE);
    originalTree->SetBranchAddress("mult",&mult);


    nentriesOriginal=originalTree->GetEntries(); // Get number of entries in the MJD Tree
    cout << "nentries " << nentriesOriginal << endl;
    
    // LOOP OVER ENTRIES IN CURRENT MJD TREE //
    for(int k=0;k<nentriesOriginal;k++)
    {
        originalTree->GetEntry(k); // Read all branches of entry and return total number of bytes read.
        chanSize = chan->size(); // Get the multiplicity of channels within the entry
        for(int j=0; j<chanSize; j++) // loop through channels in the entry
        {
            ch = chan->at(j);
            chIndex = CheckChannel(ch);
            E = shiftedE->at(j);
            mH = mult->at(j);
            trapETailMin = trapETM->at(j);
            toe = ToE->at(j);
            
            
            if(chIndex!=1000) //if(ch%2==0) // use just high gain channels
            {
                if(E<50.0 && mH==1 && trapETailMin<0.0 && toe < 1.59 && toe > 1.0)
                {
                    hLowE_Cleaned->Fill(E);
                    if (CheckEnriched(ch)==0) hLowE_Cleaned_Nat->Fill(E);
                    if (CheckEnriched(ch)==1) hLowE_Cleaned_Enr->Fill(E);
                    
                    if(promptTag->at(j)==1)
                    {
                        hLowE_Prompt_Cleaned->Fill(E);
                        if (CheckEnriched(ch)==0) hLowE_Prompt_Cleaned_Nat->Fill(E);
                        if (CheckEnriched(ch)==1) hLowE_Prompt_Cleaned_Enr->Fill(E);
                    }
                }
                
                if (delayedTag->at(j)==1 && mH==1 && trapETailMin<0.0 && toe < 1.59 && toe > 1.0)
                {
                    hDelayed_Cleaned->Fill(E);
                    if (CheckEnriched(ch)==0) hDelayed_Cleaned_Nat->Fill(E);
                    if (CheckEnriched(ch)==1) hDelayed_Cleaned_Enr->Fill(E);
                }
                
                if (mH==1 && trapETailMin<0.0 && toe < 1.59 && toe > 1.0)
                {
                    hFullE_Cleaned->Fill(E);
                    if (CheckEnriched(ch)==0) hFullE_Cleaned_Nat->Fill(E);
                    if (CheckEnriched(ch)==1) hFullE_Cleaned_Enr->Fill(E);
                    
                    if(delayedTag->at(j)!=1)
                    {
                        hFullE_SansDelayed_Cleaned->Fill(E);
                        if (CheckEnriched(ch)==0) hFullE_SansDelayed_Cleaned_Nat->Fill(E);
                        if (CheckEnriched(ch)==1) hFullE_SansDelayed_Cleaned_Enr->Fill(E);
                    }
                }
                
                
                
                
                //if (delayedOutsideTag->at(j)==1)
                //{
                //    hDelayedOutside->Fill(E);
                //    if (E < 513.0 && E > 509.0) hDelayedOutsideTimeSince->Fill(E);
                //}
            }
        } // End of loop over channels j in entry k
    } // End of loop over entries k in run i
    delete originalTree;
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Plot Plots
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    char title[200];
    int color[5] = {2,4,1,30,6};//{46,9,1,30}; // {Nat,Enr,All,Aux,Tagged}
    double natM = 4.0914;
    double enrM = 10.6960;
    double totM = 14.7874;
    double exposure = 45.0508;
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //--------Low E & Prompt ------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    cLowE_Cleaned->cd();
    cLowE_Cleaned->SetLogy();
    sprintf(title,"Energy, 'Good' Detectors, Cleaned");
    hLowE_Cleaned->SetTitle(title);
    hLowE_Cleaned->GetYaxis()->SetTitle("cnt/kg/dy");
    hLowE_Cleaned->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
    hLowE_Cleaned->SetLineWidth(2);
    hLowE_Cleaned->SetLineColor(color[2]);
    hLowE_Cleaned->Scale(1/(totM*exposure));
    hLowE_Cleaned->Draw();
    cLowE_Cleaned->Update();
    
    hLowE_Prompt_Cleaned->SetLineWidth(2);
    hLowE_Prompt_Cleaned->SetLineColor(color[4]);
    hLowE_Prompt_Cleaned->Scale(1/(totM*exposure));
    hLowE_Prompt_Cleaned->Draw("same");
    //cLowE_Prompt_Cleaned->Update();
    cLowE_Cleaned->Update();
    cLowE_Cleaned->Print("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/Plots/hLowE_Cleaned_All.png");
    cout << "Full prompt tag" << hLowE_Prompt_Cleaned->Integral(0,100) << endl;

    //-----------------------------//
    cLowE_Cleaned_Nat->cd();
    //cLowE_Cleaned_Nat->SetLogy();
    sprintf(title,"Energy, 'Good' Natural Detectors, Cleaned");
    hLowE_Cleaned_Nat->SetTitle(title);
    hLowE_Cleaned_Nat->GetYaxis()->SetTitle("cnt/kg/dy");
    hLowE_Cleaned_Nat->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
    hLowE_Cleaned_Nat->SetLineWidth(2);
    hLowE_Cleaned_Nat->SetLineColor(color[0]);
    hLowE_Cleaned_Nat->Scale(1/(exposure*natM));
    hLowE_Cleaned_Nat->Draw();
    cLowE_Cleaned_Nat->Update();
    
    hLowE_Prompt_Cleaned_Nat->SetLineWidth(2);
    hLowE_Prompt_Cleaned_Nat->SetLineColor(color[4]);
    hLowE_Prompt_Cleaned_Nat->Scale(1/(exposure*natM));
    hLowE_Prompt_Cleaned_Nat->Draw("same");
    cLowE_Cleaned_Nat->Update();
    cLowE_Cleaned_Nat->Print("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/Plots/hLowE_Cleaned_Nat.png");
    cout << "Nat prompt tag" << hLowE_Prompt_Cleaned_Nat->Integral(0,100) << endl;

    //-----------------------------//
    cLowE_Cleaned_Enr->cd();
    cLowE_Cleaned_Enr->SetLogy();
    sprintf(title,"Energy, 'Good' Enriched Detectors, Cleaned");
    hLowE_Cleaned_Enr->SetTitle(title);
    hLowE_Cleaned_Enr->GetYaxis()->SetTitle("cnt/kg/dy");
    hLowE_Cleaned_Enr->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
    hLowE_Cleaned_Enr->SetLineWidth(2);
    hLowE_Cleaned_Enr->SetLineColor(color[1]);
    hLowE_Cleaned_Enr->Scale(1/(exposure*enrM));
    hLowE_Cleaned_Enr->Draw();
    cLowE_Cleaned_Enr->Update();
    
    hLowE_Prompt_Cleaned_Enr->SetLineWidth(2);
    hLowE_Prompt_Cleaned_Enr->SetLineColor(color[4]);
    hLowE_Prompt_Cleaned_Enr->Scale(1/(exposure*enrM));
    hLowE_Prompt_Cleaned_Enr->Draw("same");
    cLowE_Cleaned_Enr->Update();
    cLowE_Cleaned_Enr->Print("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/Plots/hLowE_Cleaned_Enr.png");
    cout << "Enr prompt tag" << hLowE_Prompt_Cleaned_Enr->Integral(0,100) << endl;

    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //--------- Delayed  ----------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    cDelayed_Cleaned->cd();
    cDelayed_Cleaned->SetLogy();
    sprintf(title,"E of Delayed Events, 'Good' Detectors, Cleaned");
    hDelayed_Cleaned->SetTitle(title);
    hDelayed_Cleaned->GetYaxis()->SetTitle("cnt/kg/dy");
    hDelayed_Cleaned->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
    hDelayed_Cleaned->SetLineWidth(2);
    hDelayed_Cleaned->SetLineColor(color[2]);
    hDelayed_Cleaned->Scale(1/(totM*exposure));
    hDelayed_Cleaned->Draw();
    cDelayed_Cleaned->Update();

    cDelayed_Cleaned->Print("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/Plots/hDelayed_Cleaned_All.png");
    //-----------------------------//
    cDelayed_Cleaned_Nat->cd();
    cDelayed_Cleaned_Nat->SetLogy();
    sprintf(title,"E of Delayed Events, 'Good' Natural Detectors, Cleaned");
    hDelayed_Cleaned_Nat->SetTitle(title);
    hDelayed_Cleaned_Nat->GetYaxis()->SetTitle("cnt/kg/dy");
    hDelayed_Cleaned_Nat->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
    hDelayed_Cleaned_Nat->SetLineWidth(2);
    hDelayed_Cleaned_Nat->SetLineColor(color[0]);
    hDelayed_Cleaned_Nat->Scale(1/(natM*exposure));
    hDelayed_Cleaned_Nat->Draw();
    cDelayed_Cleaned_Nat->Update();
    
    cDelayed_Cleaned_Nat->Print("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/Plots/hDelayed_Cleaned_Nat.png");
    //-----------------------------//
    cDelayed_Cleaned_Enr->cd();
    cDelayed_Cleaned_Enr->SetLogy();
    sprintf(title,"E of Delayed Events, 'Good' Enriched Detectors, Cleaned");
    hDelayed_Cleaned_Enr->SetTitle(title);
    hDelayed_Cleaned_Enr->GetYaxis()->SetTitle("cnt/kg/dy");
    hDelayed_Cleaned_Enr->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
    hDelayed_Cleaned_Enr->SetLineWidth(2);
    hDelayed_Cleaned_Enr->SetLineColor(color[1]);
    hDelayed_Cleaned_Enr->Scale(1/(enrM*exposure));
    hDelayed_Cleaned_Enr->Draw();
    cDelayed_Cleaned_Enr->Update();
    
    cDelayed_Cleaned_Enr->Print("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/Plots/hDelayed_Cleaned_Enr.png");
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //------- Full Spectrum -------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    //-----------------------------//
    cFullE_Cleaned->cd();
    cFullE_Cleaned->SetLogy();
    sprintf(title,"Energy, 'Good' Detectors, Cleaned");
    hFullE_Cleaned->SetTitle(title);
    hFullE_Cleaned->GetYaxis()->SetTitle("cnt/kg/dy");
    hFullE_Cleaned->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
    hFullE_Cleaned->SetLineWidth(2);
    hFullE_Cleaned->SetLineColor(color[2]);
    hFullE_Cleaned->Scale(1/(totM*exposure));
    hFullE_Cleaned->Draw();
    cFullE_Cleaned->Update();

    hFullE_SansDelayed_Cleaned->SetLineWidth(2);
    hFullE_SansDelayed_Cleaned->SetLineColor(color[3]);
    hFullE_SansDelayed_Cleaned->Scale(1/(totM*exposure));
    hFullE_SansDelayed_Cleaned->Draw("same");
    cFullE_Cleaned->Update();
    cFullE_Cleaned->Print("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/Plots/hFullE_Cleaned_All.png");
    //-----------------------------//
    
    cFullE_Cleaned_Nat->cd();
    cFullE_Cleaned_Nat->SetLogy();
    sprintf(title,"Energy, 'Good' Natural Detectors, Cleaned");
    hFullE_Cleaned_Nat->SetTitle(title);
    hFullE_Cleaned_Nat->GetYaxis()->SetTitle("cnt/kg/dy");
    hFullE_Cleaned_Nat->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
    hFullE_Cleaned_Nat->SetLineWidth(2);
    hFullE_Cleaned_Nat->SetLineColor(color[0]);
    hFullE_Cleaned_Nat->Scale(1/(natM*exposure));
    hFullE_Cleaned_Nat->Draw();
    cFullE_Cleaned_Nat->Update();
    
    hFullE_SansDelayed_Cleaned_Nat->SetLineWidth(2);
    hFullE_SansDelayed_Cleaned_Nat->SetLineColor(color[3]);
    hFullE_SansDelayed_Cleaned_Nat->Scale(1/(natM*exposure));
    hFullE_SansDelayed_Cleaned_Nat->Draw("same");
    cFullE_Cleaned_Nat->Update();
    cFullE_Cleaned_Nat->Print("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/Plots/hFullE_Cleaned_Nat.png");
    //-----------------------------//
    cFullE_Cleaned_Enr->cd();
    cFullE_Cleaned_Enr->SetLogy();
    sprintf(title,"Energy, 'Good' Enriched Detectors, Cleaned");
    hFullE_Cleaned_Enr->SetTitle(title);
    hFullE_Cleaned_Enr->GetYaxis()->SetTitle("cnt/kg/dy");
    hFullE_Cleaned_Enr->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
    hFullE_Cleaned_Enr->SetLineWidth(2);
    hFullE_Cleaned_Enr->SetLineColor(color[1]);
    hFullE_Cleaned_Enr->Scale(1/(enrM*exposure));
    hFullE_Cleaned_Enr->Draw();
    cFullE_Cleaned_Enr->Update();
    
    hFullE_SansDelayed_Cleaned_Enr->SetLineWidth(2);
    hFullE_SansDelayed_Cleaned_Enr->SetLineColor(color[3]);
    hFullE_SansDelayed_Cleaned_Enr->Scale(1/(enrM*exposure));
    hFullE_SansDelayed_Cleaned_Enr->Draw("same");
    cFullE_Cleaned_Enr->Update();
    cFullE_Cleaned_Enr->Print("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/Plots/hFullE_Cleaned_Enr.png");
    //-----------------------------//
    
    
    
    
    
    
    
    
    
    cout << "F I N I S H E D" << endl;
    App->Run();
} // End of main()




/*
 cDelayed->cd();
 cDelayed->SetLogy();
 sprintf(title,"E of Delayed Events");
 hDelayed->SetTitle(title);
 hDelayed->GetYaxis()->SetTitle("Counts");
 hDelayed->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
 hDelayed->Draw();
 cDelayed->Update();
 
 cDelayedOutside->cd();
 cDelayedOutside->SetLogy();
 sprintf(title,"E of Delayed Outside Det of Prompt");
 hDelayedOutside->SetTitle(title);
 hDelayedOutside->GetYaxis()->SetTitle("Counts");
 hDelayedOutside->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
 hDelayedOutside->Draw();
 cDelayedOutside->Update();
 
 cDelayedOutsideTimeSince->cd();
 cDelayedOutsideTimeSince->SetLogy();
 sprintf(title,"Time of Delayed Escaped 511 Since Prompt");
 hDelayedOutsideTimeSince->SetTitle(title);
 hDelayedOutsideTimeSince->GetYaxis()->SetTitle("Counts");
 hDelayedOutsideTimeSince->GetXaxis()->SetTitle("'time since' (s)");
 hDelayedOutsideTimeSince->Draw();
 cDelayedOutsideTimeSince->Update();
 
 cFullE->cd();
 cFullE->SetLogy();
 sprintf(title,"Energy, All 'Good' Detectors");
 hFullE->SetTitle(title);
 hFullE->GetYaxis()->SetTitle("Counts");
 hFullE->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
 hFullE->Draw();
 cFullE->Update();
 
 cFullE_SansDelayed->cd();
 cFullE_SansDelayed->SetLogy();
 sprintf(title,"Energy, All 'Good' Detectors, No Delayed Events");
 hFullE_SansDelayed->SetTitle(title);
 hFullE_SansDelayed->GetYaxis()->SetTitle("Counts");
 hFullE_SansDelayed->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
 hFullE_SansDelayed->Draw();
 cFullE_SansDelayed->Update();
 */
