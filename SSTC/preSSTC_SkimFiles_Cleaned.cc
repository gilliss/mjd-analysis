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
/*
// FUNCTION FOR SPECIAL_CASE PLOTS //
char SpecialName(int i)
{
    int nGoodChannels = 21; // number of dets to be plotted (add 3 for "combined", "enr", and "nat" plots)
    char specialName[200];
    if(i==nGoodChannels+0){sprintf(specialName,"Natural");
    if(i==nGoodChannels+1){sprintf(specialName,"Enriched");
    if(i==nGoodChannels+2){sprintf(specialName,"AllDets");
    return specialName;
}
*/
    
int main(int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initialize Vars
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    // RUN LIST CONSIDERATIONS //

    // MJD DATA & MJD TREE BRANCHES //
    double startTime = 0;
    double stopTime = 0;
    vector<int>* channel = 0;
    vector<double>* trapECal = 0;
    vector<double>* trapETailMin = 0;
    vector<double>* toe = 0;
    vector<double>* tloc_s = 0; // time since start of run in seconds
    int mH;
    vector<bool>* isEnr = 0;
    vector<bool>* isGood = 0;
    int run = 0;
    

    // OtheR and FLOW CONTROL //
    char infile[200], infilename[200];
    int nentriesOriginal; // number of entries in the MJD tree
    int ch;
    int chIndex;
    int chanSize;
    
    double shiftE = 0;
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initialize Plots and TFile
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    //TFile *f = new TFile("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/data/preSSTC_Plots_5610-5680.root","NEW"); // Local
    TFile *f = new TFile("/global/u2/g/gilliss/dev/Analysis/Cryo1/SSTC/preSSTC_Plots_P3JDYGolden_SkimFiles_Cleaned1st.root","NEW"); // PDSF
        //http://www-numi.fnal.gov/offline_software/srt_public_context/WebDocs/Companion/root_crib/tfile.html
    
    int nGoodChannels = 21; // number of dets to be plotted (add 3 for "combined +2", "enr +1", and "nat +0" plots)
    
    int nbins_E = 1000;
    int nbins_ToE = 50;
    
    TCanvas* cToEvE;
    cToEvE = new TCanvas();
    
    TCanvas* cToE;
    cToE = new TCanvas();
    
    TCanvas* cLowE;
    cLowE = new TCanvas();
    
    TH2D *hToEvE[nGoodChannels+3];
    TH1D *hToE[nGoodChannels+3];
    TH1D *hLowE[nGoodChannels+3];
    
    char hname_ToEvE[200];
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
            hToE[i] = new TH1D(hname_ToE, hname_ToE, nbins_ToE, 0, 5);
            
            sprintf(hname_LowE,"hist_LowE_%s",specialName);
            hLowE[i] = new TH1D(hname_LowE, hname_LowE, 100, 0, 50);
            
            sprintf(hname_ToEvE,"hist_ToEvE_%s",specialName);
            hToEvE[i] = new TH2D(hname_ToEvE, hname_ToEvE, nbins_E, 0, 3000, nbins_ToE, 0, 5);//(name,title,nbinsx,xlower,xupper,nbinsy,ylower,yupper)
        }
        if (i<nGoodChannels)
        {
            sprintf(hname_ToEvE,"hist_ToEvE_%d",IndexToChannel(i));
            hToEvE[i] = new TH2D(hname_ToEvE, hname_ToEvE, nbins_E, 0, 3000, nbins_ToE, 0, 5);//(name,title,nbinsx,xlower,xupper,nbinsy,ylower,yupper)
            
            sprintf(hname_ToE,"hist_ToE_%d",IndexToChannel(i));
            hToE[i] = new TH1D(hname_ToE, hname_ToE, nbins_ToE, 0, 5);
            
            sprintf(hname_LowE,"hist_LowE_%d",IndexToChannel(i));
            hLowE[i] = new TH1D(hname_LowE, hname_LowE, 100, 0, 50);
        }
    }
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Fill Plots
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    TChain* originalTree = new TChain("skimTree", "skimTree");
    //originalTree->Add("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/data/skim_5610_5680.root"); // Local
    originalTree->Add("/global/homes/j/jasondet/myProject/Skim/skim_*.root"); // PDSF
    originalTree->SetBranchAddress("startTime",&startTime); // unix time at start of run
    originalTree->SetBranchAddress("stopTime",&stopTime); // unix time at end of run
    originalTree->SetBranchAddress("channel",&channel); // set addresses for branches in the MJD Tree
    originalTree->SetBranchAddress("trapECal",&trapECal);
    originalTree->SetBranchAddress("trapETailMin",&trapETailMin);
    originalTree->SetBranchAddress("toe",&toe);
    originalTree->SetBranchAddress("tloc_s",&tloc_s);
    originalTree->SetBranchAddress("mH",&mH);
    originalTree->SetBranchAddress("isEnr",&isEnr);
    originalTree->SetBranchAddress("isGood",&isGood);
    originalTree->SetBranchAddress("run",&run);

    nentriesOriginal=originalTree->GetEntries(); // Get number of entries in the MJD Tree
    
    // LOOP OVER ENTRIES IN CURRENT MJD TREE //
    for(int k=0;k<nentriesOriginal;k++)
    {
        originalTree->GetEntry(k); // Read all branches of entry and return total number of bytes read.
        chanSize = channel->size(); // Get the multiplicity of channels within the entry
        if(channel->at(0)!=0) // check for bogus 0-value
        {
            for(int j=0; j<chanSize; j++) // loop through channels in the entry
            {
                ch = channel->at(j);
                if(CheckChannel(ch)!=1000) //if(ch%2==0) // use just high gain channels
                {
                    chIndex = CheckChannel(ch);
                    if(trapETailMin->at(j)<0.0 && mH==1)
                    {
                        shiftE = trapECal->at(j) - EnergyShift(ch);
                        hToEvE[chIndex]->Fill(shiftE,toe->at(j));
                        hToEvE[nGoodChannels+2]->Fill(shiftE,toe->at(j));
                        if(CheckEnriched(ch)==0) hToEvE[nGoodChannels+0]->Fill(shiftE,toe->at(j));
                        if(CheckEnriched(ch)==1) hToEvE[nGoodChannels+1]->Fill(shiftE,toe->at(j));
                        
                        if (shiftE<=50.0 && toe->at(j) < 1.59 && toe->at(j) > 1.0)
                        {
                            hLowE[chIndex]->Fill(shiftE);
                            hLowE[nGoodChannels+2]->Fill(shiftE);
                            if(CheckEnriched(ch)==0) hLowE[nGoodChannels+0]->Fill(shiftE);
                            if(CheckEnriched(ch)==1) hLowE[nGoodChannels+1]->Fill(shiftE);
                        }
                    }
                }
            } // End of loop over channels j within the current entry k
        }
    } // End of loop over entries k in run i
    delete originalTree;
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Plot Plots
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    char title[200];
    for (int i=0; i<(nGoodChannels+3); i++)
    {
        if(i==nGoodChannels+0) sprintf(specialName,"Natural");
        if(i==nGoodChannels+1) sprintf(specialName,"Enriched");
        if(i==nGoodChannels+2) sprintf(specialName,"AllDets");
        
        if(i>=nGoodChannels) sprintf(hname_ToE,"hist_ToE_%s",specialName);
        else sprintf(hname_ToE,"hist_ToE_%d",IndexToChannel(i));
        
        if (i>=nGoodChannels)
        {
            cToEvE->cd();
            sprintf(title,"ToE vs E, P3JDY Golden BG Runs, %s",specialName);
            hToEvE[i]->SetTitle(title);
            hToEvE[i]->GetYaxis()->SetTitle("ToE");
            hToEvE[i]->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
            hToEvE[i]->Draw(); //Draw("colz"); if want binning/averaging ("blocks rather than points") in the 2d hist
            if(i==nGoodChannels+2) cToEvE->Update();
            
            // PROJECT ToEvE onto ToE // reference:PlotFitSaveHists_KrisRootFile.cc
            hToE[i] = hToEvE[i]->ProjectionY(hname_ToE,0,nbins_ToE);
            cToE->cd();
            cToE->SetLogy();
            sprintf(title,"ToE, P3JDY Golden BG Runs, %s",specialName);
            hToE[i]->SetTitle(title);
            hToE[i]->GetYaxis()->SetTitle("Counts");
            hToE[i]->GetYaxis()->SetTitle("ToE");
            hToE[i]->Draw();
            if(i==nGoodChannels+2) cToE->Update();
            
            cLowE->cd();
            cLowE->SetLogy();
            sprintf(title,"Energy, P3JDY Golden BG Runs, %s",specialName);
            hLowE[i]->SetTitle(title);
            hLowE[i]->GetYaxis()->SetTitle("Counts");
            hLowE[i]->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
            hLowE[i]->Draw();
            if(i==nGoodChannels+2) cLowE->Update();
        }
        if (i<nGoodChannels)
        {
            cToEvE->cd();
            sprintf(title,"ToE vs E, P3JDY Golden BG Runs, %d",IndexToChannel(i));
            hToEvE[i]->SetTitle(title);
            hToEvE[i]->GetYaxis()->SetTitle("ToE");
            hToEvE[i]->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
            hToEvE[i]->Draw(); //Draw("colz"); if want binning/averaging ("blocks rather than points") in the 2d hist
            //cToEvE->Update();
            
            // PROJECT ToEvE onto ToE // reference:PlotFitSaveHists_KrisRootFile.cc
            hToE[i] = hToEvE[i]->ProjectionY(hname_ToE,0,nbins_ToE);
            cToE->cd();
            cToE->SetLogy();
            sprintf(title,"ToE, P3JDY Golden BG Runs, %d",IndexToChannel(i));
            hToE[i]->SetTitle(title);
            hToE[i]->GetYaxis()->SetTitle("Counts");
            hToE[i]->GetYaxis()->SetTitle("ToE");
            hToE[i]->Draw();
            //cToE->Update();
            
            cLowE->cd();
            cLowE->SetLogy();
            sprintf(title,"Energy, P3JDY Golden BG Runs, %d",IndexToChannel(i));
            hLowE[i]->SetTitle(title);
            hLowE[i]->GetYaxis()->SetTitle("Counts");
            hLowE[i]->GetXaxis()->SetTitle("trapECal - 0PtCalShift (keV)");
            hLowE[i]->Draw();
            //cLowE->Update();
        }
        hToEvE[i]->Write();
        hToE[i]->Write();
        hLowE[i]->Write();
    }
    f->Close();
    
    
    cout << "F I N I S H E D" << endl;
    App->Run();
} // End of main()
