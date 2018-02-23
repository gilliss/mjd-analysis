/*************************************************
-sub fsteps.edep for vector to get "chansize" from

 *************************************************/
#include "TProof.h"
#include "TROOT.h"
#include <getopt.h>

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
//#include <vector>

using namespace std;

int main(int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);
    
    //gROOT->ProcessLine(".x $MGDODIR/Root/LoadMGDOClasses.C");


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initialize Vars
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    // MJD DATA & MJD TREE BRANCHES //

    ULong64_t fEventID = 0;
    vector<Double_t>* fT = 0;
    vector<Int_t>* fParticleID = 0;
    vector<Double_t>* fEdep = 0;
    vector<Int_t>* fSensVolID = 0;
    Double_t fTotalEnergy = 0;

    // OTHER and FLOW CONTROL //
    char infile[200], infilename[200];
    int nentriest; // number of entries in the MJD tree
    Int_t ch;
    int eventSize;
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initialize Plots
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    TH1D *hEnergy;
    int nbins = 100;
    int minbin = 0;
    int maxbin = 50;
    char hname[200];
    sprintf(hname,"hist_energy");
    hEnergy= new TH1D(hname, hname, nbins, minbin, maxbin);
    
    TCanvas* cEnergy;
    cEnergy = new TCanvas();
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Loop over data
//////////////////////////////////////////////////////////////////////////////////////////////////////////

        // LOAD & PREP THE ORIGINAL MJD TREE //
        cout << "----------RUNNING-------------" << endl;
    
        sprintf(infilename,"/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/68GeSimulation/68GeInCrystal_ElectronCaptures_10000Events.root");

        TChain *t = new TChain("fTree"); // Having this inside the for loop prevents files from being chained--the chain just restarts for each run
        t->AddFile(infilename); // Associate ROOT file with the current TChain
/*
        t->SetBranchAddress("fEventID",&fEventID);
        t->SetBranchAddress("fT",&fT);
        t->SetBranchAddress("fParticleID",&fParticleID);
        t->SetBranchAddress("fEdep",&fEdep);
        t->SetBranchAddress("fSensVolID",&fSensVolID);
        t->SetBranchAddress("fTotalEnergy",&fTotalEnergy);
*/
    //MGTMCEventSteps *eventSteps = 0;
    //t->SetBranchAddress("eventSteps", &eventSteps);
    
    //MGTMCStepData *step;

       	nentriest=t->GetEntries(); // Get number of entries in the MJD Tree
        cout << "nentriest = " << nentriest << endl;
        // LOOP OVER ENTRIES IN CURRENT MaGe TREE //
        for(int k=0;k<nentriest;k++)
        {
            t->GetEntry(k); // Read all branches of entry and return total number of bytes read.
            //if(fSensVolID->at(0)==1000401)
            //{
                cout << k << " " << fTotalEnergy << endl;
                hEnergy->Fill(fTotalEnergy);
            //}
        } // End of loop over entries k in run i
    
    cEnergy->cd();
    hEnergy->Draw();
        
    cout << "F I N I S H E D" << endl;
    App->Run();
} // End of main()
