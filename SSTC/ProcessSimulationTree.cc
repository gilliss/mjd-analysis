/*
 Need a handle on the initial decay particle energies (initial KE doesn't seem to do the trick) in order to calculate efficiency of capturing the full E of >10 keV K-shell e-captures
    -sometimes the initKineticEnergy (of step 2) is less than the total totalEnergyDep
    ^should be using highest KE value, not that of step 2, but that would take time to implement
    -could see how often a >10keV event leaves 1000401
    ^the above would be perfect, but I have only implemented a counting of how many exscaped events then find another detector. This does not include escaped events that fully escape the array of sensitive volumes. Based on what is currently implemented, the lower limit to % escape is 4/10000=.0004=0.04%
    ^there are also three events, btw 8 & 9 keV, that don't match the e-capture spectrum. Could add those in as possible escapes.
*/
 

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include "TH2D.h"
#include "TProof.h"
#include "TROOT.h"
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <getopt.h>

//#include "TBranchElement.h"
#include "MGTMCEventSteps.hh"

#include <TApplication.h>

using namespace std;

int main(int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);

    cout.precision(15);
    cout << "Start " << endl;
    cout << "-------" << endl;
    //gROOT->ProcessLine(".x $MGDODIR/Root/LoadMGDOClasses.C");
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Initialize Tree from Simulation Data
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TChain chain("fTree");
    TString path = "/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/InitialData/";
    TString file = "68GeInCrystal_ElectronCaptures_10000Events.root";
    TString filename;
    filename = path;
    filename += file;
    chain.Add(filename);
    
    cout << "Added files" << endl;
    cout << "-------" << endl;
    Long64_t nentries = (Long64_t)chain.GetEntries();
    cout << chain.GetNtrees() << " tree(s) with " << nentries << " entries " <<endl;
    cout << "-------" << endl;
    
    //Branch Structure with member vars??
    MGTMCEventSteps *eventSteps = 0; // some object loosely related to a TBranchElement
    chain.SetBranchAddress("eventSteps", &eventSteps);
    const MGTMCStepData *step;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Initialize Vars
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    Int_t ParticleID;
    Int_t sensitiveVolume;
    Double_t timestamp; //timestamp = step->GetT()*1e-9;
    Double_t stepEnergyDep; //stepEnergyDep = step->GetEdep()*1000;
    Double_t totalEnergyDep;
    Double_t initKineticEnergy;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Initialize Plots
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TH1D *hStepEnergy;
    int nbins = 120;
    double minbin = 0;
    double maxbin = 0.012*1000;
    char hname[200];
    sprintf(hname,"hist_StepEnergy");
    hStepEnergy= new TH1D(hname, hname, nbins, minbin, maxbin);
    
    TH1D *hTotalEnergy;
    sprintf(hname,"hist_TotalEnergy");
    hTotalEnergy= new TH1D(hname, hname, nbins, minbin, maxbin);
    
    TH1D *hInitKEnergy;
    sprintf(hname,"hist_InitKEnergy");
    hInitKEnergy= new TH1D(hname, hname, nbins, minbin, maxbin);
    
    TCanvas* cStepEnergy;
    cStepEnergy = new TCanvas();
    
    TCanvas* cTotalEnergy;
    cTotalEnergy = new TCanvas();
    
    TCanvas* cInitKEnergy;
    cInitKEnergy = new TCanvas();
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Fill Plots
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    int mismatchCount=0;
    int escapeCount=0;
    int lowTotalECount=0;
    for (Int_t i=0;i<nentries;i++)
    {
        chain.GetEntry(i);
        totalEnergyDep = 0.0;
        initKineticEnergy = 0.0;
        for (Int_t k=0;k<eventSteps->GetNSteps();k++)
        {
            step = eventSteps->GetStep(k);
            sensitiveVolume = step->GetSensitiveVolumeID(); // MGMCStepData.hh for all "Get" commands
            stepEnergyDep = step->GetEdep();
            
            if (sensitiveVolume == 1000401)
            {
                hStepEnergy->Fill(stepEnergyDep*1000);
                totalEnergyDep += stepEnergyDep;
                if (k==2) {
                    initKineticEnergy = step->GetKineticE();
                    hInitKEnergy->Fill(initKineticEnergy*1000);
                }
            }
            
            cout << i << " " <<  k << " " << sensitiveVolume << "  " << stepEnergyDep << "  " << initKineticEnergy << endl;
            //cout << i << " " << k << " " << initKineticEnergy << endl;
        }
        hTotalEnergy->Fill(totalEnergyDep*1000);
        cout << "   ^^" << totalEnergyDep << endl;
        if (totalEnergyDep > .01 && (initKineticEnergy/totalEnergyDep)<0.998) {
            cout << "   ====mismatch in event i = " << i << endl;
            mismatchCount++;
        }

        if(sensitiveVolume!=1000401) {
            escapeCount++;
            cout << "   ====escape in event i = " << i << ", escapeCount = " << escapeCount << endl;
        }
        
        if(totalEnergyDep < 0.009 && totalEnergyDep > 0.008) {
            lowTotalECount++;
            cout << "   ====~.008 in event i = " << i << endl;
        }
    }
    cout << "lowTotalECount = " << lowTotalECount << endl;
    cout << "mismatchCount = " << mismatchCount << endl;
    cout << "escapeCount = " << escapeCount-1 << endl; // -2 for fake first and last events
    cout << "Percentage of e-captures escaping to other detectors = ..." << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Plot Plots
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    cStepEnergy->cd();
    cStepEnergy->SetLogy();
    hStepEnergy->SetTitle("Energy Depositions of Simulation Steps");
    hStepEnergy->GetYaxis()->SetTitle("Counts");
    hStepEnergy->GetXaxis()->SetTitle("Edeps (keV)");
    hStepEnergy->Draw();
    cStepEnergy->Update();
    
    cTotalEnergy->cd();
    cTotalEnergy->SetLogy();
    hTotalEnergy->SetTitle("Total Energy Depositions of Simulation Events");
    hTotalEnergy->GetYaxis()->SetTitle("Counts");
    hTotalEnergy->GetXaxis()->SetTitle("TotalEnergy (keV)");
    hTotalEnergy->Draw();
    cTotalEnergy->Update();
        TAxis *xaxis = hTotalEnergy->GetXaxis();
        Int_t binx1_FullE = xaxis->FindBin(10);
        Int_t binx2_FullE = xaxis->FindBin(maxbin);
        Int_t binx1_LesserE = xaxis->FindBin(minbin); // maybe don't include the lowest population
        Int_t binx2_LesserE = xaxis->FindBin(10);
        Double_t integral_TotalE = hTotalEnergy->Integral(binx1_FullE,binx2_FullE);
        //Double_t integral_LesserE = hTotalEnergy->Integral(binx1_LesserE,binx2_LesserE);
        cout << "TotalE>10 Integral = " << integral_TotalE << endl; //" " << "LesserE Integral = " << integral_LesserE <<endl;
        //cout << "Efficiency of tagging 68Ge e-capture btw 10-11keV = " << integral_FullE/(integral_FullE+integral_LesserE-2) << " =? " << integral_FullE/10000 << endl; // "-2" for fake first and last events
    
    cInitKEnergy->cd();
    cInitKEnergy->SetLogy();
    hInitKEnergy->SetTitle("Initial Kinetic Energies of Simulated Decay Particles");
    hInitKEnergy->GetYaxis()->SetTitle("Counts");
    hInitKEnergy->GetXaxis()->SetTitle("Energy (keV)");
    hInitKEnergy->Draw();
    cInitKEnergy->Update();
        xaxis = hInitKEnergy->GetXaxis();
        binx1_FullE = xaxis->FindBin(10);
        binx2_FullE = xaxis->FindBin(maxbin);
        binx1_LesserE = xaxis->FindBin(minbin);
        binx2_LesserE = xaxis->FindBin(10);
        Double_t integral_KineticE = hInitKEnergy->Integral(binx1_FullE,binx2_FullE);
        cout << "KineticE>10 Integral = " << integral_KineticE << endl;

    //cout << "Efficiency of tagging 68Ge e-capture btw 10-11keV = " << integral_TotalE/integral_KineticE << endl;
    
    cout << "F I N I S H E D" << endl;
    App->Run();
}
