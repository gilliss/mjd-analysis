/*
This script reads in Ben Shanks' file_for_ben.root which contains simulated WFs. This script makes WF objects out of the simulated WFs and tries to find their t0.
For a given WF index and noise level, the script plots one waveform and marks where the threshold and estimated t0 are.
 
 Current algorithm
 1. Sample the baseline and CALCULATE the mean and stddev
 2. Find the WF maximum
 3. Set threshold to mean+n*stddev
 4. Find the last timepoint, prior to the max, that the WF crossed the threshold
*/

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include <TApplication.h>
#include <MGTEvent.hh>
#include <MGTRun.hh>
#include <MGTWaveform.hh>
#include <MGWFTrapezoidalFilter.hh>
#include <MGWFTimePointCalculator.hh>
#include <MGWFExtremumFinder.hh>
#include <MGWFBaselineRemover.hh>
#include <MGWFLinearFit.hh>

using namespace std;


double getMean(vector<double> wfVec)
{
    double mean = 0.0;
    double sum = 0.0;
    for(int i = 50; i<=700; i++) sum += wfVec.at(i);
    mean = sum/(700-50+1.0);
    return mean;
}
double getStdDev(vector<double> wfVec)
{
    double stdDev = 0.0;
    double variance = 0.0;
    double temp = 0.0;
    double mean = 0.0;
    double sum = 0.0;
    
    for(int i = 50; i<=700; i++)
    {
        sum += wfVec.at(i);
        //cout << i << " " << wfVec.at(i) << endl;
    }
    mean = sum/(700-50+1.0);
    for(int i = 50; i<=700; i++)
    {
        temp += (mean-wfVec.at(i))*(mean-wfVec.at(i));
    }
    
    variance = temp/(700-50+1-1.0);
    stdDev = sqrt(variance);
    //cout << stdDev << endl;
    return stdDev;
}

int main (int argc, char* argv[])
{
    int entryIndex = atoi(argv[1]);
    double idealSigma = atof(argv[2]); // 10*0.00125
    
    TApplication *App = new TApplication("App", 0, NULL);

    char infilepath[200];
    sprintf (infilepath,"/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/t0_Finder/SimulatedWFs_BenShanks/file_for_ben.root");

    int branchVectorSize;
    vector<double> newWFVec;
    double rGausNum;
    double newWFVal;
    
    TChain *t1 = new TChain("waveformTree");
    t1->AddFile(infilepath);
    vector<double>* waveform;
    t1->SetBranchAddress("waveform",&waveform);
    
    int nentries=t1->GetEntries();
    cout << "nentries " << nentries << endl;
    
    
    
    t1->GetEntry(entryIndex);
    branchVectorSize = waveform->size();
    TRandom* rNum = new TRandom();
    rNum->SetSeed(0); // using 0 should invoke a new UUID each time this line is called
    int m = 0;
    for (int k=0; k<2000; k++)
    {
        rGausNum = rNum->Gaus(0,idealSigma);
        if(k>=1200 && m<branchVectorSize) // the simulated t0 should be at newWFVec[1200]
        {
            newWFVal = rGausNum + waveform->at(m);
            newWFVec.push_back(newWFVal);
            m++;
        }
        else newWFVec.push_back(rGausNum);
    }
    cout << newWFVec.size() << endl;
    
    MGTWaveform* simWF = new MGTWaveform();
    simWF->SetData(newWFVec);//*waveform);
    
    TCanvas* cWF = new TCanvas();
    TH1D* hWF = new TH1D();
    simWF->LoadIntoHist(hWF); // MGTWaveform.hh
    cWF->cd();
    hWF->Draw();
    cWF->Update();
    
    
    //////////////////////////////////////////////////////////////////////
    // Implementing the t0 estimator (below). The code above can be run indpendently by commenting out everything below
    //////////////////////////////////////////////////////////////////////
    
    
    int nPlotsOnCanvas = 2; // rawWF, trapWF
    TCanvas* cDiv = new TCanvas();
    cDiv->Divide(1,nPlotsOnCanvas);
    char title[200];
    
    TH1D* hwave[2];// = new TH1D(); // hist for waveform
    for (int i=0; i<2; i++) hwave[i] = new TH1D(); // hist for waveform
    int extrmPt;
    double extrmVl;
    TLine *l;
    TLine *l2;
    TLine *lMax;
    
    double stdDev, mean;
    
    double crossing;
    double threshold;
    
    MGWFExtremumFinder* EXF = new MGWFExtremumFinder();
    EXF->SetFindMaximum(true);
    EXF->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
    EXF->Transform(simWF);
    extrmPt = EXF->GetTheExtremumPoint();
    extrmVl = EXF->GetTheExtremumValue();
    cout << "Extremum (Max) Point (in samples)" << extrmPt << endl;
    cout << "Extremum (Max) Value (in samples)" << extrmVl << endl;

    simWF->LoadIntoHist(hwave[0]); // MGTWaveform.hh
    sprintf(title,"Entry %d",entryIndex);
    hwave[0]->SetTitle(title);
    hwave[0]->SetStats(0);
    cDiv->cd(1);
    hwave[0]->Draw();
    
    lMax = new TLine(0,extrmVl,20000,extrmVl); // (x1,y1,x2,y2) // 10ns/sample
    lMax->SetLineStyle(7);
    lMax->SetLineColor(15);
    lMax->Draw();
    
    mean = getMean(newWFVec);
    stdDev = getStdDev(newWFVec);
    threshold = mean+0.5*stdDev;
    if (threshold<=0)
    {
        cout << "WARNING: threshold = " << threshold << " < 0 ... setting threshold to 0" << endl;
        cout << "This is b/c of a suspected MGWFTimePointCalculator issue, in which its methods do not like finding 0%-of-max timepoints" << endl;
        threshold = 1e-20;
    }
    cout << "mean " << mean << " stdDev " << stdDev << " threshold " << threshold << endl;
    
    MGWFTimePointCalculator* TPC = new MGWFTimePointCalculator();
    TPC->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
    TPC->Transform(simWF);
    TPC->AddPoint(threshold/extrmVl);
    TPC->FindTimePoints(*simWF);
    crossing = TPC->GetFromMaxRiseTime(0);
    cout << "threshold crossing prior to max " << crossing << " (ideal is 1200)" << endl;
    
    cDiv->cd(2);
    simWF->LoadIntoHist(hwave[1]); // MGTWaveform.hh
    hwave[1]->GetXaxis()->SetRangeUser(crossing-200,crossing+200);
    hwave[1]->GetYaxis()->SetRangeUser(threshold-3*idealSigma,threshold+4*idealSigma);
    hwave[1]->SetStats(0);
    hwave[1]->Draw();
    l = new TLine(crossing,threshold-1,crossing,threshold+1); // (x1,y1,x2,y2) // 10ns/sample
    l->SetLineStyle(7);
    l->SetLineColor(15);
    l2 = new TLine(crossing-200,threshold,crossing+200,threshold); // (x1,y1,x2,y2) // 10ns/sample
    l2->SetLineStyle(7);
    l2->SetLineColor(8);
    l->Draw();
    l2->Draw();
    
    cDiv->Update();
    
    //for (int i=0; i<nPlotsOnCanvas; i++) delete hwave[i];
    delete EXF;
    delete TPC;

    delete t1;
    delete simWF;
    delete rNum;
    delete lMax;
    
    App->Run();
}

/*cout << "looping over events ..." << endl;
 for(int i=0;i<nentries;i++)
 {
 
 }*/

/*
 for (int j=0; j<branchVectorSize; j++)
 {
 cout << waveform->at(j) << endl;
 }*/
