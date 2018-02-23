/*
This script reads in Ben Shanks' file_for_ben.root which contains simulated WFs. This script makes WF objects out of the simulated WFs and tries to find their t0.
For a given WF index and noise level, the script plots one waveform and marks where the threshold and estimated t0 are.
 
 Current algorithm
 0. Smooth the WF using Savitsky-Golay method
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
#include <TImage.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <TApplication.h>
#include <MGTEvent.hh>
#include <MGTRun.hh>
#include <MGTWaveform.hh>
#include <MGWFTrapezoidalFilter.hh>
#include <MGWFTimePointCalculator.hh>
#include <MGWFExtremumFinder.hh>
#include <MGWFBaselineRemover.hh>
#include <MGWFLinearFit.hh>
#include <MGWFSavitzkyGolaySmoother.hh>

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
    double idealSigma = atof(argv[2]);
    
    TApplication *App = new TApplication("App", 0, NULL);

    char infilepath[200];
    sprintf (infilepath,"/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/t0_Finder/SimulatedWFs_BenShanks/file_for_ben.root");
    
    /* LOOP THROUGH THE FILE, PULL OUT A WF, APPLY NOISE TO IT */

    int branchVectorSize;
    vector<double> newWFVec; // lengthened and noisified WF vector, built from the original sim'd WF
    
    TChain *t1 = new TChain("waveformTree");
    t1->AddFile(infilepath);
    vector<double>* waveform; // original sim'd WF vector
    t1->SetBranchAddress("waveform",&waveform);
    
    int nentries=t1->GetEntries();
    cout << "nentries in waveformTree = " << nentries << endl;
    
    t1->GetEntry(entryIndex);
    branchVectorSize = waveform->size();
    TRandom* rNum = new TRandom();
    rNum->SetSeed(0); // using 0 should invoke a new UUID each time this line is called
    double rGausNum;
    double newWFVal;
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
    cout << "newWFVec (simWF) size = " << newWFVec.size() << endl;
    
    MGTWaveform* simWF = new MGTWaveform();
    simWF->SetData(newWFVec);
    
    /* DRAW THE NOISIFIED WF */
    
    TCanvas* cWF = new TCanvas();
    TH1D* hWF = new TH1D();
    simWF->LoadIntoHist(hWF); // MGTWaveform.hh
    cWF->cd();
    hWF->GetYaxis()->SetRangeUser(-0.2,1.2);
    hWF->SetStats(0);
    TString hWFTitle;
    hWFTitle.Form("Sim'd WF Entry %d, GausNoiseStdDev %2.1f%% of MaxWF",entryIndex,idealSigma*100);
    hWF->SetTitle(hWFTitle);
    hWF->SetLineWidth(2);
    hWF->SetLineColor(38);
    hWF->Draw();
    cWF->Update();
    
    /*
    double ns = 1.e-9;
    double us = 1.e-6;
    double ms = 1.e-3;
    MGTWaveform* trapWF = new MGTWaveform(); // Declare WF object for TrapWF
    MGWFTrapezoidalFilter* TZF = new MGWFTrapezoidalFilter(); // declare trapezoidal filter
    TZF->SetRampTime(50);//.1*us/ns);
    TZF->SetFlatTime(50);//.010*us/ns);
    TZF->SetRestingBaseline(0.0);
    TZF->SetDecayConstant(90000);
    TZF->TransformOutOfPlace(*simWF,*trapWF); // Transform rawWF into trapWF
    TCanvas* cWF22 = new TCanvas();
    TH1D* hWF22 = new TH1D();
    trapWF->LoadIntoHist(hWF22); // MGTWaveform.hh
    cWF22->cd();
    hWF22->GetYaxis()->SetRangeUser(-0.2,1.2);
    hWF22->SetStats(0);
    hWF22->SetLineWidth(2);
    hWF22->SetLineColor(38);
    hWF22->Draw();
    cWF22->Update();
    delete trapWF;
    delete TZF;
    */
    
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///// Implementing the t0 estimator (below). The code above can be run indpendently by commenting out
///// everything below
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

    
    int nPlotsOnCanvas = 2; // plots for the noisified WF, and the same WF with t0 estimates
    TCanvas* cDiv = new TCanvas();
    cDiv->Divide(1,nPlotsOnCanvas);
    char title[200];
    
    TH1D* hwave[2];// = new TH1D(); // hist for waveform
    for (int i=0; i<2; i++) hwave[i] = new TH1D(); // hist for waveform
    TLine *lCross;
    TLine *lThresh;
    TLine *lMax;
    
    clock_t t; // http://www.cplusplus.com/reference/ctime/clock/
    t = clock();
    
    int extrmPt;
    double extrmVl;
    double stdDev, mean;
    double crossing;
    double threshold;
    
    /* SMOOTH THE WF */
    /* Both S/N ~and signal distortion~ increase as smoothing width increases and as the degree of the polynomial decreases */
    
    MGWFSavitzkyGolaySmoother* SGS = new MGWFSavitzkyGolaySmoother(24,0,5,"MGWFSavitzkyGolaySmoother"); // (smoothSize,derivativeOrder,polynomialDegree,name)
    MGTWaveform* smoothSimWF = new MGTWaveform();
    SGS->TransformOutOfPlace(*simWF,*smoothSimWF);
    TH1D* hWF2 = new TH1D();
    smoothSimWF->LoadIntoHist(hWF2); // MGTWaveform.hh
    cWF->cd();
    hWF2->GetYaxis()->SetRangeUser(-0.2,1.2);
    hWF2->SetStats(0);
    hWF2->SetLineWidth(2);
    hWF2->SetLineColor(4);
    hWF2->Draw("SAME");
    cWF->Update();
    /*
    TImage *img = TImage::Create(); // https://root.cern.ch/root/html/tutorials/image/pad2png.C.html
    img->FromPad(cWF);
    TString imageString;
    imageString.Form("SimWFvsSmooth_%2.0fPctGausNoise_Smooth4-0-2.png",idealSigma*100);
    img->WriteImage(imageString);
    delete img;
    */
    
    /* FIND THE WF MAX AND DRAW THE RESULT */
    
    MGWFExtremumFinder* EXF = new MGWFExtremumFinder();
    EXF->SetFindMaximum(true);
    EXF->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
    EXF->Transform(smoothSimWF);
    extrmPt = EXF->GetTheExtremumPoint();
    extrmVl = EXF->GetTheExtremumValue();
    cout << "Extremum (Max) Point (in samples)" << extrmPt << endl;
    cout << "Extremum (Max) Value (in samples)" << extrmVl << endl;

    smoothSimWF->LoadIntoHist(hwave[0]); // MGTWaveform.hh
    sprintf(title,"Smoothed Entry %d",entryIndex);
    hwave[0]->SetTitle(title);
    hwave[0]->SetStats(0);
    cDiv->cd(1);
    hwave[0]->Draw();
    
    lMax = new TLine(0,extrmVl,20000,extrmVl); // (x1,y1,x2,y2) // 10ns/sample
    lMax->SetLineStyle(7);
    lMax->SetLineColor(15);
    lMax->Draw();
    
    /* SET THE THRESHOLD, FIND CROSSING, DRAW RESULT */

    vector<double> smoothSimWFVector;
    smoothSimWFVector = smoothSimWF->GetVectorData(); // MGWaveform.hh
    cout << "smoothSimWFVector (smoothSimWF) size = " << smoothSimWFVector.size() << endl;
    
    mean = getMean(smoothSimWFVector);
    stdDev = getStdDev(smoothSimWFVector);
    threshold = mean+3*stdDev;
    if (threshold<=0)
    {
        cout << "WARNING: threshold = " << threshold << " < 0 ... setting threshold to 0" << endl;
        cout << "This is b/c of a suspected MGWFTimePointCalculator issue, in which its methods do not like finding 0%-of-max timepoints" << endl;
        threshold = 1e-20;
    }
    cout << "mean " << mean << " stdDev " << stdDev << " threshold " << threshold << endl;
    
    MGWFTimePointCalculator* TPC = new MGWFTimePointCalculator();
    TPC->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
    TPC->Transform(smoothSimWF);
    TPC->AddPoint(threshold/extrmVl);
    TPC->FindTimePoints(*smoothSimWF);
    crossing = TPC->GetFromMaxRiseTime(0);
    cout << "threshold crossing prior to max " << crossing << " (ideal is 1200)" << endl;
    
    t = clock() - t;
    //printf ("Program took %d clicks (%f seconds).\n",((int)t),((float)t)/CLOCKS_PER_SEC);
    cout << "Program took " << t << " clicks, which is " << ((float)t)/CLOCKS_PER_SEC << " seconds" << endl;
    
    cDiv->cd(2);
    smoothSimWF->LoadIntoHist(hwave[1]); // MGTWaveform.hh
    hwave[1]->GetXaxis()->SetRangeUser(crossing-200,crossing+200);
    hwave[1]->GetYaxis()->SetRangeUser(threshold-3*idealSigma,threshold+4*idealSigma);
    sprintf(title,"Smoothed Entry %d, Zoom on Result",entryIndex);
    hwave[1]->SetTitle(title);
    hwave[1]->SetStats(0);
    hwave[1]->Draw();
    lCross = new TLine(crossing,threshold-1,crossing,threshold+1); // (x1,y1,x2,y2) // 10ns/sample
    lCross->SetLineStyle(7);
    lCross->SetLineColor(15);
    lThresh = new TLine(crossing-200,threshold,crossing+200,threshold); // (x1,y1,x2,y2) // 10ns/sample
    lThresh->SetLineStyle(7);
    lThresh->SetLineColor(8);
    lCross->Draw();
    lThresh->Draw();
    
    cDiv->Update();
    
    /* CLOSE OUT */
    
    //for (int i=0; i<nPlotsOnCanvas; i++) delete hwave[i];
    delete SGS; // Smoothing transform
    delete EXF; // Extremum finding transform
    delete TPC; // timepoint finding transform

    delete t1;
    delete simWF;
    delete smoothSimWF;
    delete rNum;
    //delete lMax;
    //delete lCross;
    //delete lThresh;
    
    App->Run();
}

//Scratch work from figuring out MGWFSavitzkyGolaySmoother
//MGWaveformRegion* wfRegion = new MGWaveformRegion();
//wfRegion->SetBeginning(0);
//wfRegion->SetEnd(2000);
//cout << "smoothing region vs simWF size: " << wfRegion->GetLength() << " = " << simWF->GetLength() << endl;
//vector<double> smoothNewWFVec; // smoothed version of newWFVec (simWF)