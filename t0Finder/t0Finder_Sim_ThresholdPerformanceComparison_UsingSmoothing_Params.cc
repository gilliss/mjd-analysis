/*
This script reads in Ben Shanks' file_for_ben.root which contains simulated WFs. The script makes WF objects out of the simulated WFs and tries to find their t0 based on a threshold-crossing algorithm.
For a given noise level, this script plots the results of multiple threshold settings and overlays them for comparison.
Current algorithm
    0. Smooth the WF using Savitsky-Golay method
    1. Sample the baseline and CALCULATE the mean and stddev
    2. Find the WF maximum
    3. Set threshold to mean+n*stddev
    4. Find the last timepoint, prior to the max, that the WF crossed the threshold
 
 Use:  ./t0Finder_Sim_ThresholdPerformanceComparison_UsingSmoothing_Params <idealSigma> <nStdDevs>

Issues:
 -Note that for a 0% of max threshold, MGWFTimePointCalculator::GetFromMaxRiseTime returns a NaN. Or at least some bug in my code makes it seem like this is the case.
 -As the size of the smoothing window grows and, possibly as the degree of the polynomial increases, the fitting takes increasingly longer.
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
#include <TColor.h>

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
    double idealSigma = atof(argv[1]);
    int nStdDevs = atof(argv[2]);
    
    TApplication *App = new TApplication("App", 0, NULL);

    char infilepath[200];
    sprintf (infilepath,"/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/t0_Finder/SimulatedWFs_BenShanks/file_for_ben.root");

    int branchVectorSize;
    vector<double> newWFVec;
    double rGausNum;
    double newWFVal;
    int m;
    
    vector<double> smoothSimWFVector;
    
    TChain *t1 = new TChain("waveformTree");
    t1->AddFile(infilepath);
    vector<double>* waveform;
    t1->SetBranchAddress("waveform",&waveform);
    vector<double>* r;
    t1->SetBranchAddress("r",&r);
    vector<double>* z;
    t1->SetBranchAddress("z",&z);
    
    int nentries=t1->GetEntries();
    cout << "nentries in waveformTree " << nentries << endl;
    
    int extrmPt;
    double extrmVl;
    double stdDev, mean;
    double crossing;
    double threshold;
    double result;
    double distance;
    
    vector<double> resultVec;
    double resultVecMean, resultVecStdDev;
    double goodness;
    int smoothSizeMin = 0;
    int smoothSizeMax = 30; // 20
    int polyOrderMin = 0; // 6
    int polyOrderMax = 10; // 6
    TH2D* h3dParams = new TH2D("GoodnessOfSmoothingParams","Goodness of Smoothing Parameters",smoothSizeMax-smoothSizeMin,smoothSizeMin,smoothSizeMax,polyOrderMax-polyOrderMin,polyOrderMin,polyOrderMax);
    char title[200], name[200];
    sprintf(title,"Distance vs Deviation from Sim'd t0, GausNoiseSigma=%2.2f%% of WF Max",idealSigma*100);
    
    /////////////// looping //////////////////
    /////////////// looping //////////////////
    /////////////// looping //////////////////

    for (int smoothRegion = smoothSizeMin; smoothRegion < smoothSizeMax; smoothRegion++)
    {
        for (int polyOrder = polyOrderMin; polyOrder < polyOrderMax; polyOrder++)
        {
            if (2*smoothRegion >= polyOrder)
            {
                resultVec.resize(0);

                cout << "smoothSize " << smoothRegion << " polyOrder " << polyOrder << endl;
                for(int i=0;i<nentries;i++)
                {
                    /* CLEAR OBJECTS */
                    newWFVec.resize(0);
                    smoothSimWFVector.resize(0);
                    
                    /* MAKE THE SIM'D WF */
                    t1->GetEntry(i);
                    branchVectorSize = waveform->size();
                    TRandom* rNum = new TRandom();
                    rNum->SetSeed(0); // using 0 should invoke a new UUID each time this line is called
                    m = 0;
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
                    //cout << newWFVec.size() << endl;
                    
                    MGTWaveform* simWF = new MGTWaveform();
                    simWF->SetData(newWFVec);//*waveform);
                    
                    ////////////////////////////////////
                    /////// apply transforms ///////////
                    ////////////////////////////////////
                    
                    /* SMOOTH THE WF */
                    MGWFSavitzkyGolaySmoother* SGS = new MGWFSavitzkyGolaySmoother(smoothRegion,0,polyOrder,"MGWFSavitzkyGolaySmoother"); // (smoothSize,derivativeOrder,polynomialDegree,name)  // 4,0,2 seems good
                    MGTWaveform* smoothSimWF = new MGTWaveform();
                    SGS->TransformOutOfPlace(*simWF,*smoothSimWF);
                    smoothSimWFVector = smoothSimWF->GetVectorData(); // MGWaveform.hh
                    
                    /* FIND THE WF MAX */
                    MGWFExtremumFinder* EXF = new MGWFExtremumFinder();
                    EXF->SetFindMaximum(true);
                    EXF->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
                    EXF->Transform(smoothSimWF);
                    extrmPt = EXF->GetTheExtremumPoint();
                    extrmVl = EXF->GetTheExtremumValue();
                    //cout << " extrmVl: " << extrmVl << endl;
                    
                    mean = getMean(smoothSimWFVector);
                    stdDev = getStdDev(smoothSimWFVector);
                    threshold = mean+nStdDevs*stdDev;
                    if (threshold<=0)
                    {
                        cout << "WARNING: threshold = " << threshold << " < 0 ... setting threshold to 0" << endl;
                        cout << "This is b/c of a suspected MGWFTimePointCalculator issue, in which its methods do not like finding 0%-of-max timepoints" << endl;
                        threshold = 1e-20;
                    }
                    //cout << "mean " << mean << " stdDev " << stdDev << " threshold " << threshold << endl;
                    
                    MGWFTimePointCalculator* TPC = new MGWFTimePointCalculator();
                    TPC->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
                    TPC->Transform(smoothSimWF);//TransformInPlace(*simWF);
                    TPC->AddPoint(threshold/extrmVl);
                    TPC->AddPoint(extrmVl/extrmVl);
                    TPC->AddPoint(0.5);
                    TPC->FindTimePoints(*smoothSimWF);
                    crossing = TPC->GetFromMaxRiseTime(0); // finding crossing for time point index 0
                    
                    result = crossing - 1200.0;
                    //distance = sqrt(r->at(0)*r->at(0)+z->at(0)*z->at(0));
                    
                    resultVec.push_back(result);
                    
                    delete rNum;
                    delete simWF;
                    delete smoothSimWF;
                    delete SGS;
                    delete EXF;
                    delete TPC;
                } // End loop over sim'd WFs
                
                //cout << "size of resultVec = " << resultVec.size() << endl;
                resultVecMean = getMean(resultVec);
                resultVecStdDev = getStdDev(resultVec);
                //goodness = abs(resultVecMean) + resultVecStdDev;
                goodness = resultVecStdDev;
                cout << "resultVecMean " << resultVecMean << " resultVecStdDev " << resultVecStdDev << " goodness " << goodness << endl;
                
                
                h3dParams->Fill(smoothRegion,polyOrder,goodness);
            } // End of if condition 2*smoothSize+1 >= degree of polynomial
        } // End loop over polynomial orders
    } // End loop over smooth region size
    
    /////////////// Plotting //////////////////
    /////////////// Plotting //////////////////
    /////////////// Plotting //////////////////
    
    TCanvas* c3d = new TCanvas();
    c3d->cd();
    h3dParams->GetXaxis()->SetTitle("smoothSize (fit over 2*smoothSize+1)");
    h3dParams->GetYaxis()->SetTitle("order of polynomial fit");
    h3dParams->GetZaxis()->SetTitle("goodness");
    h3dParams->SetStats(0);
    
    //////////// https://root.cern.ch/root/roottalk/roottalk04/0641.html
    Int_t ncol = 100;
    Int_t colors[ncol];
    TColor *col;
    Double_t dg=1/(Double_t)ncol;
    Double_t grey=0;
    //Double_t green = 0.9;
    for (Int_t i=0; i<ncol; i++)
    {
        colors[i]= i+2000; // +2000 just gets the index out of the way of any predefined indices
        //col = gROOT->GetColor(colors[i]);
        //col->SetRGB(grey, grey, grey);
        col = new TColor(colors[i],grey,grey,0.8); // (color index,r,g,b)
        grey =grey+dg;
        //green = green + dg;
    }
    h3dParams->SetContour(ncol);
    gStyle->SetPalette(100,colors);
    ////////////
    
    //gStyle->SetPalette(1); // https://root.cern.ch/doc/v606/classTHistPainter.html#HP01d  and TStyle
    // https://root.cern.ch/doc/master/classTColor.html#C05
    //gStyle->SetNumberContours(255);
    
    h3dParams->Draw("COLZ"); // https://root.cern.ch/doc/v606/classTHistPainter.html#HP01d
    c3d->Update();
    
    delete col;

    App->Run();
}



/*
 for (int j=0; j<branchVectorSize; j++)
 {
 cout << waveform->at(j) << endl;
 }*/
