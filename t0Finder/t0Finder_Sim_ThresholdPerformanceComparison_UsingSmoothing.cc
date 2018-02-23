/*
This script reads in Ben Shanks' file_for_ben.root which contains simulated WFs. The script makes WF objects out of the simulated WFs and tries to find their t0 based on a threshold-crossing algorithm.
For a given noise level, this script plots the results of multiple threshold settings and overlays them for comparison.
Current algorithm
    0. Smooth the WF using Savitsky-Golay method
    1. Sample the baseline and CALCULATE the mean and stddev
    2. Find the WF maximum
    3. Set threshold to mean+n*stddev
    4. Find the last timepoint, prior to the max, that the WF crossed the threshold

Issues:
Note that for a 0% of max threshold, MGWFTimePointCalculator::GetFromMaxRiseTime returns a NaN. Or at least some bug in my code makes it seem like this is the case.
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
    //cout << idealSigma << endl;
    
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
    cout << "nentries " << nentries << endl;
    
    int extrmPt;
    double extrmVl;
    double stdDev, mean;
    double crossing;
    double threshold;
    double result;
    double distance;
    
    int nThreshCases = 2;
    TH2D* h2[nThreshCases];
    char title[200], name[200];
    sprintf(title,"Distance vs Deviation from Sim'd t0, GausNoiseSigma=%2.2f%% of WF Max",idealSigma*100);
    
    /////////////// looping //////////////////
    /////////////// looping //////////////////
    /////////////// looping //////////////////

    for (int threshCase = 0; threshCase < nThreshCases; threshCase++)
    {
        //cout << " threshold = mean + " << (1+threshCase*0.5) << "*stdDev" << endl;
        //sprintf(name,"ThreshMeanPlus%2.1fStdDev",1+threshCase*0.5);
        switch (threshCase)
        {
            case 0:
                cout << " threshold = mean + 1*stdDev" << endl;
                sprintf(name,"ThreshMeanPlus1stdDev");
                break;
            case 1:
                cout << " threshold = mean + 3*stdDev" << endl;
                sprintf(name,"ThreshMeanPlus3stdDev");
                break;
            default:
                cout << "Not a valid threshold case" << endl;
                break;
        }
        h2[threshCase] = new TH2D(name,title,220,-10,100,130,0,65);
        
        cout << "looping over events ..." << endl;
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
            MGWFSavitzkyGolaySmoother* SGS = new MGWFSavitzkyGolaySmoother(24,0,5,"MGWFSavitzkyGolaySmoother"); // (smoothSize,derivativeOrder,polynomialDegree,name)  // 4,0,2 seems good
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
            //threshold = mean+(1+threshCase*0.5)*stdDev;
            switch (threshCase)
            {
                case 0:
                    threshold = mean+1*stdDev;
                    break;
                case 1:
                    threshold = mean+3*stdDev;
                    break;
                default:
                    cout << "Not a valid threshold case" << endl;
                    break;
            }
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
            //cout << "   pt 100%: " <<  TPC->GetFromMaxRiseTime(1) << endl; /// for zero noise, this returns a nan
            //cout << "   pt 50%: " <<  TPC->GetFromMaxRiseTime(2) << endl; /// for zero noise, this returns a nan
            result = crossing - 1200.0;
            distance = sqrt(r->at(0)*r->at(0)+z->at(0)*z->at(0));
            
            h2[threshCase]->Fill(result,distance);
            //cout << "(" << result << "," << distance << ")" << endl;
            
            delete rNum;
            delete simWF;
            delete smoothSimWF;
            delete EXF;
            delete TPC;
        } // End loop over sim'd WFs
    } // End loop over threshold settings
    
    /////////////// Plotting //////////////////
    /////////////// Plotting //////////////////
    /////////////// Plotting //////////////////
    
    TCanvas* c2 = new TCanvas();
    c2->cd();
    int colorarray[10]={1,6};//{1,2,3,4,6,7,8,9,48}; //https://root.cern.ch/doc/master/classTAttMarker.html#M1
    int colorarray_index=0;
    TLegend* lOverlay = new TLegend(0.8,0.2,0.97,0.5); //https://root.cern.ch/doc/master/classTLegend.html
    lOverlay->SetHeader("Threshold");
    TString threshString;
    
    for (int threshCase = 0; threshCase < nThreshCases; threshCase++)
    {
        h2[threshCase]->SetStats(0);
        h2[threshCase]->SetMarkerStyle(6); // https://root.cern.ch/doc/master/classTAttMarker.html#M2
        h2[threshCase]->SetMarkerSize(1);
        h2[threshCase]->SetMarkerColor(colorarray[colorarray_index]);
        //threshString.Form("mean+%2.1f*sigma",(1+threshCase*0.5));
        switch (threshCase)
        {
            case 0:
                threshString.Form("mean+1*sigma");
                break;
            case 1:
                threshString.Form("mean+3*sigma");
                break;
            default:
                cout << "Not a valid threshold case" << endl;
                break;
        }
        lOverlay->AddEntry(h2[threshCase],threshString,"p"); // https://root.cern.ch/doc/v606/classTLegend.html
        if (threshCase == 0)
        {
            h2[threshCase]->GetXaxis()->SetTitle("Deviation from Sim'd t0 (samples)");
            h2[threshCase]->GetYaxis()->SetTitle("Distance from point contact");
            h2[threshCase]->Draw();
        }
        else
        {
            h2[threshCase]->Draw("SAME");
        }
        colorarray_index++;
    }
    lOverlay->Draw();
    c2->Update();
    
    TImage *img = TImage::Create(); // https://root.cern.ch/root/html/tutorials/image/pad2png.C.html
    img->FromPad(c2);
    TString imageString;
    imageString.Form("smoothThresholdOverlay_%2.0fPctGausNoise_24-0-5_1and3StdDev.png",idealSigma*100);
    img->WriteImage(imageString);
    
    //delete c2;
    //delete img;
    //delete lOverlay;
    //for (int threshCase = 0; threshCase < nThreshCases; threshCase++) delete h2[threshCase];

    App->Run();
}



/*
 for (int j=0; j<branchVectorSize; j++)
 {
 cout << waveform->at(j) << endl;
 }*/
