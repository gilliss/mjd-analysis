/*
 This is based off of GATStartTimeFinder.hh/.cc
 It does not work right now. It segfaults at the line reading "t1->GetEntry(entryIndex);"
 I have a feeling that the smoothing isn't being implemented as expected either.
 
 */



#include "MGWFBaselineRemover.hh"
#include "MGWFExtremumFinder.hh"
#include "MGWFTimePointCalculator.hh"
#include "MGWFSavitzkyGolaySmoother.hh"
#include <TApplication.h>

#include <TChain.h>
#include <TRandom.h>
#include <MGTWaveform.hh>

#include <math.h>

using namespace std;

////////////////////////////////////////////////////
////// Define functions for mean and stdDev ////////
////////////////////////////////////////////////////

double GetMean(vector<double> wfVec, int nBaselineSamples = 200, int iBaselineStartSample = 0) {
    double mean = 0.0;
    double sum = 0.0;
    for(int i = iBaselineStartSample; i<(iBaselineStartSample+nBaselineSamples); i++) sum += wfVec.at(i);
    mean = sum/nBaselineSamples;
    return mean;
}

double GetStdDev(vector<double> wfVec, int nBaselineSamples = 200, int iBaselineStartSample = 0) {
    double mean = 0.0;
    double sum = 0.0;
    double temp = 0.0;
    double variance = 0.0;
    double stdDev = 0.0;
    
    for(int i = iBaselineStartSample; i<(iBaselineStartSample+nBaselineSamples); i++) sum += wfVec.at(i);
    mean = sum/nBaselineSamples;
    for(int i = iBaselineStartSample; i<(iBaselineStartSample+nBaselineSamples); i++)
    {
        temp += (mean-wfVec.at(i))*(mean-wfVec.at(i));
    }
    
    variance = temp/(nBaselineSamples-1.0);
    stdDev = sqrt(variance);
    return stdDev;
}

////////////////////////////////////////////////////
////// Define classes  /////////////////////////////
////////////////////////////////////////////////////

class GATStartTimeInfo
{
public:
    GATStartTimeInfo() {}
    virtual ~GATStartTimeInfo() {}
    
    virtual size_t GetNWaveforms() { return fStartTimes.size(); }
    virtual double GetStartTime(size_t iWF) { return fStartTimes[iWF]; }
    
    virtual void Resize(size_t nWFs) { fStartTimes.resize(nWFs); }
    virtual void SetStartTime(size_t iWF, double startTime) { fStartTimes[iWF] = startTime; }
    
protected:
    vector<double> fStartTimes;

    //ClassDef(GATStartTimeInfo,1)
};


class GATStartTimeFinder// : public GATMGTEventProcBase
{
    size_t blSamps, blStartSamp, thresholdSigma;
    double extrmVl;
    vector<double> waveformVector;
    double threshold;
    
public:
    GATStartTimeFinder(size_t nBaselineSamples = 200, size_t iBaselineStartSample = 0, size_t iThresholdSigma = 3, size_t iSmoothSize = 24, size_t iPolyFitOrder = 5);
    //virtual ~GATStartTimeFinder() {}
    
public:
    void Process(int,double);
    
protected:
    MGWFBaselineRemover fBaselineRemover;
    MGWFSavitzkyGolaySmoother fSavitzkyGolaySmoother;
    MGWFExtremumFinder fExtremumFinder;
    MGWFTimePointCalculator fTimePointCalculator;
    GATStartTimeInfo fStartTimeInfo;
    
    //ClassDef(GATStartTimeFinder,1); //ROOT standalone GATStartTimeFinder
};

GATStartTimeFinder::GATStartTimeFinder(size_t nBaselineSamples, size_t iBaselineStartSample, size_t iThresholdSigma, size_t iSmoothSize, size_t iPolyFitOrder)
{
    fBaselineRemover.SetBaselineSamples(nBaselineSamples);
    fBaselineRemover.SetStartSample(iBaselineStartSample);
    fSavitzkyGolaySmoother.ResetSmootherAttributes(iSmoothSize,0,iPolyFitOrder);
    blSamps = nBaselineSamples;
    blStartSamp = iBaselineStartSample;
    thresholdSigma = iThresholdSigma;
    cout << "blSamps " << blSamps << " blStartSamp " << blStartSamp << endl;
}

void GATStartTimeFinder::Process(int entryIndex, double idealSigma)
{
    ////////////////////////////////////////////////////////
    ////// Make sim'd WF  //////////////////////////////////
    ////////////////////////////////////////////////////////
    
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
    cout << "entry index " << entryIndex << endl;
    t1->GetEntry(entryIndex);
    cout << "HERE" << endl;
    branchVectorSize = waveform->size();
    cout << "waveform size " << branchVectorSize << endl;
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
    
    ////////////////////////////////////////////////////////
    ////// Apply transforms to sim'd WF  ///////////////////
    ////////////////////////////////////////////////////////
    
    fStartTimeInfo.Resize(1);
    fBaselineRemover.TransformInPlace(*simWF);
    fSavitzkyGolaySmoother.TransformInPlace(*simWF);
    fExtremumFinder.SetFindMaximum(true);
    fExtremumFinder.SetLocalMaximumTime(20000); // argument in ns (2e4 ns WF trace)
    fExtremumFinder.Transform(simWF);
    extrmVl = fExtremumFinder.GetTheExtremumValue();
    cout << "extrmVl " << extrmVl << endl;
    waveformVector = simWF->GetVectorData();
    threshold = GetMean(waveformVector,blSamps,blStartSamp) + thresholdSigma*GetStdDev(waveformVector,blSamps,blStartSamp);
    if (threshold<=0)
    {
        cout << "WARNING: threshold = " << threshold << " < 0 ... setting threshold to 0" << endl;
        cout << "This is b/c of a suspected MGWFTimePointCalculator issue, in which its methods do not like finding 0%-of-max timepoints" << endl;
        threshold = 1e-20;
    }
    cout << "mean " << GetMean(waveformVector,blSamps,blStartSamp) << " stdDev " << GetStdDev(waveformVector,blSamps,blStartSamp) << " threshold " << threshold << endl;

    fTimePointCalculator.SetLocalMaximumTime(20000); // argument in ns (2e4 ns WF trace)
    fTimePointCalculator.Transform(simWF);
    fTimePointCalculator.AddPoint(threshold/extrmVl);
    fTimePointCalculator.FindTimePoints(*simWF);
    //fTimePointCalculator.GetFromMaxRiseTime(0);
    fStartTimeInfo.SetStartTime(1,fTimePointCalculator.GetFromMaxRiseTime(0));
    waveformVector.resize(0);
    
    cout << "StartTime " << fStartTimeInfo.GetStartTime(0) << endl;
    //return fStartTimeInfo.GetStartTime(1);
}

//ClassImp(GATStartTimeInfo);
//ClassImp(GATStartTimeFinder);

////////////////////////////////////////////////////////
////// Implement member functions from classes  ////////
////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    int entryIndex = atoi(argv[1]);
    double idealSigma = atof(argv[2]);

    TApplication *App = new TApplication("App", 0, NULL);

    GATStartTimeFinder* finder = new GATStartTimeFinder(650, 50, 3, 24, 5);
    finder->Process(entryIndex,idealSigma);
    
    App->Run();
    return 0;
}
