#include "MGStartTimeFinder.hh"

#include <MGWaveform.hh> //"/Users/tomgilliss/Dev/mjswDev/MGDO/Root/MGTWaveform.hh"
#include "MGWFBaselineRemover.hh"
#include "MGWFSavitzkyGolaySmoother.hh"
#include "MGWFExtremumFinder.hh"
#include "MGWFTimePointCalculator.hh"

using namespace std;

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

MGStartTimeFinder::MGStartTimeFinder(const string& aName) :
  MGVWaveformParameterCalculator(aName)
{
  AddParameter("wfStartTime");
}

void MGStartTimeFinder::CalculateParameters(const MGWaveform& wf)
{
    //cout << " wf length " << wf.GetLength() << endl;
    if(wf.GetLength() < 1999) {
        MGDOerr << "Waveform (wf) has length < 1999" << std::endl;
        return;
    }
    int nBaselineSamples = 200;
    int iBaselineStartSample = 0;
    double extrmVl = 0;
    vector<double> waveformVec;
    waveformVec.resize(0);
    double mean = 0;
    double stdDev = 0;
    double threshold = 0;
    
    MGWFBaselineRemover* BLR = new MGWFBaselineRemover();
    BLR->SetBaselineSamples(nBaselineSamples);
    BLR->SetStartSample(iBaselineStartSample);
    MGWaveform blrWF = (MGWaveform) wf;  //new MGWaveform();
    BLR->TransformInPlace(blrWF);
    //cout << " blrWF length " << blrWF.GetLength() << endl;
    if(blrWF.GetLength() < 1999) {
        MGDOerr << "Waveform (blrWF) has length < 1999" << " : " << blrWF.GetLength() << std::endl;
        return;
    }
    
    MGWFSavitzkyGolaySmoother* SGS = new MGWFSavitzkyGolaySmoother(24,0,5,"MGWFSavitzkyGolaySmoother");
    MGWaveform* smoothWF = new MGWaveform();
    SGS->TransformOutOfPlace(blrWF,*smoothWF);
    //cout << " smoothWF length " << smoothWF->GetLength() << endl;
    if(smoothWF->GetLength() < 1999) {
        MGDOerr << "Waveform (smoothWF) has length < 1999" << std::endl;
        return;
    }
    
    MGWFExtremumFinder* EXF = new MGWFExtremumFinder();
    EXF->SetFindMaximum(true);
    EXF->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
    EXF->Transform(smoothWF);
    extrmVl = EXF->GetTheExtremumValue();
    //cout << "extrmVl " << extrmVl << endl;
    
    waveformVec = smoothWF->GetVectorData();
    mean = GetMean(waveformVec,nBaselineSamples,iBaselineStartSample);
    stdDev = GetStdDev(waveformVec,nBaselineSamples,iBaselineStartSample);
    threshold = mean+3*stdDev;
    //cout << "threshold " << threshold << endl;
    if (threshold/extrmVl>=1)
    {
        cout << "WARNING (MGStartTimeFinder): threshold/extrmVl " << threshold/extrmVl << " > 100%" << endl;
        cout << "... setting this WF's t0 estimate to 0.0 ns" << endl;
        SetParameterValue(0,0.0);
    }
    else
    {
        MGWFTimePointCalculator* TPC = new MGWFTimePointCalculator();
        TPC->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
        TPC->Transform(smoothWF);
        TPC->AddPoint(threshold/extrmVl);
        TPC->FindTimePoints(*smoothWF);
        
        SetParameterValue(0, TPC->GetFromMaxRiseTime(0));
        //cout << "t0 " << TPC->GetFromMaxRiseTime(0) << endl;
        delete TPC;
    }
    
    delete BLR;
    delete SGS;
    delete smoothWF;
    delete EXF;
}
