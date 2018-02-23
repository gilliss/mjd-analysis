... other includes

#include <vector>
#include <complex>

#include <MGWFFastFourierTransformDefault.hh> // these guys are in MGDO/base and MGDO/Root
#include <MGWaveformFT.hh>
#include <MGTWaveformRegion.hh>

...

looping through entries and waveforms
{
    ... pull out your favorite waveform
    
    ... begin working on a fourier transform of the baseline from the start of waveform:
    
    // INVOKE FFT TRANSFORM
    MGWFFastFourierTransformDefault* FT = new MGWFFastFourierTransformDefault();
    // INVOKE OBJECT TO HOLD FFT RESULT
    MGWaveformFT* wfFT = new MGWaveformFT();
    // PREP WF
    MGWaveform* rawWF_MG = (MGWaveform*)rawWF_MGT; // convert object type (a better soln to working with MG* objects and derived MGT* objects?)
    cout<<"   rawWF_MG length (samples):"<<rawWF_MG->GetLength()<<endl;
    // SET WF REGION TO STUDY
    MGTWaveformRegion* wfRegion_MGT = new MGTWaveformRegion(startSample,endSample);
    MGWaveformRegion* wfRegion_MG = (MGTWaveformRegion*)wfRegion_MGT; // another conversion
    // PERFOM FFT
    FT->PerformFFT(rawWF_MG, wfFT, wfRegion_MG);
    // PULL OUT THE DATA AND LOOK AT IT
    vector<complex<double> > wfFT_vec = wfFT->GetVectorData();
    
    ... loop through FT result data to calculate power vs freq. Can use things like:
    
    complex<double> z = complex<double>(wfFT_vec[vec_i].real(), wfFT_vec[vec_i].imag());
    power = (z * conj(z)).real()
    
    ... can look at the data with things like:

    double samplingPeriod = wfFT->GetSamplingPeriod();
    int wfRegion_beg = wfRegion_MG->GetBeginning();
    int wfRegion_end = wfRegion_MG->GetEnd();
    int wfRegion_samples = wfRegion_end - wfRegion_beg; // see MGDO/MGWFFastFourierTransformDefault.cc:53
    
    ... plot power spectrum as a TGraph where the x and y data are vectors you feed to the TGraph
    
    ... convert the x axis into units of frequency. Should be f_n = n/(N*T). Bin n, # samples N, sampling period T. Mind the two-sided vs one-sided power spectrum concept ...
    
}
