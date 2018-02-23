///////////////////////////////
// reference for writing new tree with branches: GeneralAnalysis_SSTC_FriendTrees_Working2-23-15.cc
///////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>
#include <vector>
#include <complex>

#include <MGTEvent.hh>
#include <MGTWaveform.hh>
#include <MGWFBaselineRemover.hh>
#include <MGWFExtremumFinder.hh>
#include <MGWFFastFourierTransformDefault.hh>
#include <MGWaveformFT.hh>
#include <MGTWaveformRegion.hh>

using namespace std;

int main ()
{
    TApplication *App = new TApplication("App", 0, NULL);

    // INPUT FILE AND INPUT TREE
    TFile* f = new TFile("/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3KJR/OR_run17605.root","READ");
    TTree* t = (TTree*)f->Get("MGTree");
    
    // BRANCHES
    MGTEvent* event = 0; //new MGTEvent();
    int fEventNumber = 0;
    t->SetBranchAddress("event",&event);
    t->SetBranchAddress("fEventNumber",&fEventNumber);
    
    // FLOW VARS
    int nentries = 0, nwaveforms = 0;
    int channel = 0;
    int wfCounter = 0, wfCounterLimit = 100;
    
    // TRANSFORM VARS
    int blStatsLower = 50;
    int blStatsUpper = 700;
    double maxVl = 0., minVl = 0.;
    int maxPt = 0, minPt = 0;

    // PLOTS
    TCanvas* c = new TCanvas();
    TH1D* h[wfCounterLimit];
    for (int i=0; i<wfCounterLimit; i++) h[i] = new TH1D();
    TCanvas* cFT = new TCanvas();
    cFT->Divide(1,2);
    TGraph* gFT;
    
    // READ TREE
    nentries=t->GetEntries();
    cout<<"t nentries="<<nentries<<endl;
    cout << "looping over events ..." << endl;
    for(int i=0;i<nentries;i++)
    {
        t->LoadTree(i);
        t->GetEntry(i);
        nwaveforms = event->GetNWaveforms();
        for(int i_digit=0; i_digit<nwaveforms; i_digit++)
        {
            // READ WF
            MGTWaveform* rawWF = event->GetWaveform(i_digit); // MGTEvent.hh
            // REMOVE BASELINE
            MGWFBaselineRemover* BLR = new MGWFBaselineRemover();
            BLR->SetStartSample(blStatsLower);
            BLR->SetBaselineSamples(blStatsUpper-blStatsLower+1);
            BLR->TransformInPlace(*rawWF);
            // FIND EXTREMA
            MGWFExtremumFinder* EXF = new MGWFExtremumFinder();
            EXF->SetFindMaximum(true);
            EXF->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
            EXF->Transform(rawWF);
            maxVl = EXF->GetTheExtremumValue();
            maxPt = EXF->GetTheExtremumPoint();
            EXF->SetFindMinimum(true);
            EXF->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
            EXF->Transform(rawWF);
            minVl = EXF->GetTheExtremumValue();
            
            // WF CUTS & DRAW
            channel = rawWF->GetID(); // GetDigitizerData
            if(wfCounter < wfCounterLimit &&
               channel%2==0 &&
               maxVl > 50 && maxVl < 1000 &&
               minVl > -20)
            {
                // DRAW WFs ON SAME PLOT
                //cout<<"Here wfCounter "<<wfCounter<<" "<<maxPt<<endl;
                rawWF->LoadIntoHist(h[wfCounter]); // GET RAW WF HIST
                h[wfCounter]->SetLineColor(wfCounter);
                c->cd();
                if(wfCounter==0){
                    h[wfCounter]->Draw();
                    h[wfCounter]->SetTitle("OR_run17605");
                }
                if(wfCounter>0){
                    h[wfCounter]->Draw("SAME");
                }
                if(wfCounter==wfCounterLimit-1){
                    f->Close();
                    App->Run();
                }
                
                // CALL OUT ONE WF AND DO FFT
                if(wfCounter == 50)
                {
                    // PLOT WF
                    cFT->cd(1);
                    h[wfCounter]->Draw();
                    h[wfCounter]->SetTitle(Form("OR_run17605 Entry %d WF %d",i,i_digit));
                    
                    // INVOKE FFT TRANSFORM
                    MGWFFastFourierTransformDefault* FT = new MGWFFastFourierTransformDefault();
                    // INVOKE OBJECT TO HOLD FFT RESULT
                    MGWaveformFT* wfFT = new MGWaveformFT();
                    // PREP WF
                    MGWaveform* rawWF_MG = (MGWaveform*)rawWF;
                    cout<<"   rawWF_MG length (samples):"<<rawWF_MG->GetLength()<<endl;
                    // SET WF REGION TO STUDY
                    MGTWaveformRegion* wfRegion_MGT = new MGTWaveformRegion(0,900);
                    MGWaveformRegion* wfRegion_MG = (MGTWaveformRegion*)wfRegion_MGT;
                    // PERFOM FFT
                    FT->PerformFFT(rawWF_MG, wfFT, wfRegion_MG);
                    // PULL OUT THE DATA AND LOOK AT IT
                    vector<complex<double> > wfFT_vec = wfFT->GetVectorData();
                    vector<double> wfFT_pwr_vec;
                    vector<double> wfFT_f_vec;
                    double samplingPeriod = wfFT->GetSamplingPeriod();
                    int wfRegion_beg = wfRegion_MG->GetBeginning();
                    int wfRegion_end = wfRegion_MG->GetEnd();
                    int wfRegion_samples = wfRegion_end - wfRegion_beg; // do not include the "end" ... MGDO/MGWFFastFourierTransformDefault.cc:53 : "for(size_t i=region->GetBeginning(); i<region->GetEnd(); ++i)"
                    cout<<"    wfFT_vec.size() = "<<wfFT_vec.size()<<endl;
                    cout<<"    samplingPeriod = "<<samplingPeriod<<endl;
                    cout<<"    wfRegion_beg = "<<wfRegion_beg<<endl;
                    cout<<"    wfRegion_end = "<<wfRegion_end<<endl;
                    cout<<"    wfRegion_samples = "<<wfRegion_samples<<endl;
                    double samplingPeriod_s = samplingPeriod*1e-9; // convert from ns to s
                    
                    // PLOT FFT RESULT
                    double f_n = 0.; // f_n = n/(N*dT)
                    double f_n_MHz = 0.;
                    for(unsigned int vec_i = 0; vec_i < wfFT_vec.size(); vec_i++)
                    {
                        complex<double> z = complex<double>(wfFT_vec[vec_i].real(), wfFT_vec[vec_i].imag());
                        //cout<<"    "<<wfFT_vec[vec_i]<<", XX* = "<< z * conj(z) <<endl;
                        //cout<<"    "<<wfFT_vec[vec_i]<<" ("<<wfFT_vec[vec_i].real()<<","<<wfFT_vec[vec_i].imag()<<")"<<endl;
                        wfFT_pwr_vec.push_back((z * conj(z)).real());
                        f_n = vec_i/(wfRegion_samples*samplingPeriod_s);
                        f_n_MHz = f_n/1e6; // convert from Hz to MHz
                        wfFT_f_vec.push_back(f_n_MHz);
                    }
                    TGraph* gFT = new TGraph(wfFT_pwr_vec.size(), &wfFT_f_vec[0], &wfFT_pwr_vec[0]);
                    gFT->GetHistogram()->GetXaxis()->SetTitle("frequency (MHz)");
                    gFT->GetHistogram()->GetYaxis()->SetTitle("power (zz*, not normalized)");
                    gFT->SetTitle(Form("OR_run17605 Entry %d WF %d",i,i_digit));
                    //gFT->SetTitle("graph title;x title;y title");
                    cFT->cd(2);
                    gFT->Draw();
                    cFT->Update();
                    cout<<"    Power spectrum drawn"<<endl;
                    
                }
                wfCounter++;
            } // end cuts
        } // end wfs
    } // end event

    f->Close();
    
    cout << "F I N I S H E D" << endl;
    App->Run();
}