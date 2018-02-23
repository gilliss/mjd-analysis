/*
This script reads in waveforms simulated by pysiggen. The waveforms have varied risetime, energy, QT time constants, etc.
Once sim'd waveforms are read in, the script corrects them and makes plots to show how the correction works.
*/

#include <fstream>

#include <TApplication.h>
#include <MGTWaveform.hh>
#include <TCanvas.h>
#include <TH1D.h>
#include <TMath.h>
#include <TLegend.h>

#include <MGWFPoleZeroCorrection.hh>

using namespace std;

double ns = 1e-9;
double us = 1e-6;

int main (int argc, char* argv[])
{
    cout << "Script requires an input argument"<<endl;
    double coeff = atoi(argv[1]); // for qt=0.5 and rc=72.6, coeff = 116 seems OK
    
    TApplication *App = new TApplication("App", 0, NULL);
    
    double rcTimeConstant = 30.0*us/ns/10; // 72.6 us // evaluate to units of 10 ns to match sampling rate
    double qtTimeConstant = 50.0*us/ns/10; // 0.5 us
    double effPZTimeConstant = coeff*us/ns/10; //(rcTimeConstant*qtTimeConstant)/(qtTimeConstant-rcTimeConstant);//coeff*us/ns/10;//-rcTimeConstant+coeff*qtTimeConstant;
    double baseline = 0;
    TH1D* hWF[3][2];
    TCanvas* cWF = new TCanvas();
    cWF->cd();
    MGWFPoleZeroCorrection *PZC = new MGWFPoleZeroCorrection();
    MGTWaveform* simWF = new MGTWaveform();

    int nbins = 200;
    double new_bins[nbins+1];
    for(int i=0; i <= nbins; i++){
        new_bins[i] = i*10;
    }
    
    ////////////////////////////////
    ////////////////////////////////
    //// NO TRAPPING NO RC PULSE, FAST
    ////////////////////////////////
    ////////////////////////////////
    
    ifstream infile0("/global/u2/g/gilliss/ChargeTrapping/ExplorationOfProcedure/EffPZTrapFilter/WaveformFiles/WF_e80_r15_z15_QTrcINF_rcINF_noRing.txt");
    vector<double> simWF_vec;
    double val = 0;
    while (infile0 >> val)
    {
        simWF_vec.push_back(val);
    }
    simWF->SetData(simWF_vec);
    hWF[0][0] = new TH1D();
    simWF->LoadIntoHist(hWF[0][0]);
    cout<<"nbins "<<hWF[0][0]->GetNbinsX()<<endl;
    for(int i=0; i<20; i++){cout<<"bin "<<i<<" "<<hWF[0][0]->GetBinContent(i)<<endl;}
    hWF[0][0]->SetBins(nbins, new_bins); // (nx, xBins) xBins is supposed to be of length nx+1
    cout<<"nbins "<<hWF[0][0]->GetNbinsX()<<endl;
    for(int i=0; i<20; i++){cout<<"bin "<<i<<" "<<hWF[0][0]->GetBinContent(i)<<endl;}
    hWF[0][0]->SetTitle("");
    hWF[0][0]->GetXaxis()->SetTitle("time (ns)");
    hWF[0][0]->GetXaxis()->SetTitleSize(.04);
    hWF[0][0]->GetYaxis()->SetTitle("Charge Signal");
    hWF[0][0]->GetYaxis()->SetTitleOffset(0.5);
    hWF[0][0]->GetYaxis()->SetTitleSize(.04);
    hWF[0][0]->GetYaxis()->SetLabelSize(0.);
    hWF[0][0]->SetLineWidth(2);
    hWF[0][0]->SetLineColor(kBlack);
    hWF[0][0]->SetLineStyle(9);
    hWF[0][0]->SetStats(0);
    hWF[0][0]->Draw("C");
    cWF->Update();
    
    ////////////////////////////////
    ////////////////////////////////
    //// RC & QT PULSE, FAST
    ////////////////////////////////
    ////////////////////////////////
    
    ifstream infile2("/global/u2/g/gilliss/ChargeTrapping/ExplorationOfProcedure/EffPZTrapFilter/WaveformFiles/WF_e80_r15_z15_QTrc50_rc30_noRing.txt");
    simWF_vec.clear();
    val = 0;
    while (infile2 >> val)
    {
        simWF_vec.push_back(val);
    }
    simWF->SetData(simWF_vec);
    hWF[1][0] = new TH1D();
    simWF->LoadIntoHist(hWF[1][0]);
    hWF[1][0]->SetBins(nbins, new_bins);
    hWF[1][0]->SetTitle("");
    hWF[1][0]->GetXaxis()->SetTitle("time (10ns/sample)");
    hWF[1][0]->SetLineWidth(2);
    hWF[1][0]->SetLineColor(kBlue);
    hWF[1][0]->SetStats(0);
    hWF[1][0]->Draw("C SAME");
    cWF->Update();
    
    PZC->SetDecayConstant(effPZTimeConstant);
    PZC->SetRestingBaselineValue(baseline);
    MGTWaveform* WF11 = new MGTWaveform();
    PZC->TransformOutOfPlace(*simWF,*WF11);
    hWF[1][1] = new TH1D();
    WF11->LoadIntoHist(hWF[1][1]);
    hWF[1][1]->SetBins(nbins, new_bins);
    hWF[1][1]->SetLineWidth(2);
    hWF[1][1]->SetLineColor(kPink);
    hWF[1][1]->SetStats(0);
    hWF[1][1]->Draw("C SAME");
    
    ////////////////////////////////
    ////////////////////////////////
    //// RC & QT PULSE, SLOW
    ////////////////////////////////
    ////////////////////////////////
    
    ifstream infile3("/global/u2/g/gilliss/ChargeTrapping/ExplorationOfProcedure/EffPZTrapFilter/WaveformFiles/WF_e80_r30_z30_QTrc50_rc30_noRing.txt");
    simWF_vec.clear();
    val = 0;
    while (infile3 >> val)
    {
        simWF_vec.push_back(val);
    }
    simWF->SetData(simWF_vec);
    hWF[2][0] = new TH1D();
    simWF->LoadIntoHist(hWF[2][0]); // MGTWaveform.hh
    hWF[2][0]->SetBins(nbins, new_bins);
    hWF[2][0]->SetTitle("");
    hWF[2][0]->GetXaxis()->SetTitle("time (10ns/sample)");
    hWF[2][0]->SetLineWidth(2);
    hWF[2][0]->SetLineColor(kBlue);
    hWF[2][0]->SetStats(0);
    hWF[2][0]->Draw("C SAME");
    cWF->Update();
    
    PZC->SetDecayConstant(effPZTimeConstant);
    PZC->SetRestingBaselineValue(baseline);
    MGTWaveform* WF21 = new MGTWaveform();
    PZC->TransformOutOfPlace(*simWF,*WF21);
    hWF[2][1] = new TH1D();
    WF21->LoadIntoHist(hWF[2][1]);
    hWF[2][1]->SetBins(nbins, new_bins);
    hWF[2][1]->SetLineWidth(2);
    hWF[2][1]->SetLineColor(kPink);
    hWF[2][1]->SetStats(0);
    hWF[2][1]->Draw("C SAME");
    
    ////////////////////////////////
    ////////////////////////////////
    //// LEGEND
    ////////////////////////////////
    ////////////////////////////////
    
    TLegend *legend=new TLegend(0.6,0.62,0.8,0.85);
    //    legend->SetTextFont(72);
    legend->SetTextSize(0.03);
    legend->AddEntry(hWF[0][0],"Unmodified","l");
    legend->AddEntry(hWF[1][0],"RC & QT Modified","l");
    legend->AddEntry(hWF[1][1],"RC & QT Corrected","l");
    
    legend->Draw();
    
    App->Run();
}