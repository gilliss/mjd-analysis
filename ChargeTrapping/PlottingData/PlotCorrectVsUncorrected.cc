/*
[position] [channel] [energy parameter/decay constant] [mu] [mu unc] [sigma] [sigma unc] [fwhm] [fwhm unc] [ratio] [ratio unc]
 
http://stackoverflow.com/questions/7868936/read-file-line-by-line
https://root.cern.ch/root/html/tutorials/tree/basic.C.html
http://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
*/

#include <iostream>
#include <fstream>
#include <algorithm> // min_element, sort, binary_search
#include <vector>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMultiGraph.h>

using namespace std;

inline int GetTauIndex(int tau) {
    return (tau/10)-6; // this will give 0 for 60, 1 for 70, ..., 19 for 250
}
inline int GetTauFromIndex(int i) {
    return (i+6)*10; // this will give 0 for 60, 1 for 70, ..., 19 for 250
}
inline double GetkeV(double fwhm, double muADC, double keV) {return fwhm*(keV/muADC);}

int main (int argc, char* argv[])
{
///////////////////////////////
//// INITIALIZATIONS
///////////////////////////////
    TApplication *App = new TApplication("App", 0, NULL);

    if (argc < 2)
    {
        cout << "too few arguments " << argv[0] << endl;
        return 1;
    }
    int DS = atoi(argv[1]);
    
    // input ROOT file
    char fileName[50];
    sprintf(fileName,"TauvsChan_DS01234.root");
    TFile *f = new TFile(fileName,"READ");
    
    TH1D* hECal = (TH1D*)f->Get("FWHMHist_DS3_Uncorrected");
    hECal->SetLineWidth(2);
    TH1D* hENFCal = (TH1D*)f->Get("FWHMHist_DS3");
    hENFCal->SetLineWidth(2);
    
    TCanvas *c = new TCanvas();
    c->cd();
    hECal->Draw();
    hENFCal->Draw("SAME");
    
//    TLegend *l = new TLegend(0.6,0.62,0.75,0.85);
//    l->AddEntry(hECal,"Before Correction","l");
//    l->AddEntry(hENFCal,"After Correction","l");
//    l->Draw();
//    c->Update();
    
    TCanvas *c1 = new TCanvas();
    c1->cd();
    TGraphErrors *grECal = (TGraphErrors*)f->Get("FWHMvsChan_Uncorrected_DS3");
    TGraphErrors *grENFCal = (TGraphErrors*)f->Get("FWHMvsChan_DS3");
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(grECal);
    mg->Add(grENFCal);
    mg->SetTitle("FWHM, Each Detector at 2614 keV");
    mg->Draw("AP");
    mg->GetXaxis()->SetLabelSize(0.0);
    mg->GetXaxis()->SetTitle("Detector");
    mg->GetYaxis()->SetTitle("FWHM (keV)");
    c1->Modified();
    

    cout << "F I N I S H E D" << endl;
    App->Run();
}