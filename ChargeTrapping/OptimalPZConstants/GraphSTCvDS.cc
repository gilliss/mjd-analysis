///////////////////////////////
//
// 60us to 250 us in steps fo 10
///////////////////////////////

#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TMultiGraph.h>
#include <TList.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TVectorD.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TApplication.h>

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>

using namespace std;

inline int GetTauIndex(int tau) {
    return (tau/10)-6; // this will give 0 for 60, 1 for 70, ..., 19 for 250
}
inline double ConvertToQTTau(int tPZ) {
    double tRC = 72.6; // us
    double tQT = (tRC*tPZ)/(tPZ-tRC); // tPZ = (tRC*tQT)/(tQT-tRC)
    return tQT;
}

int main ()
{
    //////////////////////////////
    // Title, Stats, Date off by default
    gStyle->SetOptTitle(0);
    //gStyle->SetTitleSize(133, "TITLE");
    gStyle->SetTitleOffset(0.1);
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);
    
    // Axis titles
    Int_t kAxisTitleFont = 133; // 13 (Times New Roman) + 3 (size in pixels)
    Float_t kAxisTitleSize = 30;
    Float_t kAxisTitleOffset = 1.0;
    gStyle->SetTitleSize(kAxisTitleSize, "XYZ");
    gStyle->SetTitleFont(kAxisTitleFont, "XYZ");
    gStyle->SetTitleXOffset(kAxisTitleOffset);
    gStyle->SetTitleYOffset(kAxisTitleOffset);
    
    // Axis labels
    Int_t kLabelFont = 133;
    Float_t kLabelSize = 25;
    Float_t kLabelOffset = 0.006;
    gStyle->SetLabelFont(kLabelFont, "XYZ");
    gStyle->SetLabelSize(kLabelSize, "XYZ");
    gStyle->SetLabelOffset(kLabelOffset, "XY");
    gStyle->SetLabelOffset(kLabelOffset*0.5, "Z");
    
    // Other text (e.g. legends, stats)
    Int_t kTextFont = 133;
    Float_t kTextSize = 25;
    gStyle->SetTextSize(kTextSize);
    gStyle->SetTextFont(kTextFont);
    gStyle->SetTextColor(1);
    
    // Fill solid by default
    Int_t kFillStyle=1001;
    
    // No little lines at the ends of error bars.
    gStyle->SetEndErrorSize(0);
    
    // Canvas width and height: 600x800
    gStyle->SetCanvasDefH(600);
    gStyle->SetCanvasDefW(800);
    gStyle->SetCanvasBorderMode(0);
    
    // Pads and margins: use same margin on top / bottom and left / right for easy
    // alignment. User can crop if desired.
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadTopMargin(0.12);
    gStyle->SetPadBottomMargin(0.15); // 0.12
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);
    
    //Frame
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameLineWidth(2);
    
    // Legend: borderless with white background.
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    
    // Colors: use rainbow scale
    gStyle->SetPalette(1);
    //gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
    /////////////////////////////////////////
    
    TApplication *App = new TApplication("App", 0, NULL);
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Command line arguments and initializations
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    TFile *f = new TFile("STCvDS_tauQT_graphs.root","READ");
    
    vector<string> binLabelVec;
    binLabelVec.push_back("B8470");
    binLabelVec.push_back("P42665B");
    binLabelVec.push_back("B8487");
    binLabelVec.push_back("P42538A");
    binLabelVec.push_back("P42575A");
    binLabelVec.push_back("P42661C");
    binLabelVec.push_back("P42698B");
    binLabelVec.push_back("P42574A");
    binLabelVec.push_back("P42661B");
    binLabelVec.push_back("P42573B");
    binLabelVec.push_back("P42661A");
    binLabelVec.push_back("P42573A");
    binLabelVec.push_back("P42698A");
    binLabelVec.push_back("P42538B");
    
    int dsList[4] = {0,1,3,4};
    //int modList[4] = {1,1,1,2};
    int colorList[4] = {4,8,2,6};
    int markerList[4] = {20,21,22,23};
    
    TMultiGraph *mg = new TMultiGraph();
    TLegend *legend=new TLegend(0.6,0.62,0.75,0.85);
    legend->SetTextFont(72);
    legend->SetTextSize(0.04);
    
    TGraphErrors *gSTC = (TGraphErrors*)f->Get("STC_Detectors_tauQT");
    gSTC->SetMarkerColor(1);
    gSTC->SetMarkerStyle(24);
    mg->Add(gSTC);
    legend->AddEntry(gSTC,"STC","p");
    
    for(int i=0; i<4; i++)
    {
        TGraphErrors *gDS = (TGraphErrors*)f->Get(Form("DS%d_Detectors_tauQT",dsList[i]));
        gDS->SetMarkerColor(colorList[i]);
        gDS->SetMarkerStyle(markerList[i]);
        
        mg->Add(gDS);
        
        legend->AddEntry(gDS,Form("DS%d",dsList[i]),"p");
    }
    
    TCanvas *canv = new TCanvas();
    canv->SetBottomMargin(0.25);
    canv->cd();
    mg->Draw("AP");
    mg->SetTitle("Inferred #tau_{QT}, STC vs Modules");
    mg->GetYaxis()->SetTitle("mean free drift time #tau_{QT} (#mu s)");
    //mg->GetYaxis()->SetRangeUser(-0.01,0.2);
    mg->GetYaxis()->SetRangeUser(120,1040);
    //mg->GetXaxis()->SetRangeUser(0.0,18.0);
    TAxis* a = mg->GetXaxis();
    int j = 1;
    int bin_j;
    while (j<((a->GetXmax()))) // may need (a->GetXmax())-1
    {
        bin_j = a->FindBin(j);
        a->SetBinLabel(bin_j,binLabelVec[j-1].c_str());
        j++;
    }
    //TAxis* a2 = mg->GetYaxis();
    //a2->SetBinLabel(a2->FindBin(2000),"-"); // ~inf mean free drift time
    
    TLine* line[20];
    int PZi=0;
    for(int PZ = 80; PZ<=150; PZ+=10)
    {
        PZi = GetTauIndex(PZ);
        line[PZi] = new TLine(a->GetXmin(),ConvertToQTTau(PZ),a->GetXmax(),ConvertToQTTau(PZ));
        line[PZi]->SetLineColorAlpha(17, 1.0);
        line[PZi]->Draw("SAME");
    }

    legend->Draw();
    canv->Modified();
    canv->Update();

    cout << "F I N I S H E D" << endl;
    App->Run();
}


/*
 binLabelVec.push_back("B8455");
 binLabelVec.push_back("B8470");
 binLabelVec.push_back("B8465");
 binLabelVec.push_back("P42665B");
 binLabelVec.push_back("B8487");
 binLabelVec.push_back("P42538A");
 binLabelVec.push_back("P42575A");
 binLabelVec.push_back("P42661C");
 binLabelVec.push_back("P42698B");
 binLabelVec.push_back("P42662A");
 binLabelVec.push_back("P42574A");
 binLabelVec.push_back("P42661B");
 binLabelVec.push_back("P42573B");
 binLabelVec.push_back("P42661A");
 binLabelVec.push_back("P42573A");
 binLabelVec.push_back("P42698A");
 binLabelVec.push_back("P42538B");
 */


//for(int p = 1; p <= maxPos; p++)
//{
//    for(int d = 1; d <= maxDet; d++)
//    {
//        
//    } // d
//} // p
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// ...
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////