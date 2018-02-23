///////////////////////////////
//
// see ResultsFromText.cc for useful plotting reference
// see 2nuBB systematics plotting scripts for useful reference
// Quick & Dirty Method: http://mjwiki.npl.washington.edu/pub/Majorana/AnalysisReports/dirty_and_quick_upper_limit.pdf
//
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
#include <TApplication.h>

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>

#include <GATDataSet.hh>
#include <MJTChannelMap.hh>
#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201706Prod/GAT/Apps/DataSetInfo.hh"

using namespace std;

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
    // ISOTOPES TO USE
    vector <string> xNameVec;
    string xName;
    map < string , vector<double> > isotopeMap;
    string tempi;
    vector<double> tempe;
    tempi = "228Ac"; tempe.push_back(338.320); tempe.push_back(911.204); tempe.push_back(968.971); tempe.push_back(1588.2);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "108mAg"; tempe.push_back(433.937); tempe.push_back(614.276); tempe.push_back(722.907);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "214Bi"; tempe.push_back(609.312); tempe.push_back(1120.3); tempe.push_back(1764.494); tempe.push_back(2204.1);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "60Co"; tempe.push_back(1173.24); tempe.push_back(1332.5);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "40K"; tempe.push_back(1460.83);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "234mPa"; tempe.push_back(1001.03);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "210Pb"; tempe.push_back(46.539);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "214Pb"; tempe.push_back(351.9321); // progenitor of 214Bi
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "226Ra"; tempe.push_back(186.211);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "208Tl"; tempe.push_back(583.191); tempe.push_back(1593.0); tempe.push_back(2105.0); tempe.push_back(2614.53);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    int truncatedLine;
    for(std::map< string , vector<double> >::iterator map_i=isotopeMap.begin(); map_i!=isotopeMap.end(); map_i++)
    {
        for(unsigned int lin_i = 0; lin_i<map_i->second.size(); lin_i++)
        {
            cout<<map_i->first<<" "<<map_i->second[lin_i]<<endl;
            truncatedLine = map_i->second[lin_i];
            xName = map_i->first + " " + to_string(truncatedLine);
            xNameVec.push_back(xName);
        }
    }
    cout<<"------------"<<endl;
    
    TFile *f = new TFile("DS_MultiLine_Graphs.root","READ");
    
    int nDS = 9;
    int dsList[9] = {0,1,2,3,4,5,5,6,6};
    int modList[9] = {1,1,1,1,2,1,2,1,2};
    int colorList[9] = {28,30,34,58,51,98,91,209,79};//{4,8,93,2,6};
    int markerList[9] = {22,22,22,22,23,22,23,22,23};//{20,21,24,22,23};
    
    TMultiGraph *mg = new TMultiGraph();
    TLegend *legend=new TLegend(0.6,0.62,0.75,0.85);
    legend->SetTextFont(72);
    legend->SetTextSize(0.04);
    
    for(int i=0; i<nDS; i++)
    {
        TGraphErrors *g = (TGraphErrors*)f->Get(Form("DS%d_C%d_RoiCnt",dsList[i],modList[i]));
        //TGraphErrors *g = (TGraphErrors*)f->Get(Form("DS%d_C%d_PkCnt",dsList[i],modList[i]));
        g->SetMarkerColor(colorList[i]);
        g->SetMarkerStyle(markerList[i]);
        
        mg->Add(g);
        
        legend->AddEntry(g,Form("DS%d",dsList[i]),"p");
    }
    
    TCanvas *canv = new TCanvas();
    canv->SetBottomMargin(0.25);
    canv->cd();
    mg->Draw("AP");
    mg->SetTitle("Counts in various ROIs");
    mg->GetYaxis()->SetTitle("Cnt/kg/dy");
    mg->GetYaxis()->SetRangeUser(-0.01,0.2);
    TAxis* a = mg->GetXaxis();
    int j = 1;
    int bin_j;
    while (j<((a->GetXmax()))) // may need (a->GetXmax())-1
    {
        bin_j = a->FindBin(j);
        a->SetBinLabel(bin_j,xNameVec[j-1].c_str());
        j++;
    }
    legend->Draw();
    canv->Modified();
    canv->Update();

    cout << "F I N I S H E D" << endl;
    App->Run();
}
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