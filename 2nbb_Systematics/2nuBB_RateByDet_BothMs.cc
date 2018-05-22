///////////////////////////////
//
// see ResultsFromText.cc for useful plotting reference
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
#include <TLatex.h>

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>

#include <MJTChannelMap.hh>

using namespace std;

int main (int argc, char* argv[])
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Command line arguments and initializations
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 2)
    {
        cout << "Error: Too few arguments" << argv[0] << endl;
        cout << "Usage is: ./2nuBB_RateByDet_BothMs DS enrType save" << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // DS#
    // int c = atoi(argv[2]); // M#
    int enrType = atoi(argv[2]); // Nat(0) or Enr(2)
    if(enrType!=0 && enrType!=2) {cout<<"Error:  Nat(0) or Enr(2)"<<endl; return 1;}
    int save = atoi(argv[3]);
    string dataType = "open";

    TApplication *App;
    if(save == 0) {App = new TApplication("App", 0, NULL);}

    // PRINT SETTINGS
    cout<<"------------"<<endl;
    cout<<"Settings:"<<endl;
    cout<<" Plotting DS" << DS << endl;
    cout<<" dataType: " << dataType << endl;
    cout<<"------------"<<endl;

    // INPUT MJTCHANNELMAP FILE
    TFile *mapFile = new TFile("channel_map_data/DSChannelMaps.root","READ");

    // INPUT DATA FROM TEXT FILE
    char inTextFilePath[50];
    int m, p, d;
    /////int dInt;
    /////double dMkg, dLT;
    int dInt; // integral
    double dIntN, dEX, dEXU; // normed integral, exposure, exposure uncert

    // INITIALIZATIONS
    int maxMod = 2, maxPos = 7, maxDet = 5;
    char canvName[50];
    char graphName[50];
    int cutColor[4]={1,16,92,59}; // to match Lukas
    double norm_dInt = 0;
    double enorm_dInt = 0;
    double scale = 0;
    int i = 0;
    map <string, double> detCoord; // <cpd,xCoordForGraph>
    vector <string> cpdNameVec;
    string cpdName; //char cpdName[50];

    vector<double> avgnormrate;
    vector<double> avgCounts;

    vector<double> groupingFill;
    vector<string> groupingName;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Plotting and saving
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // PLOTTING AND SAVING
    sprintf(canvName,"c_DS%d_DetType%d",DS,enrType);
    sprintf(graphName,"g_DS%d_DetType%d",DS,enrType);
    TCanvas *canv = new TCanvas(canvName,canvName);
    TGraphErrors *graph[4];

    for(int cutScheme = 0; cutScheme <= 3; cutScheme++)
    {
        int chMapDS = DS; if(DS==51 || DS==52 || DS==53){chMapDS=5;}
        MJTChannelMap *chMap = (MJTChannelMap*) mapFile->Get(Form("ChannelMapDS%d",chMapDS));
        if(chMap==NULL) cout << "null chMap"<<endl;
        graph[cutScheme] = new TGraphErrors();
        graph[cutScheme]->SetTitle(Form("DS%d_DetType%d_CutScheme%d",DS,enrType,cutScheme));
        graph[cutScheme]->SetMarkerColor(cutColor[cutScheme]);
        graph[cutScheme]->SetMarkerStyle(20);
        i = 0;
        for(int c = 1; c <= 2; c++)
        {
          sprintf(inTextFilePath,"output_data/%s/DS%d_M%d_CutScheme%d.txt",dataType.c_str(),DS,c,cutScheme);
          ifstream inTextFile(inTextFilePath); // scope of loop makes it OK to redefine inTextFile
          while(inTextFile >> m >> p >> d >> dInt >> dIntN >> dEX >> dEXU)
          {
              if(chMap->GetInt(m,p,d,"kDetectorType")==enrType)
              {
                  //cout << "cpd "<<m<<p<<d<<" detType "<<chMap->GetInt(m,p,d,"kDetectorType")<< endl;
                  cpdName = "C" + to_string(m) + "P" + to_string(p) + "D" + to_string(d);
                  if(cutScheme == 0)
                  {
                      detCoord.insert(pair<string,double>(cpdName,double(i+1)));
                      cpdNameVec.push_back(cpdName);
                  }
                  scale = 1/(dEX);
                  norm_dInt = dInt*scale;
                  // enorm_dInt = sqrt(dInt)*scale; // could also fold in error on exposure dEXU
                  enorm_dInt = sqrt( ( (1/dEX) * sqrt(dInt) )*( (1/dEX) * sqrt(dInt) ) + ( (-dInt/(dEX*dEX)) * dEXU )*( (-dInt/(dEX*dEX)) * dEXU ) );
                  graph[cutScheme]->SetPoint(i,detCoord.find(cpdName)->second,norm_dInt);
                  if(cutScheme==3){graph[cutScheme]->SetPointError(i,0.0,enorm_dInt);}
                  i++;

                  if(cutScheme==3){
                      avgnormrate.push_back(norm_dInt);
                      avgCounts.push_back(dInt);
                      groupingFill.push_back(norm_dInt);
                      groupingName.push_back(cpdName);
                  }
              }
          }
        }
        delete chMap;
        chMap = 0;
    } // cutScheme

    // MULTIGRAPH
    canv->cd();
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(graph[0]);
    mg->Add(graph[1]);
    mg->Add(graph[2]);
    mg->Add(graph[3]);
    if(enrType==0)mg->SetTitle(Form("Rate 1000-1400keV, DS%d Natural Detectors",DS));
    if(enrType==2)mg->SetTitle(Form("Rate 1000-1400keV, DS%d Enriched Detectors",DS));
    mg->Draw("AP");
    mg->GetYaxis()->SetTitle("Cnt/kg/dy");
    if(enrType==0)mg->GetYaxis()->SetRangeUser(0.0,1.7);
    if(enrType==2)mg->GetYaxis()->SetRangeUser(0.0,3.1); // if(enrType==2)mg->GetYaxis()->SetRangeUser(0.0,2.7);
    TAxis* a = mg->GetXaxis();
    int j = 1;
    int bin_j;
    cout<<"groupingFill size "<<groupingFill.size()<<" a->GetXmax() "<< a->GetXmax()<<" int(a->GetXmax()) "<< int(a->GetXmax()) << endl;
    while (j<(int(a->GetXmax()))) // may need (a->GetXmax())-1
    {
        bin_j = a->FindBin(j);
        a->SetBinLabel(bin_j,cpdNameVec[j-1].c_str());
        j++;
    }
    canv->Modified();
   // TLegend *legend=new TLegend(0.6,0.62,0.75,0.85);
   // legend->SetTextFont(72);
   // legend->SetTextSize(0.04);
   // legend->AddEntry(graph[0],"Raw","p");
   // legend->AddEntry(graph[1],"Raw+DC","p");
   // legend->AddEntry(graph[2],"Raw+DC+AE","p");
   // legend->AddEntry(graph[3],"Raw+DC+AE+DCR","p");
   // legend->Draw();
    canv->Update();

    // AVERAGE RATE
    double avgrate = 0.;
    double avgCnts = 0.;
    cout << "size avgnormrate and avgCounts vec "<<avgnormrate.size()<<" "<<avgCounts.size()<<endl;
    for(unsigned int k = 0; k < avgnormrate.size(); k++)
    {
        cout<<"  k "<<k<<" "<<avgnormrate[k]<<" "<<avgCounts[k]<<endl;
        avgrate+=avgnormrate[k];
        avgCnts+=avgCounts[k];
    }
    cout <<"cut scheme 3: avg norm rate = " << avgrate/double(avgnormrate.size()) << " avg counts = " << avgCnts/double(avgCounts.size()) << endl;
    avgCnts = avgCnts/double(avgCounts.size());

    // Fit TGraphErrors // https://root.cern.ch/root/html530/TGraph.html#TGraph:Fit%1
    // Unweighted fit
    TF1  *f1 = new TF1("f1","pol0");
    graph[3]->Fit(f1,"W N 0");
    double p0 = f1->GetParameter(0);
    double p0Err = f1->GetParError(0);
    cout << "p0 " << p0 << " +/- " << p0Err << endl;
    // Weighted fit
    TF1  *f1_2 = new TF1("f1","pol0");
    graph[3]->Fit(f1_2,"N 0 EX0");
    double p0_2 = f1_2->GetParameter(0);
    double p0Err_2 = f1_2->GetParError(0);
    cout << "p0 " << p0_2 << " +/- " << p0Err_2 << endl;

    // DRAW FIT LINE ON THE FIRST CANVAS
    canv->cd();
    TLine *line = new TLine(a->GetXmin(),p0,a->GetXmax(),p0); // x1,y1,x2,y2
    line->SetLineColorAlpha(cutColor[3], 0.45);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw("same");
    canv->Update();

    canv->cd();
    TLine *line_2 = new TLine(a->GetXmin(),p0_2,a->GetXmax(),p0_2); // x1,y1,x2,y2
    line_2->SetLineColorAlpha(cutColor[3], 0.45);
    line_2->SetLineWidth(2);
    line_2->SetLineStyle(9);
    line_2->Draw("same");
    canv->Update();

    // DRAW TEXT ON FIRST CANVAS
    canv->cd();
    TLatex tlatex;
    tlatex.SetTextSize(13);
    tlatex.DrawLatex(0.78*groupingFill.size(), 2.9, Form("Avg Counts: %.2f", avgCnts));
    tlatex.DrawLatex(0.78*groupingFill.size(), 2.75, Form("Unweighted Fit: %.2f +/- %.2f", p0, p0Err));
    tlatex.DrawLatex(0.78*groupingFill.size(), 2.6, Form("Weighted Fit: %.2f +/- %.2f", p0_2, p0Err_2));
    canv->Update();

    // CLOSEOUT
    mapFile->Close();
    delete mapFile;
    mapFile = 0;

    if(save==1) {canv->Print(Form("DS%d_DetType%d_Det_Fit.pdf",DS,enrType));}

    cout << "F I N I S H E D" << endl;
    if(save == 0) {App->Run();}
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

// HIST FOR FIT
// double binLo = 1.0-0.5;
// double binHi = groupingFill.size()+0.5;
// int nBins = binHi - binLo;
// sprintf(canvName,"c2_DS%d_DetType%d",DS,enrType);
// sprintf(graphName,"h_DS%d_DetType%d",DS,enrType);
// TCanvas *canv2 = new TCanvas(canvName,canvName);
// TH1D *h = new TH1D(graphName,graphName,nBins,binLo,binHi);
// h->SetMarkerStyle(20);
// h->SetMarkerColor(cutColor[3]);
// h->Sumw2(); // the sum of squares of weights is also stored
// TAxis* a2 = h->GetXaxis();
// for(unsigned int j = 0; j<groupingFill.size(); j++)
// {
//     h->Fill(j+1,groupingFill[j]);
//     a2->SetBinLabel(j+1,groupingName[j].c_str());
// }
// canv2->cd();
// h->Draw("P HIST");
// TF1  *f1 = new TF1("f1","pol0");
// h->Fit(f1,"WL"); //h->Fit("pol0", "WL"); //“WL” Weighted log likelihood method. To be used when the histogram has been filled with weights different than 1.
// f1->Draw("SAME");
// double p0 = f1->GetParameter(0);
// double p0Err = f1->GetParError(0);
// cout << "p0 " << p0 << " +/- " << p0Err << endl;
// canv2->Update();
// delete canv2; canv2 = NULL;
