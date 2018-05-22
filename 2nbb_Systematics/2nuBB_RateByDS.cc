///////////////////////////////
//
// This script takes the individual detector results from 2nuBB_Systematics.cc, defines some way of grouping the detectors, and combines the detector results according to that grouping.
// Examples of groupings are grouping by controller card, mother board, DS, string, and-- in the trivial case-- by individual detector.
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
    //gROOT->ProcessLine(".x MJDTalkPlotStyle.C");

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
/////////////////////////////////////////

    TApplication *App = new TApplication("App", 0, NULL);
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Command line arguments and initializations
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 2)
    {
        cout << "Error: Too few arguments" << argv[0] << endl;
        cout << "Usage is: ./2nuBB_RateByGrouping <enrType> <save>" << endl;
        return 1;
    }
    //int DS = atoi(argv[1]); // DS#
    //int c = atoi(argv[2]); // M#
    int enrType = atoi(argv[1]); // Nat(0) or Enr(2)
    if(enrType!=0 && enrType!=2) {cout<<"Error:  Nat(0) or Enr(2)"<<endl; return 1;}
    int save = atoi(argv[2]);

    // PRINT SETTINGS
    cout<<"------------"<<endl;
    cout<<"Settings:"<<endl;
    cout<<" enrType = " << enrType << endl;
    cout<<" save = " << save << endl;
    cout<<"------------"<<endl;

    // INPUT MJTCHANNELMAP FILE
    TFile *mapFile = new TFile("channel_map_data/DSChannelMaps.root","READ");

    // INPUT DATA FROM TEXT FILE
    char inTextFilePath[50]; // will be filled with relative path to text file of results from 2nuBB_Systematics.cc
    int dInt; // integral counts for detector
    double dIntN, dEX, dEXU; // normed integral, exposure, exposure uncert

    // INITIALIZATIONS
    int m, p, d;
    char canvName[50];
    char graphName[50];
    int cutColor[4]={1,16,92,59}; // to match Lukas
    double norm_dInt = 0;
    double enorm_dInt = 0;
    double scale = 0;
    int i = 0;
    int nDets = 0;
    vector <string> tickLabelVec;
    vector <int> dsVec;
    double groupInt = 0;
    double norm_groupInt = 0;
    double enorm_groupInt = 0;
    double groupMass = 0;
    double groupEX = 0;

    tickLabelVec.push_back("0");
    tickLabelVec.push_back("1");
    tickLabelVec.push_back("2");
    tickLabelVec.push_back("3");
    tickLabelVec.push_back("4");
    // tickLabelVec.push_back("51");
    // tickLabelVec.push_back("52");
    tickLabelVec.push_back("53");
    tickLabelVec.push_back("6");
    dsVec.push_back(0);
    dsVec.push_back(1);
    dsVec.push_back(2);
    dsVec.push_back(3);
    dsVec.push_back(4);
    // dsVec.push_back(51);
    // dsVec.push_back(52);
    dsVec.push_back(53);
    dsVec.push_back(6);

    string tickLabel; //char cpdName[50];
    vector<double> groupingFill;
    vector<string> groupingName;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Filling plots according to grouping
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // PLOTTING AND SAVING
    sprintf(canvName,"cPlot");
    sprintf(graphName,"gPlot");
    TCanvas *canv = new TCanvas(canvName,canvName);
    TGraphErrors *graph[4]; // 4 is the number of cut schemes
    canv->cd();

    for(int cutScheme = 0; cutScheme <= 3; cutScheme++)
    {
        cout << "cutScheme" << cutScheme << endl;
        graph[cutScheme] = new TGraphErrors();
        graph[cutScheme]->SetTitle(Form("DSx_DetType%d_CutScheme%d",enrType,cutScheme));
        graph[cutScheme]->SetMarkerColor(cutColor[cutScheme]);
        graph[cutScheme]->SetMarkerStyle(20);
        i = 0;
        for(unsigned int ds_i = 0; ds_i < dsVec.size(); ds_i++)
        {
          int dsTemp = dsVec[ds_i];
          int chMapDS = dsTemp; if(dsTemp==51 || dsTemp==52 || dsTemp==53){chMapDS=5;}
          MJTChannelMap *chMap = (MJTChannelMap*) mapFile->Get(Form("ChannelMapDS%d",chMapDS));
          if(chMap==NULL) cout << "null chMap"<<endl;
          for (int M = 1; M <= 2; M++)
          {
              sprintf(inTextFilePath,"output_data/open/DS%d_M%d_CutScheme%d.txt",dsTemp,M,cutScheme);
              ifstream inTextFile(inTextFilePath); // scope of loop makes it OK to redefine inTextFile
              cout << dsTemp << " " << M << " " << inTextFilePath << endl;
              while(inTextFile >> m >> p >> d >> dInt >> dIntN >> dEX >> dEXU) // m,p,d,integral,integralNormed,exposure,exposureUncert
              {
                  if(chMap->GetInt(m,p,d,"kDetectorType")==enrType)
                  {
                      if(m == M)
                      {
                          groupInt+=dInt;
                          groupEX+=dEX;
                          nDets++;
                      }
                  }
              }
              inTextFile.close();
          } // end loop over M
          if(nDets!=0){
              norm_groupInt=groupInt/groupEX;
              enorm_groupInt=sqrt(groupInt)/groupEX;
          }
          else{ // nDets == 0
              norm_groupInt=0.;
              enorm_groupInt=0.;
          }
          graph[cutScheme]->SetPoint(i,double(i+1),norm_groupInt);
          if(cutScheme==3){graph[cutScheme]->SetPointError(i,0.0,enorm_groupInt);} // graph[cutScheme]->SetPointError(i,0.0,enorm_groupInt);
          i++;

          // GATHER INPUTS TO THE FIT TO THE cutScheme 3 DATA
          if(cutScheme==3 && nDets!=0){
              tickLabel = to_string(dsTemp);
              groupingFill.push_back(norm_groupInt);
              groupingName.push_back(tickLabel);
          }

          groupInt = 0.;
          norm_groupInt = 0.;
          enorm_groupInt = 0.;
          groupMass = 0.;
          groupEX = 0.;
          nDets = 0;

          delete chMap;
          chMap = NULL;
        } // end loop over DSs
    } // cutScheme

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Plotting and fitting
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // MULTIGRAPH
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(graph[0]);
    mg->Add(graph[1]);
    mg->Add(graph[2]);
    mg->Add(graph[3]);
    if(enrType==0)mg->SetTitle(Form("Rate 1000-1400keV, DSx Natural Detectors by Grouping"));
    if(enrType==2)mg->SetTitle(Form("Rate 1000-1400keV, DSx Enriched Detectors by Grouping"));
    mg->Draw("AP");
    mg->GetYaxis()->SetTitle("Cnt/kg/dy");
    if(enrType==0)mg->GetYaxis()->SetRangeUser(0.002,1.2);
    if(enrType==2)mg->GetYaxis()->SetRangeUser(1.0,2.5);
    mg->GetXaxis()->SetTitle("DS");
    mg->GetXaxis()->SetTitleOffset(1.35);
    TAxis* a = mg->GetXaxis();
    int j = 1;
    int bin_j;
    while (j<((a->GetXmax()))) // may need (a->GetXmax())-1
    {
        bin_j = a->FindBin(j);
        a->SetBinLabel(bin_j, tickLabelVec[j-1].c_str());
        j++;
    }
    canv->Modified();
    //    TLegend *legend=new TLegend(0.6,0.62,0.75,0.85);
    //    legend->SetTextFont(72);
    //    legend->SetTextSize(0.04);
    //    legend->AddEntry(graph[0],"Raw","p");
    //    legend->AddEntry(graph[1],"Raw+DC","p");
    //    legend->AddEntry(graph[2],"Raw+DC+AE","p");
    //    legend->AddEntry(graph[3],"Raw+DC+AE+DCR","p");
    //    legend->Draw();
    canv->Update();

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
    double ytext = 0.;
    if(enrType == 0) {ytext = 1.1;}
    if(enrType == 2) {ytext = 2.4;}
    canv->cd();
    TLatex tlatex;
    tlatex.SetTextSize(13);
    tlatex.DrawLatex(0.78*groupingFill.size(), ytext, Form("Unweighted Fit: %.2f +/- %.2f", p0, p0Err));
    tlatex.DrawLatex(0.78*groupingFill.size(), ytext-.15, Form("Weighted Fit: %.2f +/- %.2f", p0_2, p0Err_2));
    canv->Update();

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Closeout
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    mapFile->Close();
    delete mapFile;
    mapFile = NULL;

    if(save==1) canv->Print(Form("DSx_DetType%d_Fit.pdf",enrType));

    cout << "F I N I S H E D" << endl;
    App->Run();
}
