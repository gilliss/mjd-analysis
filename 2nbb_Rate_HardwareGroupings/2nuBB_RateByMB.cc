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

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>
#include <algorithm>    // std::binary_search

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
        cout << "Usage is: ./2nuBB_RateByDet DS C enrType" << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // DS#
    int c = atoi(argv[2]); // M#
    int enrType = atoi(argv[3]); // Nat(0) or Enr(2)
    if(enrType!=0 && enrType!=2) {cout<<"Error:  Nat(0) or Enr(2)"<<endl; return 1;}
    int save = atoi(argv[4]);
    
    // PRINT SETTINGS
    cout<<"------------"<<endl;
    cout<<"Settings:"<<endl;
    cout<<" Plotting DS" << DS << " for M"<<c<< endl;
    cout<<"------------"<<endl;
    
    // INPUT MJTCHANNELMAP FILE
    TFile *mapFile = new TFile("/global/u2/g/gilliss/2nuBB_Systematics/Update_October2017/channel_map_data/DSChannelMaps.root","READ");
    
    // INPUT DATA FROM TEXT FILE
    char inTextFilePath[50];
    int m, p, d;
    /////int dInt;
    /////double dMkg, dLT;
    int dInt; // integral
    double dIntN, dEX, dEXU; // normed integral, exposure, exposure uncert
    
    // INITIALIZATIONS
    int maxMod = 2, maxPos = 7, maxDet = 5, maxCC = 4, maxMBperCC = 2;
    char canvName[50];
    char graphName[50];
    int cutColor[4]={1,16,92,59}; // to match Lukas
    double norm_dInt = 0;
    double enorm_dInt = 0;
    double scale = 0;
    int i = 0;
    int nDets = 0;
    map <string, double> ccCoord; // <cpd,xCoordForGraph>
    vector <string> ccNameVec;
    double mbInt = 0;
    double norm_mbInt = 0;
    double enorm_mbInt = 0;
    /////double mbMass = 0;
    int ccList[4];
    int chMapDS = DS; if(DS==51 || DS==52){chMapDS=5;} // probably need the MB mappings for other DSs
    if(chMapDS==5 && c==1) {ccList[0] = 6; ccList[1] = 7; ccList[2] = 9; ccList[3] = 10;}
    if(chMapDS==5 && c==2) {ccList[0] = 5; ccList[1] = 9; ccList[2] = 13; ccList[3] = 17;}
    if(chMapDS==5 && c==1) {
        ccNameVec.push_back("6,0-5");
        ccNameVec.push_back("6,8-13");
        ccNameVec.push_back("7,0-5");
        ccNameVec.push_back("7,8-13");
        ccNameVec.push_back("9,0-5");
        ccNameVec.push_back("9,8-13");
        ccNameVec.push_back("10,0-5");
        ccNameVec.push_back("10,8-13");
    }
    if(chMapDS==5 && c==2){
        ccNameVec.push_back("5,0-5");
        ccNameVec.push_back("5,8-13");
        ccNameVec.push_back("9,0-5");
        ccNameVec.push_back("9,8-13");
        ccNameVec.push_back("13,0-5");
        ccNameVec.push_back("13,8-13");
        ccNameVec.push_back("17,0-5");
        ccNameVec.push_back("17,8-13");
    }
    vector< vector<int> > mbList(2, vector<int>(6)); // these vectors must be ordered for binary_search
    mbList[0][0] = 0; mbList[0][1] = 1; mbList[0][2] = 2; mbList[0][3] = 3; mbList[0][4] = 4; mbList[0][5] = 5;
    mbList[1][0] = 8; mbList[1][1] = 9; mbList[1][2] = 10; mbList[1][3] = 11; mbList[1][4] = 12; mbList[1][5] = 13;

    string ccName; //char cpdName[50];
    vector<double> groupingFill;
    vector<string> groupingName;
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Plotting and saving
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // PLOTTING AND SAVING
    sprintf(canvName,"c_DS%d_M%d_DetType%d",DS,c,enrType);
    sprintf(graphName,"g_DS%d_M%d_DetType%d",DS,c,enrType);
    TCanvas *canv = new TCanvas(canvName,canvName);
    TGraphErrors *graph[4];
    canv->cd();

    for(int cutScheme = 0; cutScheme <= 3; cutScheme++)
    {
        cout << "cutScheme" << cutScheme << endl;
        MJTChannelMap *chMap = (MJTChannelMap*) mapFile->Get(Form("ChannelMapDS%d",chMapDS));
        if(chMap==NULL) cout << "null chMap"<<endl;
        graph[cutScheme] = new TGraphErrors();
        graph[cutScheme]->SetTitle(Form("DS%d_C%d_DetType%d_CutScheme%d",DS,c,enrType,cutScheme));
        graph[cutScheme]->SetMarkerColor(cutColor[cutScheme]);
        graph[cutScheme]->SetMarkerStyle(20);
        i = 0;
        for (int CC = 0; CC <= (maxCC-1); CC++)
        {
            //cout << "CC "<<ccList[CC]<<endl;
            for(int MBperCC = 0; MBperCC <= (maxMBperCC-1); MBperCC++)
            {
                //cout << "  MB Chan Bank "<<MBperCC<<endl;
                sprintf(inTextFilePath,"output_data/DS%d_M%d_CutScheme%d.txt",DS,c,cutScheme);
                ifstream inTextFile(inTextFilePath); // scope of loop makes it OK to redefine inTextFile
                while(inTextFile >> m >> p >> d >> dInt >> dIntN >> dEX >> dEXU)
                {
                    if(chMap->GetInt(m,p,d,"kDetectorType")==enrType)
                    {
                        if(chMap->GetInt(m,p,d,"kPreAmpDigitizer") == ccList[CC])
                        {
                            //cout << "    CC " << chMap->GetInt(m,p,d,"kPreAmpDigitizer") << endl;
                            //cout << "    MB " << chMap->GetInt(m,p,d,"kPreAmpChan") << endl;

                            if(binary_search(mbList[MBperCC].begin(), mbList[MBperCC].end(), chMap->GetInt(m,p,d,"kPreAmpChan"))) // see if channel is in MB list
                            //if(chMap->GetInt(m,p,d,"kPreAmpChan") > mbList[MBperCC][0] && chMap->GetInt(m,p,d,"kPreAmpChan") < mbList[MBperCC][5])
                            {
                                //cout << "      AccMB " << chMap->GetInt(m,p,d,"kPreAmpChan") << endl;
                                scale = 1/(dEX);
                                norm_dInt = dInt*scale;
                                norm_mbInt+=norm_dInt;
                                enorm_mbInt+=sqrt(dInt)*scale;
                                mbInt+=dInt;
                                /////mbMass+=dMkg;
                                nDets++;
                            }
                        }
                    }
                }
                if(nDets!=0){
                    norm_mbInt=norm_mbInt/nDets;
                    enorm_mbInt=enorm_mbInt/nDets;
                }
                else{ // nDets == 0
                    norm_mbInt=0.;
                    enorm_mbInt=0.;
                }
                graph[cutScheme]->SetPoint(i,double(i+1),norm_mbInt);
                graph[cutScheme]->SetPointError(i,0.0,enorm_mbInt);
                i++;
                
                if(cutScheme==3 && nDets!=0){
                    ccName = to_string(ccList[CC]);
                    groupingFill.push_back(norm_mbInt);
                    groupingName.push_back(ccName);
                }
                
                mbInt = 0.;
                norm_mbInt = 0.;
                enorm_mbInt = 0.;
                /////mbMass = 0.;
                nDets = 0;
            }
        }
        delete chMap;
        chMap = 0;
    } // cutScheme
    
    // MULTIGRAPH
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(graph[0]);
    mg->Add(graph[1]);
    mg->Add(graph[2]);
    mg->Add(graph[3]);
    if(enrType==0)mg->SetTitle(Form("Rate 1000-1400keV, DS%d M%d Natural Detectors by MB",DS,c));
    if(enrType==2)mg->SetTitle(Form("Rate 1000-1400keV, DS%d M%d Enriched Detectors by MB",DS,c));
    mg->Draw("AP");
    //mg->SetTitleOffset(0.5);
    mg->GetYaxis()->SetTitle("Cnt/kg/dy");
    if(enrType==0)mg->GetYaxis()->SetRangeUser(0.005,1.7);
    if(enrType==2)mg->GetYaxis()->SetRangeUser(0.5,2.6);
    mg->GetXaxis()->SetTitle("CC,CC_Chan (kPreAmpDigitizer,kPreAmpChan)");
    mg->GetXaxis()->SetTitleOffset(2.0);
    mg->GetXaxis()->SetLabelSize(25); // 0.04
    TAxis* a = mg->GetXaxis();
    int j = 1;
    int bin_j;
    while (j<((a->GetXmax()))) // may need (a->GetXmax())-1
    {
        bin_j = a->FindBin(j);
        a->SetBinLabel(bin_j,ccNameVec[j-1].c_str());
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
    canv->SetBottomMargin(0.2);
    canv->Update();
    

    // HIST FOR FIT
    double binLo = 1.0-0.5;
    double binHi = groupingFill.size()+0.5;
    int nBins = binHi - binLo;
    sprintf(canvName,"c2_DS%d_M%d_DetType%d",DS,c,enrType);
    sprintf(graphName,"h_DS%d_M%d_DetType%d",DS,c,enrType);
    TCanvas *canv2 = new TCanvas(canvName,canvName);
    TH1D *h = new TH1D(graphName,graphName,nBins,binLo,binHi);
    h->SetMarkerStyle(20);
    h->SetMarkerColor(cutColor[3]);
    h->Sumw2(); // the sum of squares of weights is also stored
    TAxis* a2 = h->GetXaxis();
    for(int j = 0; j<groupingFill.size(); j++)
    {
        h->Fill(j+1,groupingFill[j]);
        a2->SetBinLabel(j+1,groupingName[j].c_str());
        cout << "plot name fill "<<groupingFill[j]<<" "<<groupingName[j]<<endl;
    }
    canv2->cd();
    h->Draw("P HIST");
    TF1  *f1 = new TF1("f1","pol0");
    h->Fit(f1,"WL"); //h->Fit("pol0", "WL"); //“WL” Weighted log likelihood method. To be used when the histogram has been filled with weights different than 1.
    f1->Draw("SAME");
    double p0 = f1->GetParameter(0);
    cout << "p0 " << p0 << endl;
    canv2->Update();
    canv->cd();
    TLine *line = new TLine(a->GetXmin(),p0,a->GetXmax(),p0); // x1,y1,x2,y2
    line->SetLineColorAlpha(cutColor[3], 0.45);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw("same");
    canv->Update();


    // CLOSEOUT
    mapFile->Close();
    delete mapFile;
    mapFile = 0;
    
    if(save==1) canv->Print(Form("DS%d_M%d_DetType%d_MB_Fit.pdf",DS,c,enrType));
    
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