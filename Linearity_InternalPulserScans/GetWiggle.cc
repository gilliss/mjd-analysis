/***************************************************
July 2017
****************************************************/

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include "TF1.h"
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TApplication.h>

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>
#include <math.h>

#include <GATDataSet.hh>
#include <MJTChannelMap.hh>

using namespace std;

void RetrievePlots(int M, int chan, char* eChar, char* aChar)
{
    /* TFILE TO OPEN */
    char resultsPath[200], inFile[100], inFilePath[300];
    sprintf (resultsPath,"/global/u2/g/gilliss/PulserLinearity/Update_July2017_DS34/results/");
    sprintf (inFile,"PulserLin_M%d_%s_%s_10-11-2016.root",M,aChar,eChar);
    sprintf (inFilePath,"%s%s",resultsPath,inFile);
    TFile *f = new TFile(inFilePath); // Open an existing file for reading (default) ... READ
    TCanvas *c = new TCanvas("c","c",1300,400);
    c->Divide(2,1);
    
//    if(f != 0 &&
//       f->Get(Form("EvR_%d_%s_%s",chan,aChar,eChar)) != 0 &&
//       f->Get(Form("Fit_%d_%s_%s",chan,aChar,eChar)) != 0 &&
//       f->Get(Form("LinRes_%d_%s_%s",chan,aChar,eChar)) != 0) {cout<< chan <<" ...All files exist"<< endl;}
    
    if(f == 0 ||
       f->Get(Form("EvR_%d_%s_%s",chan,aChar,eChar)) == 0 ||
       f->Get(Form("Fit_%d_%s_%s",chan,aChar,eChar)) == 0 ||
       f->Get(Form("LinRes_%d_%s_%s",chan,aChar,eChar)) == 0) {cout<< chan <<" ...No good!!"<< endl;}
    
    TProfile *prof = (TProfile*) f->Get(Form("EvR_%d_%s_%s",chan,aChar,eChar));
        gStyle->SetErrorX(0.);
    TF1 *fit;
    TGraph *graph;
    if(f->Get(Form("Fit_%d_%s_%s",chan,aChar,eChar))!=0)
    {
        fit = (TF1*) f->Get(Form("Fit_%d_%s_%s",chan,aChar,eChar));
        graph = (TGraph*) f->Get(Form("LinRes_%d_%s_%s",chan,aChar,eChar));
    }
    c->cd(1);
    prof->Draw();
    if(f->Get(Form("Fit_%d_%s_%s",chan,aChar,eChar))!=0)
    {
        fit->Draw("SAME");
        c->Update();
        c->cd(2);
        graph->GetYaxis()->SetRangeUser(-2.0,2.0);//(-4.0,4.0);
        graph->Draw("AP");
    }
    
    
    //if(chan==1302)
    //{
        double* xdata = graph->GetX();
        double* ydata = graph->GetY();
        // Compute Std Dev in the style of ROOT's GetRMS()
        double sumx = 0, sumx2 = 0;
        int fNpoints = graph->GetN();
        for(int i_data = 0; i_data < fNpoints; i_data++)
        {
            if(xdata[i_data]>=100 && xdata[i_data]<=275)
            {
                //cout<<xdata[i_data]<<" "<<ydata[i_data]<<endl;
                sumx += ydata[i_data];
                sumx2 += ydata[i_data] * ydata[i_data];
            }
        }
        Double_t x = sumx / fNpoints;
        Double_t rms2 = abs(sumx2 / fNpoints - x * x);
        if(fNpoints>1 && sqrt(rms2)>0.0){cout<<sqrt(rms2)<<endl;}
    //}
    
    f->Close();
    delete c;
}

int main (int argc, char* argv[])
{
    //TApplication *App = new TApplication("App", 0, NULL);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Read in command line arguments, make TChain
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    /* READ IN ARGUMENTS/SETTINGS FROM COMMAND */
    if (argc < 3)
    {
        cout << "too few arguments" << argv[0] << endl;
        return 1;
    }
    int module = atoi(argv[1]);
    int e_i = atoi(argv[2]);
    int att_i = atoi(argv[3]); // 0=Un, 1=Inter, 2=Final
    int gain = atoi(argv[4]); // 0=HG, 1=LG
    char eEstimator[10];
        if (e_i == 0) {sprintf(eEstimator,"trapECal");}
        if (e_i == 1) {sprintf(eEstimator,"trapENFCal");}
        if (e_i == 2) {sprintf(eEstimator,"trapEN2FCal");}
        if (e_i == 3) {sprintf(eEstimator,"trapE");}
        if (e_i == 4) {sprintf(eEstimator,"trapENF");}
        if (e_i == 5) {sprintf(eEstimator,"trapEN2F");}
    char attChar[10];
        if (att_i == 0) {sprintf(attChar,"Unatten");}
        if (att_i == 1) {sprintf(attChar,"InterAtten");}
        if (att_i == 2) {sprintf(attChar,"FinalAtten");}

    /* DS5 P3LQK Channel Lists */
    int M1ChanArray[] = {/*584,*/585,582,583,/*580,*/581,578,579,/*692,*/693,648,649,640,641,/*642,*/643,/*616,*/617,610,611,608,609,664,665,624,625,/*628,*/629,688,689,694,695,/*614,*/615,/*680,*/681,678,679,672,673,/*696,*/697,632,633,/*630,*/631,626,627,690,691,600,601,598,599,594,595,592,593};
    int M2ChanArray[] = {1140,1141,1142,1143,1110,1111,1204,1205,/*1174,*/1175,/*1144,*/1145,/*1106,*/1107,/*1108,*/1109,/*1138,*/1139,1176,1177,1172,1173,/*1202,*/1203,1170,1171,/*1208,*/1209,/*1206,*/1207,1136,1137,/*1168,*/1169,1330,1331,/*1304,*/1305,/*1332,*/1333,/*1302,*/1303,/*1296,*/1297,1298,1299,/*1328,*/1329,/*1234,*/1235,/*1268,*/1269,/*1238,*/1239,/*1236,*/1237,1232,1233};
    
    if(module==1){
        for(unsigned int i = 0; i<sizeof(M1ChanArray)/sizeof(*M1ChanArray); i++)
        {
            if(gain==0 && M1ChanArray[i]%2==0) RetrievePlots(module, M1ChanArray[i], eEstimator, attChar);
            if(gain==1 && M1ChanArray[i]%2!=0) RetrievePlots(module, M1ChanArray[i], eEstimator, attChar);
        }
    }
    if(module==2){
        for(unsigned int i = 0; i<sizeof(M2ChanArray)/sizeof(*M2ChanArray); i++)
        {
            if(gain==0 && M2ChanArray[i]%2==0) RetrievePlots(module, M2ChanArray[i], eEstimator, attChar);
            if(gain==1 && M2ChanArray[i]%2!=0) RetrievePlots(module, M2ChanArray[i], eEstimator, attChar);
        }
    }
    
    cout << "F I N I S H E D" << endl;
    //App->Run();
}

//    double offset = fit->GetParameter(0);
//    double slope = fit->GetParameter(1);
//    int xmin = fit->GetXmin();
//    int xmax = fit->GetXmin();
//    TF1 *fitNew = new TF1("fitNew",Form("%f*x+%f",slope,offset),xmin,xmax);
//fitNew->SetLineColor(kRed);
//fitNew->SetLineWidth(2);
//fitNew->Draw("SAME");