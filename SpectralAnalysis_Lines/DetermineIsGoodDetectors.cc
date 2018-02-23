/*

*/

#include <TApplication.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TFile.h>

#include <MJTChannelMap.hh>

#include <iostream>
#include <string>
#include <stdlib.h>     /* exit(), EXIT_FAILURE */

using namespace std;

void Error(string errormessage) {
    // print message and exit program
    cout<<"Error: "<<errormessage<<endl;
    exit(EXIT_FAILURE);
}

int main(int argc, char* argv[]) {
    TApplication *App = new TApplication("App", 0, NULL);
    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 2)
    {
        cout << "Error: Too few arguments" << argv[0] << endl;
        cout << "Usage is: ./DetermineIsGoodDetectors DS C" << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // DS#
    int c = atoi(argv[2]); // M#
    
    char skim[50];

    // READ IN DATA
    TChain *ch;
    if(DS==0) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS0/GAT-v01-06-134-g3f44fab/skimDS0*.root");
        sprintf(skim,"GAT-v01-06-134-g3f44fab");
    }
    if(DS==1) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS1/GAT-v01-06-134-g3f44fab/skimDS1*.root");
        sprintf(skim,"GAT-v01-06-134-g3f44fab");
    }
    if(DS==2) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS2/GAT-v01-06-134-g3f44fab/skimDS2*.root");
        sprintf(skim,"GAT-v01-06-134-g3f44fab");
    }
    if(DS==3) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS3/GAT-v01-06-134-g3f44fab/skimDS3*.root");
        sprintf(skim,"GAT-v01-06-134-g3f44fab");
    }
    if(DS==4) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS4/GAT-v01-06-134-g3f44fab/skimDS4*.root");
        sprintf(skim,"GAT-v01-06-134-g3f44fab");
    }
    if(DS==5) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS5/GAT-v01-06-134-g3f44fab/skimDS5*.root");
        sprintf(skim,"GAT-v01-06-134-g3f44fab");
    }
    if(DS==6) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS6/GAT-v01-06-22-g5c7b61c/skimDS6*.root");
        sprintf(skim,"GAT-v01-06-134-g3f44fab");
    }
    
    char title[50];
    
    // PLOT 2D HIST
    TCanvas* canv = new TCanvas();
//    canv->Divide(2,1);
    canv->cd(1);
    sprintf(title,"DS%d M%d isGood==1 HG Skim:%s",DS,c,skim);
    TH2D *h = new TH2D("h",title,7,0.5,7.5,5,0.5,5.5);
    h->SetStats(0);
    ch->Draw("D:P>>h",Form("C==%d && isGood==1 && channel%2==0",c),"COLZ");
    
//    canv->cd(2);
    TCanvas* canv2 = new TCanvas();
    canv2->cd(1);
    sprintf(title,"DS%d M%d isGood==0 HG Skim:%s",DS,c,skim);
    TH2D *h2 = new TH2D("h2",title,7,0.5,7.5,5,0.5,5.5);
    h2->SetStats(0);
    ch->Draw("D:P>>h2",Form("C==%d && isGood==0 && channel%2==0",c),"COLZ");
    
    // LOOP THROUGH HIST TO MAKE USEFUL OUTPUT
    TFile *mapFile = new TFile("./ChannelMaps/DSChannelMaps.root","READ");
    MJTChannelMap *chMap = (MJTChannelMap*) mapFile->Get(Form("ChannelMapDS%d",DS));
    if(chMap==NULL) cout << "null chMap"<<endl;
    int iIsNotGood = 0, iIsGood = 1, iDNE = -1;
    int detStatus = iDNE;
    char funcOutput[50];
    for(int p = 1; p <= 7; p++)
    {
        for(int d = 1; d <= 5; d++)
        {
            //cout<<"P"<<p<<"D"<<d;
            detStatus = iDNE;
            if(chMap->GetDetectorName(c,p,d)!="") // avoid nonexistent detectors
            {
                detStatus = iIsNotGood;
                //cout<<"  "<<h->GetBinContent((h->GetBin(p,d)));
                if(h->GetBinContent((h->GetBin(p,d)))>0) {detStatus = iIsGood;}
            }
            //cout<<"  "<<detStatus<<endl;
            sprintf(funcOutput,"if(DS==%d && C==%d && P==%d && D==%d) {detStatus=%d;}",DS,c,p,d,detStatus);
            cout<<funcOutput<<endl;
        }
    }
    
    App->Run();
}