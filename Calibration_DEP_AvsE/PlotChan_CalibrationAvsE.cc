///////////////
//
// ref: /Users/tomgilliss/Desktop/UNC/ENAP/Analysis/Livetime/ChannelPlots.cc
//
// http://mjwiki.npl.washington.edu/bin/view/Majorana/DS4Cal
//
///////////////

// C++
#include <string>
#include <iostream> // cout
#include <algorithm>    // std::binary_search, sort, find
#include <vector>
#include <map>
//#include <iomanip>      // std::setprecision
//#include <stdlib.h> // exit(), EXIT_FAILURE
//#include <math.h> // ceil

// ROOT
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TApplication.h>

// MJD ... MJSW201705PROD environment
//#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201705Prod/GAT/BaseClasses/GATDataSet.hh"
//#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201705Prod/GAT/Apps/DataSetInfo.hh"
//#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201705Prod/MGDO/Majorana/MJTChannelMap.hh"

using namespace std;

int main (int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Command line arguments and initializations
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 4)
    {
        cout << "Error: Too few arguments" << argv[0] << endl;
        cout << "Usage is: ..." << endl;
        return 1;
    }
    int channel = atoi(argv[1]);
    int startRun = atoi(argv[2]);
    int endRun = atoi(argv[3]);
    
    int UnitsPerBin = 1;
    double xlow = -200., xup = 100.;
    int nbinsx = (xup-xlow)/UnitsPerBin;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Open input file
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    char fName[200];
    sprintf(fName,Form("hists_%d_%d_calAvsE.root",startRun,endRun));
    
    TFile* f = new TFile(fName,"READ");
    f->cd();
    
    
    
    char histName[200];
    char histTitle[200];
    TH1D* h[4];
    TCanvas* c = new TCanvas();
    c->cd();
    TLegend *legend=new TLegend(0.740688,0.527311,0.891117,0.756303);
    int colorList[4] = {1,2,3,4}; // black=FEP,red=SEP,green=DEP,blue=ROI
    char* eRegion[] = {"FEP","SEP","DEP","ROI"};
    
    for(int hist_i = 0; hist_i < 4; hist_i++)
    {
        sprintf(histName,Form("calAvsE_%d_%d_%d_%d",startRun,endRun,channel,hist_i));
        if(f != 0 && f->Get(histName) != 0)
        {
            cout<<"File and hist_i "<<hist_i<<" are good."<<endl;
            h[hist_i] = (TH1D*)f->Get(histName);
            h[hist_i]->SetLineColor(colorList[hist_i]);
            if(hist_i==0){
                sprintf(histTitle,Form("calAvsE, Run %d-%d, Chan %d",startRun,endRun,channel));
                h[hist_i]->SetTitle(histTitle);
                h[hist_i]->SetStats(0);
                h[hist_i]->Draw("HIST");
            }
            else {
                h[hist_i]->Draw("HIST SAME");
            }
            legend->AddEntry(h[hist_i],eRegion[hist_i],"l");
        }
        if(f == 0 && f->Get(histName) == 0) {cout<<"File or hist no good. hist_i "<<hist_i<< endl;}
    }
    legend->Draw();
    c->Update();
    


    
        
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Closeout
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"FINISHED"<<endl;
    
    App->Run();

}