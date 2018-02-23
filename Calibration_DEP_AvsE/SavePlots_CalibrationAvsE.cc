///////////////
//
// ref: /Users/tomgilliss/Desktop/UNC/ENAP/Analysis/Livetime/ChannelPlots.cc
//
// http://mjwiki.npl.washington.edu/bin/view/Majorana/DS4Cal
//
//
// ./SavePlots_CalibrationAvsE 4 60001442 60001461
// ./SavePlots_CalibrationAvsE 4 60001508 60001521
// ./SavePlots_CalibrationAvsE 4 60001578 60001592
// ./SavePlots_CalibrationAvsE 4 60001720 60001732
// ./SavePlots_CalibrationAvsE 4 60001855 60001868
// ./SavePlots_CalibrationAvsE 4 60001914 60001926
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
#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201705Prod/MGDO/Majorana/MJTChannelMap.hh"

using namespace std;

int main (int argc, char* argv[])
{
    //TApplication *App = new TApplication("App", 0, NULL);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Command line arguments and initializations
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 3)
    {
        cout << "Error: Too few arguments" << argv[0] << endl;
        cout << "Usage is: ..." << endl;
        return 1;
    }
    int DS = atoi(argv[1]);
    int startRun = atoi(argv[2]);
    int endRun = atoi(argv[3]);
    cout<<"DS "<<DS<<" "<<startRun<<"-"<<endRun<<endl;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Make channel list
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TFile* ChannelMapFile = new TFile("./DSChannelMaps.root");
    MJTChannelMap *chmap = (MJTChannelMap*)ChannelMapFile->Get(Form("ChannelMapDS%d",DS));
    std::map<int,string> channelmap; // this will map between channel # and previous-event timestamp
    vector<int> channellist;
    int tempchan = 0;
    string tempstring = "";
    std::ostringstream tempstringStream;
    for(int c = 1; c <= 2; c++)
    {
        for(int p = 1; p <= 7; p++)
        {
            for(int d = 1; d <= 5; d++)
            {
                if(chmap->GetDetectorName(c,p,d)!="") // avoid nonexistent detectors
                {
                    tempstringStream << "C" << c << "P" << p << "D" << d;
                    tempstring = tempstringStream.str();
                    
                    tempchan = chmap->GetInt(chmap->GetDetectorName(c,p,d),"kIDHi"); // load HG chan #s
                    channellist.push_back(tempchan);
                    channelmap.insert(std::pair<int,string>(tempchan,tempstring));
                    
                    // Do not include LG channels
                    
                    tempstringStream.str(""); // reset the stringstream
                }
            }
        }
    }
    sort(channellist.begin(), channellist.end()); // list must be sorted for binary_search
    //for(unsigned int list_i = 0; list_i < channellist.size(); list_i++) {
    //    cout<<channellist[list_i]<<" "<<channelmap.at(channellist[list_i])<<endl;
    //}

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
    int channel = 0;
    int colorList[4] = {1,2,3,4}; // black=FEP,red=SEP,green=DEP,blue=ROI
    char* eRegion[] = {"FEP","SEP","DEP","ROI"};
    
    // ONLY PLOT DEP HIST FOR EACH CHANNEL
    TH1D* h;
    TCanvas* c = new TCanvas();
    c->cd();
    TLegend *legend=new TLegend(0.740688,0.527311,0.891117,0.756303);
    int hist_i = 2; // hard code the DEP hist index
    for(unsigned int list_i = 0; list_i < channellist.size(); list_i++)
    {
        channel = channellist[list_i];
        sprintf(histName,Form("calAvsE_%d_%d_%d_%d",startRun,endRun,channel,hist_i));
        
        if(f != 0 && f->Get(histName) != 0 && channel!=1144)
        {
            if(list_i==0){
                h = (TH1D*)f->Get(histName);
                h->SetLineColor(colorList[hist_i]);
                sprintf(histTitle,Form("calAvsE, Run %d-%d, Good HG Chans",startRun,endRun));
                h->SetTitle(histTitle);
                h->SetStats(0);
                h->Draw("HIST GOFF");
                legend->AddEntry(h,eRegion[hist_i],"l");
            }
            if(list_i!=0){
                TH1D* hOther;
                hOther = (TH1D*)f->Get(histName);
                hOther->SetLineColor(colorList[hist_i]);
                hOther->SetStats(0);
                h->Add(hOther);
                // CLEAR MEMORY FOR NEXT ITERATION OF LOOP
                delete hOther;
            }
        }
        if(f == 0 || f->Get(histName) == 0 || channel==1144) {cout<<"File or hist hist_i "<<hist_i<<" no good. Or channel "<<channel<<" excluded."<< endl;}
    }
    c->SetLogy();
    h->GetXaxis()->SetRangeUser(-500.,500.);
    legend->Draw();
    c->Update();
    c->Print(Form("PDFs/AvsE_DEP_%d_%d.pdf",startRun,endRun));
    


    
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Closeout
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"FINISHED"<<endl;
    
    //App->Run();

}