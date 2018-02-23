///////////////
//
// ref: /Users/tomgilliss/Desktop/UNC/ENAP/Analysis/Livetime/ChannelPlots.cc
//
// http://mjwiki.npl.washington.edu/bin/view/Majorana/DS4Cal
// x x 60001442 60001461
// 60001442 60001507 <--- covered region
// x x 60001508 60001521 <--- last "good" cal of DS
// x x 60001578 60001592
// 60001578 60001719 <--- covered region
// x x 60001720 60001732
// 60001720 60001854 <--- covered region
// x x 60001855 60001868
// 60001855 60001913 <--- covered region
// x x 60001914 60001926
// 60001914 60002394 <--- covered region
//
//
//
//AvsE distributions, from the 208Tl FEP/SEP/DEP/0nbbROI, for each channel in the following DS4 calibration sets:
//60001442 60001461
//60001508 60001521 <--- last "good" cal of DS4
//60001578 60001592
//60001720 60001732
//60001855 60001868
//60001914 60001926
//
//
// ./CalibrationAvsE 4 60001442 60001461
// ./CalibrationAvsE 4 60001508 60001521
// ./CalibrationAvsE 4 60001578 60001592
// ./CalibrationAvsE 4 60001720 60001732
// ./CalibrationAvsE 4 60001855 60001868
// ./CalibrationAvsE 4 60001914 60001926
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
#include <TApplication.h>

// MJD ... MJSW201705PROD environment
#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201705Prod/GAT/BaseClasses/GATDataSet.hh"
#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201705Prod/GAT/Apps/DataSetInfo.hh"
#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201705Prod/MGDO/Majorana/MJTChannelMap.hh"

using namespace std;

double GetHalfROI(double E) {
    /* compute HWHM based on Pinghan's DS3 resolution data */
    double p0 = 0.164225116870199;
    double p1 = 0.0171317589551445;
    double p2 = 0.000305455041318674;
    double sigma = sqrt(p0*p0 + p1*p1*E + p2*p2*E*E);
    double HalfROI = 0.5*(5*sigma);
    return HalfROI;
}

int main (int argc, char* argv[])
{
    //TApplication *App = new TApplication("App", 0, NULL);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Command line arguments and initializations
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 4)
    {
        cout << "Error: Too few arguments" << argv[0] << endl;
        cout << "Usage is: ./CalibrationAvsE DS startRun endRun" << endl;
        return 1;
    }
    int DS = atoi(argv[1]);
    int startRun = atoi(argv[2]);
    int endRun = atoi(argv[3]);
    cout<<"DS "<<DS<<" "<<startRun<<"-"<<endRun<<endl;
    
    double FEP = 2614.533; int iFEP = 0;
    double SEP = FEP - 511.; int iSEP = 1;
    double DEP = SEP - 511.; int iDEP = 2;
    double ROI = 2039.; int iROI = 3;
    cout<<"Example: FEP window is "<<FEP - GetHalfROI(FEP)<<"-"<<FEP + GetHalfROI(FEP)<<endl;
    
    int UnitsPerBin = 1;
    double xlow = -3000., xup = 1000.;
    int nbinsx = (xup-xlow)/UnitsPerBin;

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
    /////////// Make channel-to-TH1D_array map
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::map<int, vector<TH1D*> > chanHistMap;
    char histName[200];
    
    // INITIALIZE THE MAP OF VECTORS TO HAVE BLANK HISTS
    for(unsigned int list_i = 0; list_i < channellist.size(); list_i++)
    {
        chanHistMap.insert(std::pair<int, vector<TH1D*> >(channellist[list_i],vector<TH1D*>()));
        for(int hist_i = 0; hist_i < 4; hist_i++)
        {
            sprintf(histName,Form("calAvsE_%d_%d_%d_%d",startRun,endRun,channellist[list_i],hist_i));
            chanHistMap.at(channellist[list_i]).push_back(new TH1D(histName,histName,nbinsx,xlow,xup));
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Make chain
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TChain* ch = new TChain("mjdTree", "mjdTree");
    for(int i_run = startRun; i_run<=endRun; i_run++)
    {
        ch->Add(Form("/global/project/projectdirs/majorana/data/mjd/surfmjd/data/gatified/P3LQG/mjd_run%d.root",i_run));
    }
    
    vector<double>* channel = NULL;
    vector<double>* TSCurrent50nsMax = NULL;
    vector<double>* TSCurrent100nsMax = NULL;
    vector<double>* TSCurrent200nsMax = NULL;
    vector<double>* trapENF = NULL;
    vector<double>* trapENFCal = NULL;
    int run = 0;
    
    double avse = 0.;
    
    ch->SetBranchAddress("channel",&channel);
    ch->SetBranchAddress("TSCurrent50nsMax",&TSCurrent50nsMax);
    ch->SetBranchAddress("TSCurrent100nsMax",&TSCurrent100nsMax);
    ch->SetBranchAddress("TSCurrent200nsMax",&TSCurrent200nsMax);
    ch->SetBranchAddress("trapENF",&trapENF);
    ch->SetBranchAddress("trapENFCal",&trapENFCal);
    ch->SetBranchAddress("run",&run);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Fill hists
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    int nentries = ch->GetEntries();
    cout <<"nentries "<<nentries<<endl;
    for(int i_iter=0; i_iter < nentries; i_iter++)
    {
        ch->GetEntry(i_iter);
        if(i_iter%500000==0) {cout<<"at entry "<<i_iter<<"/"<<nentries<<endl;}
        for (unsigned int j_iter=0; j_iter<channel->size(); j_iter++)
        {
            //cout<<channel->size()<<" "<<trapENFCal->size()<<endl;
            tempchan = channel->at(j_iter);
            if(tempchan%2==0)
            {
                // FEP
                if(trapENFCal->at(j_iter) < FEP + GetHalfROI(FEP) &&
                   trapENFCal->at(j_iter) > FEP - GetHalfROI(FEP))
                {
                    avse = GetAvsE(channel->at(j_iter),TSCurrent50nsMax->at(j_iter),TSCurrent100nsMax->at(j_iter),TSCurrent200nsMax->at(j_iter),trapENF->at(j_iter),trapENFCal->at(j_iter),DS,run);
                    chanHistMap.at(channel->at(j_iter)).at(iFEP)->Fill(avse);
                }
                //SEP
                if(trapENFCal->at(j_iter) < SEP + GetHalfROI(SEP) &&
                   trapENFCal->at(j_iter) > SEP - GetHalfROI(SEP))
                {
                    avse = GetAvsE(channel->at(j_iter),TSCurrent50nsMax->at(j_iter),TSCurrent100nsMax->at(j_iter),TSCurrent200nsMax->at(j_iter),trapENF->at(j_iter),trapENFCal->at(j_iter),DS,run);
                    chanHistMap.at(channel->at(j_iter)).at(iSEP)->Fill(avse);
                }
                //DEP
                if(trapENFCal->at(j_iter) < DEP + GetHalfROI(DEP) &&
                   trapENFCal->at(j_iter) > DEP - GetHalfROI(DEP))
                {
                    avse = GetAvsE(channel->at(j_iter),TSCurrent50nsMax->at(j_iter),TSCurrent100nsMax->at(j_iter),TSCurrent200nsMax->at(j_iter),trapENF->at(j_iter),trapENFCal->at(j_iter),DS,run);
                    chanHistMap.at(channel->at(j_iter)).at(iDEP)->Fill(avse);
                }
                //ROI
                if(trapENFCal->at(j_iter) < ROI + GetHalfROI(ROI) &&
                   trapENFCal->at(j_iter) > ROI - GetHalfROI(ROI))
                {
                    avse = GetAvsE(channel->at(j_iter),TSCurrent50nsMax->at(j_iter),TSCurrent100nsMax->at(j_iter),TSCurrent200nsMax->at(j_iter),trapENF->at(j_iter),trapENFCal->at(j_iter),DS,run);
                    chanHistMap.at(channel->at(j_iter)).at(iROI)->Fill(avse);
                }
            }
        }
    }
    
    // EXAMPLE PLOT
    /*
     TCanvas* c = new TCanvas();
     c->cd();
     chanHistMap.at(1298).at(iFEP)->Draw("HIST");
     */
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Write to file
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    char fName[200];
    sprintf(fName,Form("hists_%d_%d_calAvsE.root",startRun,endRun));
    
    TFile* f = new TFile(fName,"CREATE");
    f->cd();
    
    // LOOP MAP AND VECTORS OF HISTS
    int histEntries = 0;
    int bin_0 = 0;
    int content_0 = 0;
    for(unsigned int list_i = 0; list_i < channellist.size(); list_i++)
    {
        chanHistMap.insert(std::pair<int, vector<TH1D*> >(channellist[list_i],vector<TH1D*>()));
        for(int hist_i = 0; hist_i < 4; hist_i++)
        {
            // CHECK THAT HIST IS NOT EMPTY and that all contents are not in the AvsE=0 bin
            histEntries = chanHistMap.at(channellist[list_i]).at(hist_i)->GetEntries();
            TAxis* a = chanHistMap.at(channellist[list_i]).at(hist_i)->GetXaxis();
            bin_0 = a->FindBin(0.);
            content_0 = chanHistMap.at(channellist[list_i]).at(hist_i)->GetBinContent(bin_0);
            if(histEntries != 0 && content_0 != histEntries)
            {
                chanHistMap.at(channellist[list_i]).at(hist_i)->GetYaxis()->SetTitle("counts");
                chanHistMap.at(channellist[list_i]).at(hist_i)->GetXaxis()->SetTitle("avse");
                chanHistMap.at(channellist[list_i]).at(hist_i)->Write();
            }
            else {cout<<"Excluded channel: "<<channellist[list_i]<<endl;}
        }
    }
    f->Write();
    f->Close();
    //delete f;
    //f = 0;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Closeout
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"FINISHED"<<endl;

    //App->Run();

}