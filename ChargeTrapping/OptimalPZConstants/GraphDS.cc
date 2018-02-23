/*
A script to read in the optimal tau/FWHM for each detector in the STCs and make a graph.

The format of the .txt file is the following:
[chan] [optTau] [eoptTau] [optFWHM] [eoptFWHM] [qtTau]
 
Example use: ./GraphSTC
*/

#include <iostream>
#include <fstream>
#include <algorithm> // min_element, sort, binary_search
#include <vector>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TFile.h>
#include <TLatex.h>

#include <MJTChannelMap.hh> // https://github.com/mppmu/MGDO/blob/master/Majorana/MJTChannelMap.hh

using namespace std;

//string GetDetSN(int DS, int channel)
//{
//    //within this routine, open chanmap root file and find the det SN based on the channel //
//    string detSN="";
//
//    return detSN;
//}

int main (int argc, char* argv[])
{
///////////////////////////////
//// INITIALIZATIONS
///////////////////////////////
    TApplication *App = new TApplication("App", 0, NULL);
    
    if (argc < 2)
    {
        cout << "too few arguments " << argv[0] << endl;
        return 1;
    }
    int DS = atoi(argv[1]);
    
    // IN FILE PARAMETERS //
    string infilepath = "/global/u2/g/gilliss/ChargeTrapping/Update_March2017/";
    string infilename = "";
    int channel=0;
    double optTau=0, eoptTau=0, optFWHM=0, eoptFWHM=0, qtTau=0;
    
    // INPUT MJTCHANNELMAP FILE
    TFile *mapFile = new TFile("/global/u2/g/gilliss/2nuBB_Systematics/data/DSChannelMaps.root","READ");
    MJTChannelMap *chMap = (MJTChannelMap*) mapFile->Get(Form("ChannelMapDS%d",DS));
    if(chMap==NULL) cout << "null chMap"<<endl;
    
    // SETUP TGRAPHERRORS //
    TGraphErrors* graph = new TGraphErrors();
    map <string, double> detCoord; // <SN,xCoordForGraph>
    detCoord.insert(pair<string,double>("B8470",1.));
    detCoord.insert(pair<string,double>("P42665B",2.));
    detCoord.insert(pair<string,double>("B8487",3.));
    detCoord.insert(pair<string,double>("P42538A",4.));
    detCoord.insert(pair<string,double>("P42575A",5.));
    detCoord.insert(pair<string,double>("P42661C",6.));
    detCoord.insert(pair<string,double>("P42698B",7.));
    detCoord.insert(pair<string,double>("P42574A",8.));
    detCoord.insert(pair<string,double>("P42661B",9.));
    detCoord.insert(pair<string,double>("P42573B",10.));
    detCoord.insert(pair<string,double>("P42661A",11.));
    detCoord.insert(pair<string,double>("P42573A",12.));
    detCoord.insert(pair<string,double>("P42698A",13.));
    detCoord.insert(pair<string,double>("P42538B",14.));
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
    //for(unsigned int i = 0; i < binLabelVec.size(); i++)
    //{
    //    graph->SetPoint(i,detCoord.find(chMap->GetDetectorName(channel))->second,qtTau);
    //    graph->SetPointError(point_i,0.0,0.0);
    //}

    
    // READ FILE AND FILL GRAPH //
    int point_i = 0;
    graph->SetPoint(point_i,1.,0.); graph->SetPointError(point_i,0.0,0.0); point_i++; // ensure full range
    graph->SetPoint(point_i,14.,0.); graph->SetPointError(point_i,0.0,0.0); point_i++; // ensure full range
    //int lowestpoint = binLabelVec.size();
    ifstream infile;
    infilename = infilepath + "OptimalFWHMvsTau_DS"+to_string(DS)+".txt";
    infile.open(infilename);
    while (infile >> channel >> optTau >> eoptTau >> optFWHM >> eoptFWHM >> qtTau)
    {
        //cout<<channel<<" "<<optTau<<" "<<eoptTau<<" "<<optFWHM<<" "<<eoptFWHM<<" "<<qtTau<<endl;
        if(channel%2==0)
        {
            // set the point's x-value according to the detSN in order to match STC plot...
            if(detCoord.find(chMap->GetDetectorName(channel)) != detCoord.end())
            {
                cout<<"   "<<chMap->GetDetectorName(channel)<<" "<<detCoord.find(chMap->GetDetectorName(channel))->second<<" point_i="<<point_i<<endl;
                graph->SetPoint(point_i,detCoord.find(chMap->GetDetectorName(channel))->second,qtTau);
                graph->SetPointError(point_i,0.0,0.0);
                point_i++;
                //if(detCoord.find(chMap->GetDetectorName(channel))->second < lowestpoint){lowestpoint = detCoord.find(chMap->GetDetectorName(channel))->second;}
            }
        }
    }
    infile.close();
    
    // DRAW GRAPH
    TCanvas* canv = new TCanvas();
    canv->cd();
    graph->Draw("AP");
    graph->SetTitle("STC Detectors Inferred #tau_{QT}");
    graph->GetYaxis()->SetTitle("Cnt/kg/dy");
    graph->SetMarkerColor(kBlack);
    graph->SetMarkerStyle(20);
    TAxis* a = new TAxis();
    a = graph->GetXaxis();
    int j = 1; // =lowestpoint; // this prevents SetBinLabel() from trying to set a label on an out-of-range point whose bin is the underflow bin
    int bin_j = 0;
    //cout<<"TEST binLabelVec.size() = "<<binLabelVec.size()<<endl;
    //cout<<"TEST n points in graph = "<<point_i<<endl;
    //cout<<"TEST a->GetXmax() = "<<a->GetXmax()<<" ... a->GetXmax()-1 = "<<a->GetXmax()-1<<endl;
    while (j<((a->GetXmax())-1))
    {
        //cout<<"TEST j = "<<j<<endl;
        bin_j = a->FindBin(j);
        //cout<<"TEST bin_j = "<<bin_j<<endl;
        //cout<<"TEST j-1 = "<<j-1<<" ... < "<<binLabelVec.size()<<" ?"<<endl;
        a->SetBinLabel(bin_j,binLabelVec[j-1].c_str());
        j++;
    }
    canv->Modified();
    canv->Update();
    
    // WRITE GRAPH TO FILE //
    TFile* f = new TFile("STCvDS_tauQT_graphs.root","UPDATE");
    f->cd();
    graph->Write(Form("DS%d_Detectors_tauQT",DS));
    f->Close();
    

    

///////////////////////////////
//// Close out
///////////////////////////////

    cout << "F I N I S H E D" << endl;
    App->Run();
}



/*
// SETUP TGRAPHERRORS //
TGraphErrors* graph = new TGraphErrors();
map <string, double> detCoord; // <SN,xCoordForGraph>
detCoord.insert(pair<string,double>("B8455",1.));
detCoord.insert(pair<string,double>("B8470",2.));
detCoord.insert(pair<string,double>("B8465",3.));
detCoord.insert(pair<string,double>("P42665B",4.));
detCoord.insert(pair<string,double>("B8487",5.));
detCoord.insert(pair<string,double>("P42538A",6.));
detCoord.insert(pair<string,double>("P42575A",7.));
detCoord.insert(pair<string,double>("P42661C",8.));
detCoord.insert(pair<string,double>("P42698B",9.));
detCoord.insert(pair<string,double>("P42662A",10.));
detCoord.insert(pair<string,double>("P42574A",11.));
detCoord.insert(pair<string,double>("P42661B",12.));
detCoord.insert(pair<string,double>("P42573B",13.));
detCoord.insert(pair<string,double>("P42661A",14.));
detCoord.insert(pair<string,double>("P42573A",15.));
detCoord.insert(pair<string,double>("P42698A",16.));
detCoord.insert(pair<string,double>("P42538B",17.));
vector<string> binLabelVec;
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
 
 graph->SetPoint(point_i,1.,0.); graph->SetPointError(point_i,0.0,0.0); point_i++; // ensure full range
 graph->SetPoint(point_i,17.,0.); graph->SetPointError(point_i,0.0,0.0); point_i++; // ensure full range
*/



/*
 int j = 1;
 int bin_j = 0;
 //cout<<"TEST graphNameVec.size() = "<<graphNameVec.size()<<endl;
 //cout<<"TEST n points in graph = "<<point_i + 1<<endl;
 //cout<<"TEST a->GetXmax() = "<<a->GetXmax()<<" ... a->GetXmax()-1 = "<<a->GetXmax()-1<<endl;
 while (j<((a->GetXmax())-1))
 {
 //cout<<"TEST j = "<<j<<endl;
 bin_j = a->FindBin(j);
 //cout<<"TEST bin_j = "<<bin_j<<endl;
 //cout<<"TEST j-1 = "<<j-1<<" ... < "<<graphNameVec.size()<<" ?"<<endl;
 a->SetBinLabel(bin_j,binLabelVec[j-1].c_str());
 j++;
 }
*/