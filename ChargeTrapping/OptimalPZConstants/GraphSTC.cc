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

using namespace std;

string GetDetSN(string stcName, int channel)
{
    //within this routine, open HW file and find the det SN based on the channel //
    string detSN="";
    ifstream stcHWMap;
    string dirPath = "/global/u2/g/gilliss/ChargeTrapping/Update_March2017/";
    string stcPath = "STCDataSets/";
    string stcHWMapName = dirPath + stcPath +stcName+ "_HWMap.txt";
    int chan=0;
    string SN="";
    double bias=0.;
    stcHWMap.open(stcHWMapName);
    while (stcHWMap >> chan >> SN >> bias)
    {
        if(chan==channel){detSN = SN; cout<<"detSNVec.push_back("<<detSN<<");"<<endl;}
    }
    stcHWMap.close();
    return detSN;
}

int main ()
{
///////////////////////////////
//// INITIALIZATIONS
///////////////////////////////
    TApplication *App = new TApplication("App", 0, NULL);
    // DEFINE THE STC NAMES IN SOME ORDER //
    vector<string> stcNameVec;
    stcNameVec.push_back("P3CLR");stcNameVec.push_back("P3DCR");stcNameVec.push_back("P3DNQ");stcNameVec.push_back("P3EVV");stcNameVec.push_back("P3FPV");stcNameVec.push_back("P3FPW");
    //for(unsigned int i = 0; i<stcNameVec.size(); i++) {cout<<stcNameVec[i]<<endl;}
    
    // IN FILE PARAMETERS //
    string infilepath = "/global/u2/g/gilliss/ChargeTrapping/Update_March2017/";
    string infilename = "";
    int channel=0;
    double optTau=0, eoptTau=0, optFWHM=0, eoptFWHM=0, qtTau=0;
    
    // SETUP TGRAPHERRORS //
    TGraphErrors* graph = new TGraphErrors();
    vector<string> binLabelVec;
    
    // READ FILE AND FILL GRAPH //
    int point_i = 0;
    for(unsigned int i_stc = 0; i_stc<stcNameVec.size(); i_stc++)
    {
        ifstream infile;
        infilename = infilepath + "OptimalFWHMvsTau_Stc"+stcNameVec[i_stc]+".txt";
        infile.open(infilename);
        //cout<<"Reading "<<stcNameVec[i_stc]<<endl;
        while (infile >> channel >> optTau >> eoptTau >> optFWHM >> eoptFWHM >> qtTau)
        {
            //cout<<channel<<" "<<optTau<<" "<<eoptTau<<" "<<optFWHM<<" "<<eoptFWHM<<" "<<qtTau<<endl;
            if(channel%2==0 && GetDetSN(stcNameVec[i_stc],channel)!="B8455" && GetDetSN(stcNameVec[i_stc],channel)!="B8465" && GetDetSN(stcNameVec[i_stc],channel)!="P42662A")
            {
                graph->SetPoint(point_i,double(point_i+1),qtTau);
                graph->SetPointError(point_i,0.0,0.0);
                binLabelVec.push_back(GetDetSN(stcNameVec[i_stc],channel)); //.push_back(to_string(channel));
                point_i++;
            }
        }
        infile.close();
    }
    
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
    canv->Modified();
    canv->Update();
    
    // WRITE GRAPH TO FILE //
    TFile* f = new TFile("STCvDS_tauQT_graphs.root","UPDATE");
    f->cd();
    graph->Write("STC_Detectors_tauQT");
    f->Close();
    

    

///////////////////////////////
//// Close out
///////////////////////////////

    cout << "F I N I S H E D" << endl;
    App->Run();
}