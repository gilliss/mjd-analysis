///////////////////////////////
//July 2017 Update
//
// Script that makes plots for each detector. Calculates count rate in a region. Saves plots to file. Saves info
// and rates to file. This script will pull detector-by-detector exposures from the livetime unidoc.
//
// NOTE: This script is not robust against good detectors that may somehow show no rate. Here we are looking at the 2vBB region, so there should always be counts, but if we were looking at a sparse or narrow region this script would discount detectors that have zero counts. ... BUT I can avoid this susceptibility by relying on the livetime unidoc's list of channels. That is, I can include/exclude detectors based on the channels included/excluded in the livetime unidoc's list.
//
// STEPS IN BUILDING SCRIPT:
//      1. xGather unidoc source and DS/detector exposure, mass, livetime values
//          https://mjdoc.npl.washington.edu/record/1863/files/livetime_DS0-5_unidoc_2017_10_18b.pdf
//      2. xRead in those exposure values given DS number
//      3. xstore exposures in STL map, set precision (need 4 decimal places)
//      4. Work out routine for looping through channels (need channel map) and testing if in STL map
//      5. Figure out DS5a and DS5b run lists, check DSInfo
//      6. #include most recent mjsw
///////////////////////////////

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>    // std::binary_search, sort
#include <iomanip>      // std::setprecision


#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1D.h>

#include <MJTChannelMap.hh>

using namespace std;

struct channelExposureStruct {
    // STRUCT TO HOLD ACTIVE MASS, LIVETIME, EXPOSURE FOR EACH CHANNEL
    double AM, LT, EX, EXU;
};

int main (int argc, char* argv[])
{
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Command line arguments and initializations
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 4)
    {
        cout << "Error: Too few arguments " << argv[0] << endl;
        cout << "Usage is: ./dev_2nuBB_Systematics DS c cutScheme" << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // Dataset
    int c = atoi(argv[2]); // Module
    int cutScheme = atoi(argv[3]); // [raw,raw+DC,raw+DC+AE,raw+DC+AE+DCR]
    int gain = 0; //[0=HG,1=LG]
    cout<<"DS"<<DS<<", c"<<c<<endl;
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Setup MJD data
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    TChain* ch;
    if(DS==0) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS0/GAT-v01-07/skimDS0*.root");
    }
    if(DS==1) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS1/GAT-v01-07/skimDS1*.root");
    }
    if(DS==2) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS2/GAT-v01-07/skimDS2*.root");
    }
    if(DS==3) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS3/GAT-v01-07/skimDS3*.root");
    }
    if(DS==4) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS4/GAT-v01-07/skimDS4*.root");
    }
    if(DS==5) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS5/GAT-v01-07/skimDS5*.root");
    }
    if(DS==51 /*5a*/) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS5/GAT-v01-07/skimDS5*.root");
    }
    if(DS==52 /*5b*/) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS5/GAT-v01-07/skimDS5*.root");
    }
    
    // MAKE FULL CUT [raw,raw+DC,raw+DC+AE,raw+DC+AE+DCR] cuts could be made individ for dets, or for DSs above
    string eEst = "trapENFCal"; // energy estimator
    char rawCut[1000];
    if(gain==0){sprintf(rawCut,"channel%%2==%d && mH==1",gain);}
    if(gain==1){sprintf(rawCut,"channel%%2==%d && mL==1",gain);}
    char dcCut[1000]; sprintf(dcCut,"isGood && !(isLNFill1 && C==1) && !(isLNFill2 && C==2) && !muVeto && !wfDCBits");
    char aeCut[1000]; sprintf(aeCut,"avse>-1");
    char dcrCut[1000]; sprintf(dcrCut,"dcr99<0");
    char ds5a[1000]; sprintf(ds5a,"run >= 18623 && run <= 22392");
    char ds5b[1000]; sprintf(ds5b,"run >= 22393 && run <= 23958");
    char tempCut[1000];
    if(cutScheme==0) sprintf(tempCut,"%s",rawCut);
    if(cutScheme==1) sprintf(tempCut,"%s && %s",rawCut,dcCut);
    if(cutScheme==2) sprintf(tempCut,"%s && %s && %s",rawCut,dcCut,aeCut);
    if(cutScheme==3) sprintf(tempCut,"%s && %s && %s && %s",rawCut,dcCut,aeCut,dcrCut);
    if(cutScheme==0 && DS==51) sprintf(tempCut,"%s && %s",rawCut,ds5a);
    if(cutScheme==1 && DS==51) sprintf(tempCut,"%s && %s && %s",rawCut,dcCut,ds5a);
    if(cutScheme==2 && DS==51) sprintf(tempCut,"%s && %s && %s && %s",rawCut,dcCut,aeCut,ds5a);
    if(cutScheme==3 && DS==51) sprintf(tempCut,"%s && %s && %s && %s && %s",rawCut,dcCut,aeCut,dcrCut,ds5a);
    if(cutScheme==0 && DS==52) sprintf(tempCut,"%s && %s",rawCut,ds5b);
    if(cutScheme==1 && DS==52) sprintf(tempCut,"%s && %s && %s",rawCut,dcCut,ds5b);
    if(cutScheme==2 && DS==52) sprintf(tempCut,"%s && %s && %s && %s",rawCut,dcCut,aeCut,ds5b);
    if(cutScheme==3 && DS==52) sprintf(tempCut,"%s && %s && %s && %s && %s",rawCut,dcCut,aeCut,dcrCut,ds5b);
    char finalCut[1000]; // for when C,P,D get appended

    // OUTPUT ROOT FILE
    char outFileName[50];
    sprintf(outFileName,"DS%d_M%d_CutScheme%d.root",DS,c,cutScheme);
    TFile *outFile = new TFile(outFileName,"RECREATE");
    
    //OUTPUT TEXT FILE
    char outTextFileName[50];
    sprintf(outTextFileName,"DS%d_M%d_CutScheme%d.txt",DS,c,cutScheme);
    ofstream outTextFile;
    outTextFile.open (outTextFileName);
    
    // PRINT SETTINGS
    cout<<"------------"<<endl;
    cout<<"Settings:"<<endl;
    cout<<" Plotting with "<< eEst << " over DS" << DS << " for M"<<c<< endl;
    cout<<" Gain: " << gain << " (0=Hi,1=Lo)" << endl;
    cout<<" Cuts: " << tempCut << " (and C,P,D)" << endl;
    cout<<" outFiles: " << outFileName << " " << outTextFileName << endl;
    cout<<"------------"<<endl;
    
    // PLOT SETTINGS
    int binLo = 1000, binHi = 1400;
    double keVPerBin = 0.5;
    int nBins = (binHi-binLo)/keVPerBin;
    char canvName[50];
    char histName[50];
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Read and Map Data
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // PREP INPUT TEXT FILE AND STL MAP
    char inTextFilePath[50];
    sprintf(inTextFilePath,"exposure_data/DS%d_exposure_data.txt",DS);
    ifstream inTextFile(inTextFilePath);
    int inChan, inDetID, inNRuns;
    double inAM, inRuntime, inLivetime, inLTExpo, inLTExpoUnc, inAvgLTFrac, inLTUnc;
    map<int, channelExposureStruct> channelExposureMap;
    vector<int> chanList; // can use this list to check for included/excluded channels
    
    // LOOP INPUT TEXT FILE AND FILL STL MAP
    while(inTextFile >> inChan >> inDetID >> setprecision(10)>>inAM >> inRuntime >> setprecision(10)>>inLivetime >> setprecision(10)>>inLTExpo >> setprecision(10)>>inLTExpoUnc >> inAvgLTFrac >> inLTUnc >> inNRuns)
    {
        //cout << inChan <<" "<< inDetID <<" "<< inAM <<" "<< inRuntime <<" "<< inLivetime <<" "<< inLTExpo <<" "<< inLTExpoUnc <<" "<< inAvgLTFrac <<" "<< inLTUnc <<" "<< inNRuns << endl;
        channelExposureStruct tmpStruct = {inAM,inLivetime,inLTExpo,inLTExpoUnc};
        channelExposureMap.insert(std::pair<int,channelExposureStruct>(inChan,tmpStruct));
        chanList.push_back(inChan);
    }
    sort(chanList.begin(), chanList.end()); // list must be sorted for binary_search
    
    /* 
    // CHECK STL MAP
    for(unsigned int i = 0; i<chanList.size(); i++)
    {
        int chtmp = chanList[i];
        cout<< chtmp <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).AM <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).LT <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EX <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EXU <<endl;
    }
    */
     
    // PREP CHANNELMAP
    int chMapDS = DS; if(DS==51 || DS==52){chMapDS=5;}
    TFile *mapFile = new TFile("channel_map_data/DSChannelMaps.root","READ");
    MJTChannelMap *chMap = (MJTChannelMap*) mapFile->Get(Form("ChannelMapDS%d",chMapDS));
    
    // LOOP CHANNELS AND MAP DATA
    int maxMod = 2, maxPos = 7, maxDet = 5;
    int chtmp = 0;
    int hIntegral = 0;
    double hIntegralNormed = 0;
    if(c>maxMod) {cout<<"Error: Bad module number"<<endl; return 1;}
    for(int p = 1; p <= maxPos; p++)
    {
        for(int d = 1; d <= maxDet; d++)
        {
            if(chMap->GetDetectorName(c,p,d)!="")
            {
                chtmp = chMap->GetInt(chMap->GetDetectorName(c,p,d),"kIDHi");
                if(binary_search(chanList.begin(), chanList.end(), chtmp))
                {
                    //cout<<"Including "<<chtmp<<endl;
                    //cout<<"    "<< chtmp <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).AM <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).LT <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EX <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EXU <<endl;
                    
                    sprintf(histName,"h_C%dP%dD%d_CutScheme%d",c,p,d,cutScheme);
                    sprintf(canvName,"c_C%dP%dD%d_CutScheme%d",c,p,d,cutScheme);
                    
                    // DRAW HIST
                    TCanvas *cE = new TCanvas(canvName,canvName);
                    cE->cd();
                    TH1D *hE = new TH1D(histName,histName,nBins,binLo,binHi);
                    hE->GetYaxis()->SetTitle(Form("Cnt/kg/dy/%.2f keV",keVPerBin));
                    hE->GetXaxis()->SetTitle(Form("%s (keV)",eEst.c_str()));
                    hE->SetTitle(histName);
                    hE->SetName(histName);
                    hE->SetMarkerColor(1);
                    hE->SetStats(0);
                    sprintf(finalCut,"C==%d && P==%d && D==%d && %s",c,p,d,tempCut);
                    ch->Draw(Form("%s >> %s",eEst.c_str(),histName),finalCut,"P E0"); // plot points, and errors are sqrt(cnts)
                    if(hE->GetEntries()>0) // this constraint weeds out the !isGood detector hists for which there are no entries; an isGood detector hist with no entries would be missed
                    {
                        hIntegral = hE->Integral(hE->GetXaxis()->FindBin(binLo),hE->GetXaxis()->FindBin(binHi)-1); // -1 b/c FindBin adds 1
                        hE->Scale(1/channelExposureMap.at(chtmp).EX); // kg*dy
                        hIntegralNormed = hE->Integral(hE->GetXaxis()->FindBin(binLo),hE->GetXaxis()->FindBin(binHi)-1);  // -1 b/c FindBin adds 1
                        outFile->cd();
                        hE->Write();
                        // write other useful stuff to a file here
                        outTextFile << c << " " << p << " " << d << " " << hIntegral << " " << hIntegralNormed << " " << channelExposureMap.at(chtmp).EX << " " << channelExposureMap.at(chtmp).EXU << endl;
                    }
                    else {cout<<"!!!No Entries: c"<<c<<"p"<<p<<"d"<<d<<endl;}
                    
                    // CLEAR DYNAMICALLY ALLOCATED MEMORY
                    delete hE;
                    hE = 0;
                    delete cE;
                    cE = 0;
                }
                else {cout<<"Channel not included: c"<<c<<"p"<<p<<"d"<<d<<", "<<chtmp<<endl;}
            } // end check kDetectorName string
            else {cout<<"Not in channel map: c"<<c<<"p"<<p<<"d"<<d<<endl;}
        } // d
    } // p
    
    outFile->Close();
    outTextFile.close();
    delete mapFile; mapFile = 0;
    delete ch; ch = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Closeout
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "F I N I S H E D" << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// ...
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////