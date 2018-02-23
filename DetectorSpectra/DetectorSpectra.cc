///////////////////////////////
// Makes individual detector spectra for a variety of cuts.
// Input desired DS, module, cut scheme. Corresponding spectra are produced for each detector, scaled by exposure, and saved in a ROOT file.
// Makes use of MJD skim data, MJTChannelMap objects from MJD data files, and exposure data from the Livetime Unidoc.
//
// T. Gilliss 2017-11-15
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
        cout << "Usage is: ./DetectorSpectra DS c cutScheme" << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // Dataset [0,1,2,3,4,5,51(5a),52(5b),...]
    int c = atoi(argv[2]); // Module [1,2]
    int cutScheme = atoi(argv[3]); // Cut Scheme [raw,raw+DC,raw+DC+AE,raw+DC+AE+DCR]
    
    // HARD CODED (FOR NOW) (NEED TO RECOMPILE IF CHANGED)
    int gain = 0; //[0=HG,1=LG]
    string skimVersion = "GAT-v01-07";
    string inFileBasePath = "/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/";
    int baseDS = DS; // a hack for sub-datasets
    if(DS==51 || DS==52) baseDS = 5; // the DS5a,b run ranges are hard coded as cuts below
    int binLo = 0, binHi = 3000; // keV plot bounds
    double keVPerBin = 0.5;
    
    // PRINT SETTINGS
    cout<<"DS"<<DS<<", c"<<c<<endl;
    cout<<"cutScheme "<<cutScheme<<", gain "<<gain<<endl;
    cout<<"skimVersion "<<skimVersion<<endl;
    
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Setup MJD data and cuts
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    TChain* ch;
    ch = new TChain("skimTree", "skimTree");
    ch->Add(Form("%s/DS%d/%s/skimDS%d*.root",inFileBasePath.c_str(),baseDS,skimVersion.c_str(),baseDS));
    
    // PREPARE CUTS
    string eEst = "trapENFCal"; // energy estimator
    char rawCut[1000];
    if(gain==0){sprintf(rawCut,"channel%%2==%d && mH==1",gain);}
    if(gain==1){sprintf(rawCut,"channel%%2==%d && mL==1",gain);}
    char dcCut[1000]; sprintf(dcCut,"isGood && !(isLNFill1 && C==1) && !(isLNFill2 && C==2) && !muVeto && !wfDCBits");
    char aeCut[1000]; sprintf(aeCut,"avse>-1");
    char dcrCut[1000]; sprintf(dcrCut,"dcr99<0");
    char ds5a[1000]; sprintf(ds5a,"run >= 18623 && run <= 22392");
    char ds5b[1000]; sprintf(ds5b,"run >= 22393 && run <= 23958");
    // ASSEMBLE CUTS (cutSchemes [0,1,2,3] = [raw,raw+DC,raw+DC+AE,raw+DC+AE+DCR])
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
    char finalCut[1000]; // used later when C,P,D get appended

    // OUTPUT ROOT FILE
    char outFileName[50];
    sprintf(outFileName,"output_data/DS%d_M%d_CutScheme%d.root",DS,c,cutScheme);
    TFile *outFile = new TFile(outFileName,"RECREATE");
    
    // PRINT SETTINGS
    cout<<"------------"<<endl;
    cout<<"Settings:"<<endl;
    cout<<" Plotting with "<< eEst << " over DS" << DS << " for M"<<c<< endl;
    cout<<" Gain: " << gain << " (0=Hi,1=Lo)" << endl;
    cout<<" Cuts: " << tempCut << " (and C,P,D)" << endl;
    cout<<" Output File: " << outFileName << endl;
    cout<<"------------"<<endl;
    
    // PLOT SETTINGS
    int nBins = (binHi-binLo)/keVPerBin;
    char canvName[50];
    char histName[50];
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Read in mapping info
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
        /* // debugging
        cout << inChan <<" "<< inDetID <<" "<< inAM <<" "<< inRuntime <<" "<< inLivetime <<" "<< inLTExpo <<" "<< inLTExpoUnc <<" "<< inAvgLTFrac <<" "<< inLTUnc <<" "<< inNRuns << endl;
        */
        channelExposureStruct tmpStruct = {inAM,inLivetime,inLTExpo,inLTExpoUnc};
        channelExposureMap.insert(std::pair<int,channelExposureStruct>(inChan,tmpStruct));
        chanList.push_back(inChan);
    }
    inTextFile.close();
    sort(chanList.begin(), chanList.end()); // list must be sorted for binary_search
    
    /* // debugging
    // CHECK STL MAP
    for(unsigned int i = 0; i<chanList.size(); i++)
    {
        int chtmp = chanList[i];
        cout<< chtmp <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).AM <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).LT <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EX <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EXU <<endl;
    }
    */
     
    // PREP CHANNELMAP
    TFile *mapFile = new TFile("channel_map_data/DSChannelMaps.root","READ");
    MJTChannelMap *chMap = (MJTChannelMap*) mapFile->Get(Form("ChannelMapDS%d",baseDS));
    
    // LOOP CHANNELS AND MAP DATA
    int maxMod = 2, maxPos = 7, maxDet = 5;
    int chtmp = 0;
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
                    /* // debugging
                    cout<<"Including "<<chtmp<<endl;
                    cout<<"    "<< chtmp <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).AM <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).LT <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EX <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EXU <<endl;
                    */
                    
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
                        hE->Scale(1/channelExposureMap.at(chtmp).EX); // kg*dy
                        outFile->cd();
                        hE->Write();
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
    delete mapFile; mapFile = 0;
    delete ch; ch = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Closeout
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "F I N I S H E D" << endl;
}