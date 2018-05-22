///////////////////////////////
//
// Script that makes plots for each detector. Calculates count rate in a region. Saves plots to file. Saves info
// and rates to file. This script will pull detector-by-detector exposures from the livetime unidoc.
//
// NOTE: This script relies on the livetime unidoc's list of channels.
// That is, it includes/excludes detectors based on the channels included/excluded in the livetime unidoc's list.
//
// Refs:
//  Active masses
//  -https://mjwiki.npl.washington.edu/pub/Majorana/AnalysisReports/ActiveMassCalc.pdf
//  DS0-5b channel selection and exposure info
//  -https://mjdoc.npl.washington.edu/record/1863/files/livetime_DS0-5_unidoc_2017_11_5.pdf
//  DS5c-6 channel selection and exposure info
//  -https://mjcalendar.npl.washington.edu/indico/event/2660/session/7/contribution/30/material/slides/0.pdf
//  -
//
///////////////////////////////

#include <iostream> // i/o stream: cout
#include <fstream> // file i/o: ifstream, ofstream
#include <cstdio> // sprintf
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <algorithm> // binary_search, sort
#include <iomanip> // setprecision

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

string GetCut(int ds, int cutScheme, int c = 0, int p = 0, int d = 0)
{
    string cutEnergy = "trapENFCal>0.0 && trapENFCal<10000.0"; // 1000-1400 is 2nbb ROI for this analysis
    string cutDC = " && isGood && !(isLNFill1 && C==1) && !(isLNFill2 && C==2) && !muVeto && !wfDCBits";
    string cutHG = " && gain==0";
    string cutEnr = " && isEnr";
    string cutNat = " && isNat";
    string cutMult = " && mH==1"; //string cutMult_inv = " && mH>1";
    string cutAvsE = " && avse>-1"; //string cutAvsE_inv = " && avse<=-1";
    string cutDCR = " && dcr99<0"; //string cutDCR_inv = " && dcr99>=0";

    map<int,string> cutSchemeMap = { // { cutScheme, cutString } // cutEnergy needs to be first; see above
        { 0, cutEnergy + cutHG + cutMult }, // RAW
        { 1, cutEnergy + cutHG + cutMult + cutDC }, // RAW+DC
        { 2, cutEnergy + cutHG + cutMult + cutDC + cutAvsE }, // RAW+DC+AE
        { 3, cutEnergy + cutHG + cutMult + cutDC + cutAvsE + cutDCR } // RAW+DC+AE+DCR
    };

    string finalCut = cutSchemeMap[cutScheme];
    if(ds == 51) {finalCut = finalCut + " && run >= 18623 && run <= 22392";} // DS5a
    if(ds == 52) {finalCut = finalCut + " && run >= 22393 && run <= 23958";} // DS5b

    if(c != 0 && p != 0 && d != 0)
    {
      finalCut = finalCut + " && C==" + to_string(c) + " && P==" + to_string(p) + " && D==" + to_string(d);
    }
    return finalCut;
}

int GetTempDS(int ds)
{
    int dsTemp;
    if(ds == 51 || ds == 52 || ds == 53) {dsTemp = 5;}
    else {dsTemp = ds;}
    return dsTemp;
}

string GetTempDSStr(int ds)
{
    int dsTemp = GetTempDS(ds);
    string dsTempStr;
    if(ds == 53) {dsTempStr = "5c";}
    else {dsTempStr = to_string(dsTemp);}
    return dsTempStr;
}

string GetGatRevStr(int ds, string fileType)
{
    string gatRevStr;
    if(fileType == "open")
    {
      if(ds == 0) {gatRevStr = "GAT-v01-07";} // updated 2017-09
      if(ds == 1) {gatRevStr = "GAT-v01-07";} // updated 2017-09
      if(ds == 2) {gatRevStr = "GAT-v01-07";} // updated 2017-09
      if(ds == 3) {gatRevStr = "GAT-v01-07";} // updated 2017-09
      if(ds == 4) {gatRevStr = "GAT-v01-07";} // updated 2017-09
      if(ds == 51) {gatRevStr = "GAT-v01-07";} // updated 2017-09
      if(ds == 52) {gatRevStr = "GAT-v01-07";} // updated 2017-09
      if(ds == 53) {gatRevStr = "GAT-v02-00-51-g69c5025";} // updated 2018-05-16
      if(ds == 6) {gatRevStr = "GAT-v02-00-66-gf078278";} // updated 2018-05-16
    }
    if(fileType == "blind")
    {
      if(ds == 1) {gatRevStr = "GAT-v02-00-72-gc94f3a1";} // updated 2018-05-17 1421 ET
      if(ds == 2) {gatRevStr = "GAT-v02-00-66-gf078278";} // updated 2018-05-16
      if(ds == 53) {gatRevStr = "GAT-v02-00-66-gf078278";} // updated 2018-05-17 1421 ET
      if(ds == 6) {gatRevStr = "GAT-v02-00-74-g794b9f8";} // updated 2018-05017 1914 ET
    }
    return gatRevStr;
}

string GetSkimPath(int ds, string fileType)
{
    string skimPath;
    if(fileType == "open")
    {
      skimPath = "/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS" + GetTempDSStr(ds) + "/" + GetGatRevStr(ds,fileType) + "/skimDS" + to_string(GetTempDS(ds)) + "*.root";
    }
    if(fileType == "blind")
    {
      skimPath = "/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS" + GetTempDSStr(ds) + "_blind/" + GetGatRevStr(ds,fileType) + "/skimDS" + to_string(GetTempDS(ds)) + "*.root";
    }
    return skimPath;
}

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
        cout << "Usage is: ./2nuBB_Systematics DS c cutScheme" << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // Dataset [0,1,2,3,4,51,52,53,6]
    int c = atoi(argv[2]); // Module [1,2]
    int cutScheme = atoi(argv[3]); // [0,1,2,3] :: [raw,raw+DC,raw+DC+AE,raw+DC+AE+DCR]
    string fileType = "open";

    // DEBUG MODE
    bool debug = false;

    // PLOT SETTINGS
    int xlow = 0, xup = 10000; // int xlow = 1000, xup = 1400;
    double keVPerBin = 1.0; // 0.5;
    int nbinsx = (xup-xlow)/keVPerBin;
    string eEst = "trapENFCal";
    int xlowROI = 1000, xupROI = 1400;

    // OUTPUT ROOT FILE
    char outFileName[50];
    sprintf(outFileName,"DS%d_M%d_CutScheme%d.root",DS,c,cutScheme);
    TFile *outFile = new TFile(outFileName,"RECREATE");

    // INPUT TEXT FILE // These files contain only good HG channels
    char inTextFilePath[50];
    if(fileType == "open")
    {
      sprintf(inTextFilePath,"exposure_data/DS%d_exposure_data.txt",DS);
    }
    if(fileType == "blind")
    {
      sprintf(inTextFilePath,"exposure_data/DS%d_blind_exposure_data.txt",DS);
    }
    ifstream inTextFile(inTextFilePath);

    //OUTPUT TEXT FILE
    char outTextFileName[50];
    sprintf(outTextFileName,"DS%d_M%d_CutScheme%d.txt",DS,c,cutScheme);
    ofstream outTextFile;
    outTextFile.open (outTextFileName);

    // PRINT SETTINGS
    cout<<"------------"<<endl;
    cout<<"Settings:"<<endl;
    cout<<" Plotting DS" << DS << " for M"<<c<< endl;
    cout<<" Gain: gain==0 (0=Hi,1=Lo)" << endl;
    cout<<" Cuts: " << GetCut(DS, cutScheme) << " (and C,P,D)" << endl;
    cout<<" skim files: " << GetSkimPath(DS,fileType) << endl;
    cout<<" outFiles: " << outFileName << " " << outTextFileName << endl;
    cout<<" inTextFile: " << inTextFilePath << endl;
    cout<<" fileType: " << fileType << endl;
    cout<<"------------"<<endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Setup MJD data
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // PREP TCHAIN of SKIM FILES
    TChain* ch;
    ch = new TChain("skimTree", "skimTree");
    ch->Add(GetSkimPath(DS,fileType).c_str());

    // PREP CHANNELMAP
    TFile *chMapFile = new TFile("channel_map_data/DSChannelMaps.root","READ");
    MJTChannelMap *chMap = (MJTChannelMap*) chMapFile->Get(Form("ChannelMapDS%d",GetTempDS(DS)));

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Read and Map Data
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // PREP VARS FROM INPUT TEXT FILE AND STL MAP
    int inChan, inDetID, inNRuns;
    double inAM, inRuntime, inLivetime, inLTExpo, inLTExpoUnc, inAvgLTFrac, inLTUnc;
    map<int, channelExposureStruct> channelExposureMap;
    vector<int> chanList; // can use this list to check for included/excluded channels

    // LOOP INPUT TEXT FILE AND FILL STL MAP
    while(inTextFile >> inChan >> inDetID >> setprecision(10)>>inAM >> inRuntime >> setprecision(10)>>inLivetime
     >> setprecision(10)>>inLTExpo >> setprecision(10)>>inLTExpoUnc >> inAvgLTFrac >> inLTUnc >> inNRuns)
    {
        if(debug == true)
        {
          cout << inChan <<" "<< inDetID <<" "<< inAM <<" "<< inRuntime <<" "<< inLivetime <<" "<< inLTExpo <<" "<< inLTExpoUnc <<" "<< inAvgLTFrac <<" "<< inLTUnc <<" "<< inNRuns << endl;
        }
        channelExposureStruct tmpStruct = {inAM,inLivetime,inLTExpo,inLTExpoUnc};
        channelExposureMap.insert(std::pair<int,channelExposureStruct>(inChan,tmpStruct));
        chanList.push_back(inChan);
    }
    inTextFile.close();
    sort(chanList.begin(), chanList.end()); // list must be sorted for binary_search

    // CHECK STL MAP
    if(debug == true)
    {
      for(unsigned int i = 0; i<chanList.size(); i++)
      {
        int chtmp = chanList[i];
        cout<< chtmp <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).AM <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).LT <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EX <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EXU <<endl;
      }
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Do the work
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // LOOP OVER ALL DETECTORS, OPERATE ON THE RELEVANT ONES
    int maxMod = 2, maxPos = 7, maxDet = 5;
    int chtmp = 0;
    int hIntegral = 0;
    double hIntegralNormed = 0;
    int detType = 1000;
    char canvName[50];
    char histName[50];
    string finalCut = "";
    if(c>maxMod) {cout<<"Error: Bad module number"<<endl; return 1;}

    // Setup hTotEnr
    sprintf(canvName,"c_Enr_CutScheme%d", cutScheme);
    sprintf(histName,"h_Enr_CutScheme%d", cutScheme);
    TH1D *hTotEnr = new TH1D(histName,histName,nbinsx,xlow,xup);
    hTotEnr->Sumw2();
    hTotEnr->GetYaxis()->SetTitle(Form("Cnt/%.2f keV",keVPerBin)); // Form("Cnt/kg/dy/%.2f keV",keVPerBin)
    hTotEnr->GetXaxis()->SetTitle(Form("%s (keV)",eEst.c_str()));
    hTotEnr->SetMarkerColor(1);

    // Main loop
    for(int p = 1; p <= maxPos; p++)
    {
        for(int d = 1; d <= maxDet; d++)
        {
            if(chMap->GetDetectorName(c,p,d)!="")
            {
                // GET CHANNEL NUMBER AND CHECK IF PRESENT IN LIVETIME UNIDOC's GOOD LIST. IF SO, PROCEED
                chtmp = chMap->GetInt(chMap->GetDetectorName(c,p,d),"kIDHi");
                detType = chMap->GetInt(chMap->GetDetectorName(c,p,d),"kDetectorType");
                if(binary_search(chanList.begin(), chanList.end(), chtmp))
                {
                    if(debug == true)
                    {
                      cout<<"Including "<<chtmp<<endl;
                      cout<<"    "<< chtmp <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).AM <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).LT <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EX <<" "<< setprecision(10)<<channelExposureMap.at(chtmp).EXU <<endl;
                    }

                    // PREP HIST
                    sprintf(canvName,"c_C%dP%dD%d_CutScheme%d", c, p, d, cutScheme);
                    sprintf(histName,"h_C%dP%dD%d_CutScheme%d", c, p, d, cutScheme);
                    TCanvas *cE = new TCanvas(canvName,canvName);
                    cE->cd();
                    TH1D *hE = new TH1D(histName,histName,nbinsx,xlow,xup);
                    hE->Sumw2();
                    hE->GetYaxis()->SetTitle(Form("Cnt/%.2f keV",keVPerBin)); // Form("Cnt/kg/dy/%.2f keV",keVPerBin)
                    hE->GetXaxis()->SetTitle(Form("%s (keV)",eEst.c_str()));
                    hE->SetMarkerColor(1);
                    hE->SetStats(0);
                    finalCut = GetCut(DS, cutScheme, c, p, d);

                    // DRAW HIST
                    ch->Draw(Form("%s >> %s",eEst.c_str(),histName),finalCut.c_str(),"P E0"); // plot points, and errors are sqrt(cnts)
                    if(hE->GetEntries() == 0) {cout<<"hE->GetEntries() == 0 ... "<<"c"<<c<<"p"<<p<<"d"<<d<<endl;} // may indicate !isGood detector?
                    hIntegral = hE->Integral(hE->GetXaxis()->FindBin(xlowROI),hE->GetXaxis()->FindBin(xupROI)-1); // -1 b/c upper-edge xup is EXCLUDED (FindBin(xup) returns overflow bin)
                    // hE->Scale(1/channelExposureMap.at(chtmp).EX); // kg*dy
                    // hIntegralNormed = hE->Integral(hE->GetXaxis()->FindBin(xlow),hE->GetXaxis()->FindBin(xup)-1); // -1 b/c upper-edge xup is EXCLUDED (FindBin(xup) returns overflow bin)
                    hIntegralNormed = hIntegral/channelExposureMap.at(chtmp).EX;

                    // ADD hE, for enr dets, into total hTotEnr
                    if(detType == 2) {hTotEnr->Add(hE);}

                    // WRITE USEFUL STUFF TO OUT FILES
                    outFile->cd();
                    hE->Write();
                    outTextFile << c << " " << p << " " << d << " " << hIntegral << " " << hIntegralNormed << " " << channelExposureMap.at(chtmp).EX << " " << channelExposureMap.at(chtmp).EXU << endl;

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

    // WRITE hTotEnr to root file
    hTotEnr->Write();

    // CLOSE ROOT AND TEXT FILES AND DELETE DYNAMICALLY ALLOCATED MEMORY
    outFile->Close();
    outTextFile.close();
    delete chMapFile; chMapFile = 0;
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
