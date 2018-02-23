//
//  CreateRootFile_ResolutionData_AllDets.cc  Created by Tom Gilliss on 1/18/16.
//  
//
//  Create a root tree and file from resolution data in a text file.
//  Based on basic.c: https://root.cern.ch/root/html/tutorials/tree/basic.C.html
//

#include <TApplication.h>  //This class creates the ROOT Application Environment that interfaces to the windowing system eventloop and eventhandlers
#include "Riostream.h"
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>

#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

int main(int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);
    
    cout << "Running..." << endl;
    
    ifstream in;
    in.open(Form("/Users/tomgilliss/Desktop/STCvM1_EnergyRes/ResolutionStudy_AllDets_STC/main/Data_STCDets_SpaceSeparated.txt"));
    
    Char_t detectorSN[20]; // https://root.cern.ch/doc/master/classTString.html
    Int_t systemBOOL, stringNumber, detectorNumber, channel, fixedPtBOOL, tauIndex, tau, peak;
    Float_t fwhm;
    
    Int_t nlines = 0;
    TFile *f = new TFile("ResolutionData_AllDets.root","RECREATE");
    TNtuple *ntuple = new TNtuple("ntuple","data from ascii file",
        "systemBOOL:stringNumber:detectorNumber:channel:fixedPtBOOL:tauIndex:tau:peak:fwhm"); // TNtuple does not support detectorSN; will need to make a mapping from S#D# to detector
    
    while (1) {
        in >> detectorSN >> systemBOOL >> stringNumber >> detectorNumber >> channel >> fixedPtBOOL >> tauIndex >> tau >> peak >> fwhm;
        if (!in.good()) break;
        if (nlines < 5) printf("%s S%dD%d,%d,%d,%d,fwhm%5f\n",detectorSN,stringNumber,detectorNumber,tauIndex,tau,peak,fwhm);
        ntuple->Fill(systemBOOL,stringNumber,detectorNumber,channel,fixedPtBOOL,tauIndex,tau,peak,fwhm); // TNtuple does not support char* detectorSN; will need to make a mapping from S#D# to detector
        nlines++;
    }
    printf(" found %d points\n",nlines);
    
    in.close();
    
    f->Write();
    
    App->Run();
    
} // End of main()