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
    
    TFile *f = new TFile("OptimalResolutionData.root","RECREATE");
    TTree *t = new TTree("ntuple","data from ascii file");
    t->ReadFile("/Users/tomgilliss/Desktop/STCvM1_EnergyRes/ResolutionStudy_AllDets_STC/main/OptimalData_AllDets_SpaceSeparated.txt","detectorSN:systemBOOL:stringNumber:detectorNumber:tau:fwhm");
    t->Write();
    f->Close();
    
    cout << "Finished..." << endl;
    
    App->Run();
    
} // End of main()