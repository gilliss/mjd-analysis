/*
SpecAnMultiLine.hh
A class to measure activities of various spectral lines.
Particulars for each DS are hard coded in SpecAnMultiLine(), Cut(), Fill()
 
Notes: to simplify DS-specific code, like varied energy estimators, could use a preprocessor replace command. i.e. could write in the code DS1eEst and then have the preprocessor replace that with trapENFCal
Notes: it is awkward how the cuts are distributed amongst LoopData(), Cut(), and Fill(). They should be more localized for ease of use and understanding.
        Cut() has the main cuts. LoopData() has the HG cut (and optional mH cut). Fill() has the energy range cut.
DS1 Burst Cut: http://mjwiki.npl.washington.edu/pub/Majorana/AnalysisReports/Burst_cut_160630.pdf
*/

#ifndef SPECANMULTILINE_H // preprocessor: don't redefine header info if it's already been defined
#define SPECANMULTILINE_H

#include <string>
#include <iostream> // cout
#include <vector>
#include <stdlib.h>     /* exit(), EXIT_FAILURE */
#include <stdio.h> // FILE
#include <fstream> // file i/o
#include <map>

#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>
#include <TF1.h>
#include <TVectorD.h>

#include "/project/projectdirs/majorana/software/sl64/mjsw/mjsw201706Prod/GAT/Apps/DataSetInfo.hh"
//#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201706Prod/GAT/MJDAnalysis/GATChannelSelectionInfo.hh"
#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201706Prod/GAT/MGTEventProcessing/GATDetInfoProcessor.hh"
#include "/project/projectdirs/majorana/software/sl64/mjsw/mjsw201706Prod/MGDO/Majorana/MJTChannelMap.hh"

using namespace std;

/////////////////////////////////
// OTHER FUNCTIONS
/////////////////////////////////
int Get_DetectorStatus(int DS, int C, int P, int D)
{
    int detStatus=100; // int iIsNotGood = 0, iIsGood = 1, iDNE = -1;
    
    if(DS==0 && C==1 && P==1 && D==1) {detStatus=1;}
    if(DS==0 && C==1 && P==1 && D==2) {detStatus=1;}
    if(DS==0 && C==1 && P==1 && D==3) {detStatus=1;}
    if(DS==0 && C==1 && P==1 && D==4) {detStatus=1;}
    if(DS==0 && C==1 && P==1 && D==5) {detStatus=-1;}
    if(DS==0 && C==1 && P==2 && D==1) {detStatus=0;}
    if(DS==0 && C==1 && P==2 && D==2) {detStatus=1;}
    if(DS==0 && C==1 && P==2 && D==3) {detStatus=1;}
    if(DS==0 && C==1 && P==2 && D==4) {detStatus=0;}
    if(DS==0 && C==1 && P==2 && D==5) {detStatus=-1;}
    if(DS==0 && C==1 && P==3 && D==1) {detStatus=0;}
    if(DS==0 && C==1 && P==3 && D==2) {detStatus=0;}
    if(DS==0 && C==1 && P==3 && D==3) {detStatus=0;}
    if(DS==0 && C==1 && P==3 && D==4) {detStatus=1;}
    if(DS==0 && C==1 && P==3 && D==5) {detStatus=-1;}
    if(DS==0 && C==1 && P==4 && D==1) {detStatus=1;}
    if(DS==0 && C==1 && P==4 && D==2) {detStatus=1;}
    if(DS==0 && C==1 && P==4 && D==3) {detStatus=1;}
    if(DS==0 && C==1 && P==4 && D==4) {detStatus=1;}
    if(DS==0 && C==1 && P==4 && D==5) {detStatus=1;}
    if(DS==0 && C==1 && P==5 && D==1) {detStatus=1;}
    if(DS==0 && C==1 && P==5 && D==2) {detStatus=1;}
    if(DS==0 && C==1 && P==5 && D==3) {detStatus=1;}
    if(DS==0 && C==1 && P==5 && D==4) {detStatus=1;}
    if(DS==0 && C==1 && P==5 && D==5) {detStatus=-1;}
    if(DS==0 && C==1 && P==6 && D==1) {detStatus=0;}
    if(DS==0 && C==1 && P==6 && D==2) {detStatus=1;}
    if(DS==0 && C==1 && P==6 && D==3) {detStatus=1;}
    if(DS==0 && C==1 && P==6 && D==4) {detStatus=0;}
    if(DS==0 && C==1 && P==6 && D==5) {detStatus=-1;}
    if(DS==0 && C==1 && P==7 && D==1) {detStatus=1;}
    if(DS==0 && C==1 && P==7 && D==2) {detStatus=1;}
    if(DS==0 && C==1 && P==7 && D==3) {detStatus=1;}
    if(DS==0 && C==1 && P==7 && D==4) {detStatus=0;}
    if(DS==0 && C==1 && P==7 && D==5) {detStatus=-1;}
    
    if(DS==1 && C==1 && P==1 && D==1) {detStatus=0;}
    if(DS==1 && C==1 && P==1 && D==2) {detStatus=1;}
    if(DS==1 && C==1 && P==1 && D==3) {detStatus=1;}
    if(DS==1 && C==1 && P==1 && D==4) {detStatus=1;}
    if(DS==1 && C==1 && P==1 && D==5) {detStatus=-1;}
    if(DS==1 && C==1 && P==2 && D==1) {detStatus=1;}
    if(DS==1 && C==1 && P==2 && D==2) {detStatus=1;}
    if(DS==1 && C==1 && P==2 && D==3) {detStatus=1;}
    if(DS==1 && C==1 && P==2 && D==4) {detStatus=0;}
    if(DS==1 && C==1 && P==2 && D==5) {detStatus=-1;}
    if(DS==1 && C==1 && P==3 && D==1) {detStatus=0;}
    if(DS==1 && C==1 && P==3 && D==2) {detStatus=1;}
    if(DS==1 && C==1 && P==3 && D==3) {detStatus=1;}
    if(DS==1 && C==1 && P==3 && D==4) {detStatus=1;}
    if(DS==1 && C==1 && P==3 && D==5) {detStatus=-1;}
    if(DS==1 && C==1 && P==4 && D==1) {detStatus=0;}
    if(DS==1 && C==1 && P==4 && D==2) {detStatus=0;}
    if(DS==1 && C==1 && P==4 && D==3) {detStatus=0;}
    if(DS==1 && C==1 && P==4 && D==4) {detStatus=0;}
    if(DS==1 && C==1 && P==4 && D==5) {detStatus=0;}
    if(DS==1 && C==1 && P==5 && D==1) {detStatus=0;}
    if(DS==1 && C==1 && P==5 && D==2) {detStatus=0;}
    if(DS==1 && C==1 && P==5 && D==3) {detStatus=1;} // non-constant isGood status
    if(DS==1 && C==1 && P==5 && D==4) {detStatus=0;}
    if(DS==1 && C==1 && P==5 && D==5) {detStatus=-1;}
    if(DS==1 && C==1 && P==6 && D==1) {detStatus=1;}
    if(DS==1 && C==1 && P==6 && D==2) {detStatus=0;}
    if(DS==1 && C==1 && P==6 && D==3) {detStatus=1;}
    if(DS==1 && C==1 && P==6 && D==4) {detStatus=1;}
    if(DS==1 && C==1 && P==6 && D==5) {detStatus=-1;}
    if(DS==1 && C==1 && P==7 && D==1) {detStatus=1;}
    if(DS==1 && C==1 && P==7 && D==2) {detStatus=1;}
    if(DS==1 && C==1 && P==7 && D==3) {detStatus=1;}
    if(DS==1 && C==1 && P==7 && D==4) {detStatus=1;}
    if(DS==1 && C==1 && P==7 && D==5) {detStatus=-1;}
    
    if(DS==2 && C==1 && P==1 && D==1) {detStatus=0;}
    if(DS==2 && C==1 && P==1 && D==2) {detStatus=1;}
    if(DS==2 && C==1 && P==1 && D==3) {detStatus=1;}
    if(DS==2 && C==1 && P==1 && D==4) {detStatus=1;}
    if(DS==2 && C==1 && P==1 && D==5) {detStatus=-1;}
    if(DS==2 && C==1 && P==2 && D==1) {detStatus=1;}
    if(DS==2 && C==1 && P==2 && D==2) {detStatus=1;}
    if(DS==2 && C==1 && P==2 && D==3) {detStatus=1;}
    if(DS==2 && C==1 && P==2 && D==4) {detStatus=0;}
    if(DS==2 && C==1 && P==2 && D==5) {detStatus=-1;}
    if(DS==2 && C==1 && P==3 && D==1) {detStatus=0;}
    if(DS==2 && C==1 && P==3 && D==2) {detStatus=1;}
    if(DS==2 && C==1 && P==3 && D==3) {detStatus=1;}
    if(DS==2 && C==1 && P==3 && D==4) {detStatus=1;}
    if(DS==2 && C==1 && P==3 && D==5) {detStatus=-1;}
    if(DS==2 && C==1 && P==4 && D==1) {detStatus=0;}
    if(DS==2 && C==1 && P==4 && D==2) {detStatus=0;}
    if(DS==2 && C==1 && P==4 && D==3) {detStatus=0;}
    if(DS==2 && C==1 && P==4 && D==4) {detStatus=0;}
    if(DS==2 && C==1 && P==4 && D==5) {detStatus=0;}
    if(DS==2 && C==1 && P==5 && D==1) {detStatus=0;}
    if(DS==2 && C==1 && P==5 && D==2) {detStatus=0;}
    if(DS==2 && C==1 && P==5 && D==3) {detStatus=1;}
    if(DS==2 && C==1 && P==5 && D==4) {detStatus=0;}
    if(DS==2 && C==1 && P==5 && D==5) {detStatus=-1;}
    if(DS==2 && C==1 && P==6 && D==1) {detStatus=1;}
    if(DS==2 && C==1 && P==6 && D==2) {detStatus=0;}
    if(DS==2 && C==1 && P==6 && D==3) {detStatus=1;}
    if(DS==2 && C==1 && P==6 && D==4) {detStatus=1;}
    if(DS==2 && C==1 && P==6 && D==5) {detStatus=-1;}
    if(DS==2 && C==1 && P==7 && D==1) {detStatus=1;}
    if(DS==2 && C==1 && P==7 && D==2) {detStatus=1;}
    if(DS==2 && C==1 && P==7 && D==3) {detStatus=0;} // changed from 1 to 0 since last update
    if(DS==2 && C==1 && P==7 && D==4) {detStatus=1;}
    if(DS==2 && C==1 && P==7 && D==5) {detStatus=-1;}
    
    if(DS==3 && C==1 && P==1 && D==1) {detStatus=0;}
    if(DS==3 && C==1 && P==1 && D==2) {detStatus=1;}
    if(DS==3 && C==1 && P==1 && D==3) {detStatus=1;}
    if(DS==3 && C==1 && P==1 && D==4) {detStatus=1;}
    if(DS==3 && C==1 && P==1 && D==5) {detStatus=-1;}
    if(DS==3 && C==1 && P==2 && D==1) {detStatus=1;}
    if(DS==3 && C==1 && P==2 && D==2) {detStatus=1;}
    if(DS==3 && C==1 && P==2 && D==3) {detStatus=1;}
    if(DS==3 && C==1 && P==2 && D==4) {detStatus=0;}
    if(DS==3 && C==1 && P==2 && D==5) {detStatus=-1;}
    if(DS==3 && C==1 && P==3 && D==1) {detStatus=0;}
    if(DS==3 && C==1 && P==3 && D==2) {detStatus=1;}
    if(DS==3 && C==1 && P==3 && D==3) {detStatus=1;}
    if(DS==3 && C==1 && P==3 && D==4) {detStatus=1;}
    if(DS==3 && C==1 && P==3 && D==5) {detStatus=-1;}
    if(DS==3 && C==1 && P==4 && D==1) {detStatus=1;}
    if(DS==3 && C==1 && P==4 && D==2) {detStatus=0;}
    if(DS==3 && C==1 && P==4 && D==3) {detStatus=0;}
    if(DS==3 && C==1 && P==4 && D==4) {detStatus=1;}
    if(DS==3 && C==1 && P==4 && D==5) {detStatus=1;}
    if(DS==3 && C==1 && P==5 && D==1) {detStatus=0;}
    if(DS==3 && C==1 && P==5 && D==2) {detStatus=1;}
    if(DS==3 && C==1 && P==5 && D==3) {detStatus=1;}
    if(DS==3 && C==1 && P==5 && D==4) {detStatus=0;}
    if(DS==3 && C==1 && P==5 && D==5) {detStatus=-1;}
    if(DS==3 && C==1 && P==6 && D==1) {detStatus=1;}
    if(DS==3 && C==1 && P==6 && D==2) {detStatus=0;}
    if(DS==3 && C==1 && P==6 && D==3) {detStatus=1;}
    if(DS==3 && C==1 && P==6 && D==4) {detStatus=1;}
    if(DS==3 && C==1 && P==6 && D==5) {detStatus=-1;}
    if(DS==3 && C==1 && P==7 && D==1) {detStatus=1;}
    if(DS==3 && C==1 && P==7 && D==2) {detStatus=1;}
    if(DS==3 && C==1 && P==7 && D==3) {detStatus=1;}
    if(DS==3 && C==1 && P==7 && D==4) {detStatus=1;}
    if(DS==3 && C==1 && P==7 && D==5) {detStatus=-1;}
    
    if(DS==4 && C==2 && P==1 && D==1) {detStatus=0;}
    if(DS==4 && C==2 && P==1 && D==2) {detStatus=0;}
    if(DS==4 && C==2 && P==1 && D==3) {detStatus=0;}
    if(DS==4 && C==2 && P==1 && D==4) {detStatus=1;}
    if(DS==4 && C==2 && P==1 && D==5) {detStatus=-1;}
    if(DS==4 && C==2 && P==2 && D==1) {detStatus=1;}
    if(DS==4 && C==2 && P==2 && D==2) {detStatus=1;}
    if(DS==4 && C==2 && P==2 && D==3) {detStatus=1;}
    if(DS==4 && C==2 && P==2 && D==4) {detStatus=0;}
    if(DS==4 && C==2 && P==2 && D==5) {detStatus=0;}
    if(DS==4 && C==2 && P==3 && D==1) {detStatus=1;}
    if(DS==4 && C==2 && P==3 && D==2) {detStatus=1;}
    if(DS==4 && C==2 && P==3 && D==3) {detStatus=0;}
    if(DS==4 && C==2 && P==3 && D==4) {detStatus=-1;}
    if(DS==4 && C==2 && P==3 && D==5) {detStatus=-1;}
    if(DS==4 && C==2 && P==4 && D==1) {detStatus=1;}
    if(DS==4 && C==2 && P==4 && D==2) {detStatus=0;}
    if(DS==4 && C==2 && P==4 && D==3) {detStatus=0;}
    if(DS==4 && C==2 && P==4 && D==4) {detStatus=1;}
    if(DS==4 && C==2 && P==4 && D==5) {detStatus=0;}
    if(DS==4 && C==2 && P==5 && D==1) {detStatus=1;}
    if(DS==4 && C==2 && P==5 && D==2) {detStatus=0;}
    if(DS==4 && C==2 && P==5 && D==3) {detStatus=1;}
    if(DS==4 && C==2 && P==5 && D==4) {detStatus=0;}
    if(DS==4 && C==2 && P==5 && D==5) {detStatus=-1;}
    if(DS==4 && C==2 && P==6 && D==1) {detStatus=1;}
    if(DS==4 && C==2 && P==6 && D==2) {detStatus=1;}
    if(DS==4 && C==2 && P==6 && D==3) {detStatus=0;}
    if(DS==4 && C==2 && P==6 && D==4) {detStatus=0;}
    if(DS==4 && C==2 && P==6 && D==5) {detStatus=-1;}
    if(DS==4 && C==2 && P==7 && D==1) {detStatus=0;}
    if(DS==4 && C==2 && P==7 && D==2) {detStatus=0;}
    if(DS==4 && C==2 && P==7 && D==3) {detStatus=1;}
    if(DS==4 && C==2 && P==7 && D==4) {detStatus=1;}
    if(DS==4 && C==2 && P==7 && D==5) {detStatus=-1;}
    
    if(DS==5 && C==1 && P==1 && D==1) {detStatus=0;}
    if(DS==5 && C==1 && P==1 && D==2) {detStatus=1;}
    if(DS==5 && C==1 && P==1 && D==3) {detStatus=1;}
    if(DS==5 && C==1 && P==1 && D==4) {detStatus=1;}
    if(DS==5 && C==1 && P==1 && D==5) {detStatus=-1;}
    if(DS==5 && C==1 && P==2 && D==1) {detStatus=1;} // change from 0 to 1 since last update
    if(DS==5 && C==1 && P==2 && D==2) {detStatus=1;}
    if(DS==5 && C==1 && P==2 && D==3) {detStatus=1;}
    if(DS==5 && C==1 && P==2 && D==4) {detStatus=0;}
    if(DS==5 && C==1 && P==2 && D==5) {detStatus=-1;}
    if(DS==5 && C==1 && P==3 && D==1) {detStatus=0;}
    if(DS==5 && C==1 && P==3 && D==2) {detStatus=1;}
    if(DS==5 && C==1 && P==3 && D==3) {detStatus=1;}
    if(DS==5 && C==1 && P==3 && D==4) {detStatus=1;}
    if(DS==5 && C==1 && P==3 && D==5) {detStatus=-1;}
    if(DS==5 && C==1 && P==4 && D==1) {detStatus=1;}
    if(DS==5 && C==1 && P==4 && D==2) {detStatus=1;}
    if(DS==5 && C==1 && P==4 && D==3) {detStatus=1;}
    if(DS==5 && C==1 && P==4 && D==4) {detStatus=1;}
    if(DS==5 && C==1 && P==4 && D==5) {detStatus=1;}
    if(DS==5 && C==1 && P==5 && D==1) {detStatus=0;}
    if(DS==5 && C==1 && P==5 && D==2) {detStatus=1;}
    if(DS==5 && C==1 && P==5 && D==3) {detStatus=1;}
    if(DS==5 && C==1 && P==5 && D==4) {detStatus=0;}
    if(DS==5 && C==1 && P==5 && D==5) {detStatus=-1;}
    if(DS==5 && C==1 && P==6 && D==1) {detStatus=1;}
    if(DS==5 && C==1 && P==6 && D==2) {detStatus=0;}
    if(DS==5 && C==1 && P==6 && D==3) {detStatus=1;}
    if(DS==5 && C==1 && P==6 && D==4) {detStatus=1;}
    if(DS==5 && C==1 && P==6 && D==5) {detStatus=-1;}
    if(DS==5 && C==1 && P==7 && D==1) {detStatus=1;}
    if(DS==5 && C==1 && P==7 && D==2) {detStatus=1;}
    if(DS==5 && C==1 && P==7 && D==3) {detStatus=1;}
    if(DS==5 && C==1 && P==7 && D==4) {detStatus=1;}
    if(DS==5 && C==1 && P==7 && D==5) {detStatus=-1;}
    
    if(DS==5 && C==2 && P==1 && D==1) {detStatus=1;}
    if(DS==5 && C==2 && P==1 && D==2) {detStatus=0;}
    if(DS==5 && C==2 && P==1 && D==3) {detStatus=0;}
    if(DS==5 && C==2 && P==1 && D==4) {detStatus=1;}
    if(DS==5 && C==2 && P==1 && D==5) {detStatus=-1;}
    if(DS==5 && C==2 && P==2 && D==1) {detStatus=1;}
    if(DS==5 && C==2 && P==2 && D==2) {detStatus=1;}
    if(DS==5 && C==2 && P==2 && D==3) {detStatus=1;}
    if(DS==5 && C==2 && P==2 && D==4) {detStatus=0;}
    if(DS==5 && C==2 && P==2 && D==5) {detStatus=0;}
    if(DS==5 && C==2 && P==3 && D==1) {detStatus=1;}
    if(DS==5 && C==2 && P==3 && D==2) {detStatus=1;}
    if(DS==5 && C==2 && P==3 && D==3) {detStatus=0;}
    if(DS==5 && C==2 && P==3 && D==4) {detStatus=-1;}
    if(DS==5 && C==2 && P==3 && D==5) {detStatus=-1;}
    if(DS==5 && C==2 && P==4 && D==1) {detStatus=1;}
    if(DS==5 && C==2 && P==4 && D==2) {detStatus=1;}
    if(DS==5 && C==2 && P==4 && D==3) {detStatus=0;}
    if(DS==5 && C==2 && P==4 && D==4) {detStatus=1;}
    if(DS==5 && C==2 && P==4 && D==5) {detStatus=0;}
    if(DS==5 && C==2 && P==5 && D==1) {detStatus=1;}
    if(DS==5 && C==2 && P==5 && D==2) {detStatus=0;}
    if(DS==5 && C==2 && P==5 && D==3) {detStatus=1;}
    if(DS==5 && C==2 && P==5 && D==4) {detStatus=1;}
    if(DS==5 && C==2 && P==5 && D==5) {detStatus=-1;}
    if(DS==5 && C==2 && P==6 && D==1) {detStatus=0;}
    if(DS==5 && C==2 && P==6 && D==2) {detStatus=1;}
    if(DS==5 && C==2 && P==6 && D==3) {detStatus=0;}
    if(DS==5 && C==2 && P==6 && D==4) {detStatus=0;}
    if(DS==5 && C==2 && P==6 && D==5) {detStatus=-1;}
    if(DS==5 && C==2 && P==7 && D==1) {detStatus=0;}
    if(DS==5 && C==2 && P==7 && D==2) {detStatus=0;}
    if(DS==5 && C==2 && P==7 && D==3) {detStatus=1;}
    if(DS==5 && C==2 && P==7 && D==4) {detStatus=1;}
    if(DS==5 && C==2 && P==7 && D==5) {detStatus=-1;}
    
    if(DS==6 && C==1 && P==1 && D==1) {detStatus=0;}
    if(DS==6 && C==1 && P==1 && D==2) {detStatus=1;}
    if(DS==6 && C==1 && P==1 && D==3) {detStatus=1;}
    if(DS==6 && C==1 && P==1 && D==4) {detStatus=1;}
    if(DS==6 && C==1 && P==1 && D==5) {detStatus=-1;}
    if(DS==6 && C==1 && P==2 && D==1) {detStatus=0;}
    if(DS==6 && C==1 && P==2 && D==2) {detStatus=1;}
    if(DS==6 && C==1 && P==2 && D==3) {detStatus=1;}
    if(DS==6 && C==1 && P==2 && D==4) {detStatus=0;}
    if(DS==6 && C==1 && P==2 && D==5) {detStatus=-1;}
    if(DS==6 && C==1 && P==3 && D==1) {detStatus=0;}
    if(DS==6 && C==1 && P==3 && D==2) {detStatus=1;}
    if(DS==6 && C==1 && P==3 && D==3) {detStatus=1;}
    if(DS==6 && C==1 && P==3 && D==4) {detStatus=1;}
    if(DS==6 && C==1 && P==3 && D==5) {detStatus=-1;}
    if(DS==6 && C==1 && P==4 && D==1) {detStatus=1;}
    if(DS==6 && C==1 && P==4 && D==2) {detStatus=1;}
    if(DS==6 && C==1 && P==4 && D==3) {detStatus=1;}
    if(DS==6 && C==1 && P==4 && D==4) {detStatus=1;}
    if(DS==6 && C==1 && P==4 && D==5) {detStatus=1;}
    if(DS==6 && C==1 && P==5 && D==1) {detStatus=0;}
    if(DS==6 && C==1 && P==5 && D==2) {detStatus=1;}
    if(DS==6 && C==1 && P==5 && D==3) {detStatus=1;}
    if(DS==6 && C==1 && P==5 && D==4) {detStatus=0;}
    if(DS==6 && C==1 && P==5 && D==5) {detStatus=-1;}
    if(DS==6 && C==1 && P==6 && D==1) {detStatus=1;}
    if(DS==6 && C==1 && P==6 && D==2) {detStatus=0;}
    if(DS==6 && C==1 && P==6 && D==3) {detStatus=1;}
    if(DS==6 && C==1 && P==6 && D==4) {detStatus=1;}
    if(DS==6 && C==1 && P==6 && D==5) {detStatus=-1;}
    if(DS==6 && C==1 && P==7 && D==1) {detStatus=1;}
    if(DS==6 && C==1 && P==7 && D==2) {detStatus=1;}
    if(DS==6 && C==1 && P==7 && D==3) {detStatus=1;}
    if(DS==6 && C==1 && P==7 && D==4) {detStatus=1;}
    if(DS==6 && C==1 && P==7 && D==5) {detStatus=-1;}
    
    if(DS==6 && C==2 && P==1 && D==1) {detStatus=1;}
    if(DS==6 && C==2 && P==1 && D==2) {detStatus=0;}
    if(DS==6 && C==2 && P==1 && D==3) {detStatus=0;}
    if(DS==6 && C==2 && P==1 && D==4) {detStatus=1;}
    if(DS==6 && C==2 && P==1 && D==5) {detStatus=-1;}
    if(DS==6 && C==2 && P==2 && D==1) {detStatus=1;}
    if(DS==6 && C==2 && P==2 && D==2) {detStatus=1;}
    if(DS==6 && C==2 && P==2 && D==3) {detStatus=1;}
    if(DS==6 && C==2 && P==2 && D==4) {detStatus=0;}
    if(DS==6 && C==2 && P==2 && D==5) {detStatus=0;}
    if(DS==6 && C==2 && P==3 && D==1) {detStatus=1;}
    if(DS==6 && C==2 && P==3 && D==2) {detStatus=1;}
    if(DS==6 && C==2 && P==3 && D==3) {detStatus=0;}
    if(DS==6 && C==2 && P==3 && D==4) {detStatus=-1;}
    if(DS==6 && C==2 && P==3 && D==5) {detStatus=-1;}
    if(DS==6 && C==2 && P==4 && D==1) {detStatus=1;}
    if(DS==6 && C==2 && P==4 && D==2) {detStatus=1;}
    if(DS==6 && C==2 && P==4 && D==3) {detStatus=0;}
    if(DS==6 && C==2 && P==4 && D==4) {detStatus=1;}
    if(DS==6 && C==2 && P==4 && D==5) {detStatus=0;}
    if(DS==6 && C==2 && P==5 && D==1) {detStatus=1;}
    if(DS==6 && C==2 && P==5 && D==2) {detStatus=0;}
    if(DS==6 && C==2 && P==5 && D==3) {detStatus=1;}
    if(DS==6 && C==2 && P==5 && D==4) {detStatus=1;}
    if(DS==6 && C==2 && P==5 && D==5) {detStatus=-1;}
    if(DS==6 && C==2 && P==6 && D==1) {detStatus=0;}
    if(DS==6 && C==2 && P==6 && D==2) {detStatus=1;}
    if(DS==6 && C==2 && P==6 && D==3) {detStatus=0;}
    if(DS==6 && C==2 && P==6 && D==4) {detStatus=0;}
    if(DS==6 && C==2 && P==6 && D==5) {detStatus=-1;}
    if(DS==6 && C==2 && P==7 && D==1) {detStatus=0;}
    if(DS==6 && C==2 && P==7 && D==2) {detStatus=0;}
    if(DS==6 && C==2 && P==7 && D==3) {detStatus=1;}
    if(DS==6 && C==2 && P==7 && D==4) {detStatus=1;}
    if(DS==6 && C==2 && P==7 && D==5) {detStatus=-1;}
    
    if(detStatus!=-1 && detStatus!=0 && detStatus!=1)
    {
        cout<<"Error: bad detector status"<<endl;
        cout<<"  DS"<<DS<<" C"<<C<<"P"<<P<<"D"<<D<<"   detStatus="<<detStatus<<endl;
    }
    
    return detStatus;
}

// PREP FOR MASSES and stuff FROM DATASETINFO and ChannelSelectionInfo
vector <map<int,double> > vecActiveMassMapDS;
//int dsExampleRuns[7] = {2690,9639,15070,16932,60001034,23970,25800}; // from DataSetInfo run lists
//string chSelPath;
//vector<GATChannelSelectionInfo> ch_select;
GATDetInfoProcessor gp;
void PrepForDSIandCSI()
{
    for(int i_map_ds = 0; i_map_ds < 7; i_map_ds++)
    {
        vecActiveMassMapDS.push_back(LoadActiveMasses(i_map_ds));
        
//        if (FILE *file = fopen(chSelPath.c_str(), "r")) // check if file exists
//        {
//            fclose(file); // close FILE used to check if file existst
//            chSelPath = GetChannelSelectionPath(i_map_ds);
//            ch_select.push_back(GATChannelSelectionInfo(chSelPath, dsExampleRuns[i_map_ds]));
//        }
//        else {cout<<"DS no good for ChannelSelectionInfo ... DS"<<i_map_ds<<endl;}
    }
}

double Get_mAct_kg(int DS, int C, int P, int D)
{
    double mAct_g = 0.;
    TFile *mapFile = new TFile("./ChannelMaps/DSChannelMaps.root","READ");
    MJTChannelMap *chMap = (MJTChannelMap*) mapFile->Get(Form("ChannelMapDS%d",DS));
    int tempID = gp.GetDetIDFromName( chMap->GetString(C,P,D,"kDetectorName") );
    //int tempChan = ch_select[DS].GetChannelFromCPDG(C,P,D,0); // 0 for HG channel
    //int tempID = ch_select[DS].GetDetIDFromChannel(tempChan);
    mAct_g = vecActiveMassMapDS.at(DS)[tempID];
    return mAct_g/1000.;
}

double Get_liveTime_dy(int DS)
{
    double runTime_s = 0., startTime0 = 0.;
    GetDSRunAndStartTimes(DS,runTime_s,startTime0);
    double liveTime_dy = runTime_s/(60.*60.*24.);
    return liveTime_dy;
}

/////////////////////////////////
// CLASS
/////////////////////////////////
class SpecAnMultiLine {
private:
    /* isotopes, lines, results */
    struct speclines { // data structure to hold info for each isotope
        string isotope;
        vector<int> success;
        vector<TH1D*> henergy;
        /* sbLo...roiLo...energy...roiHi...sbHi */
        vector<double> energy;
        vector<double> roiLo; // peak integration bounds
        vector<double> roiHi;
        vector<double> sbLo; // sideband fit bounds
        vector<double> sbHi;
        vector<double> binLo; // henergy hist bounds
        vector<double> binHi;
        /* upper limit calculation */
        vector<int> roiCnt; // counts in ROI after cuts; from unbinned integral
        vector<double> bgCnt; // counts in BG sidebands from unbinned integral (used to be integral of BG TF1)
        vector<double> pkCnt; // BG subtracted roiCnt
        vector<double> ePkCnt; // 1-sigma error on BG subtracted roiCnt
    };
    vector<speclines> sl;
    /* input files, cuts */
    TChain *ch = NULL;
    int DS;
    int c;
    int p;
    int d;
    /* analysis parameters */
    vector<double>* trapENFCal = NULL;
    vector<bool>* isGood = NULL;
    bool* isLNFill1 = NULL; //vector<bool>*
    bool* isLNFill2 = NULL; //vector<bool>*
    bool* muVeto = NULL; //vector<bool>*
    vector<unsigned int>* wfDCBits = NULL;
    //vector<double>* time_s = NULL;
    int* mH = 0;
    //vector<double>* dcr90 = NULL;
    vector<int>* channel = NULL;
    vector<int>* C = NULL;
    vector<int>* P = NULL;
    vector<int>* D = NULL;
    /* plot settings */
    string eEst = "";
    double keVPerBin = 0.5;
    int nBins = 0;
    /* output file and results */
    TFile *f = NULL;
    char* fname;
    int doSave = 0;
    /* fit */
    static Bool_t reject; // status to set rejection of ROI from BG fit
    static double roiLoVal;
    static double roiHiVal;
    /* upper limit calculation */
    
public:
    /* set */
    SpecAnMultiLine(int arg_DS, int arg_c, int arg_p, int arg_d);
    void AddIsotope(string arg_isotope, vector<double> arg_energy); // add isotopes and lines to speclines struct
    void DoSave(char* arg_fname) {doSave = 1; fname = arg_fname;} // save histograms to output file
    void CloseFile(); // close file
    /* get */
    void DumpIsotopes(); // dump speclines struct
    int GetNIsotopes() {return sl.size();}
    int GetNLines(int arg_i) {return sl[arg_i].energy.size();}
    double GetRoiCnt(int arg_i, int arg_j) {return sl[arg_i].roiCnt[arg_j];}
    double GetBgCnt(int arg_i, int arg_j) {return sl[arg_i].bgCnt[arg_j];}
    double GetPkCnt(int arg_i, int arg_j) {return sl[arg_i].pkCnt[arg_j];}
    double GetePkCnt(int arg_i, int arg_j) {return sl[arg_i].ePkCnt[arg_j];}
    string GetIsotope(int arg_i) {return sl[arg_i].isotope;}
    double GetEnergy(int arg_i, int arg_j) {return sl[arg_i].energy[arg_j];}
    /* other */
    void WriteResultsToTextFile(ofstream &outTextFile);
    void LoopIsotopes(); // generic looping routine over speclines struct
    void LoopData(); // loop over data, contains hist-filling and unbinned integration
    int Cut(int j); // set the cuts for each DS
    void Fill(int j); // fill histograms
    void DoWrite(); // write histograms to output file
    void UnbinnedIntegration(int j); // tallies up counts in ROIs
    double HalfWidthHalfMax(double E); // determine ROI
    void DoFit(); // BG fit
    static Double_t SidebandFit(Double_t *x, Double_t *par); // BG fit TF1 that excludes ROI
    void CalculateUpperLimit();
    void Error(string errormessage); // prints error and ends program
    //void SetHistList(vector<TH1D*> arg_histList) {histList = arg_histList;}     //vector<TH1D*> histList;
};
Bool_t SpecAnMultiLine::reject = kTRUE;
double SpecAnMultiLine::roiLoVal = 0.0;
double SpecAnMultiLine::roiHiVal = 0.0;

SpecAnMultiLine::SpecAnMultiLine(int arg_DS, int arg_c, int arg_p, int arg_d) {
    /* set vars and check for valid values */
    DS = arg_DS;
    c = arg_c;
    p = arg_p;
    d = arg_d;
    if(DS<0 || DS>6) {Error("Unknown DS");}
    if(c!=1 && c!=2) {Error("Incorrect module number c");}
    if(p<1 || p>7) {Error("Incorrect string/position number p");}
    if(d<1 || d>5) {Error("Incorrect detector number d");}
    
    /* set DS-specific files and e-estimators */
    if(DS==0) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS0/GAT-v01-06-134-g3f44fab/skimDS0*.root");
        eEst = "trapENFCal";
    }
    if(DS==1) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS1/GAT-v01-06-134-g3f44fab/skimDS1*.root");
        eEst = "trapENFCal";
    }
    if(DS==2) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS2/GAT-v01-06-134-g3f44fab/skimDS2*.root");
        eEst = "trapENFCal";
    }
    if(DS==3) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS3/GAT-v01-06-134-g3f44fab/skimDS3*.root");
        eEst = "trapENFCal";
    }
    if(DS==4) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS4/GAT-v01-06-134-g3f44fab/skimDS4*.root");
        eEst = "trapENFCal";
    }
    if(DS==5) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS5/GAT-v01-06-134-g3f44fab/skimDS5*.root");
        eEst = "trapENFCal";
    }
    if(DS==6) {
        ch = new TChain("skimTree", "skimTree");
        ch->Add("/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS6/GAT-v01-06-22-g5c7b61c/skimDS6*.root");
        eEst = "trapENFCal";
    }
    
    /* SET DS-SPECIFIC Branches */ // subject to version of reprocessing
    ch->SetBranchAddress("isLNFill1",&isLNFill1);
    ch->SetBranchAddress("isLNFill2",&isLNFill2);
    ch->SetBranchAddress("C",&C);
    ch->SetBranchAddress("trapENFCal",&trapENFCal);
    ch->SetBranchAddress("isGood",&isGood);
    ch->SetBranchAddress("muVeto",&muVeto);
    ch->SetBranchAddress("wfDCBits",&wfDCBits);
    //ch->SetBranchAddress("time_s",&time_s);
    ch->SetBranchAddress("mH",&mH);
    //ch->SetBranchAddress("dcr90",&dcr90);
    ch->SetBranchAddress("channel",&channel);
    ch->SetBranchAddress("P",&P);
    ch->SetBranchAddress("D",&D);
    
    /* output */
    cout<<"-------------------------------------"<<endl;
    cout<<"DS"<<DS<<" "<<"C"<<c<<"P"<<p<<"D"<<d<<endl;
    /*
    cout<<"Run List: set in SpecAnMultiLine::SpecAnMultiLine()"<<endl;
    cout<<"Energy Estimator: set in SpecAnMultiLine::Fill()"<<endl;
    cout<<"Cuts: set in SpecAnMultiLine::Cut()"<<endl;
    cout<<"Plot Settings: set in SpecAnMultiLine::AddIsotope()"<<endl;
    */
    cout<<"-------------------------------------"<<endl;
}

int SpecAnMultiLine::Cut(int j) {
    /* define cuts for each DS */
    int cutresult = 0;
    if (DS==0) {
        if(C->at(j)==c && P->at(j)==p && D->at(j)==d && isGood->at(j) && !(isLNFill1 && C->at(j)==1) && !muVeto && !wfDCBits->at(j)) {cutresult = 1;}
    }
    if (DS==1) {
        if(C->at(j)==c && P->at(j)==p && D->at(j)==d && isGood->at(j) && !(isLNFill1 && C->at(j)==1) && !muVeto && !wfDCBits->at(j) /*&& (!(time_s->at(j) > 2192e3 && time_s->at(j) < 2195e3) && !(time_s->at(j) > 7370e3 && time_s->at(j) < 7371e3) && !(time_s->at(j) > 7840e3 && time_s->at(j) < 7860e3) && !(time_s->at(j) > 8384e3 && time_s->at(j) < 8387e3) && !(time_s->at(j) > 8984e3 && time_s->at(j) < 8985e3) && !(time_s->at(j) > 9002e3 && time_s->at(j) < 9005e3))*/) {cutresult = 1;}
    }
    if (DS==2) {
        if(C->at(j)==c && P->at(j)==p && D->at(j)==d && isGood->at(j) && !(isLNFill1 && C->at(j)==1) && !muVeto && !wfDCBits->at(j)) {cutresult = 1;}
    }
    if (DS==3) {
        if(C->at(j)==c && P->at(j)==p && D->at(j)==d && isGood->at(j) && !(isLNFill1 && C->at(j)==1) && !muVeto && !wfDCBits->at(j)) {cutresult = 1;}
    }
    if (DS==4) {
        if(C->at(j)==c && P->at(j)==p && D->at(j)==d && isGood->at(j) && !(isLNFill2 && C->at(j)==2) && !muVeto && !wfDCBits->at(j)) {cutresult = 1;}
    }
    if (DS==5) {
        if(C->at(j)==c && P->at(j)==p && D->at(j)==d && isGood->at(j) && !(isLNFill1 && C->at(j)==1) && !(isLNFill2 && C->at(j)==2) && !muVeto && !wfDCBits->at(j)) {cutresult = 1;}
    }
    if (DS==6) {
        if(C->at(j)==c && P->at(j)==p && D->at(j)==d && isGood->at(j) && !(isLNFill1 && C->at(j)==1) && !(isLNFill2 && C->at(j)==2) && !muVeto && !wfDCBits->at(j)) {cutresult = 1;}
    }
    return cutresult;
}

void SpecAnMultiLine::Fill(int j) {
    /* loop specline struct */
    for (unsigned int i_fill = 0; i_fill < sl.size(); i_fill++)
    {
        for (unsigned int j_fill = 0; j_fill < sl[i_fill].energy.size(); j_fill++)
        {
                if(trapENFCal->at(j)>sl[i_fill].binLo[j_fill] && trapENFCal->at(j)<sl[i_fill].binHi[j_fill])
                {
                    sl[i_fill].henergy[j_fill]->Fill(trapENFCal->at(j));
                }
        }
    }
}

void SpecAnMultiLine::UnbinnedIntegration(int j) {
    /* loop specline struct */
    for (unsigned int i_fill = 0; i_fill < sl.size(); i_fill++)
    {
        for (unsigned int j_fill = 0; j_fill < sl[i_fill].energy.size(); j_fill++)
        {
                if(trapENFCal->at(j)>sl[i_fill].roiLo[j_fill] && trapENFCal->at(j)<sl[i_fill].roiHi[j_fill])
                {
                    sl[i_fill].roiCnt[j_fill]++;
                }
                if((trapENFCal->at(j)>sl[i_fill].sbLo[j_fill] && trapENFCal->at(j)<=sl[i_fill].roiLo[j_fill]) || (trapENFCal->at(j)>=sl[i_fill].roiHi[j_fill] && trapENFCal->at(j)<sl[i_fill].sbHi[j_fill]))
                {
                    /* If DoFit() is invoked, the value of bgCnt will be overwritten */
                    sl[i_fill].bgCnt[j_fill]++;
                }
        }
    }
}

double SpecAnMultiLine::HalfWidthHalfMax(double E) {
    /* compute HWHM based on Pinghan's DS3 resolution data */
    double p0 = 0.164225116870199;
    double p1 = 0.0171317589551445;
    double p2 = 0.000305455041318674;
    double sigma = sqrt(p0*p0 + p1*p1*E + p2*p2*E*E);
    double HalfWHM = 0.5*(5*sigma);
    return HalfWHM;
}

void SpecAnMultiLine::AddIsotope(string arg_isotope, vector<double> arg_energy) {
    /* add arguments into the speclines struct */
    char hnametitle[50];
    speclines tempstruc;
    tempstruc.isotope = arg_isotope;
    tempstruc.energy = arg_energy;
    for(unsigned int i = 0; i < tempstruc.energy.size(); i++) {
        /* fill zeroes */
        tempstruc.success.push_back(0);
        tempstruc.roiCnt.push_back(0);
        tempstruc.bgCnt.push_back(0);
        tempstruc.pkCnt.push_back(0);
        tempstruc.ePkCnt.push_back(0);
        
        /* set up hist for each energy of each isotope */
        tempstruc.binLo.push_back(tempstruc.energy[i] - 0.035*tempstruc.energy[i]);
        tempstruc.binHi.push_back(tempstruc.energy[i] + 0.035*tempstruc.energy[i]);
        nBins = (tempstruc.binHi[i]-tempstruc.binLo[i])/keVPerBin;
        sprintf(hnametitle,"h_%s_%.0f_DS%d_C%dP%dD%d",tempstruc.isotope.c_str(),tempstruc.energy[i],DS,c,p,d);
        tempstruc.henergy.push_back(new TH1D(hnametitle,hnametitle,nBins,tempstruc.binLo[i],tempstruc.binHi[i]));
        tempstruc.henergy[i]->GetXaxis()->SetTitle(Form("%s (keV)",eEst.c_str()));
        tempstruc.henergy[i]->GetYaxis()->SetTitle(Form("Counts per %.2f keV",keVPerBin));
        /* set bounds for fit (sb) and roi */
        tempstruc.sbLo.push_back(tempstruc.energy[i] - 0.025*tempstruc.energy[i]);
        tempstruc.sbHi.push_back(tempstruc.energy[i] + 0.025*tempstruc.energy[i]);
        tempstruc.roiLo.push_back(tempstruc.energy[i] - HalfWidthHalfMax(tempstruc.energy[i]));
        tempstruc.roiHi.push_back(tempstruc.energy[i] + HalfWidthHalfMax(tempstruc.energy[i]));
    }
    sl.push_back(tempstruc);
}

void SpecAnMultiLine::LoopData() {
    /* loop over the input ROOT files */
    cout<<"-------------------------------------"<<endl;
    cout << "Begin looping data" << endl;
    if(doSave==0) {cout<<"WARNING: No output file is set. Not saving."<<endl;} //if(f==NULL)
    int nentries = ch->GetEntries();
    cout <<"nentries "<<nentries<<endl;
    for(int i_iter=0; i_iter < nentries; i_iter++)
    {
        ch->GetEntry(i_iter);
        if(i_iter%500000==0) {cout<<"at entry "<<i_iter<<"/"<<nentries<<endl;}
        for (unsigned int j_iter=0; j_iter<channel->size(); j_iter++)
        {
            if(Cut(j_iter))
            {
                //cout<<"PassedCut "<<i_iter<<" "<<C->at(j_iter)<<P->at(j_iter)<<D->at(j_iter)<<endl;
                /* do stuff */
                if(channel->at(j_iter)%2==0) {Fill(j_iter); UnbinnedIntegration(j_iter);} // mH->at(j_iter)==1
            }
        }
    }
    //if(f!=NULL) {DoWrite(); /*cout<<"Closing output file"<<endl;*/}
    if(doSave == 1) {DoWrite();}
    cout<<"-------------------------------------"<<endl;
}

void SpecAnMultiLine::DoWrite() {
    /* write hists to file */
    cout<<"-------------------------------------"<<endl;
    cout<<"Opening output file "<<fname<<endl;
    f = new TFile(fname,"UPDATE");
    f->cd();
    cout<<"Writing hists to ouptut file"<<endl;
    for (unsigned int i = 0; i < sl.size(); i++) {
        for (unsigned int j = 0; j < sl[i].energy.size(); j++)
        {
            sl[i].henergy[j]->Write();
            //cout<<sl[i].isotope<<endl;
            //cout<<"    "<<sl[i].energy[j]<<" roiCnt="<<sl[i].roiCnt[j]<<" "<<sl[i].henergy[j]->GetName()<<endl;
        }
    }
    f->Write(); // needed to fully write/flush tree to file
    cout<<"-------------------------------------"<<endl;
}

void SpecAnMultiLine::DoFit() {
    /* do sideband fit for each hist */
    cout<<"-------------------------------------"<<endl;
    cout<<"Performing sideband BG fits"<<endl;
    if(doSave == 1) {f->cd();}
    for (unsigned int i = 0; i < sl.size(); i++)
    {
        for (unsigned int j = 0; j < sl[i].energy.size(); j++)
        {
            roiLoVal = sl[i].roiLo[j];
            roiHiVal = sl[i].roiHi[j];
            
            /* perform fit */
            TF1 *sbFit = new TF1("fl",SpecAnMultiLine::SidebandFit,sl[i].sbLo[j],sl[i].sbHi[j],2);
            sbFit->SetParameters(0.0,0.0);
            reject = kTRUE;
            if(sl[i].henergy[j]->GetEntries()>0) {sl[i].henergy[j]->Fit(sbFit,"0 Q");} // 0=don't draw fit
            reject = kFALSE;
            
            /* record integral of BG fit */
            /* The value for bgCnt set by this DoFit routine supercedes that set by UnbinnedIntegration() */
            sl[i].bgCnt[j] = sbFit->Integral(sl[i].sbLo[j],sl[i].roiLo[j]) + sbFit->Integral(sl[i].roiHi[j],sl[i].sbHi[j]);
            
            /* save fit pars as vector in file */
            /*
            TVectorD par(2);
            for (int par_i=0; par_i<2; par_i++) {
                par[par_i] = sbFit->GetParameter(par_i);
            }
            if(doSave == 1) {par.Write(Form("sbFitPars_%s_%.0f_DS%d_C%dP%dD%d",sl[i].isotope.c_str(),sl[i].energy[j],DS,c,p,d));}
             */
            
            delete sbFit;
            sbFit = 0;
        }
    }
    if(doSave == 1) {f->Write();} // needed to fully write/flush tree to file
    cout<<"-------------------------------------"<<endl;
}

Double_t SpecAnMultiLine::SidebandFit(Double_t *x, Double_t *par) {
    if (reject && x[0] > roiLoVal && x[0] < roiHiVal) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0];
}

void SpecAnMultiLine::CalculateUpperLimit() {
    /* calculate peak upper limit according to Wenqin's Quick & Dirty prescription */
    cout<<"-------------------------------------"<<endl;
    cout<<"Calculating peak upper limits "<<endl;
    for (unsigned int i = 0; i < sl.size(); i++) {
        cout<<sl[i].isotope<<endl;
        for (unsigned int j = 0; j < sl[i].energy.size(); j++) {
            double S = sl[i].roiCnt[j];
            double B = sl[i].bgCnt[j];
            double dB = (sl[i].roiLo[j] - sl[i].sbLo[j]) + (sl[i].sbHi[j] - sl[i].roiHi[j]);
            double tau = dB/(sl[i].roiHi[j] - sl[i].roiLo[j]); // background normalization coeff
            double pCnt = S - (B/tau); // BG subtracted roiCnt
            sl[i].pkCnt[j] = pCnt;
            sl[i].ePkCnt[j] = sqrt(S + (B/(tau*tau))); // one-sigma pk error
            cout<<"    "<<sl[i].energy[j]<<" S="<<S<<" B="<<B<<" dB="<<dB<<" tau="<<tau<<" pkCnt="<<pCnt<<" ePkCnt="<<sl[i].ePkCnt[j]<<endl;
        }
    }
    cout<<"-------------------------------------"<<endl;
}

void SpecAnMultiLine::DumpIsotopes() {
    /* print the contents of the speclines struct */
    cout<<"-------------------------------------"<<endl;
    cout<<"Dumping list of isotopes"<<endl;
    cout<<"energy (analysis success) hist binLo|sbLo|roiLo|roiHi|sbHi|binHi:"<<endl;
    for (unsigned int i = 0; i < sl.size(); i++) {
        cout<<sl[i].isotope<<endl;
        if (sl[i].energy.size() == sl[i].success.size() && sl[i].energy.size() == sl[i].henergy.size()) {
            for (unsigned int j = 0; j < sl[i].energy.size(); j++) {
                cout<<sl[i].energy[j]<<" ("<<sl[i].success[j]<<") ";
                cout<<sl[i].henergy[j]->GetName()<<" ";
                cout<<sl[i].binLo[j]<<"|"<<sl[i].sbLo[j]<<"|"<<sl[i].roiLo[j]<<"|"<<sl[i].roiHi[j]<<"|"<<sl[i].sbHi[j]<<"|"<<sl[i].binHi[j]<<endl;
            }
        }
        else {Error("mismatch vector sizes");}
        cout<<endl;
    }
    cout<<"-------------------------------------"<<endl;
}

void SpecAnMultiLine::WriteResultsToTextFile(ofstream &outTextFile) {
    /* generic loop over the speclines struct */
    //cout<<"-------------------------------------"<<endl;
    for (unsigned int i = 0; i < sl.size(); i++) {
        //cout<<sl[i].isotope<<" (GetRoiCnt,GetBgCount)"<<endl;
        if (sl[i].energy.size() == sl[i].success.size() && sl[i].energy.size() == sl[i].henergy.size()) {
            for (unsigned int j = 0; j < sl[i].energy.size(); j++) {
                /* do stuff */
                //cout<<" ("<<GetRoiCnt(i,j)<<","<<GetBgCnt(i,j)<<")"<<endl;
                outTextFile << c << " " << p << " " << d << " " << sl[i].isotope << " " << sl[i].energy[j] << " " << GetRoiCnt(i,j) << " " << GetBgCnt(i,j) << " " << sl[i].pkCnt[j] << " " << sl[i].ePkCnt[j] << " " << Get_DetectorStatus(DS,c,p,d) << " " << Get_mAct_kg(DS,c,p,d) << " " << Get_liveTime_dy(DS) << endl;
            }
        }
        else {Error("mismatch vector sizes");}
        //cout<<endl;
    }
    //cout<<"-------------------------------------"<<endl;
}

void SpecAnMultiLine::LoopIsotopes() {
    //generic loop over the speclines struct
    cout<<"-------------------------------------"<<endl;
    for (unsigned int i = 0; i < sl.size(); i++) {
        cout<<sl[i].isotope<<endl;
        if (sl[i].energy.size() == sl[i].success.size() && sl[i].energy.size() == sl[i].henergy.size()) {
            for (unsigned int j = 0; j < sl[i].energy.size(); j++) {
                //do stuff
            }
        }
        else {Error("mismatch vector sizes");}
        cout<<endl;
    }
    cout<<"-------------------------------------"<<endl;
}

void SpecAnMultiLine::CloseFile() {
    /* close the file */
    if(doSave==1) {
        cout<<"Closing file"<<fname<<endl;
        f->cd();
        f->Close();
        doSave = 0;
        delete f;
        f = NULL;
    }
}

void SpecAnMultiLine::Error(string errormessage) {
    /* print message and exit program */
    cout<<"Error: "<<errormessage<<endl;
    exit(EXIT_FAILURE);
}

#endif /* SPECANMULTILINE_H */