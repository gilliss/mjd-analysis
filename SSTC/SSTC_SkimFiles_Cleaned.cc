/*************************************************
Passes through a list of gatified runs, one run at a time. Tags prompt events and stores their global timestamps in a vector. Checks all events for correlation (of chosen type) with prompt events. If correlation is found, event is tagged as delayed event.
 
Data of the tagged prompt and delayed events are stored in a skim tree that can later be friended to the original MJD tree from which the skim tree was made. Once friended, the MJD data can be analyzed using the tags and data of the skim tree.
 
The most complicated thing here is probably the fact that this script makes its own tree and branches in order to hold the data (tags and timestamps) that it creates. The trees are structured such that each entry contains the WFs from each detector that was hit at particular time (i.e. the n WFs that determine multiplicity).
 
 Notes/Issues:
 x-If I cut out ghost/bad channels, than the skim tree may end up a different size as the real tree. As a soln, could skip ghost/bad channels only when it comes to the prompt and delayed tagging methods.
 x-Something funny in skim tree such that alot of events have 0 globalTimestamp and 0 timeSince, as opposed to -1 timeSince, since first prompt event has not yet occured. SOLN: THIS WAS RELATED TO THE POINT ABOVE, i.e. LG CHANNELS WERE EXCLUDED AND THUS GETTING 0-VALUES WRITTEN IN FOR THEM
 -Sometimes an HG-/LG-channel pair do not fire together, you might see one or the other, but not both

 *************************************************/

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TEventList.h>
#include <TEntryList.h>
#include <TCut.h>
#include <TMath.h>
#include <TApplication.h>  //This class creates the ROOT Application Environment that interfaces to the windowing system eventloop and eventhandlers

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h> 	// atof, atoi
#include <iomanip>      // std::setprecision
#include <utility>      // pair<Type1,Type2>
//#include <vector>

//#include <GATBaseClassesDICT.h>   //#include <GATDataSet.hh>   //#include <GATPeakShape.hh>
//#include <MJTChannelSettings.hh>
//#include <MJTChannelMap.hh>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initializing
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// A FUNCTION TO CHECK FOR GOOD CHANNELS, RETURNS 1000 IF NOT A CHANNEL TO BE INCLUDED //
int CheckChannel(int ch)
{
    int cheasy = 1000;
    if(ch == 646) cheasy = 0; //HG channels, * * * BASED ON P3JDY * * *
    if(ch == 644) cheasy = 1;
    if(ch == 642) cheasy = 2;
    if(ch == 626) cheasy = 3;
    if(ch == 624) cheasy = 4;
    if(ch == 674) cheasy = 5;
    if(ch == 576) cheasy = 6;
    if(ch == 692) cheasy = 7;
    if(ch == 690) cheasy = 8;
    if(ch == 688) cheasy = 9;
    if(ch == 640) cheasy = 10;
    if(ch == 610) cheasy = 11;
    if(ch == 664) cheasy = 12;
    if(ch == 662) cheasy = 13;
    if(ch == 656) cheasy = 14;
    if(ch == 696) cheasy = 15;
    if(ch == 608) cheasy = 16;
    if(ch == 598) cheasy = 17;
    if(ch == 600) cheasy = 18;
    if(ch == 594) cheasy = 19;
    if(ch == 592) cheasy = 20;
    /*
     if(ch == 584) cheasy = 21; //Bad HV connection
     if(ch == 680) cheasy = 22; //High leakage current1
     if(ch == 676) cheasy = 23; //Bad signal connection
     if(ch == 616) cheasy = 24; //Bad signal connection
     if(ch == 614) cheasy = 25; //Gain oscillation
     if(ch == 628) cheasy = 26; //Sparking
     if(ch == 632) cheasy = 27; //High leakage current
     if(ch == 630) cheasy = 28; //Bad signal connection
     */
    return cheasy;
}
// A FUNCTION TO SHIFT ENERGIES BASED ON ZERO-POINT CALIBRATION //
double EnergyShift(int ch)
{
    double shift = 0;
    if(ch == 646) shift = 0.022; // oldE - shift = correctedE, * * * BASED ON P3JDY * * *
    if(ch == 644) shift = 0.091;
    if(ch == 642) shift = 0.009;
    if(ch == 626) shift = -0.181;
    if(ch == 624) shift = -0.153;
    if(ch == 674) shift = -0.022;
    if(ch == 576) shift = 0.075;
    if(ch == 692) shift = 0.318;
    if(ch == 690) shift = 0.124;
    if(ch == 688) shift = 0.469;
    if(ch == 640) shift = 0.432;
    if(ch == 610) shift = -0.013;
    if(ch == 664) shift = 0.474;
    if(ch == 662) shift = 0.159;
    if(ch == 656) shift = 0.196;
    if(ch == 696) shift = 0.449;
    if(ch == 608) shift = 0.381;
    if(ch == 598) shift = 0.322;
    if(ch == 600) shift = -0.209;
    if(ch == 594) shift = 0.617;
    if(ch == 592) shift = 0.288;
    if(ch == 584) shift = 0.000; //Bad HV connection
    if(ch == 680) shift = 0.000; //High leakage current1
    if(ch == 676) shift = 0.000; //Bad signal connection
    if(ch == 616) shift = 0.000; //Bad signal connection
    if(ch == 614) shift = 0.448; //Gain oscillation
    if(ch == 628) shift = 0.587; //Sparking
    if(ch == 632) shift = 0.000; //High leakage current
    if(ch == 630) shift = 0.000; //Bad signal connection
    return shift;
}
// A FUNCTION TO CHECK FOR Enr VS Nat DETs //
int CheckEnriched(int ch)
{
    int cheasy = 1000;
    if(ch == 646) cheasy = 0; //Enr=1, Nat=0 * * * BASED ON P3JDY * * *
    if(ch == 644) cheasy = 1;
    if(ch == 642) cheasy = 1;
    if(ch == 626) cheasy = 1;
    if(ch == 624) cheasy = 1;
    if(ch == 674) cheasy = 1;
    if(ch == 576) cheasy = 1;
    if(ch == 692) cheasy = 1;
    if(ch == 690) cheasy = 1;
    if(ch == 688) cheasy = 1;
    if(ch == 640) cheasy = 1;
    if(ch == 610) cheasy = 1;
    if(ch == 664) cheasy = 0;
    if(ch == 662) cheasy = 1;
    if(ch == 656) cheasy = 1;
    if(ch == 696) cheasy = 1;
    if(ch == 608) cheasy = 0;
    if(ch == 598) cheasy = 0;
    if(ch == 600) cheasy = 0;
    if(ch == 594) cheasy = 0;
    if(ch == 592) cheasy = 0;
    if(ch == 584) cheasy = 0; //Bad HV connection
    if(ch == 680) cheasy = 1; //High leakage current1
    if(ch == 676) cheasy = 0; //Bad signal connection
    if(ch == 616) cheasy = 1; //Bad signal connection
    if(ch == 614) cheasy = 1; //Gain oscillation
    if(ch == 628) cheasy = 1; //Sparking
    if(ch == 632) cheasy = 1; //High leakage current
    if(ch == 630) cheasy = 1; //Bad signal connection
    return cheasy;
}

int main(int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initializing
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TIME CONSIDERATIONS //
    double timeCorrelationWindow = 1*(67.71)*60; // duration of prompt-delayed time window in seconds
    //double startrunStartTime, endrunStopTime; // start of first run, end of last run
    //double totalTime = 0; // The duration of the run. Also used to calculate K-shell event rate
    //double clockFreq = 100000000.0; // hardware clock cycle rate
    
    // ENERGY ROI //
    double meanROI = 10.3; // energy of prompt event
    double fwhmROI = 2.0; // width of prompt event ROI

    // MJD DATA & MJD TREE BRANCHES //
    double startTime = 0;
    double stopTime = 0;
    vector<int>* channel = 0;
    vector<double>* trapECal = 0;
    vector<double>* trapETailMin = 0;
    vector<double>* toe = 0;
    vector<double>* tloc_s = 0; // time since start of run in seconds
    int mH;
    vector<bool>* isEnr = 0;
    vector<bool>* isGood = 0;
    int run = 0;

    double toeVal;
    double enrVal;
    
    // PROMPT-DELAYED DATA & SKIM TREE BRANCHES //
    double globalTimestampValue; // The values filling globalTimestamp
    vector<double> globalTimestamp; // vector containing globalTimestampValues of each waveform of an entry
    vector<int> promptTag; // a vector containing 0 for non-K-shell WFs, 1 for K-shell WFs
    vector<int> delayedTag; // 0 for out of window, 1 for in window
    vector<double> timeSincePrompt; // Time between a prompt (K-shell) and delayed event (in the SSTC window)
    vector<double> shiftedE;
    vector<int> delayedOutsideTag; // 0 for out of window, 1 for in window
    vector<int> chan;
    vector<double> ToE;
    vector<double> trapETM;
    vector<int> mult;


    
    vector<pair<int,double> > promptVector; // vector of pairs:(ROOT channel,globalTimestampValue). The "> >" space is important.
    int promptCounter = 0;
    int promptVectorSize = 0;

    // OTHER and FLOW CONTROL //
    char sstcFileName[200]; //, skimfile[200]
    int nentriesOriginal, nentriesSSTC; // number of entries in the MJD and SSTC trees
    int ch;
    int chanSize;
    int chIndex;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Tagging Bounds
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int nGoodChannels = 21; // number of dets to be plotted (add 3 for "combined +2", "enr +1", and "nat +0" plots)
    vector<pair<double,double> > energyBounds; // vector of pairs:(ROOT mean,sigma). The "> >" space is important.
    vector<pair<double,double> > toeBounds;
    energyBounds.push_back(std::make_pair(10.479,(3*0.511))); // Natural
    energyBounds.push_back(std::make_pair(10.535,(3*0.92))); // Enriched
    //energyBounds.push_back(std::make_pair(10.0899,(4*0.122843))); // Natural
    //energyBounds.push_back(std::make_pair(10.2127,(4*0.0898953))); // Enriched
    cout << "Natural: mean, sigma: " << energyBounds[0].first << ", " << energyBounds[0].second << endl;
    cout << "Enriched: mean, sigma: " << energyBounds[1].first << ", " << energyBounds[1].second << endl;
    double meanToE, sigmaToE, meanLowE, sigmaLowE;
    int nsigma = 3; // how many sigma on either side of mean
    /*
     for (int i=0; i<(nGoodChannels+3); i++)
     {
         if(i==0) {meanToE=1.26068; sigmaToE=0.0652452; meanLowE=10.2517; sigmaLowE=0.00764049;}
         if(i==1) {meanToE=1.33909; sigmaToE=0.0480142; meanLowE=(10.535-.091); sigmaLowE=0.92;} //// use values from spec analysis Enr (-.091)
         if(i==2) {meanToE=1.3479; sigmaToE=0.0474187; meanLowE=(10.535-.009); sigmaLowE=0.92;} //// use values from spec analysis Enr (-.009)
         if(i==3) {meanToE=1.34849; sigmaToE=0.0497867; meanLowE=10.1034; sigmaLowE=0.149821;}
         if(i==4) {meanToE=1.2795; sigmaToE=0.0578463; meanLowE=10.0397; sigmaLowE=0.737556;}
         if(i==5) {meanToE=1.3935; sigmaToE=0.0775471; meanLowE=10.1333; sigmaLowE=0.314887;}
         if(i==6) {meanToE=1.43737; sigmaToE=0.0493958; meanLowE=10.4964; sigmaLowE=0.332499;}
         if(i==7) {meanToE=1.35406; sigmaToE=0.047295; meanLowE=(10.535-.318); sigmaLowE=0.92;} //// use values from spec analysis Enr (-.318)
         if(i==8) {meanToE=1.36925; sigmaToE=0.0535161; meanLowE=10.1114; sigmaLowE=0.205546;}
         if(i==9) {meanToE=1.36014; sigmaToE=0.0549928; meanLowE=10.0879; sigmaLowE=0.504255;}
         if(i==10) {meanToE=1.35726; sigmaToE=0.047857; meanLowE=(10.535-.432); sigmaLowE=0.92;} //// use values from spec analysis Enr (-.432)
         if(i==11) {meanToE=1.32464; sigmaToE=0.0704632; meanLowE=9.73321; sigmaLowE=9.53437;} //// use values from spec analysis Enr (+.013)
         if(i==12) {meanToE=1.30797; sigmaToE=0.10652; meanLowE=10.0996; sigmaLowE=0.120488;}
         if(i==13) {meanToE=1.34843; sigmaToE=0.050752; meanLowE=10.0636; sigmaLowE=0.142456;}
         if(i==14) {meanToE=1.30276; sigmaToE=0.198618; meanLowE=10.2522; sigmaLowE=0.00351379;}
         if(i==15) {meanToE=1.31056; sigmaToE=0.0696599; meanLowE=(10.535-.449); sigmaLowE=0.92;} //// use values from spec analysis Enr (-.449)
         if(i==16) {meanToE=1.54966; sigmaToE=0.0688828; meanLowE=10.0761; sigmaLowE=0.129506;}
         if(i==17) {meanToE=1.37; sigmaToE=0.0863083; meanLowE=10.1038; sigmaLowE=0.275539;}
         if(i==18) {meanToE=1.47053; sigmaToE=0.0665682; meanLowE=10.0506; sigmaLowE=0.122577;}
         if(i==19) {meanToE=1.44507; sigmaToE=0.0708437; meanLowE=10.4445; sigmaLowE=0.126637;}
         if(i==20) {meanToE=1.47079; sigmaToE=0.0801819; meanLowE=10.4403; sigmaLowE=0.133696;}
         if(i==21) {meanToE=1.42404; sigmaToE=0.162012; meanLowE=10.0899; sigmaLowE=0.122843;} // this is lower than the Nat spec an 10.479,.511
         if(i==22) {meanToE=1.29952; sigmaToE=0.211273; meanLowE=(10.535); sigmaLowE=0.92;} // this is lower than the Enr spec an 10.535,.92
         if(i==23) {meanToE=1.30869; sigmaToE=0.216365; meanLowE=10.09; sigmaLowE=0.12084;}
         
         energyBounds.push_back(std::make_pair(meanLowE,nsigma*sigmaLowE));
         toeBounds.push_back(std::make_pair(meanToE,nsigma*sigmaToE));
     }
    */
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Load and Prepare Trees
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    //GATDataSet ds(/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/data/skim_*.root); // Local
    //GATDataSet ds("/global/homes/j/jasondet/myProject/Skim/skim_5610_5680.root"); // PDSF
    
    TChain* originalTree = new TChain("skimTree", "skimTree");
    originalTree->Add("/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/data/skim2/skim_*.root"); // Local
    //originalTree->Add("/global/homes/j/jasondet/myProject/Skim/skim_*.root"); // PDSF
    originalTree->SetBranchAddress("startTime",&startTime); // unix time at start of run
    originalTree->SetBranchAddress("stopTime",&stopTime); // unix time at end of run
    originalTree->SetBranchAddress("channel",&channel); // set addresses for branches in the MJD Tree
    originalTree->SetBranchAddress("trapECal",&trapECal);
    originalTree->SetBranchAddress("trapETailMin",&trapETailMin);
    originalTree->SetBranchAddress("toe",&toe);
    originalTree->SetBranchAddress("tloc_s",&tloc_s);
    originalTree->SetBranchAddress("mH",&mH);
    originalTree->SetBranchAddress("isEnr",&isEnr);
    originalTree->SetBranchAddress("isGood",&isGood);
    originalTree->SetBranchAddress("run",&run);

    // CREATE SKIM TREE WITH BRANCHES, MIRRORING MJD TREE'S HIERARCHY/STRUCTURE //
    TTree *sstcTree = new TTree("sstcTree","sstcTree"); // Declare new tree (name, title)
    sstcTree->Branch("globalTimestamp",&globalTimestamp); // Declare branch for storing globalTimestamps
    sstcTree->Branch("promptTag",&promptTag); // Declare branch for storing K-shell tag
    sstcTree->Branch("delayedTag",&delayedTag); // Declare branch for storing SSTC tag
    sstcTree->Branch("timeSincePrompt",&timeSincePrompt); // Declare branch for storing timeSincePrompt
    sstcTree->Branch("shiftedE",&shiftedE);
    sstcTree->Branch("delayedOutsideTag",&delayedOutsideTag); // Declare branch for storing SSTC tag
    sstcTree->Branch("chan",&chan);
    sstcTree->Branch("trapETM",&trapETM);
    sstcTree->Branch("ToE",&ToE);
    sstcTree->Branch("mult",&mult);

    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////// Tagging Method
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // LOOP OVER ENTRIES IN CURRENT MJD TREE //
        nentriesOriginal=originalTree->GetEntries(); // Get number of entries in the MJD Tree
        for(int k=0;k<nentriesOriginal;k++)
        {
            originalTree->GetEntry(k); // Read all branches of entry and return total number of bytes read.
            chanSize = channel->size(); // Get the multiplicity of channels within the entry
		
            // CREATE SKIM TREE //
            globalTimestamp.resize(chanSize); // tailor vector to suit the size of the current entry
            promptTag.resize(chanSize);
            delayedTag.resize(chanSize);
            timeSincePrompt.resize(chanSize);
            shiftedE.resize(chanSize);
            delayedOutsideTag.resize(chanSize);
            chan.resize(chanSize);
            trapETM.resize(chanSize);
            ToE.resize(chanSize);
            mult.resize(chanSize);


            for(int j=0; j<chanSize; j++) // loop through channels in the entry
        	{
                ch = channel->at(j);
                
                    globalTimestampValue = startTime + tloc_s->at(j);
                    globalTimestamp[j] = globalTimestampValue;
                
                    promptTag[j] = 0; // Tag a non-prompt WF with the "false" bit
                    delayedTag[j]=0;
                    timeSincePrompt[j] = -1.0; // Using -1 as a bogus value/placeholder
                    delayedOutsideTag[j] = 0;
                
                    shiftedE[j] = trapECal->at(j) - EnergyShift(ch);
                
                    chan[j] = ch;
                    ToE[j] = toe->at(j);
                    trapETM[j] = trapETailMin->at(j);
                    mult[j] = mH;

                    chIndex = CheckChannel(ch);
                
                    toeVal = toe->at(j);
                
                    enrVal = CheckEnriched(ch);

                    ////{ put loop here that enforces trapEMinTail<0 && toe[channel] bounds
                    /////////////////////////////////////////////
                    // FIND PROMPT EVENTS, TAG & ADD TO A LIST //
                    /////////////////////////////////////////////
                    //if(chIndex!=1000 && shiftedE[j] > (energyBounds[chIndex].first-(energyBounds[chIndex].second)) && shiftedE[j] < (energyBounds[chIndex].first+(energyBounds[chIndex].second)) && toeVal > (toeBounds[chIndex].first-(toeBounds[chIndex].second)) && toeVal < (toeBounds[chIndex].first+(toeBounds[chIndex].second))) // ChanByChan3sigma
                    if(chIndex!=1000 && shiftedE[j] > (energyBounds[enrVal].first-(energyBounds[enrVal].second)) && shiftedE[j] < (energyBounds[enrVal].first+(energyBounds[enrVal].second)) && toeVal > 1.0 && toeVal < 1.59 && trapETailMin->at(j) < 0.0) // fixed bounds
                    {
                        promptCounter++;
                        promptTag[j] = 1; // Tag a prompt WF with the "true" bit

                        promptVector.push_back(std::make_pair(ch,globalTimestampValue));
                        
                        //cout << "Prompt event at entry " << k << " and chan " << ch << endl;
                        //cout << "   promptVector size = " << promptVector.size() << endl;
                    }


                    ///////////////////////////////
                    // FIND DELAYED EVENTS & TAG //
                    ///////////////////////////////
                    promptVectorSize = int(promptVector.size());
                    if(chIndex!=1000 && promptVectorSize > 0)
                    {
                        // COMPARE CURRENT EVENT AGAINST ALL PROMPT EVENTS SO FAR //
                        for(int m=0; m<promptVectorSize; m++)
                        {
                            // ENFORCE YOUR DELAYED TAGGING CRITERIA //
                            //if((globalTimestampValue-promptVector[m].second)>0 && (globalTimestampValue-promptVector[m].second)<timeCorrelationWindow && channel->at(j)==promptVector[m].first) // ChanByChan3sigma
                            if((globalTimestampValue-promptVector[m].second)>0 && (globalTimestampValue-promptVector[m].second)<timeCorrelationWindow && channel->at(j)==promptVector[m].first && toeVal > 1.0 && toeVal < 1.59 && trapETailMin->at(j) < 0.0) // fixed bounds
                            {
                                delayedTag[j]=1;
                                timeSincePrompt[j] = globalTimestampValue-promptVector[m].second;
                            }
                            //if((globalTimestampValue-promptVector[m].second)>0 && (globalTimestampValue-promptVector[m].second)<timeCorrelationWindow && channel->at(j)!=promptVector[m].first) // ChanByChan3sigma
                            if((globalTimestampValue-promptVector[m].second)>0 && (globalTimestampValue-promptVector[m].second)<timeCorrelationWindow && channel->at(j)!=promptVector[m].first && toeVal > 1.0 && toeVal < 1.59 && trapETailMin->at(j) < 0.0) // fixed bounds
                            {
                                delayedOutsideTag[j]=1;
                                timeSincePrompt[j] = globalTimestampValue-promptVector[m].second;
                            }
                        }
                    }
                    ////}
            } // End of loop over channels j within the current entry
            
            sstcTree->Fill(); // Fill the skim tree with the new branches/vectors-- one vector for each entry
        } // End of loop over entries k
	
        // SANITY CHECK OUTPUTS //
        nentriesSSTC=sstcTree->GetEntries();
        cout << "Num Entries Original Tree = " << nentriesOriginal << endl;
        cout << "Num Entries SkimTree = " << nentriesSSTC << endl;
        cout << "promptVector size = " << promptVector.size() << endl;

        // SAVE SKIM TREE TO TFile //
        //sprintf(sstcFileName,"sstc_skimfile_5610-5680.root");
        sprintf(sstcFileName,"SSTC_SkimFiles_P3JDYGolden_CleanedFixedBounds3Sigma.root");
        TFile *sstcFile = new TFile(sstcFileName,"recreate");
        sstcTree->Write(); // Write the tree to the current directory
        sstcFile->Close();

        delete originalTree, sstcTree, sstcFile;

    cout << "F I N I S H E D" << endl;
    App->Run();
} // End of main()
