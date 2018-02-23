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

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initializing
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// A FUNCTION TO CHECK FOR GOOD CHANNELS, RETURNS 1000 IF NOT A CHANNEL TO BE INCLUDED //
int CheckChannel(int ch)
{
    int cheasy = 1000;
    if(ch == 646) cheasy = 0; //HG channels
    if(ch == 644) cheasy = 1;
    if(ch == 642) cheasy = 2;
    if(ch == 630) cheasy = 3;
    if(ch == 628) cheasy = 4;
    if(ch == 626) cheasy = 5;
    if(ch == 624) cheasy = 6;
    if(ch == 632) cheasy = 7;
    if(ch == 584) cheasy = 8;
    if(ch == 674) cheasy = 9;
    if(ch == 576) cheasy = 10;
    if(ch == 680) cheasy = 11;
    if(ch == 692) cheasy = 12;
    if(ch == 690) cheasy = 13;
    if(ch == 688) cheasy = 14;
    if(ch == 640) cheasy = 15;
    if(ch == 676) cheasy = 16;
    if(ch == 616) cheasy = 17;
    if(ch == 614) cheasy = 18;
    if(ch == 610) cheasy = 19;
    if(ch == 664) cheasy = 20;
    if(ch == 662) cheasy = 21;
    if(ch == 656) cheasy = 22;
    if(ch == 696) cheasy = 23;
    if(ch == 608) cheasy = 24;
    if(ch == 598) cheasy = 25;
    if(ch == 600) cheasy = 26;
    if(ch == 594) cheasy = 27;
    if(ch == 592) cheasy = 28;
    return cheasy;
}

int main(int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initializing
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    // RUN LIST CONSIDERATIONS //
    int startrun = 4854;
    int endrun = 4854;
    double run = 0;

    // CHANNEL CONSIDERATIONS //

    // TIME CONSIDERATIONS //
    double timeCorrelationWindow = 1*(67.71)*60; // duration of prompt-delayed time window in seconds
    double startrunStartTime, endrunStopTime; // start of first run, end of last run
    double totalTime = 0; // The duration of the run. Also used to calculate K-shell event rate
    double clockFreq = 100000000.0; // hardware clock cycle rate
    int corruptTimestampBuffer = 0; // buffer to skip n possibly-corrupt entries at the start of a run
    
    // ENERGY ROI //
    double meanROI = 11.0216; // energy of prompt event
    double fwhmROI = 4.99634; // width of prompt event ROI

    // MJD DATA & MJD TREE BRANCHES //
    double startTime = 0;
    double stopTime = 0;
    vector<double>* channel = 0;
    vector<double>* energyCal = 0;
    vector<double>* timestamp = 0;

    // PROMPT-DELAYED DATA & SKIM TREE BRANCHES //
    vector<double> globalTimestamp; // vector containing globalTimestampValues of each waveform of an entry
    double globalTimestampValue; // The values filling globalTimestamp
    double initialTimestamp = 0; // initial timestamp from gatified data (in clock cycles)
    vector<int> promptTag; // a vector containing 0 for non-K-shell WFs, 1 for K-shell WFs
    vector<int> delayedTag; // 0 for out of window, 1 for in window
    vector<double> timeSincePrompt; // Time between a prompt (K-shell) and delayed event (in the SSTC window)
    
    vector<pair<int,double> > promptVector; // vector of pairs:(ROOT channel,globalTimestampValue). The "> >" space is important.
    int promptCounter = 0;
    int promptVectorSize = 0;

    // OTHER and FLOW CONTROL //
    char infile[200], infilename[200], skimfilename[200]; //, skimfile[200]
    int nentriest, nentriestf; // number of entries in the MJD and SSTC trees
    int initTimestampTag; // This ensures only one timestamp is saved as the initial timestamp of a run
    int ch;
    int chanSize;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Tagging Method
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(int i=startrun;i<endrun+1;i++)
    {
        // LOAD & PREP THE ORIGINAL MJD TREE //
        cout << "----------RUN " << i << "-------------------" << endl;
        sprintf(infile,"mjd_run%d",i);
        //sprintf(infilename,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/gatified/P3JDY/%s.root",infile);
        sprintf(infilename,"/Users/tomgilliss/Desktop/UNC/ENAP/Analysis/SSTC/data/%s.root",infile);

        TChain *t = new TChain("mjdTree"); // Having this inside the for loop prevents files from being chained--the chain just restarts for each run
        t->AddFile(infilename); // Associate ROOT file with the current TChain
        t->SetBranchAddress("channel",&channel); // set addresses for branches in the MJD Tree
        t->SetBranchAddress("timestamp",&timestamp); // clock cycle at which WF recorded
        t->SetBranchAddress("run",&run);
        t->SetBranchAddress("startTime",&startTime); // unix time at start of run
        t->SetBranchAddress("stopTime",&stopTime); // unix time at end of run
        t->SetBranchAddress("energyCal",&energyCal);

       	nentriest=t->GetEntries(); // Get number of entries in the MJD Tree

        // CREATE SKIM TREE WITH BRANCHES, MIRRORING MJD TREE'S HIERARCHY/STRUCTURE //
        TTree *tf = new TTree("tf","tf"); // Declare new tree (name, title)
        tf->Branch("globalTimestamp",&globalTimestamp); // Declare branch for storing globalTimestamps
        tf->Branch("promptTag",&promptTag); // Declare branch for storing K-shell tag
        tf->Branch("delayedTag",&delayedTag); // Declare branch for storing SSTC tag
        tf->Branch("timeSincePrompt",&timeSincePrompt); // Declare branch for storing timeSincePrompt
        
        // LOOP OVER ENTRIES IN CURRENT MJD TREE //
        initTimestampTag = 0;
        for(int k=0;k<nentriest;k++)
        {
            t->GetEntry(k); // Read all branches of entry and return total number of bytes read.
            chanSize = channel->size(); // Get the multiplicity of channels within the entry
            

            // GET TIMESTAMPS //
            if(channel->at(0)!=0) // check for bogus 0-value
            {
                for(int j=0; j<chanSize; j++) // loop through channels in the entry
                {
                    ch = channel->at(j); // Get the value of the ROOT channel
                    if(CheckChannel(ch)!=1000) //if(ch%2==0) // use just high gain channels
                    {
                        if(k>corruptTimestampBuffer && initTimestampTag==0)
                        {
                            initialTimestamp = timestamp->at(j); // Save init stamp for use in calculating globalTimestampValue
                            initTimestampTag++;
                            //if(i==startrun) startrunStartTime = initialTimestamp; // This avoids the use of startTime which may be corrupt for some runs
                        }
                        if(i==startrun) startrunStartTime = startTime;
                        if(i==endrun) endrunStopTime = stopTime;
                    }
                }
            }
		
            // CREATE SKIM TREE //
            globalTimestamp.resize(chanSize); // tailor vector to suit the size of the current entry
            promptTag.resize(chanSize);
            delayedTag.resize(chanSize);
            timeSincePrompt.resize(chanSize);
            for(int j=0; j<chanSize; j++) // loop through channels in the entry
        	{
                ch = channel->at(j);
                //if(CheckChannel(ch)!=1000) //if(ch%2==0) // use just high gain channels
                //{
                    // cout << "Current channel = " << ch << " energy = " << energyCal->at(j) << endl;
                    globalTimestampValue = startTime + ((timestamp->at(j)-initialTimestamp)/clockFreq);
                    globalTimestamp[j] = globalTimestampValue;

                    /////////////////////////////////////////////
                    // FIND PROMPT EVENTS, TAG & ADD TO A LIST //
                    /////////////////////////////////////////////
                    if(CheckChannel(ch)!=1000 && energyCal->at(j) > (meanROI-(fwhmROI/2)) && energyCal->at(j) < (meanROI+(fwhmROI/2)))
                    {
                        promptCounter++;
                        promptTag[j] = 1; // Tag a prompt WF with the "true" bit

                        promptVector.push_back(std::make_pair(ch,globalTimestampValue));
                        
                        cout << "Prompt event at entry " << k << " and chan " << ch << endl;
                        cout << "   promptVector size = " << promptVector.size() << endl;
                    }
                    else promptTag[j] = 0; // Tag a non-prompt WF with the "false" bit

                    ///////////////////////////////
                    // FIND DELAYED EVENTS & TAG //
                    ///////////////////////////////
                    delayedTag[j]=0;
                    promptVectorSize = int(promptVector.size());
                    timeSincePrompt[j] = -1.0; // Using -1 as a bogus value/placeholder
                    if(CheckChannel(ch)!=1000 && promptVectorSize > 0)
                    {
                        // COMPARE CURRENT EVENT AGAINST ALL PROMPT EVENTS SO FAR //
                        for(int m=0; m<promptVectorSize; m++)
                        {
                            // ENFORCE YOUR DELAYED TAGGING CRITERIA //
                            if((globalTimestampValue-promptVector[m].second)>0 && (globalTimestampValue-promptVector[m].second)<timeCorrelationWindow && channel->at(j)==promptVector[m].first)
                            {
                                delayedTag[j]=1;
                                timeSincePrompt[j] = globalTimestampValue-promptVector[m].second;
                            }
                        }
                    }
                //}
            } // End of loop over channels j within the current entry
            
            tf->Fill(); // Fill the skim tree with the new branches/vectors-- one vector for each entry
        } // End of loop over entries k
	
        // SANITY CHECK OUTPUTS //
        cout << "- - - RUN " << i << " Summary - - -" << endl;
        nentriestf=tf->GetEntries();
        cout << "Num Entries Run " << i << " = " << nentriest << endl;
        cout << "Num Entries SkimTree " << i << " = " << nentriestf << endl;
        cout << "promptVector size = " << promptVector.size() << endl;
        cout << "startTime "<<std::setprecision(12)<<startTime<<" initialTimestamp "<<std::setprecision(12)<<initialTimestamp<<" intitialTimestamp in sec "<<initialTimestamp/clockFreq<<endl;
        totalTime = endrunStopTime - startrunStartTime;
        cout << "totalTime run " << i << " = " << totalTime << " s = " << totalTime/3600.0 << " hr "<< endl;

        // SAVE SKIM TREE TO TFile, LET LOOP TAKE YOU TO NEXT RUN //
        //sprintf(skimfile,"sstc_skimfile_run%d",i);
        //sprintf(skimfilename,"%s.root",skimfile);
        sprintf(skimfilename,"sstc_skimfile_run%d.root",i);
        TFile *skimFile = new TFile(skimfilename,"recreate");
        tf->Write(); // Write the tree to the current directory
        skimFile->Close();

        delete t, tf;
    }
    cout << "F I N I S H E D" << endl;
    App->Run();
} // End of main()
