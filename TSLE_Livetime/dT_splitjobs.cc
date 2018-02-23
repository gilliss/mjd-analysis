///////////////
//
// for DS5 might want to split this up into multiple run ranges. Make dT_splitjobs
//
///////////////

#include "dT.hh"

int main (int argc, char* argv[])
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Command line arguments and initializations
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 2)
    {
        cout << "Error: Too few arguments" << argv[0] << endl;
        cout << "Usage is: ..." << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // DS#
    int dataType = atoi(argv[2]); // [0,1]::[BG,Calib]
    int subsetLo = atoi(argv[3]);
    int subsetHi = atoi(argv[4]);
    
    // GET INPUT GATIFIED DATA
    cout << "DS" << DS << endl;
    GATDataSet ds;
    TChain* ch = new TChain();
    int dsSubsetMax = 0;
    if(dataType == 0 /* BG data */){
        cout << "Running on background data" << endl;
        dsSubsetMax = GetDataSetSequences(DS);
        cout << "Max Subset: " << dsSubsetMax << endl;
        cout << "SubsetLo - SubsetHi: " << subsetLo << " - " << subsetHi << endl;
        for(int subset_i = subsetLo; subset_i <= subsetHi; subset_i++) {LoadDataSet(ds, DS, subset_i);}
        ch = ds.GetGatifiedChain();
        cout<< " Number of runs in mjdTree chain " << ds.GetNRuns() << endl;
        cout<< " Number of events in mjdTree chain " << ch->GetEntries() << endl;
    }
    if(dataType == 1 /* Calib data */){
        cout << "Running on calibration data" << endl;
        if(DS==2){
            ds.AddRunRange(14508, 14515);
            ds.AddRunRange(14568, 14575);
            ds.AddRunRange(14699, 14706);
            ds.AddRunRange(14855, 14861);
            ds.AddRunRange(14927, 14933);
        }
        if(DS==5){
            ds.AddRunRange(18740, 18760);
            ds.AddRunRange(19055, 19069);
            ds.AddRunRange(19094, 19109);
            ds.AddRunRange(19243, 19258);
            ds.AddRunRange(19520, 19542);
        }
        ch = ds.GetGatifiedChain();
    }
    
    // SETUP THE INPUT BRANCHES
    double run = 0.;
    double startClockTime = 0.;
    double startTime = 0.;
    vector<double>* timestamp = 0;
    vector<double>* channel = NULL;
    vector<double>* trapENFCal = 0;
    ch->SetBranchAddress("run", &run);
    ch->SetBranchAddress("startTime", &startTime);
    ch->SetBranchAddress("startClockTime", &startClockTime);
    ch->SetBranchAddress("timestamp", &timestamp);
    ch->SetBranchAddress("channel", &channel);
    ch->SetBranchAddress("trapENFCal", &trapENFCal);
    
    // SETUP THE OUTPUT FILE
    TFile* f = new TFile(Form("DS%d_dT_%d_sub%d_%d.root",DS,dataType,subsetLo,subsetHi),"RECREATE");
    
    // SETUP THE OUTPUT TREE
    TTree* auxTree = new TTree("auxTree","auxTree");
    
    // SETUP THE OUTPUT BRANCHES
    int channel_aux = 0.;
    int run_aux = 0.;
    double dT = 0.;
    double globalTime = 0.;
    //double trapENFCal_aux = 0.;
    //double timestamp_aux = 0.;
    auxTree->Branch("channel",&channel_aux);
    auxTree->Branch("run",&run_aux);
    auxTree->Branch("dT",&dT);
    auxTree->Branch("globalTime",&globalTime);
    //auxTree->Branch("trapENFCal",&trapENFCal_aux);
    //auxTree->Branch("timestamp",&timestamp_aux);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Make channel map b/c what we want is dT within each channel
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    MJTChannelMap *chmap = ds.GetChannelMap();
    std::map<int,double> channelmap; // this will map between channel # and previous-event timestamp
    vector<int> channellist;
    int tempchan = 0;
    for(int c = 1; c <= 2; c++)
    {
        for(int p = 1; p <= 7; p++)
        {
            for(int d = 1; d <= 5; d++)
            {
                if(chmap->GetDetectorName(c,p,d)!="") // avoid nonexistent detectors
                {
                    tempchan = chmap->GetInt(chmap->GetDetectorName(c,p,d),"kIDHi"); // load HG chan #s
                    channelmap.insert(std::pair<int,double>(tempchan,0.)); // load with zeros to start
                    channellist.push_back(tempchan);
                    
                    tempchan = chmap->GetInt(chmap->GetDetectorName(c,p,d),"kIDLo"); // load LG chan #s
                    channelmap.insert(std::pair<int,double>(tempchan,0.)); // load with zeros to start
                    channellist.push_back(tempchan);
                }
            }
        }
    }
    sort(channellist.begin(), channellist.end()); // list must be sorted for binary_search
    for(unsigned int list_i = 0; list_i < channellist.size(); list_i++) {
        cout<<channellist[list_i]<<" "<<channelmap.at(channellist[list_i])<<endl;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Do the work
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    // LOOP DATA
    int nEvents = ch->GetEntries();
    int numHits = 0;
    double t1 = 0;
    double t2 = 0;
    int run1 = 0;
    int run2 = 0;
    for(Long64_t iEvent=0; iEvent < nEvents; iEvent++)
    {
        if(iEvent%1000000==0) {cout << "  iEvent: " << iEvent << std::endl;}
        ch->GetEntry(iEvent);
        
        // RESET THE AUX BRANCH VALS JUST FOR KICKS
        dT = 0.;
        channel_aux = 0;
        globalTime = 0.;
        
        // LOOP OVER WFs IN EVENT
        run1 = run; // to check for run boundary
        numHits = channel->size();
        for(int ich = 0; ich<numHits; ich++)
        {
            channel_aux = channel->at(ich);
            //cout <<"channel_aux "<<channel_aux<<" "<<chmap->GetDetectorName(channel_aux)<<"   "<<binary_search(channellist.begin(), channellist.end(), channel_aux)<<endl;
            if(binary_search(channellist.begin(), channellist.end(), channel_aux) /*&& trapENFCal->at(ich)>1.0*/) // cuts on channel and noise
            {
                if(run2!=run1) { // if first event of a run, just set every chan back to 0
                    //cout<<"run2!=run1"<<endl;
                    for(unsigned int list_i = 0; list_i < channellist.size(); list_i++) {
                        channelmap.at(channellist[list_i]) = 0.;
                    }
                }
                
                if(channelmap.at(channel_aux)==0) { // if first event for this channel, just set dT = 0
                    dT = 0.0;
                    channelmap.at(channel_aux) = timestamp->at(ich);
                    
                    run_aux = run;
                    globalTime = timestamp->at(ich)/1e8 - startClockTime/1e9 + startTime;
                    
                    //cout <<run<<" "<<iEvent<<" "<<numHits<<" "<<dT<< endl;
                }
                else {
                    t2 = timestamp->at(ich);
                    t1 = channelmap.at(channel_aux);
                    if(dataType==0) {dT = (t2-t1)*1e-8;} // time since last event in s
                    if(dataType==1) {dT = (t2-t1)*1e-2;} // time since last event in us
                    channelmap.at(channel_aux) = timestamp->at(ich);
                    
                    run_aux = run;
                    globalTime = timestamp->at(ich)/1e8 - startClockTime/1e9 + startTime;
                    
                    //cout <<run<<" "<<iEvent<<" "<<numHits<<" "<<dT<< endl;
                }
                
                auxTree->Fill(); // fill a row of the output tree
            }
        }
        run2 = run; // to check for run boundary
    }
    
    // WRITE OUTPUT TREE TO OUTPUT FILE
    f->cd();
    auxTree->Write();
    f->Close();

    cout << "FINISHED" << endl;
}



// TEST CHANNEL MAP AND STL MAP
/*
cout<<channelmap.at(648)<<endl;
channelmap.at(648)=10.;
cout<<channelmap.at(648)<<endl; // only need to make a channel check and replace t1. t2 can basically stay as is
cout<<chmap->GetDetectorName(648)<<endl;
*/
