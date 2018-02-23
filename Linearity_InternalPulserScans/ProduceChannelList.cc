/***************************************************
could also put together lists of kPreampDigitizers ...
****************************************************/

#include <TFile.h>

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>
#include <vector>

#include <MJTChannelMap.hh>

using namespace std;

int main (int argc, char* argv[])
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Read in command line arguments, make TChain
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    /* READ IN ARGUMENTS/SETTINGS FROM COMMAND */
    if (argc < 2)
    {
        cout << "too few arguments" << argv[0] << endl;
        return 1;
    }
    int DS = atoi(argv[1]);

    TFile* inFile = new TFile("./DSChannelMaps.root");
    MJTChannelMap *chmap = (MJTChannelMap*)inFile->Get(Form("ChannelMapDS%d",DS));

    int channel = 0;
    for(int c = 1; c <= 2; c++)
    {
        cout<<"Module "<<c<<endl;
        cout<<"[";
        for(int p = 1; p <= 7; p++)
        {
            for(int d = 1; d <= 5; d++)
            {
                if(chmap->GetDetectorName(c,p,d)!="") // avoid nonexistent detectors
                {
                    channel = chmap->GetInt(chmap->GetDetectorName(c,p,d),"kIDHi");
                    cout << channel;
                    cout << ",";
                    channel = chmap->GetInt(chmap->GetDetectorName(c,p,d),"kIDLo");
                    cout << channel;
                    cout << ",";
                }
            }
        }
        cout<<"]"<<endl;
    }
    
    
    
}