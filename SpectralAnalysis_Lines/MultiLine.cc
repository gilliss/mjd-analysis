///////////////
//
// only uses high gain right now
//
///////////////

#include "MultiLine.hh"

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
        cout << "Usage is: ./MultiLine DS C CutScheme gain" << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // DS#
    int c = atoi(argv[2]); // M#
    //int cutScheme = atoi(argv[3]); // 0-3 [raw,raw+DC,raw+DC+AE,raw+DC+AE+DCR]
    //int gain =  atoi(argv[4]); // 0-1 [HG,LG]

    char outTextFileName[50];
    sprintf(outTextFileName,Form("./DS%d_M%d_hDetLines.txt",DS,c));
    ofstream outTextFile;
    outTextFile.open (outTextFileName);
    
    PrepForDSIandCSI(); // prep for masses and stuff from DataSetInfo and ChannelSelectionInfo

    for(int p = 1; p <= 7; p++)
    {
        for(int d = 1; d <= 5; d++)
        {
            // CREATE INSTANCE OF CLASS
            SpecAnMultiLine SAML(DS,c,p,d);
            // DEFINE ISOTOPES AND LINES OF INTEREST
            string tempi;
            vector<double> tempe;
            tempi = "228Ac"; tempe.push_back(338.320); tempe.push_back(911.204); tempe.push_back(968.971); tempe.push_back(1588.2);
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            tempi = "108mAg"; tempe.push_back(433.937); tempe.push_back(614.276); tempe.push_back(722.907);
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            tempi = "214Bi"; tempe.push_back(609.312); tempe.push_back(1120.3); tempe.push_back(1764.494); tempe.push_back(2204.1);
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            tempi = "60Co"; tempe.push_back(1173.237); tempe.push_back(1332.501);
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            tempi = "40K"; tempe.push_back(1460.830);
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            tempi = "234mPa"; tempe.push_back(1001.03);
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            tempi = "210Pb"; tempe.push_back(46.539);
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            tempi = "214Pb"; tempe.push_back(351.9321); // progenitor of 214Bi
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            tempi = "226Ra"; tempe.push_back(186.211);
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            tempi = "208Tl"; tempe.push_back(583.191); tempe.push_back(1593.0); tempe.push_back(2105.0); tempe.push_back(2614.533);
            SAML.AddIsotope(tempi,tempe);
            tempe.clear();
            if(p==1 && d==1) {SAML.DumpIsotopes();}

            SAML.DoSave(Form("./DS%d_M%d_hDetLines.root",DS,c)); // sets doSave=1
            SAML.LoopData();
            SAML.CloseFile(); // closes the root file and sets doSave=0
            SAML.CalculateUpperLimit();
            SAML.WriteResultsToTextFile(outTextFile);


        } // end p
    } // end d
    
    outTextFile.close();

}