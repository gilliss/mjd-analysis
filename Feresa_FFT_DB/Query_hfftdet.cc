///////////////////////////////
// Script to get integrated PSD, within a few frequency bands, for past runs. Write results to JSON files to be
// put on CouchDB.
//
// Tom Gilliss, gilliss@unc.edu
// v1 2017-12-16
//
// Refs and notes:
// -http://feresa.physics.unc.edu/mjdMonitor/script/fftMap.js
// -Tom C's run_plots.cc with functions hist_to_json and write_table
// -To get Json library working, need to add the library path and flag for libJSON to your LD_LIBRARY_PATH
// -Good c++ JSON ref: https://www.codeproject.com/Articles/1102603/Accessing-JSON-Data-with-Cplusplus
///////////////////////////////

#include <string> // string, find, substr
#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdlib.h>
#include <sstream> // peak, ignore, istringstream
#include <map>

#include "./json_library/json_value.hh"
#include "./json_library/json_writer.hh"
#include "./json_library/json_reader.hh"

using namespace std;

// STRUCT TO HOLD FREQ RANGES
struct FreqRangeStruct
{
    double minF, maxF;
    int minBin, maxBin;
};
int nRanges = 5; // fill channelDataMap with 5 zeros
map<int,FreqRangeStruct> freqRangeMap = {
    { 0, {0., 2., 0, 0} }, // { ds, dsExposureStruct (in kg-dy) }
    { 1, {2., 6., 0, 0} },
    { 2, {6., 9., 0, 0} },
    { 3, {9., 15., 0, 0} },
    { 4, {36., 48., 0, 0} }
    };

// FUNCTION TO FIND BIN
int GetBin(double nbins, double minval, double maxval, double val)
{
    int bin = 0;
    bin = 1 + int( nbins*(val-minval)/(maxval-minval) );
    return bin;
}

/////////////////////////
// MAIN
/////////////////////////
int main (int argc, char* argv[])
{
    // COMMAND LINE INPUTS
    int run = atoi(argv[1]);
    bool verbose = false;
    cout << "S T A R T E D" << " run " << run << endl;

    // CURL COMMNAND W/ OUTPUT TO TEXT FILE
    string curl = "curl -X GET ";
    string url = "\"http://feresa.physics.unc.edu/~feresa/mjdMonitor/php/designDocQuery.php?db=history_onsite_analysis_mac_pro&design=history&view=hfftdet&start=RunNumber," + to_string(run) + "&end=RunNumber," + to_string(run) + "\" ";
    string redirection = "> ";
    string queryFile = "hfftdet" + to_string(run) + ".txt";

    curl += url;
    curl += redirection;
    curl += queryFile;
    system (curl.c_str()); //cout<<"curl command:\n    "<<curl<<endl;

    // READ CURLED TEXT FILE INTO STRING
    ifstream inTextFile(queryFile.c_str());
    string content;
    content.assign(
                   (std::istreambuf_iterator<char>(inTextFile)),
                   (std::istreambuf_iterator<char>())
                  ); // https://stackoverflow.com/questions/2912520/read-file-contents-into-a-string-in-c

    // READ STRING INTO JSON VALUE
    Json::Reader reader;
    Json::Value root;
    if(reader.parse(content, root))
    {
        cout << "Parse succeeded: root[\"rows\"][0][\"key\"] = " << root["rows"][0]["key"] << endl; // have to use index [0] b/c the json file has this part in an array
    }
    else{cout<<"Warning: Parse failed"<<endl;}

    // THINGS TO TAKE FROM INPUT FILE
    vector<double> bincontents;
    int runNumber;
    int runType;
    string time;
    int timestamp;
    int xbins;
    double xmax, xmin;
    int ybins;
    double ymax, ymin;
    vector<string> ylabels;

    for(unsigned int i = 0; i < root["rows"][0]["value"]["bincontents"].size(); i++)
    {
        bincontents.push_back(root["rows"][0]["value"]["bincontents"][i].asDouble());
    }
    runNumber = root["rows"][0]["value"]["runNumber"].asInt();
    runType = root["rows"][0]["value"]["runType"].asInt();
    time = root["rows"][0]["value"]["time"].asString();
    timestamp = root["rows"][0]["value"]["timestamp"].asInt();
    xbins = root["rows"][0]["value"]["xbins"].asInt();
    xmax = root["rows"][0]["value"]["xmax"].asDouble();
    xmin = root["rows"][0]["value"]["xmin"].asDouble();
    ybins = root["rows"][0]["value"]["ybins"].asInt();
    ymax = root["rows"][0]["value"]["ymax"].asDouble();
    ymin = root["rows"][0]["value"]["ymin"].asDouble();
    for(unsigned int i = 0; i < root["rows"][0]["value"]["ylabels"].size(); i++)
    {
        ylabels.push_back(root["rows"][0]["value"]["ylabels"][i].asString());
    }

    // PRINT THINGS TAKEN FROM INPUT FILE
    if(verbose)
    {
        for(unsigned int vec_i = 0; vec_i < 5; vec_i++) {cout<<bincontents[vec_i]<<endl;} //bincontents.size()
        cout << runNumber << endl;
        cout << runType << endl;
        cout << time << endl;
        cout << timestamp << endl;
        cout << xbins << endl;
        cout << xmax << endl;
        cout << xmin << endl;
        cout << ybins << endl;
        cout << ymax << endl;
        cout << ymin << endl;
        for(unsigned int vec_i = 0; vec_i < 5; vec_i++) {cout<<ylabels[vec_i]<<endl;} //ylabels.size()
    }

    // MAKE MAP OF FREQ RANGES FOR NEW PLOT
    for(int map_i = 0; map_i<nRanges; map_i++)
    {
        freqRangeMap[map_i].minBin = GetBin(xbins,xmin,xmax,freqRangeMap[map_i].minF);
        freqRangeMap[map_i].maxBin = GetBin(xbins,xmin,xmax,freqRangeMap[map_i].maxF);
    }
    if(verbose)
    {
        cout<<"FFT Bands (index, freq band, bin range):"<<endl;
        for(int map_i = 0; map_i<nRanges; map_i++)
        {
            cout<<map_i<<endl;
            cout<<"    "<<freqRangeMap[map_i].minF<<" - "<<freqRangeMap[map_i].maxF<<endl;
            cout<<"    "<<freqRangeMap[map_i].minBin<<" - "<<freqRangeMap[map_i].maxBin<<endl;
            //cout<<"    "<<GetBin(xbins,xmin,xmax,freqRangeMap[map_i].minF)<<" - "<<GetBin(xbins,xmin,xmax,freqRangeMap[map_i].maxF)<<endl;
        }
    }

    // MAKE MAP TO HOLD INTEGRATED POWER FOR NEW PLOT
    map<int,vector<double> > channelDataMap;
    for(unsigned int vec_i = 0; vec_i < ylabels.size(); vec_i++)
    {
        channelDataMap.insert( std::pair<int,vector<double> >(vec_i, {0.,0.,0.,0.,0.}) );
    }
    if(ylabels.size() != ybins)
    {
        cout<<"Warning: ylabels.size() != ybins"<<endl;
    }

    // DRAW ORIGINAL PLOT AND INTEGRATE POWER FOR NEW PLOT
    int bincontents_length = bincontents.size();
    int d = 0;
    for(int y = 0; y < ybins; y++)
    {
        for(int x = 1; x <= xbins; x++)
        {
            //cout << x << " " << y << " " << d << endl;
            if(d<bincontents_length)
            {
                // INTEGRATE POWER
                for(int map_i = 0; map_i<nRanges; map_i++)
                {
                    if(x >= freqRangeMap[map_i].minBin && x <= freqRangeMap[map_i].maxBin)
                    {
                        channelDataMap[y][map_i]+=bincontents[d];
                    }
                }
                d++;
            }
        }
    }

    // NOTES ON OUTPUT FORMAT
    /*
     {
     "fft_int_0" : [], // compute from bincontents[] from input file
     "fft_int_1" : [],
     "fft_int_2" : [],
     "fft_int_3" : [],
     "fft_int_4" : [],
     "fft_int_label_0" : "FFT Integral 0-2 ADC^{2}-MHz", // create for output file
     "fft_int_label_1" : "FFT Integral 2-6 ADC^{2}-MHz",
     "fft_int_label_2" : "FFT Integral 6-9 ADC^{2}-MHz",
     "fft_int_label_3" : "FFT Integral 9-15 ADC^{2}-MHz",
     "fft_int_label_4" : "FFT Integral 36-48 ADC^{2}-MHz",
     "name" : "fft_integrals", // create for output file
     "runNumber" : 32148, // take from input file
     "runType" : 258, // take from input file
     "time" : "2017/11/21 03:12:49", // take from input file
     "timestamp" : 1511233969, // take from input file
     "ylabels" : [] // take from input file
     }
     */

    // BUILD UP JSON INFO
    Json::Value table;
    string key = "";

    for(int map_i = 0; map_i<nRanges; map_i++)
    {
        key = "fft_int_" + to_string(map_i);
        for(int y = 0; y < ybins; y++)
        {
            table[key.c_str()][y] = channelDataMap[y][map_i];
        }
    }
    for(int map_i = 0; map_i<nRanges; map_i++)
    {
        key = "fft_int_label_" + to_string(map_i);
        string tmpstr = "FFT Integral " + to_string(int(freqRangeMap[map_i].minF)) + "-" + to_string(int(freqRangeMap[map_i].maxF)) + " ADC^{2}-MHz";
        table[key.c_str()] = tmpstr.c_str();
    }
    table["name"] = "fft_integrals";
    table["runNumber"] = runNumber;
    table["runType"] = runType;
    table["time"] = time.c_str();
    table["timestamp"] = timestamp;
    for(int y = 0; y < ybins; y++)
    {
        table["ylabels"][y] = ylabels[y].c_str();
    }

    // OUTPUT JSON FILE
    string fname = "testfile.txt";
    Json::StyledStreamWriter writer;
    filebuf fb;
    fb.open(fname.c_str(), ios::out);
    if(fb.is_open()){
        ostream os(&fb);
        writer.write(os, table);
    }
    fb.close();

    cout << "F I N I S H E D" << " run " << run << endl;
}
