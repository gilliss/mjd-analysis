/*
A script to read Pinghan's data and deterimine the optimal tau/FWHM for each detector in a DS or STC.

Pinghan's .txt data is in the following format for DSs:
[position] [channel] [energy parameter/decay constant] [mu] [mu unc] [sigma] [sigma unc] [fwhm] [fwhm unc] [ratio] [ratio unc]
 
Pinghan's .txt data is in the following format for STCs:
[channel] [energy parameter/decay constant] [mu] [mu unc] [sigma] [sigma unc] [fwhm] [fwhm unc] [ratio] [ratio unc]
 
Output .txt file is in the following format:
[chan] [optTau] [eoptTau] [optFWHM] [eoptFWHM] [qtTau]
 
Example use: ./FindOptimal <DS#/STCserial>
*/

#include <iostream>
#include <fstream>
#include <algorithm> // min_element, sort, binary_search
#include <vector>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TFile.h>
#include <TLatex.h>

using namespace std;

inline int GetTauIndex(int tau) {
    return (tau/10)-6; // this will give 0 for 60, 1 for 70, ..., 19 for 250
}
inline int GetTauFromIndex(int i) {
    return (i+6)*10; // this will give 0 for 60, 1 for 70, ..., 19 for 250
}
inline double GetkeV(double fwhmADC, double muADC, double keV) {return fwhmADC*(keV/muADC);}
inline double ConvertToQTTau(int tPZ) {
    double tRC = 72.6; // us
    double tQT = (tRC*tPZ)/(tPZ-tRC); // tPZ = (tRC*tQT)/(tQT-tRC)
    return tQT;
}

int main (int argc, char* argv[])
{
///////////////////////////////
//// INITIALIZATIONS
///////////////////////////////
    //TApplication *App = new TApplication("App", 0, NULL);
    if (argc < 2)
    {
        cout << "too few arguments " << argv[0] << endl;
        return 1;
    }
    string argStr(argv[1]);
    
    // output ROOT file //
    char fileName[50];
    
    // input file and its data fields //
    string infilepath = "/global/u2/g/gilliss/ChargeTrapping/Update_March2017/";
    string modulepath = "ModuleDataSets/";
    string stcpath = "STCDataSets/";
    string infilename = "";
    int DS = 100;
    string STC = "None";
    if(argStr == "0") {
        infilename = infilepath + modulepath + "DS0_auxpeakfit_5064_5089.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argStr == "1") {
        infilename = infilepath + modulepath + "DS1_auxpeakfit_10529_10545.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argStr == "2") {
        infilename = infilepath + modulepath + "DS2_auxpeakfit_14568_14575.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argStr == "3") {
        infilename = infilepath + modulepath + "DS3_auxpeakfit_16836_16854.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argStr == "4") {
        infilename = infilepath + modulepath + "DS4_auxpeakfit_60001014_60001031.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argStr == "P3CLR") {
        infilename = infilepath + stcpath + "P3CLR_auxpeakfit_30000392_30000492.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argStr == "P3DCR") {
        infilename = infilepath + stcpath + "P3DCR_auxpeakfit_30001706_30001806.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argStr == "P3DNQ") {
        infilename = infilepath + stcpath + "P3DNQ_auxpeakfit_30002032_30002132.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argStr == "P3EVV") {
        infilename = infilepath + stcpath + "P3EVV_auxpeakfit_30004106_30004206.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argStr == "P3FPV") {
        infilename = infilepath + stcpath + "P3FPV_auxpeakfit_30005405_30005413.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argStr == "P3FPW") {
        infilename = infilepath + stcpath + "P3FPW_auxpeakfit_30004968_30005020.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    // input file from Pinghan //
    ifstream infile;
    int position, channel;
    string stau;
    double mu, emu, sigma, esigma, fwhm, efwhm, ratio, eratio;
    double keV = 2614.533; // the peak that Pinghan fit was the 208Tl 2615 peak
    
    // output file for chan : optFWHM : optTau : qtTau //
    ofstream outfile;
    string outfilename;
    if (DS != 100) {outfilename = "OptimalFWHMvsTau_DS"+to_string(DS)+".txt";}
    if (STC != "None") {outfilename = "OptimalFWHMvsTau_Stc"+STC+".txt";}
    //if (DS != 100) {string outfilename(Form("OptimalFWHMvsTau_DS%d.txt",DS));}
    //if (STC != "None") {string outfilename(Form("OptimalFWHMvsTau_Stc%s.txt",STC.c_str()));}
    outfile.open(outfilename);
    if (outfile.is_open()) {cout<<"  outfile is open"<<endl;}
    
    // output ROOT file //
    TFile *f = new TFile(fileName,"RECREATE");
    
    // tau handling //
    int ntaus = 20;
    if (DS != 100) {ntaus = 20;}
    if (STC != "None") {ntaus = 10;}
    double tau_array[ntaus];
    double etau_array[ntaus]; for(int i=0; i>ntaus; i++) {etau_array[i]=0.;}
    double temptau = 0;
    int tauIndex = 0;
    int taupos = 7; // for reading in tau from string
    string delim = "us"; // for reading in tau from string
    int delimpos = 0; // for reading in tau from string
    // fwhm handling //
    double fwhm_array[ntaus];
    double efwhm_array[ntaus];
    
///////////////////////////////
//// READ IN DATA to BUILD ORDERED CHANNEL LIST
///////////////////////////////
    
    vector<int> channel_vector;
    
    cout<<"  infilename = "<<infilename<<endl;
    
    infile.open(infilename);
    if (infile.is_open()) {cout<<"  infile is open"<<endl;}
    else {cout<<"Error: infile did not open"<<endl; return 1;}
    // READ IN DS DATA //
    if(DS != 100)
    {
    while (infile >> position >> channel >> stau >> mu >> emu >> sigma >> esigma >> fwhm >> efwhm >> ratio >> eratio)
    {
        sort(channel_vector.begin(),channel_vector.end()); // sort the list into chronological order
        if(!binary_search(channel_vector.begin(),channel_vector.end(),channel)) // if channel not yet in vector
        {
            channel_vector.push_back(channel);
        }
    }
    }
    // READ IN STC DATA //
    if(STC != "None")
    {
        while (infile >> channel >> stau >> mu >> emu >> sigma >> esigma >> fwhm >> efwhm >> ratio >> eratio)
        {
            sort(channel_vector.begin(),channel_vector.end()); // sort the list into chronological order
            if(!binary_search(channel_vector.begin(),channel_vector.end(),channel)) // if channel not yet in vector
            {
                channel_vector.push_back(channel);
            }
        }
    }
    infile.close();
    if (!infile.is_open()) {cout<<"  infile is closed"<<endl;}
    else {cout<<"Error: infile did not close"<<endl; return 1;}
    //for(int i=0;i<channel_vector.size();i++) {cout<<channel_vector[i]<<endl;}

    
///////////////////////////////
//// LOOP OVER CHANNEL LIST
///////////////////////////////
    int chan = 0;
    for(unsigned int chan_i = 0; chan_i < channel_vector.size(); chan_i++)
    {
        chan = channel_vector[chan_i];
        
        ///////////////////////////////
        //// READ IN DATA
        ///////////////////////////////
        cout<<"--------- Chan "<<chan<<" ---------"<<endl;
        
        infile.open(infilename);
        if (infile.is_open()) {cout<<"  infile is open"<<endl;}
        else {cout<<"Error: infile did not open"<<endl; return 1;}
        
        // READ IN DS DATA //
        if(DS != 100)
        {
        while (infile >> position >> channel >> stau >> mu >> emu >> sigma >> esigma >> fwhm >> efwhm >> ratio >> eratio)
        {
            if(channel == chan)
            {
                // read tau from string in .txt data file //
                delimpos = stau.find(delim);
                temptau = stod(stau.substr(taupos,delimpos-taupos));
                // put tau in correctly ordered spot in array //
                tauIndex = GetTauIndex(temptau);
                tau_array[tauIndex] = temptau; // note that etau_array[] is full of zeros
                // put fwhm in array in the corresponding spot //
                fwhm_array[tauIndex] = fwhm;//GetkeV(fwhm,mu,keV);
                efwhm_array[tauIndex] = efwhm;//GetkeV(efwhm,mu,keV);
            }
        }
        }
        if(STC != "None")
        {
            while (infile >> channel >> stau >> mu >> emu >> sigma >> esigma >> fwhm >> efwhm >> ratio >> eratio)
            {
                if(channel == chan)
                {
                    // read tau from string in .txt data file //
                    delimpos = stau.find(delim);
                    temptau = stod(stau.substr(taupos,delimpos-taupos));
                    // put tau in correctly ordered spot in array //
                    tauIndex = GetTauIndex(temptau);
                    tau_array[tauIndex] = temptau; // note that etau_array[] is full of zeros
                    // put fwhm in array in the corresponding spot //
                    fwhm_array[tauIndex] = fwhm;//GetkeV(fwhm,mu,keV);
                    efwhm_array[tauIndex] = efwhm;//GetkeV(efwhm,mu,keV);
                }
            }
        }
        infile.close();
        if (!infile.is_open()) {cout<<"  infile is closed"<<endl;}
        else {cout<<"Error: infile did not close"<<endl; return 1;}
        
        ///////////////////////////////
        //// RESULTS AND GRAPHS
        ///////////////////////////////
        // Find optimal tau and fwhm ///
        //for(int k=0; k<ntaus; k++){cout<<"fwhm_array.push_back("<<fwhm_array[k]<<")"<<endl;}
        double optFWHM = *min_element(fwhm_array,fwhm_array+ntaus);
        double eoptFWHM = 0.;
        cout<<"  Best fwhm = "<<optFWHM<<endl;
        double optTau = 0.;
        double eoptTau = 0.; // might need to devise some method of finding this error (it'll be <= +/- 10us)
        for(int i = 0; i<ntaus; i++) // search for which tau had the optFWHM
        {
            if(fwhm_array[i]==optFWHM)
            {
                optTau = GetTauFromIndex(i);
                eoptFWHM = efwhm_array[i];
            }
        }
        cout<<"  Best tau = "<<optTau<<endl;
        
        // Calculate qtTau by inverting the optimal effPZ tau //
        double qtTau = ConvertToQTTau(optTau);
        if(qtTau<0.){qtTau = 1000.;} // when approximately optimal tauPZ is < tauRC
        
        // Write to output file //
        if (outfile.is_open())
        {
            outfile<<chan<<" "<<optTau<<" "<<eoptTau<<" "<<optFWHM<<" "<<eoptFWHM<<" "<<qtTau<<endl;
        }
        else {cout<<"Error: outfile is not open"<<endl; return 1;}
        
        // Plot FWHM vs tau //
        //TCanvas *c = new TCanvas();
        //c->cd();
        TGraphErrors *gr = new TGraphErrors(ntaus,tau_array,fwhm_array,etau_array,efwhm_array);
        gr->SetName(Form("FWHMvsTau_%d",chan));
        gr->SetTitle(Form("FWHM vs Tau at 2614 keV"));
        gr->GetYaxis()->SetTitle("FWHM (keV)");
        gr->GetXaxis()->SetTitle("Effective PZ Time Constant (#mu s)");
        gr->SetMarkerColor(4);
        gr->SetMarkerStyle(21);
        //gr->Draw("AP");
        //c->Update();
        gr->Write();
        delete gr; gr = 0;
    } // end channel list loop

///////////////////////////////
//// Close out
///////////////////////////////
    outfile.close();
    if (!outfile.is_open()) {cout<<"  outfile is closed"<<endl;}
    else {cout<<"Error: outfile did not close"<<endl; return 1;}
    
    f->Close();
    delete f; f = 0;
    cout << "F I N I S H E D" << endl;
    //App->Run();
}

/*
References:
 http://stackoverflow.com/questions/7868936/read-file-line-by-line
 https://root.cern.ch/root/html/tutorials/tree/basic.C.html
 http://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
*/