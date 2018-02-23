/*
[position] [channel] [energy parameter/decay constant] [mu] [mu unc] [sigma] [sigma unc] [fwhm] [fwhm unc] [ratio] [ratio unc]
 
http://stackoverflow.com/questions/7868936/read-file-line-by-line
https://root.cern.ch/root/html/tutorials/tree/basic.C.html
http://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
*/

#include <iostream>
#include <fstream>
#include <algorithm> // min_element, sort, binary_search
#include <vector>

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
inline double GetkeV(double fwhm, double muADC, double keV) {return fwhm*(keV/muADC);}

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
    
    // output ROOT file
    char fileName[50];
    // input file and its data fields
    string infilepath = "/global/u2/g/gilliss/ChargeTrapping/PlottingData/";
    string modulepath = "ModuleDataSets/";
    string stcpath = "STCDataSets/";"
    string infilename = "";
    int DS = 100;
    string STC = "None";
    if(argv[1] == "0") {
        infilename = infilepath + modulepath + "DS0_auxpeakfit_5064_5089.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argv[1] == "1") {
        infilename = infilepath + modulepath + "DS1_auxpeakfit_10529_10545.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argv[1] == "2") {
        infilename = infilepath + modulepath + "DS2_auxpeakfit_14568_14575.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argv[1] == "3") {
        infilename = infilepath + modulepath + "DS3_auxpeakfit_16836_16854.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argv[1] == "4") {
        infilename = infilepath + modulepath + "DS4_auxpeakfit_60001014_60001031.txt";
        DS = atoi(argv[1]);
        sprintf(fileName,"FWHMvsTau_DS%d.root",DS);
    }
    if(argv[1] == "P3CLR") {
        infilename = infilepath + stcpath + "P3CLR_auxpeakfit_30000392_30000492.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argv[1] == "P3DCR") {
        infilename = infilepath + stcpath + "P3DCR_auxpeakfit_30001706_30001806.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argv[1] == "P3DNQ") {
        infilename = infilepath + stcpath + "P3DNQ_auxpeakfit_30002032_30002132.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argv[1] == "P3EVV") {
        infilename = infilepath + stcpath + "P3EVV_auxpeakfit_30004106_30004206.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argv[1] == "P3FPV") {
        infilename = infilepath + stcpath + "P3FPV_auxpeakfit_30005405_30005413.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }
    if(argv[1] == "P3FPW") {
        infilename = infilepath + stcpath + "P3FPW_auxpeakfit_30004968_30005020.txt";
        STC = argv[1];
        sprintf(fileName,"FWHMvsTau_Stc%s.root",argv[1]);
    }

    ifstream infile;
    int position, channel;
    string stau;
    double mu, emu, sigma, esigma, fwhm, efwhm, ratio, eratio;
    double keV = 2614.533;
    // output file for chan : optFWHM : optTau
    ofstream outfile;
    if (DS != 100) {string outfilename(Form("OptimalFWHMvsTau_DS%d.txt",DS));}
    if (STC != "None") {string outfilename(Form("OptimalFWHMvsTau_Stc%s.txt",STC.c_str()));}
    outfile.open(outfilename);
    // output ROOT file
    TFile *f = new TFile(fileName,"RECREATE");
    // tau
    int ntaus = 20;
    double tau_array[ntaus];
    double etau_array[ntaus]; for(int i=0; i>ntaus; i++) {etau_array[i]=0.;}
    double temptau = 0;
    int tauIndex = 0;
    int taupos = 7; // for reading in tau from string
    string delim = "us"; // for reading in tau from string
    int delimpos = 0; // for reading in tau from string
    // fwhm
    double fwhm_array[ntaus];
    double efwhm_array[ntaus];
    
///////////////////////////////
//// READ IN DATA to BUILD CHANNEL LIST
///////////////////////////////
    
    vector<int> channel_vector;
    
    infile.open(infilename);
    if (infile.is_open()) {cout<<"  File is open"<<endl;}
    else {cout<<"Error: File did not open"<<endl; return 1;}
    
    while (infile >> position >> channel >> stau >> mu >> emu >> sigma >> esigma >> fwhm >> efwhm >> ratio >> eratio)
    {
        sort(channel_vector.begin(),channel_vector.end());
        if(!binary_search(channel_vector.begin(),channel_vector.end(),channel))
        {
            channel_vector.push_back(channel);
        }
    }
    infile.close();
    if (!infile.is_open()) {cout<<"  File is closed"<<endl;}
    else {cout<<"Error: File did not closed"<<endl; return 1;}

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
        
        while (infile >> position >> channel >> stau >> mu >> emu >> sigma >> esigma >> fwhm >> efwhm >> ratio >> eratio)
        {
            if(channel == chan)
            {
                // read tau
                delimpos = stau.find(delim);
                temptau = stod(stau.substr(taupos,delimpos-taupos));
                // put tau in array
                tauIndex = GetTauIndex(temptau);
                tau_array[tauIndex] = temptau; // etau_array[] is full of zeros
                // put fwhm in array
                fwhm_array[tauIndex] = GetkeV(fwhm,mu,keV);
                efwhm_array[tauIndex] = GetkeV(efwhm,mu,keV);
            }
        }
        infile.close();
        if (!infile.is_open()) {cout<<"  infile is closed"<<endl;}
        else {cout<<"Error: infile did not close"<<endl; return 1;}
        
        ///////////////////////////////
        //// RESULTS AND GRAPHS
        ///////////////////////////////
        double optFWHM = *min_element(fwhm_array,fwhm_array+ntaus);
        double eoptFWHM = 0.;
        cout<<"  Best fwhm = "<<optFWHM<<endl;
        double optTau = 0.;
        double eoptTau = 0.;
        for(int i = 0; i<ntaus; i++)
        {
            if(fwhm_array[i]==optFWHM)
            {
                optTau = GetTauFromIndex(i);
                eoptFWHM = efwhm_array[i];
            }
        }
        cout<<"  Best tau = "<<optTau<<endl;
        
        if (outfile.is_open())
        {
            outfile<<chan<<" "<<optTau<<" "<<eoptTau<<" "<<optFWHM<<" "<<eoptFWHM<<endl;
        }
        else {cout<<"Error: outfile is not open"<<endl; return 1;}
        
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