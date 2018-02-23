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
#include <TH1D.h>
#include <TLatex.h>

using namespace std;

inline int GetTauIndex(int tau) {
    return (tau/10)-6; // this will give 0 for 60, 1 for 70, ..., 19 for 250
}
inline int GetTauFromIndex(int i) {
    return (i+6)*10; // this will give 0 for 60, 1 for 70, ..., 19 for 250
}
inline double GetkeV(double fwhm, double muADC, double keV) {return fwhm*(keV/muADC);}

double getMean(double array[], int size)
{
    double mean = 0.0;
    double sum = 0.0;
    for(int i = 0; i<size; i++)
    {
        sum += array[i];
        //cout << i << " " << array[i] << " " << sum << endl;
    }
    mean = sum/size;
    return mean;
}
double getStdDev(double array[], int size)
{
    double stdDev = 0.0;
    double variance = 0.0;
    double temp = 0.0;
    double mean = 0.0;
    
    mean = getMean(array, size);
    
    for(int i = 0; i<size; i++)
    {
        temp += (mean-array[i])*(mean-array[i]);
        cout<<array[i]<<endl;
    }
    
    variance = temp/(size-1.0);
    stdDev = sqrt(variance);
    //cout << stdDev << endl;
    return stdDev;
}

int main (int argc, char* argv[])
{
///////////////////////////////
//// INITIALIZATIONS
///////////////////////////////
    TApplication *App = new TApplication("App", 0, NULL);

    if (argc < 2)
    {
        cout << "too few arguments " << argv[0] << endl;
        return 1;
    }
    int DS = atoi(argv[1]);
    
    // input file and its data fields
    string infilepath = "/global/u2/g/gilliss/ChargeTrapping/PlottingData/";
    string infilename = "";
    if(DS == 0) {infilename = infilepath + "OptimalFWHMvsTau_DS0.txt";}
    if(DS == 1) {infilename = infilepath + "OptimalFWHMvsTau_DS1.txt";}
    if(DS == 2) {infilename = infilepath + "OptimalFWHMvsTau_DS2.txt";}
    if(DS == 3) {infilename = infilepath + "OptimalFWHMvsTau_DS3.txt";}
    if(DS == 4) {infilename = infilepath + "OptimalFWHMvsTau_DS4.txt";}
    ifstream infile;
    
    // output ROOT file
    char fileName[50];
    sprintf(fileName,"TauvsChan_DS01234.root");
    TFile *f = new TFile(fileName,"UPDATE");
    
///////////////////////////////
//// READ IN DATA to BUILD CHANNEL LIST
///////////////////////////////
    
    vector<int> channel_vector;
    
    infile.open(infilename);
    if (infile.is_open()) {cout<<"  File is open"<<endl;}
    else {cout<<"Error: File did not open"<<endl; return 1;}
    
    int channel;
    double optTau, eoptTau, optFWHM, eoptFWHM;
    
    while (infile >> channel >> optTau >> eoptTau >> optFWHM >> eoptFWHM)
    {
        sort(channel_vector.begin(),channel_vector.end());
        if(!binary_search(channel_vector.begin(),channel_vector.end(),channel))
        {
            if(channel%2==0) {channel_vector.push_back(channel);} // only HG channels
        }
    }
    infile.close();
    if (!infile.is_open()) {cout<<"  File is closed"<<endl;}
    else {cout<<"Error: File did not closed"<<endl; return 1;}

    //for(int i=0;i<channel_vector.size();i++) {cout<<channel_vector[i]<<endl;}

    
///////////////////////////////
//// LOOP OVER CHANNEL LIST
///////////////////////////////
    int nChans = channel_vector.size();
    double channel_array[nChans];
    double echannel_array[nChans]; for(int i=0; i>nChans; i++) {echannel_array[i]=0.;}
    double optTau_array[nChans];
    double eoptTau_array[nChans];
    double optFWHM_array[nChans];
    double eoptFWHM_array[nChans];
    
    int chan = 0;
    for(unsigned int chan_i = 0; chan_i < channel_vector.size(); chan_i++)
    {
        chan = channel_vector[chan_i];
        
        ///////////////////////////////
        //// READ IN DATA
        ///////////////////////////////
        cout<<"--------- Chan "<<chan<<" ---------"<<endl;
        
        infile.open(infilename);
        if (infile.is_open()) {/*cout<<"  infile is open"<<endl;*/}
        else {cout<<"Error: infile did not open"<<endl; return 1;}
        while (infile >> channel >> optTau >> eoptTau >> optFWHM >> eoptFWHM)
        {
            if(channel == chan)
            {
                // put chan and fwhm in array
                channel_array[chan_i] = channel;
                optTau_array[chan_i] = optTau;
                eoptTau_array[chan_i] = eoptTau;
                optFWHM_array[chan_i] = optFWHM;
                eoptFWHM_array[chan_i] = eoptFWHM;
            }
        }
        infile.close();
        if (!infile.is_open()) {/*cout<<"  infile is closed"<<endl;*/}
        else {cout<<"Error: infile did not close"<<endl; return 1;}
    } // end channel list loop
    
///////////////////////////////
//// RESULTS AND GRAPHS
///////////////////////////////
    
    TCanvas *c = new TCanvas();
    c->cd();
    TGraphErrors *gr = new TGraphErrors(nChans,channel_array,optTau_array,echannel_array,eoptTau_array);
    gr->SetName(Form("TauvsChan_DS%d",DS));
    gr->SetTitle(Form("Optimal Tau, Each Detector at 2614 keV"));
    gr->GetYaxis()->SetTitle("Effective PZ Time Constant (#mu s)");
    gr->GetXaxis()->SetTitle("Detector");
    gr->GetXaxis()->SetLabelSize(0.0);
    gr->GetYaxis()->SetTitleOffset(1.2);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.7);
    gr->Draw("AP");
    c->Update();
    
    TCanvas *c1 = new TCanvas();
    c1->cd();
    TH1D *h = new TH1D(Form("TauHist_DS%d",DS),Form("TauHist_DS%d",DS),12,60,180);
    h->SetTitle("Optimal Time Constants");
    h->GetXaxis()->SetTitle("Effective PZ Time Constant (#mu s)");
    h->GetYaxis()->SetTitle("Number of Detectors");
    h->SetLineColor(kBlue);
    h->SetStats(0);
    for(int i = 0; i<nChans; i++) {h->Fill(optTau_array[i]);}
    h->Draw();
    
    TCanvas *c2 = new TCanvas();
    c2->cd();
    TGraphErrors *gr2 = new TGraphErrors(nChans,channel_array,optFWHM_array,echannel_array,eoptFWHM_array);
    gr2->SetName(Form("FWHMvsChan_DS%d",DS));
    gr2->SetTitle(Form("FWHM for Optimal Tau, Each Detector at 2614 keV"));
    gr2->GetYaxis()->SetTitle("FWHM (keV)");
    gr2->GetXaxis()->SetTitle("Detector");
    gr2->GetXaxis()->SetLabelSize(0.0);
    gr2->GetYaxis()->SetTitleOffset(1.2);
    gr2->SetMarkerColor(kBlue);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerSize(0.7);
    gr2->Draw("AP");
    c2->Update();
    
    TCanvas *c3 = new TCanvas();
    c3->cd();
    TH1D *h2 = new TH1D(Form("FWHMHist_DS%d",DS),Form("FWHMHist_DS%d",DS),60,2,8);
    h2->SetTitle("FWHM for Optimal Time Constants");
    h2->GetXaxis()->SetTitle("FWHM (keV)");
    h2->GetYaxis()->SetTitle("Number of Detectors");
    h2->SetLineColor(kBlue);
    h2->SetStats(0);
    for(int i = 0; i<nChans; i++) {h2->Fill(optFWHM_array[i]);}
    h2->Draw();

//    gr->Write();
//    h->Write();
//    gr2->Write();
//    h2->Write();
    
    //delete gr; gr = 0;
    
///////////////////////////////
//// MEAN AND STDDEV
///////////////////////////////
    
    cout<<"Mean = "<<getMean(optFWHM_array,nChans)<<endl;
    cout<<"StdDev = "<<getStdDev(optFWHM_array,nChans)<<endl;

///////////////////////////////
//// Close out
///////////////////////////////
    
    f->Close();

    cout << "F I N I S H E D" << endl;
    App->Run();
}