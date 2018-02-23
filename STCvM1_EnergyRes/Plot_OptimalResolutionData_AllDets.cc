/* To add the corresponding tau values to each point on the TGraph, see https://root.cern.ch/phpBB3/viewtopic.php?t=4083
*/
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>

#include <iostream>
#include <fstream>
#include <stdio.h>
//#include <array>

#include <TApplication.h>

using namespace std;

char* detectorSNs(int a)
{
    char *detName = new char;

        if(a==0) detName="B8455";
        if(a==1) detName="B8463";
        if(a==2) detName="B8465";
        if(a==3) detName="B8469";
        if(a==4) detName="B8470";
        if(a==5) detName="B8474";
        if(a==6) detName="B8477";
        if(a==7) detName="B8480";
        if(a==8) detName="B8482";
        if(a==9) detName="B8487";
        if(a==10) detName="P42537A";
        if(a==11) detName="P42538A";
        if(a==12) detName="P42538B";
        if(a==13) detName="P42573A";
        if(a==14) detName="P42573B";
        if(a==15) detName="P42574A";
        if(a==16) detName="P42574B";
        if(a==17) detName="P42574C";
        if(a==18) detName="P42575A";
        if(a==19) detName="P42575B";
        if(a==20) detName="P42661A";
        if(a==21) detName="P42661B";
        if(a==22) detName="P42661C";
        if(a==23) detName="P42662A";
        if(a==24) detName="P42662B";
        if(a==25) detName="P42662C";
        if(a==26) detName="P42664A";
        if(a==27) detName="P42665A";
        if(a==28) detName="P42665B";
        if(a==29) detName="P42665C";
        if(a==30) detName="P42698A";
        if(a==31) detName="P42698B";
        if(a==32) detName="P42712B";

    return detName;
}

int main (int argc, char* argv[])
{
	TApplication *App = new TApplication("App", 0, NULL);
    
    ////////////////////////
    // ROOT FILE and TREE //
    ////////////////////////

    TFile *rootfile = TFile::Open("/Users/tomgilliss/Desktop/STCvM1_EnergyRes/ResolutionStudy_AllDets_STC/main/OptimalResolutionData.root");
    TTree *t = (TTree*)rootfile->Get("ntuple");
    Float_t detectorSN, systemBOOL, stringNumber, detectorNumber, tau, fwhm; // The TNtuple saved these ints as floats
    t->SetBranchAddress("detectorSN",&detectorSN);
    t->SetBranchAddress("systemBOOL",&systemBOOL);
    t->SetBranchAddress("stringNumber",&stringNumber);
    t->SetBranchAddress("detectorNumber",&detectorNumber);
    t->SetBranchAddress("tau",&tau);
    t->SetBranchAddress("fwhm",&fwhm);

    int nentries = t->GetEntries();
    //cout << "nentries " << nentries << endl;

    /////////////////////////
    // GRAPHING PARAMETERS //
    /////////////////////////
    int npoints=33;
    if (nentries!=(2*npoints)) {
        cout<<"Error: Size of tree does not match number of detectors to be plotted."<<endl;
        App->Run();
    }
    Double_t xData[npoints];
    Double_t STCData[npoints];
    Double_t M1Data[npoints];
    
    ///////////////////////////
    // LOAD DATA INTO ARRAYS //
    ///////////////////////////
    int isystemBOOL,istringNumber,idetectorNumber,itau;
    int STCCounter = 0;
    int M1Counter = 0;
    Double_t STCAverage = 0.0;
    Double_t M1Average = 0.0;
    Float_t fwhmProxy;
    for(int i=0;i<nentries;i++)
    {
        t->GetEntry(i);
        if (systemBOOL==0) {
            STCData[STCCounter]=fwhm;
            //cout<<i<<" "<<STCCounter<<" "<<fwhm<<" "<<STCData[STCCounter]<<endl;
            STCCounter++;
            fwhmProxy = fwhm;
            if (fwhmProxy>0) STCAverage = STCAverage + fwhm;
        }
        else {
            M1Data[M1Counter]=fwhm;
            //cout<<i<<" "<<M1Counter<<" "<<fwhm<<" "<<M1Data[M1Counter]<<endl;
            M1Counter++;
            fwhmProxy = fwhm;
            if (fwhmProxy>0) M1Average = M1Average + fwhm;
        }
    }
    for (int i = 0; i < npoints; i++) {
        xData[i]=i+1;
    }
    cout<<"This is not right: STC Avg FWHM = "<<STCAverage/33.0<<" M1 Avg FWHM = "<<M1Average/33.0<<endl;

    //cout<<"size xData,STCData,M1Data: "<<sizeof(xData)/sizeof(xData[0])<<", "<<sizeof(STCData)/sizeof(STCData[0])<<" "<<sizeof(M1Data)/sizeof(M1Data[0])<< endl;
    
    //////////////
    // GRAPHING //
    //////////////

    TCanvas *c = new TCanvas("myCanvas","myCanvas",1000,600);
    c->cd();

    TGraph *gSTC = new TGraph(npoints,xData,STCData);
    gSTC->SetTitle("STC vs M1 Energy Resolution");
    gStyle->SetTitleY(0.96);
    gSTC->GetYaxis()->SetTitle("FWHM (keV)");
    //gSTC->GetXaxis()->SetTitle("Detector SN");
    gSTC->GetYaxis()->SetRangeUser(2,8);
    gSTC->GetXaxis()->SetRangeUser(0,34);
    gSTC->GetXaxis()->SetTickLength(0.);
    gSTC->SetMarkerStyle(21);
    cout << "Made it here 8" << endl;
    gSTC->SetMarkerColor(38);
    gSTC->SetLineColor(0);
    
    TGraph *gM1 = new TGraph(npoints,xData,M1Data);
    
    //TAxis *a = gSTC->GetHistogram()->GetXaxis();
    //Double_t x1 = a->GetBinLowEdge(1);
    //Double_t x2 = a->GetBinUpEdge(a->GetNbins());
    //gSTC->GetHistogram()->GetXaxis()->Set(npoints,x1,x2);
    TAxis *xaxis = gSTC->GetXaxis();
    xaxis->CenterLabels();
    int bins[npoints];
    for (int i = 1; i<34; i++) {
        bins[i-1]=xaxis->FindBin(i);
        //cout <<xaxis->FindBin(i)<<" "<<bins[i-1]<<endl;
    }
    for (int i = 0; i < npoints; i++) {
        gSTC->GetHistogram()->GetXaxis()->SetBinLabel(bins[i],detectorSNs(i));//https://root.cern.ch/phpBB3/viewtopic.php?t=8924
    }
    
    gSTC->Draw("ALP");
    
    gM1->SetMarkerStyle(22);
    gM1->SetMarkerColor(46);
    gM1->SetLineColor(0);
    gM1->Draw("CP");
    
    TLegend* l = new TLegend(0.77,0.6,0.87,0.7); // x1,y1,x2,y2
    l->SetBorderSize(0);
    l->AddEntry(gSTC,"STC","p");
    l->AddEntry(gM1,"M1","p");
    l->Draw();
    
    c->Update();

    char titleForFile[200];
    sprintf(titleForFile,"M1vsSTC_EnergyResolutions.pdf");
    //c->Print(titleForFile);
    
    /////////////////
    // CLOSING OUT //
    /////////////////
    rootfile->Close();
    
	cout << "D O N E" << endl;
	App->Run(); 
}


////cout << "Made it here 1" << endl;
////TGraph *gSTC = new TGraph(npoints,xData,xData);
////gSTC->Draw("ALP");
////c->Update();
////for (int i = 0; i < npoints; i++) {
//cout<<xData[i]<<" "<<yData[i]<<endl;
//cout<<xData[i]<<" "<<STCData[i]<<endl;
//cout<<xData[i]<<" "<<M1Data[i]<<endl;
////}