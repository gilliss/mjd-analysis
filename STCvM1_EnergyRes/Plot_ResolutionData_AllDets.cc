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

char* detName(int a,int b,int c)
{
    char *detName = new char;
    if(a==0 && b==0 && c==1) sprintf(detName,"B8455");
    if(a==0 && b==0 && c==2) sprintf(detName,"B8470");
    if(a==0 && b==0 && c==3) sprintf(detName,"B8463");
    if(a==0 && b==0 && c==4) sprintf(detName,"B8465");
    if(a==0 && b==0 && c==5) sprintf(detName,"B8469");
    if(a==0 && b==1 && c==1) detName="P42538A";
    if(a==0 && b==1 && c==2) detName="P42575A";
    if(a==0 && b==1 && c==3) detName="P42661C";
    if(a==0 && b==1 && c==4) detName="P42698B";
    if(a==0 && b==2 && c==1) sprintf(detName,"P42662C");
    if(a==0 && b==2 && c==2) sprintf(detName,"P42665A");
    if(a==0 && b==2 && c==3) sprintf(detName,"P42664A");
    if(a==0 && b==2 && c==4) sprintf(detName,"B8474");
    if(a==0 && b==3 && c==1) sprintf(detName,"P42573B");
    if(a==0 && b==3 && c==2) sprintf(detName,"P42661A");
    if(a==0 && b==3 && c==3) sprintf(detName,"P42575B");
    if(a==0 && b==3 && c==4) sprintf(detName,"B8477");
    if(a==0 && b==4 && c==1) sprintf(detName,"P42573A");
    if(a==0 && b==4 && c==2) sprintf(detName,"P42698A");
    if(a==0 && b==4 && c==3) sprintf(detName,"B8480");
    if(a==0 && b==4 && c==4) sprintf(detName,"P42538B");
    if(a==0 && b==5 && c==1) sprintf(detName,"P42662A");
    if(a==0 && b==5 && c==2) sprintf(detName,"P42574A");
    if(a==0 && b==5 && c==3) sprintf(detName,"P42661B");
    if(a==0 && b==5 && c==4) sprintf(detName,"P42574C");
    if(a==0 && b==6 && c==1) sprintf(detName,"P42537A");
    if(a==0 && b==6 && c==2) sprintf(detName,"P42662B");
    if(a==0 && b==6 && c==3) sprintf(detName,"P42574B");
    if(a==0 && b==6 && c==4) sprintf(detName,"B8482");
    if(a==0 && b==7 && c==1) sprintf(detName,"P42665C");
    if(a==0 && b==7 && c==2) sprintf(detName,"P42712B");
    if(a==0 && b==7 && c==3) sprintf(detName,"P42665B");
    if(a==0 && b==7 && c==4) sprintf(detName,"B8487");
    //if(a==0 && b== && c==) sprintf(detName,"");
    return detName;
}

int main (int argc, char* argv[])
{
	TApplication *App = new TApplication("App", 0, NULL);
    
    ////////////////////////
    // ROOT FILE and TREE //
    ////////////////////////

    TFile *rootfile = TFile::Open("/Users/tomgilliss/Desktop/STCvM1_EnergyRes/ResolutionStudy_AllDets_STC/main/ResolutionData_AllDets.root");
    TTree *t = (TTree*)rootfile->Get("ntuple");
    Float_t systemBOOL, stringNumber, detectorNumber, channel, fixedPtBOOL, tauIndex, tau, peak; // The TNtuple saved these ints as floats
    Float_t fwhm;
    t->SetBranchAddress("systemBOOL",&systemBOOL);
    t->SetBranchAddress("stringNumber",&stringNumber);
    t->SetBranchAddress("detectorNumber",&detectorNumber);
    t->SetBranchAddress("channel",&channel);
    t->SetBranchAddress("fixedPtBOOL",&fixedPtBOOL);
    t->SetBranchAddress("tauIndex",&tauIndex);
    t->SetBranchAddress("tau",&tau);
    t->SetBranchAddress("peak",&peak);
    t->SetBranchAddress("fwhm",&fwhm);

    int nentries = t->GetEntries();
    cout << "nentries " << nentries << endl;

    /////////////////////////
    // GRAPHING PARAMETERS //
    /////////////////////////
    Int_t nsystems = 2; // number of systems (STC=0,M1=1)
    Int_t nstrings = 8; // number of possible strings in each system
    Int_t ndetectors = 5; // number of possible detectors in each string
    Int_t npeaks = 3; // number of possible peaks studied for a detector
    Int_t ntaus = 4; // number of possible tau values studied for each peak
    Int_t ndetectorsSTC=33;
    
    Float_t x[nsystems][nstrings][ndetectors][npeaks][3];
    Float_t y[nsystems][nstrings][ndetectors][npeaks][3];
    /// Float_t xMWF[nsystems][nstrings][ndetectors][npeaks];
    /// Float_t yMWF[nsystems][nstrings][ndetectors][npeaks];
    Double_t yMWF[nstrings][ndetectors][npeaks];
    
    for (int nsys = 0; nsys < nsystems ; nsys++) {
        for (int nstr = 0; nstr < nstrings; nstr++) {
            for (int ndet = 0; ndet < ndetectors; ndet++) {
                for (int npks = 0; npks < npeaks; npks++) {
                    /// xMWF[nsys][nstr][ndet][npks]=0;
                    /// yMWF[nsys][nstr][ndet][npks]=0;
                    for (int ntau = 0; ntau < 3; ntau++) {
                        x[nsys][nstr][ndet][npks][ntau]=0;
                        y[nsys][nstr][ndet][npks][ntau]=0;
                    }
                }
            }
        }
    }
    
    TGraph* g[nsystems][nstrings][ndetectors][npeaks];
    TGraph* gMWF[nsystems][nstrings][ndetectors][npeaks];
    TCanvas* c = new TCanvas("myCanvas","myCanvas",1000,400); //(name,title,width,height)
    TLegend* l = new TLegend(0.55,0.6,0.86,0.86); // x1,y1,x2,y2  https://root.cern.ch/doc/master/classTLegend.html
    l->SetBorderSize(0); // http://www.desy.de/~gbrandt/root/mkraemer.pdf
    
    ///////////////////////////
    // LOAD DATA INTO ARRAYS //
    ///////////////////////////
    int isystemBOOL,istringNumber,idetectorNumber,ipeak,itau;
    for(int i=0;i<nentries;i++)
    {
        t->GetEntry(i);
        
        isystemBOOL=(int) systemBOOL;
        istringNumber=(int) stringNumber;
        idetectorNumber=(int) (detectorNumber-1); // note index and value NOT equivalent
        ipeak=(int) peak;
        itau=(int) tauIndex;
        
        if (tauIndex==0||tauIndex==1||tauIndex==3) {
            if (tauIndex==3) itau = 2;
            x[isystemBOOL][istringNumber][idetectorNumber][ipeak][itau]=tau;
            y[isystemBOOL][istringNumber][idetectorNumber][ipeak][itau]=fwhm;
        }
        else if (tauIndex==2) {
            /// xMWF[isystemBOOL][istringNumber][idetectorNumber][ipeak]=tau;
            /// yMWF[isystemBOOL][istringNumber][idetectorNumber][ipeak]=fwhm;
            yMWF[istringNumber][idetectorNumber][ipeak]=fwhm;
        }
    }
    
    //cout<<"size x,y "<<sizeof(x)/sizeof(x[0][0][0][0][0])<<","<<sizeof(y)/sizeof(y[0][0][0][0][0])<<" size xMWF,yMWF "<<sizeof(xMWF)/sizeof(xMWF[0][0][0][0])<<","<<sizeof(yMWF)/sizeof(yMWF[0][0][0][0])<<endl;
    //cout<<"size x[0][0][0][0] "<<sizeof(x[0][0][0][0])/sizeof(x[0][0][0][0][0])<<endl;
    //cout<<"size xMWF[0][0][0] "<<sizeof(xMWF[0][0][0])/sizeof(xMWF[0][0][0][0])<<endl;

    
    /*
    /// Test the values in the arrays ///
    for (int nsys = 0; nsys < nsystems ; nsys++) {
        for (int nstr = 0; nstr < nstrings; nstr++) {
            for (int ndet = 0; ndet < ndetectors; ndet++) {
                for (int npks = 0; npks < npeaks; npks++) {
                    //cout<<"System "<<nsys<<" String "<<nstr<<" Detector "<<ndet+1<<" Tau "<<xMWF[nsys][nstr][ndet][npks]<<" Peak "<<npks<<" FWHM "<<yMWF[nsys][nstr][ndet][npks]<<endl;
                    for (int ntau = 0; ntau < 3; ntau++) {
                        cout<<"System "<<nsys<<" String "<<nstr<<" Detector "<<ndet+1<<" Tau "<<x[nsys][nstr][ndet][npks][ntau]<<" Peak "<<npks<<" FWHM "<<y[nsys][nstr][ndet][npks][ntau]<<endl;
                    }
                }
            }
        }
    }
    */
    
    //////////////
    // GRAPHING //
    //////////////
    int ndetForGraph = 4;
    int sysForGraph = 0;
    int strForGraph = 7;
    char titleForGraph[200], peaktitleForGraph[20];
    int colorarray[5]={1,2,3,4,6}; // https://root.cern.ch/doc/master/classTAttMarker.html

    c->Divide(npeaks,1,0.001,0.001,0); //(cols,rows,,,)
    int divideIndex = 0;
    for (int i = 0; i < ndetForGraph; i++) {
        for (int j = 0; j < npeaks; j++) {
            divideIndex++;
            c->cd(divideIndex);
            
            g[sysForGraph][strForGraph][i][j] = new TGraph(3,x[sysForGraph][strForGraph][i][j],y[sysForGraph][strForGraph][i][j]);

            /// //cout<<"xMWF,yMWF "<<xMWF[sysForGraph][strForGraph][i][j]<<","<<yMWF[sysForGraph][strForGraph][i][j]<<endl;
            gMWF[sysForGraph][strForGraph][i][j] = new TGraph(1);
            /// //Double_t dxMWF = xMWF[sysForGraph][strForGraph][i][j];
            /// //Double_t dyMWF = yMWF[sysForGraph][strForGraph][i][j];
            /// gMWF[sysForGraph][strForGraph][i][j]->SetPoint(0,xMWF[sysForGraph][strForGraph][i][j],yMWF[sysForGraph][strForGraph][i][j]);
            gMWF[sysForGraph][strForGraph][i][j]->SetPoint(0,70000,yMWF[strForGraph][i][j]);
            
            if (i==0) {
                if (j==0) sprintf(peaktitleForGraph,"^{214}Pb");
                if (j==1) sprintf(peaktitleForGraph,"^{40}K");
                if (j==2) sprintf(peaktitleForGraph,"^{208}Tl");
                sprintf(titleForGraph,"STC S%d %s",strForGraph,peaktitleForGraph);
                g[sysForGraph][strForGraph][i][j]->SetTitle(titleForGraph);
                gStyle->SetTitleY(0.95); // adjust main title vertical position thru gStyle pointer https://root.cern.ch/doc/master/classTStyle.html
                g[sysForGraph][strForGraph][i][j]->GetYaxis()->SetTitle("FWHM (keV)");
                g[sysForGraph][strForGraph][i][j]->GetXaxis()->SetTitle("Tau (ns)");
                g[sysForGraph][strForGraph][i][j]->GetYaxis()->SetTitleSize(.04);
                g[sysForGraph][strForGraph][i][j]->GetXaxis()->SetTitleSize(.04);
                g[sysForGraph][strForGraph][i][j]->GetYaxis()->SetLabelSize(.04);
                g[sysForGraph][strForGraph][i][j]->GetXaxis()->SetLabelSize(.04);
                g[sysForGraph][strForGraph][i][j]->GetXaxis()->SetNdivisions(6);
                g[sysForGraph][strForGraph][i][j]->GetYaxis()->SetRangeUser(1,7.5); //(1,8)
                g[sysForGraph][strForGraph][i][j]->GetXaxis()->SetRangeUser(40000,100000);
            }
            
            g[sysForGraph][strForGraph][i][j]->SetMarkerStyle(4);
            g[sysForGraph][strForGraph][i][j]->SetMarkerColor(colorarray[i]);
            g[sysForGraph][strForGraph][i][j]->SetLineColor(colorarray[i]);
            //g[sysForGraph][strForGraph][i][j]->SetLineWidth(2);
            
            gMWF[sysForGraph][strForGraph][i][j]->SetMarkerColor(colorarray[i]);
            gMWF[sysForGraph][strForGraph][i][j]->SetMarkerStyle(3);
            
            if(j==0) l->AddEntry(g[sysForGraph][strForGraph][i][j],detName(sysForGraph,strForGraph,i+1),"l");
            
            if (i==0) {
                g[sysForGraph][strForGraph][i][j]->Draw("ALP");
                gMWF[sysForGraph][strForGraph][i][j]->Draw("CP");
            }
            else {
                g[sysForGraph][strForGraph][i][j]->Draw("CP");
                gMWF[sysForGraph][strForGraph][i][j]->Draw("CP");
            }
            
        }
        divideIndex = 0;
    }
    l->AddEntry(gMWF[sysForGraph][strForGraph][0][0],"Using MaxWF","p");
    l->Draw();
    c->Update();
    
    char titleForFile[200];
    sprintf(titleForFile,"STC_S%d_EnergyResolutions.pdf",strForGraph);
    ////TPad *myPad = (TPad*) c->GetPad(1);
    c->Print(titleForFile);
    
    
    /*
    c->cd();
    g[0][0][0][0] = new TGraph(ntaus,x[0][0][0][0],y[0][0][0][0]);
    g[0][0][0][0]->SetTitle("FWHM vs Tau, STC S0D1 Pk0");
    g[0][0][0][0]->GetYaxis()->SetTitle("FWHM (keV)");
    g[0][0][0][0]->GetXaxis()->SetTitle("Tau (Index 50=0,70=1,70MWF=2,90=3)");
    g[0][0][0][0]->GetYaxis()->SetRangeUser(1,4);
    g[0][0][0][0]->GetXaxis()->SetRangeUser(-1,4);
    g[0][0][0][0]->SetMarkerStyle(21);
    g[0][0][0][0]->SetMarkerColor(38);
    g[0][0][0][0]->SetLineColor(38);
    g[0][0][0][0]->Draw("ALP");
    
    g[0][0][0][1] = new TGraph(ntaus,x[0][0][0][1],y[0][0][0][1]);
    g[0][0][0][1]->SetMarkerStyle(21);
    g[0][0][0][1]->SetMarkerColor(4);
    g[0][0][0][1]->SetLineColor(4);
    g[0][0][0][1]->Draw("CP");
    
    g[0][0][0][2] = new TGraph(ntaus,x[0][0][0][2],y[0][0][0][2]);
    g[0][0][0][2]->SetMarkerStyle(21);
    g[0][0][0][2]->SetMarkerColor(46);
    g[0][0][0][2]->SetLineColor(46);
    g[0][0][0][2]->Draw("CP");
    c->Update();
    */
     
    /////////////////
    // CLOSING OUT //
    /////////////////
    rootfile->Close();
    
	cout << "D O N E" << endl;
	App->Run(); 
}