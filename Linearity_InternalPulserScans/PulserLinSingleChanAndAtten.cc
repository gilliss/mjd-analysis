/***************************************************
 Use this scrip with ./PulserLinSingleChanAndAtten <M> <startRun> <endRun> <chan#> <0-3 for trapECal/trapENFCal> <0/1/2 for unatt/intatt/finatt> <doSave 0/1>
 
 Example Usage: ./PulserLinSingleChanAndAtten 1 18435 18498 632 1 1 0
 
 Note: This seems to require the 201610 software environment
****************************************************/

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include "TF1.h"
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TApplication.h>

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>
#include <vector>

#include <GATDataSet.hh>
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
    if (argc < 8)
    {
        cout << "too few arguments" << argv[0] << endl;
        return 1;
    }
    int M = atoi(argv[1]);
    int startRun = atoi(argv[2]);
    int endRun = atoi(argv[3]);
    int chan = atoi(argv[4]);
    int e_i = atoi(argv[5]);
    int attenuated = atoi(argv[6]);
    int doSave = atoi(argv[7]);
    char eEstimator[10];
        if (e_i == 0) {sprintf(eEstimator,"trapECal");}
        if (e_i == 1) {sprintf(eEstimator,"trapENFCal");}
        if (e_i == 2) {sprintf(eEstimator,"trapEN2FCal");}
        if (e_i == 3) {sprintf(eEstimator,"trapE");}
        if (e_i == 4) {sprintf(eEstimator,"trapENF");}
        if (e_i == 5) {sprintf(eEstimator,"trapEN2F");}
    char attChar[10];
        if (attenuated == 0) {sprintf(attChar,"Unatten");}
        if (attenuated == 1) {sprintf(attChar,"InterAtten");}
        if (attenuated == 2) {sprintf(attChar,"FinalAtten");}
    
    TApplication *App;
    if(doSave==0) {App = new TApplication("App", 0, NULL);}

    /* DATASET, TCHAIN, MJTCHANNELMAP */
    GATDataSet ds;
    for (int i_run = startRun; i_run <= endRun; i_run++)
    {
        if(i_run>startRun+2)
        {
            if(M==1){ds.AddRunNumber(i_run);}
            if(M==2){ds.AddRunNumber(i_run);}
        }
    }
    TChain* ch = ds.GetGatifiedChain();

    TFile* inFile;
    if(M==1)
    {
        inFile = new TFile(Form("/global/project/projectdirs/majorana/data/mjd/surfmjd/data/gatified/P3KJR/mjd_run%d.root",startRun)); // DS3
    }
    if(M==2)
    {
        inFile = new TFile(Form("/global/project/projectdirs/majorana/data/mjd/surfmjd/data/gatified/P3LQG/mjd_run%d.root",startRun)); // DS4
    }
    MJTChannelMap *chmap = (MJTChannelMap*)inFile->Get("ChannelMap");
    cout << "- - - - " << chmap->GetDetectorName(chan) << " " << chan << " - - - - " << endl;
    cout << "   startRun " << startRun << " endRun " << endRun << endl;
    
    /* Out TFILE */
    char resultsPath[200], outFile[100], outFilePath[300];
    sprintf (resultsPath,"/global/u2/g/gilliss/PulserLinearity/Update_July2017_DS34/results/");
    sprintf (outFile,"PulserLin_M%d_%s_%s_10-11-2016.root",M,attChar,eEstimator);
    sprintf (outFilePath,"%s%s",resultsPath,outFile);
    TFile *f;
    if(doSave==1) {f = new TFile(outFilePath,"UPDATE"); f->cd();} // open file; create & open if doesn't already exist
    
    /* RESIDUAL FILE */
    char residFile[100], residFilePath[300];
    sprintf (residFile,"Residuals_%d_%s_%s_10-11-2016.txt",chan,attChar,eEstimator);
    sprintf (residFilePath,"%s%s",resultsPath,residFile);
    ofstream myfile;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////
    /////////// Plots
    ///////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /* SETTINGS */
    double runBinLo = startRun-0.5;
    double runBinHi = endRun+0.5;
    int nRunBins = runBinHi - runBinLo;
    int runFitLo = startRun+(endRun-startRun)*0.35; // unatten'd vals, overwritten for atten'd
    int runFitHi = startRun+(endRun-startRun)*0.65; // unatten'd vals, overwritten for atten'd
    cout << "   runFitLo " << runFitLo << " runFitHi " << runFitHi << " (initial)" << endl;
    double maxEDisplay = 5000.0;
    int nEBins = maxEDisplay/100.0;
    
    /* ENERGY v RUN HIST */
    TCanvas *cEvR;
    if(doSave==0) {
        cEvR = new TCanvas("cEvR","cEvR");
        cEvR->cd();
        //ch->Draw(Form("%s:run >> hE2D(%d,%f,%f,%d,0.0,%f)",eEstimator,nRunBins,runBinLo,runBinHi,nEBins,maxEDisplay),Form("%s>25.0 && LocalEntry$>1000 && channel==%d",eEstimator,chan),"");
        //cout << "EvR plot: "<<Form("%s>25.0 && LocalEntry$>1000 && channel==%d",eEstimator,chan)<<endl;
        ch->Draw(Form("%s:run >> hE2D(%d,%f,%f,%d,0.0,%f)",eEstimator,nRunBins,runBinLo,runBinHi,nEBins,maxEDisplay),Form("channel==%d",chan),"");
        cout << "EvR plot: "<<Form("channel==%d",chan)<<endl;
    }
    
    /* TPROFILE */
    TCanvas *cProf = new TCanvas("cProf","cProf");
    cProf->cd();
    TProfile *profx = new TProfile(Form("EvR_%d_%s_%s",chan,attChar,eEstimator),Form("EvR_%d_%s_%s",chan,attChar,eEstimator),nRunBins,runBinLo,runBinHi,0.0,maxEDisplay); //Form("EvR_%d_%s_%s",chan,attChar,eEstimator)
    profx->GetYaxis()->SetTitle(Form("Processed pulser energy (%s)",eEstimator));
    profx->GetXaxis()->SetTitle("Run number");
    profx->SetTitle(Form("TProfile EvsRun, %s, %s %d",attChar,chmap->GetDetectorName(chan).c_str(),chan));
    profx->SetMarkerStyle(kFullDotMedium);
    gStyle->SetErrorX(0.); // suppress x-errors
    profx->SetErrorOption("i"); // y-errors = s.e.m. of y-vals from which avg was computed: sqrt(avg(y^2)-avg(y)^2)/sqrt(N)
    profx->SetStats(0);
    ch->Draw(Form("%s:run >> EvR_%d_%s_%s",eEstimator,chan,attChar,eEstimator),Form("%s>55.0 && channel==%d",eEstimator,chan),"PROF");
    cout << "Prof plot: "<<Form("%s>55.0 && channel==%d",eEstimator,chan)<<endl;
    if(doSave==1) {profx->Write();}

    /* ATTENUATED: ENSURE FIT DOES NOT INCLUDE 100-200KEV REGION WHERE KINK IS SUSPECTED */
    if(attenuated==1 || attenuated == 2)
    {
        int bin_i = 0;
        int binFound = 0;
        while (bin_i<=profx->GetNbinsX() && binFound==0)
        {
            //cout << "bin content " << profx->GetBinContent(bin_i) << endl;
            if(chan==658)
            {
                if(profx->GetBinContent(bin_i)>200.0 && bin_i!=9)
                {
                    cout << "   Found bin " << bin_i << " w/ profx->GetBinContent(bin_i) = " << profx->GetBinContent(bin_i) << " >200" << endl;
                    runFitLo = profx->GetBinCenter(bin_i);
                    runFitHi = profx->GetBinCenter((profx->GetNbinsX()) - 13);
                    cout << "   runFitLo " << runFitLo << " runFitHi " << runFitHi << " (adjusted rel to 200keV)" << endl;
                    binFound++;
                }
            }
            else if(profx->GetBinContent(bin_i)>200.0) // && bin_i!=9
            {
                cout << "   Found bin " << bin_i << " w/ profx->GetBinContent(bin_i) = " << profx->GetBinContent(bin_i) << " >200" << endl;
                runFitLo = profx->GetBinCenter(bin_i);
                runFitHi = profx->GetBinCenter((profx->GetNbinsX()) - 13);
                cout << "   runFitLo " << runFitLo << " runFitHi " << runFitHi << " (adjusted rel to 200keV)" << endl;
                binFound++;
            }
            if(runFitLo >= runFitHi)
            {
                cout << "   do chan " << chan << " by hand: runFitLo >= runFitHi" << endl;
                return 0;
            }
            bin_i++;
        }
        if(binFound==0) cout << "   This chan is overattenuated (does not exceed 200keV)" << endl;
    }
    
    /* FIT TPROFILE */
    cout << "   runFitLo " << runFitLo << " runFitHi " << runFitHi << " (input for the fit)" << endl;
    TF1 *myfit = new TF1(Form("Fit_%d_%s_%s",chan,attChar,eEstimator),"pol1",runFitLo,runFitHi);
    profx->Fit(myfit,"I","",runFitLo,runFitHi); // polN : [0]+[1]*x+[2]*x**2+...+[N]*x**N // "W" Set all weights to 1 for non empty bins; ignore error bars // "I" Use integral of function in bin, normalized by the bin volume, instead of value at bin center // https://root.cern.ch/doc/master/classTH1.html#a63eb028df86bc86c8e20c989eb23fb2a
    cProf->Update();
    Double_t p0 = myfit->GetParameter(0);
    Double_t p1 = myfit->GetParameter(1);
    if(doSave==1) {myfit->Write();}
    
    /* PRINT FIT RESULTS */
    cout << "   runFitLo " << myfit->GetXmin() << " runFitHi " << myfit->GetXmax() << " (min max pulled from fit)" << endl;
    cout << "   fit from: " << myfit->Eval(runFitLo) << " to " << myfit->Eval(runFitHi) << " keV" << endl;
    cout << "   myfit: " << p0 << " + " << p1 << "x" << endl;
    
    /* PLOT RESIDUALS AND WRITE RESID FILE */
    TCanvas *cGraph = new TCanvas("cGraph","cGraph");
    cGraph->cd();
    const int n = profx->GetNbinsX();
    vector<double> x; vector<double> ex; vector<double> residual; vector<double> eresidual;
    int i_pt = 0;
    if(doSave==1) {myfile.open(residFilePath);}
    for (int i=1;i<=n;i++)
    {
        if(profx->GetBinEntries(i)>0)
        {
            x.push_back(myfit->Eval(profx->GetBinCenter(i)));
            ex.push_back(0.0);
            residual.push_back(profx->GetBinContent(i) - myfit->Eval(profx->GetBinCenter(i)));
            eresidual.push_back(profx->GetBinError(i));
            if(residual[i-1]<-5 || residual[i-1]>5)
            {
                cout << " - - - - " << endl;
                cout << "i_pt="<<i_pt<<" should = vector.size() "<<x.size()<<endl;
                cout<<"bin center = "<< x.at(i_pt) << " ?= " << myfit->Eval(profx->GetBinCenter(i)) << endl; // int(x[i-1] + startRun)
                cout << "bin entries = " << profx->GetBinEntries(i) << endl;
                cout << "bin value = " << profx->GetBinContent(i) << endl;
                cout << "residual = " << residual.at(i_pt) << " ?= " << profx->GetBinContent(i) - myfit->Eval(profx->GetBinCenter(i)) << endl;
            }
            if(doSave==1) {myfile << x[i_pt] << " " << residual[i_pt] << "\n";}
            i_pt++;
        }
    }
    if(doSave==1) {myfile.close();}
    TGraphErrors *gr = new TGraphErrors(x.size(),&(x[0]),&(residual[0]),&(ex[0]),&(eresidual[0])); // https://root-forum.cern.ch/t/vector-and-tgraph/10781/8
    gr->GetYaxis()->SetTitle(Form("Deviation of energy from linear fit (%s)",eEstimator));
    gr->GetXaxis()->SetTitle(Form("Expected energy along linear fit (%s)",eEstimator));
    gr->SetMarkerStyle(kFullDotMedium);
    gStyle->SetErrorX(0.); // suppress x-errors
    gr->Draw("AP");
    gr->SetName(Form("LinRes_%d_%s_%s",chan,attChar,eEstimator));
    gr->SetTitle(Form("Pulser Linearity, %s, %s %d",attChar,chmap->GetDetectorName(chan).c_str(),chan));
    if(doSave==1) {gr->Write();}
    //gr->GetYaxis()->SetRangeUser(-2.0,1.0);
    cGraph->Update();
    //cGraph->Print(Form("results/Lin_%s_%s_%s_%d.pdf",attChar,eEstimator,chmap->GetDetectorName(chan).c_str(),chan));

    /* CLOSE OUT */
    if(doSave==1) {f->Close();}
    cout << "F I N I S H E D" << endl;
    if(doSave==0) {App->Run();}
}








/* read back in and recreate a TF1
root [0] TFile *f = new TFile("PulserLin_648_Unatten_trapENFCal_10-11-2016.root")
(TFile *) 0x21a34f0
root [1] Fit_648_Unatten_trapENFCal->GetParameter(0)
(Double_t) -529026.
root [2] Fit_648_Unatten_trapENFCal->GetParameter(1)
(Double_t) 28.7968
root [3] Fit_648_Unatten_trapENFCal->GetXmin()
(Double_t) 18393.0
root [4] Fit_648_Unatten_trapENFCal->GetXmax()
(Double_t) 18411.0
root [10] TF1 *f1 = new TF1("f1",Form("%f*x+%f",slope,offset),xmin,xmax)
(TF1 *) 0x3d459f0
root [11] f1->Draw("SAME")*/

/*
 const int n = profx->GetNbinsX();
 double x[n], ex[n], residual[n], eresidual[n];
 if(doSave==1) {myfile.open(residFilePath);}
 for (int i=1;i<=n;i++)
 {
 x[i-1] = myfit->Eval(profx->GetBinCenter(i));
 ex[i-1] = 0.0;
 residual[i-1] = profx->GetBinContent(i) - myfit->Eval(profx->GetBinCenter(i));
 eresidual[i-1] = profx->GetBinError(i);
 if(residual[i-1]<-5 || eresidual[i-1]>5)
 {
 cout << " - - - - " << endl;
 cout<<"bin center = "<< x[i-1] << " ~perhaps~ run " << int(x[i-1] + startRun) << endl;
 cout << "bin entries = " << profx->GetBinEntries(i) << endl;
 cout << "bin value = " << profx->GetBinContent(i) << endl;
 cout << "residual = " << residual[i-1] << endl;
 }
 if(doSave==1) {myfile << x[i-1] << " " << residual[i-1] << "\n";}
 }
 if(doSave==1) {myfile.close();}
 TGraphErrors *gr = new TGraphErrors(n,x,residual,ex,eresidual);
*/

/* ENERGYCAL v ENERGY HIST */
/*
TCanvas *cEvE;
if(doSave==0) {
    cEvE = new TCanvas("cEvE","cEvE");
    cEvE->cd();
    ch->Draw(Form("trapENFCal:trapENF >> hEvE2D(%d,%f,%f,%d,0.0,%f)",eEstimator,nRunBins,runBinLo,runBinHi,nEBins,maxEDisplay),Form("%s>25.0 && LocalEntry$>1000 && channel==%d",eEstimator,chan),"");
    cout << "EvE plot: "<<Form("%s>25.0 && LocalEntry$>1000 && channel==%d",eEstimator,chan)<<endl;
}
*/