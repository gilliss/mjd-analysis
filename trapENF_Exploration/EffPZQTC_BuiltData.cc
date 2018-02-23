/*
Script to get a hold of one waveform and run the EffPZ QT correction on it. Visualize each step of the way for pedagogical purposes.
 
To find reasonable waveforms, follow procedure at http://mjwiki.npl.washington.edu/bin/view/Majorana/DaqShifting#A_Few_Commands_for_Checking_Data

Examples need to be physics events so that they experience charge trapping. Pulser events do not.
Ex: ./EffPZQTC_BuiltData 2162 11
 
Note: Can time align the trap WFs' t0s in the figure for visual clarity, but this obscures what the algorithm is actually doing. Their is a deltat that accounts for the NON-alignment of the trapWFs' t0s, but one cannot show the effect of this deltat if the trap WFs are drawn as aligned.
 
Un-comment lines like
 //TTF->SetDoNormalize();
 //TTF->SetDoNotResize(false,true);
 //TTF->SetFitAndSubtractFlatBaseline(false);
in order to un-time align the WFs
*/

#include <TApplication.h>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TLatex.h>

#include <MGTEvent.hh>
#include <MGWFBaselineRemover.hh>
#include <MGWFTrapezoidalFilter.hh>
#include <MGWFStartTimeFinder.hh>
#include <GATWFFixedTimePickoffProcessor.hh>
//#include <MGVWaveformParameterCalculator.hh>


using namespace std;

int main (int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);
    
///////////////////////////////////
// INITIALIZATIONS ////////////////
///////////////////////////////////
    if (argc < 2)
    {
        cout << "too few arguments" << argv[0] << endl;
        return 1;
    }
    int eventNum = atoi(argv[1]);
    int wfNum = atoi(argv[2]);
    
    // FILE & CHAIN //
    char infilepath[200];
    sprintf (infilepath,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3JDY/OR_run5573.root");
    TChain *ch = new TChain("MGTree");
    ch->AddFile(infilepath);
    
    // EVENT //
    MGTEvent* event = 0;
    ch->SetBranchAddress("event",&event);
    int channel = 0;
    
    // WAVEFORMS & TRAP SETTINGS //
    int nWFs = 3;
    MGTWaveform* rawWF = new MGTWaveform();
    MGTWaveform* triggerTrapWF = new MGTWaveform();
    MGTWaveform* optTrapWF = new MGTWaveform();
    double ttRT = 1.0, ttFT = 1.5; // trigger trap RiseTime and FlatTopTime
    double otRT = 4.0, otFT = 2.5; // optimal trap RiseTime and FlatTopTime
    double fromFE = 0.5; // distance to walk back from start of falling edge (for pickoff)
    double decayConst = 72.0;
    int blStatsLower = 50;
    int blStatsUpper = 700;
    
    
    // CANVAS //
    int nDiv = 2; // raw WF, trigger trap WF, optimal trap WF
    TCanvas* cDiv = new TCanvas(/*"c","c",400,800*/);
    cDiv->Divide(1,nDiv);
    
    // PLOTS //
    TH1D* hWF[nWFs];
    for (int i=0; i<nWFs; i++) hWF[i] = new TH1D();
    char title[50];
    double fontSize = .09;
    double titleOffsetY = 0.5, titleOffsetX = 1.0;
    double scale = 0;
    double xrangemin = 0, xrangemax = 20000;
    
    // UNITS //
    double us = 1.e-6;
    double ns = 1.e-9;
    
///////////////////////////////////
// LOOP TO GET WF /////////////////
///////////////////////////////////
    cout << "looping over events ..." << endl;
    
    //ch->LoadTree(eventNum);
    if(ch->GetEntry(eventNum)) // test if exists
    {
                channel = event->GetWaveform(wfNum)->GetID();
                cout<<"Found event "<<eventNum<<", WF "<<wfNum<<" (0-indexed), channel "<<channel<<endl;
                rawWF = event->GetWaveform(wfNum);
    }
    
///////////////////////////////////
// TRANSFORM AND PLOT /////////////
///////////////////////////////////
    
    // RAW WF //
    MGWFBaselineRemover* BLR = new MGWFBaselineRemover();
    BLR->SetStartSample(blStatsLower);
    BLR->SetBaselineSamples(blStatsUpper-blStatsLower+1);
    BLR->TransformInPlace(*rawWF);
    rawWF->LoadIntoHist(hWF[0]);
    //sprintf(title,"Original WF");// (Chan %d, Entry %d, WF %d)",channel,eventNum,wfNum);
    hWF[0]->SetTitle("");
    //hWF[0]->SetTitleSize(0.);
    hWF[0]->GetYaxis()->SetTitle("Pulse Height");
    hWF[0]->GetYaxis()->SetTitleSize(fontSize);
    hWF[0]->GetYaxis()->SetTitleOffset(titleOffsetY);
    hWF[0]->GetYaxis()->SetLabelSize(fontSize);
    hWF[0]->GetXaxis()->SetTitle("ns (10ns/sample)");
    hWF[0]->GetXaxis()->SetTitleSize(0.);
    hWF[0]->GetXaxis()->SetLabelSize(0.);
    hWF[0]->GetXaxis()->SetRangeUser(xrangemin,xrangemax);
    hWF[0]->SetLineColor(kBlack);
    hWF[0]->SetLineWidth(4);
    hWF[0]->SetStats(0);
    cDiv->cd(1);
    TPad* pad1 = (TPad*)cDiv->GetPad(1);
    pad1->SetTopMargin(0.07);
    pad1->SetBottomMargin(0.01);//(0.01);
    scale = hWF[0]->GetMaximum();
    hWF[0]->Scale(1/scale);
    hWF[0]->Draw("");
    
    // TRIGGER TRAP WF //
    MGWFTrapezoidalFilter* TTF = new MGWFTrapezoidalFilter();
    //TTF->SetDoNormalize();
    //TTF->SetDoNotResize(false,true);
    //TTF->SetFitAndSubtractFlatBaseline(false);
    TTF->SetRampTime(ttRT*us/ns); //1.*us/ns = 1000(= 1 us); in otherwords, the argument evaluates to units of ns
    TTF->SetFlatTime(ttFT*us/ns);
    TTF->SetDecayConstant(decayConst*us/ns);
    TTF->TransformOutOfPlace(*rawWF,*triggerTrapWF); // Transform rawWF into trapWF
    triggerTrapWF->LoadIntoHist(hWF[1]);
    //sprintf(title,"Trapezoidal Filtered WFs");
    hWF[1]->SetTitle("");
    //hWF[1]->SetTitleSize(0.);
    hWF[1]->GetYaxis()->SetTitle("Pulse Height");
    hWF[1]->GetYaxis()->SetTitleSize(fontSize);
    hWF[1]->GetYaxis()->SetTitleOffset(titleOffsetY);
    hWF[1]->GetYaxis()->SetLabelSize(fontSize);
    hWF[1]->GetXaxis()->SetTitle("ns (10ns/sample)");
    hWF[1]->GetXaxis()->SetTitleSize(fontSize);
    hWF[1]->GetXaxis()->SetTitleOffset(titleOffsetX);
    hWF[1]->GetXaxis()->SetLabelSize(fontSize);
    hWF[1]->GetXaxis()->SetRangeUser(xrangemin,xrangemax);
    hWF[1]->SetLineColor(kBlue);
    hWF[1]->SetLineWidth(4);
    hWF[1]->SetStats(0);
    cDiv->cd(2);
    TPad* pad2 = (TPad*)cDiv->GetPad(2);
    pad2->SetTopMargin(0.);//(0.07);
    pad2->SetBottomMargin(0.2);
    hWF[1]->Scale(1/scale);
    hWF[1]->Draw("");
    
    // START TIME FINDER //
    MGWFStartTimeFinder* STF = new MGWFStartTimeFinder();
    STF->SetupBLR(false);
    STF->GetExtremumFinder().SetLocalMinimumTime(0.*us/ns);
    STF->GetExtremumFinder().SetLocalMaximumTime(10.*us/ns);
    STF->RestrictToRegion(0, 1000);
    STF->SetThreshold(2.);
    STF->CalculateParameters(*triggerTrapWF);
    double t0 = STF->GetParameterValue(0); // in ns
    cout<<"t0 = "<<t0<<" ns"<<endl;
    TLine* l_t0 = new TLine(t0,hWF[1]->GetMinimum(),t0,hWF[1]->GetMaximum()); // (x1,y1,x2,y2)
    l_t0->SetLineStyle(7);
    l_t0->SetLineWidth(2);
    l_t0->SetLineColor(6);
    l_t0->Draw();
    
    // OPTIMAL TRAP WF
    MGWFTrapezoidalFilter* OTF = new MGWFTrapezoidalFilter();
    //OTF->SetDoNormalize();
    //OTF->SetDoNotResize(false,true);
    //OTF->SetFitAndSubtractFlatBaseline(false);
    OTF->SetRampTime(otRT*us/ns); // 1.*us/ns = 1 us ; in otherwords, the argument evaluates to units of ns
    OTF->SetFlatTime(otFT*us/ns);
    OTF->SetDecayConstant(decayConst*us/ns);
    OTF->TransformOutOfPlace(*rawWF,*optTrapWF); // Transform rawWF into trapWF
    optTrapWF->LoadIntoHist(hWF[2]);
    hWF[2]->GetXaxis()->SetRangeUser(xrangemin,xrangemax);
    hWF[2]->SetLineColor(kOrange-3);
    hWF[2]->SetLineWidth(4);
    hWF[2]->SetStats(0);
    cDiv->cd(2);
    hWF[2]->Scale(1/scale);
    hWF[2]->Draw("SAME");
    
    // PICKOFF TIME //
    // IF Trap WFs are not time aligned wrt their t0's
    //double deltat = (2*ttRT + ttFT) - (2*otRT + otFT); // evaluates to us
    //double offset = (deltat + otRT + (otFT-fromFE))*us/ns; // evaluates to ns
    //double pickoff = offset + t0;
    // IF Trap WFs are time aligned regarding their t0's
    double deltat = 0.; // evaluates to us
    double offset = (deltat + otRT + (otFT-fromFE))*us/ns; // evaluates to ns
    double pickoff = offset + t0;
    //
    cout<<"deltat = "<<deltat<<" us"<<endl;
    cout<<"offset = "<<(deltat + otRT + (otFT-fromFE))<<" us = "<<offset<<" ns"<<endl;
    cout<<"pickoff time = offset + t0 = "<<pickoff<<" ns"<<endl;
    TLine* l_pick = new TLine(pickoff,hWF[2]->GetMinimum(),pickoff,hWF[2]->GetMaximum()); // (x1,y1,x2,y2)
    l_pick->SetLineStyle(2);
    l_pick->SetLineWidth(2);
    l_pick->SetLineColor(2);
    l_pick->Draw();
    
    // RAW WF //
    cDiv->cd(1);
    hWF[0]->GetXaxis()->SetRangeUser(hWF[1]->GetXaxis()->GetXmin(),hWF[1]->GetXaxis()->GetXmax());
    
    // LEGEND //
    cDiv->cd(2);
    TLegend* legend = new TLegend(0.6,0.62,0.75,0.85);
    legend->AddEntry(hWF[0],Form("Original WF"),"l");
    legend->AddEntry(hWF[1],Form("%.1f-%.1f-%.1f #mu s",ttRT,ttFT,ttRT),"l");
    legend->AddEntry(hWF[2],Form("%.1f-%.1f-%.1f #mu s",otRT,otFT,otRT),"l");
    legend->AddEntry(l_t0,"t0","l");
    legend->AddEntry(l_pick,"pickoff","l");
    legend->SetTextSize(fontSize);
    legend->Draw();
    cDiv->Update();
    
    App->Run();
}

/*
 Regarding GAT/Apps/process_mjd_gat.cc:416
 
 `GATWFFixedTimePickoffProcessor trapENFDBProc(-7.0*us+4.0*us+2.0*us, trapStartTimeProc.GetNameOfPostedVector(), "dbTrap", "trapENF");`
 
 The first argument, `-7.0*us+4.0*us+2.0*us`, comes from `(t0_dbTrap - t0_triggerTrap) + (riseTime_dbTrap) + (flatTime_dbTrap - 0.5us)`
 
 And `t0_dbTrap - t0_triggerTrap` can perhaps be re-expressed as `(2*riseTime_triggerTrap + flatTime_triggerTrap) - (2*riseTime_dbTrap + flatTime_dbTrap)`.
*/