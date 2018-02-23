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

#include <iostream>
#include <fstream>
#include <stdio.h>

#include <TApplication.h>
#include <MGWFTrapezoidalFilter.hh>
#include <MGTEvent.hh>
#include <MGTRun.hh>
#include <MGTWaveform.hh>
#include <MGWFLinearFit.hh>
#include <MGWFBaselineRemover.hh>
#include <MGWFSavitzkyGolaySmoother.hh>
#include <MGWFExtremumFinder.hh>
#include <MGWFTimePointCalculator.hh>

using namespace std;

double getMean(vector<double> wfVec)
{
    double mean = 0.0;
    double sum = 0.0;
    for(int i = 50; i<=300; i++) sum += wfVec.at(i);
    mean = sum/(300-50+1.0);
    return mean;
}
double getStdDev(vector<double> wfVec)
{
    double stdDev = 0.0;
    double variance = 0.0;
    double temp = 0.0;
    double mean = 0.0;
    double sum = 0.0;
    
    for(int i = 50; i<=300; i++)
    {
        sum += wfVec.at(i);
        //cout << i << " " << wfVec.at(i) << endl;
    }
    mean = sum/(300-50+1.0);
    for(int i = 50; i<=300; i++)
    {
        temp += (mean-wfVec.at(i))*(mean-wfVec.at(i));
    }
    
    variance = temp/(300-50+1-1.0);
    stdDev = sqrt(variance);
    //cout << stdDev << endl;
    return stdDev;
}

int main (int argc, char* argv[])
{
    if (argc < 11) 
    {
        cout << "too few arguments" << argv[0] << endl;
        return 1;
    }
    
    
    char path_ch[100];
    sprintf (path_ch,"./input/%s",argv[1]); // hardware map   
    bool usetauexternal = atoi(argv[2]); // BOOL to use external tau set below 
    double extau = atof(argv[3]); // use a separate tau defined in .py file
    char path[100];
    sprintf (path,"/global/project/projectdirs/majorana/data/mjd/%s/data/built/%s",argv[4],argv[5]);
    int startrun = atoi(argv[6]);
    int endrun = atoi(argv[7]);
    bool datasaved = atoi(argv[8]);
    bool firstrun = atoi(argv[9]);
    bool lookatwf = atoi(argv[10]);
    double risetime = atof(argv[11]);
    double flattop = atof(argv[12]);
    bool fixtrap = atoi(argv[13]);
    double pickuptime = atof(argv[14]);
    string ext = argv[15];
    double ADC = atof(argv[16]);
    int flatfitsample = atoi(argv[17]);
    int topfitsample = atoi(argv[18]);
    int nchannels = atoi(argv[19]);
    
    
    cout << "These are your settings: " << endl;
    cout << "   system = " << argv[4] << ", partnumber " << argv[5] << endl;
    cout << "   runs = " << startrun << " - " << endrun << endl; 
    cout << "   use external tau = " << usetauexternal << ", tau = " << extau << "ns" << endl;
    cout << "   risetime = " << risetime << "us, flattop " << flattop << "us" << endl;
    cout << "   pick up time = " << pickuptime << endl;
    cout << "   using fixed pos of trap filter: " << fixtrap << endl;
    cout << "   first time running = " << firstrun << ", data already saved " << datasaved << ", look at wf " << lookatwf << endl;
    cout << "   filenameextension = " << ext << endl;
    cout << "   ADC Guess = " << ADC << endl;
    cout << "   nchannels = " << nchannels << endl;
    cout << "   flatfitsample = " << flatfitsample << endl;
    cout << "   topfitsample = " << topfitsample << endl;
    
    
    
    
    
    
    TApplication *App = new TApplication("App", 0, NULL);    
    char hname[200],rootname[200],infilename[200],infile[200], hnameCalArray[200];   
    
    //int nchannels=4; // number of detectors... should match number of lines in HW map
    int npeaksCal = 2; // 3
    int npeaks= 2; // 3
    string npeakNames[npeaks];
    TH1D* histCal[nchannels];
    TH1D *histCalArray[nchannels][npeaksCal]; // histCalArray , this is used for fitting the fancy fit and being able to draw all fits at once.
    double peak_HG[npeaks];
    double peak_LG[npeaks];
    double peakpos[npeaks];
    double energy[npeaks];    
    double fitpos[npeaks];
    double fitpos_corr[npeaks];

    double ns = 1.e-9;
    double us = 1.e-6;
    double ms = 1.e-3;
    int nbins = 10000;
    int minbin = 0;
    int maxbin_HG = 10000;
    double slope[nchannels];
    double offset[nchannels];
    
    double peak_HG_Cal[npeaksCal];
    double peak_LG_Cal[npeaksCal];
    double energy_Cal[npeaks];
    
    double baseline[nchannels];
    double decaytime[nchannels];
    //read in hardware map:    The HW map gets used [here]...
    ifstream hwmap;
    hwmap.open(path_ch);
    cout << "opening file : " << path_ch << endl;
    std::map<int,int> chmap;
    std::string names[nchannels];
    int i = 0;
    int channelHG, channelLG;
    double bl, dc;
    string name, stringname;
    while(hwmap >> name >> stringname >> channelLG >> channelHG >> bl >> dc)
    {
        if(dc>0) // this only creates the chmap if the dc values being read in are > 0
        {
            names[i] = name;
            chmap[channelHG] = i+1;  // I think this makes a "pair" <channelHG (key, first),i+1 ()>. Maps all the channels into an indexed list
            baseline[i] = bl; // should this be the baseline calculated from data?
            decaytime[i] = dc; // should this be the decay time calculated from data?
            cout << name << " " << channelHG << " " << bl << " " << dc << endl;
            i++;
        }
    }
    
    string system = argv[4];
    cout << system << endl;
    if(system == "surfmjd")
    {
        //double ADC;
        //if(fixtrap==true) ADC = 6350;
        //else ADC = 6855;
        
        double ADC2614keV = ADC;
        double Energy2614keV = 2614.55;
        double ratio = ADC2614keV/Energy2614keV;
        
        //peaks for fitting
        energy[0]=238.5;
        peak_HG[0]=ratio*energy[0];
        peak_LG[0]=1/3.*ratio*energy[0];  
        npeakNames[0]="238.5 keV";
        
        energy[1]=Energy2614keV;
        peak_HG[1]=ADC2614keV; 
        peak_LG[1]=1/3.*ADC2614keV; 
        npeakNames[1]="208Tl";
        
        //peaks for calibration
        energy_Cal[0]=238.5;
        peak_HG_Cal[0]=ratio*energy_Cal[0];
        peak_LG_Cal[0]=1/3.*ratio*energy_Cal[0]; 
        
        energy_Cal[1]=Energy2614keV;
        peak_HG_Cal[1]=ADC2614keV;
        peak_LG_Cal[1]=1/3.*ADC2614keV;
        
    }
    else if(system == "surfst")
    {
        //double ADC;
        //if(fixtrap==true) ADC = 3600; //user input
        //else ADC = 3600; //user input

        double ADC1460keV = ADC; //user input 
        double Energy1460keV = 1460.859;
        double ratio = ADC1460keV/Energy1460keV;
        
        energy[0]=351.932; // This seems equivalent to energy_Cal
        peak_HG[0]=ratio*energy[0]; // This seems equivalent to peak_HG_Cal
        peak_LG[0]=1/3.*ratio*energy[0];   // don't think this gets used
        npeakNames[0]="214Pb";       
 
        energy[1]=Energy1460keV;
        peak_HG[1]=ADC1460keV;
        peak_LG[1]=1/3.*ADC1460keV;     // don't think this gets used
        
        energy[2]=2614.55;
        peak_HG[2]=ratio*energy[2];
        peak_LG[2]=1/3.*ratio*energy[2];   // don't think this gets used

        //peaks for calibration
        energy_Cal[0]=351.932;
        peak_HG_Cal[0]=ratio*energy_Cal[0];
        peak_LG_Cal[0]=1/3.*ratio*energy_Cal[0];  // don't think this gets used
        
        energy_Cal[1]=Energy1460keV;
        peak_HG_Cal[1]=ADC1460keV;
        peak_LG_Cal[1]=1/3.*ADC1460keV;  // don't think this gets used
        npeakNames[1]="40K";

        energy_Cal[2]=2614.55;
        peak_HG_Cal[2]=ratio*energy_Cal[2]; 
        peak_LG_Cal[2]=1/3.*ratio*energy_Cal[2];  // don't think this gets used
        npeakNames[2]="208Tl";
    }
    else 
    {
        cout << "There is an issue" << endl;
        exit(-1);
    }
    
    if(datasaved == false)
    {
        cout << "peak positions used for calibration" << endl;
        for(int i = 0; i<npeaksCal; i++)
        {
            cout << "   keV: " << energy_Cal[i] << " ADC HG: " << peak_HG_Cal[i] << endl;
        }
        cout << endl;  
    }
    else 
    {
        cout << "will try to fit these peaks  " << endl;
        for(int i = 0; i<npeaks; i++)
        {
            cout << "   keV: " << energy[i] << " ADC HG: " << peak_HG[i] << endl;
        }
        cout << endl;
    }

    if(datasaved==false)
    {    
        TGraphErrors *graphCal[nchannels]; // Calibration map , to be linear fit
        TH1D *histEnergyArray[nchannels]; // Uncal energy hist
        TH1D *histEnergyArrayROI[nchannels][npeaksCal]; // Uncal energy hist zoom-in , to be fit by simple gaus
        TH1D* hbase[nchannels]; // hist of baseline values
        TH1D* htau[nchannels]; // hist of tau values
        TH1D* hwave = new TH1D(); // hist for waveform
        MGWFTrapezoidalFilter* trap[nchannels];   //This is a trapezoidal filter for each channel

        for(int j=0;j<nchannels;j++) // Declare the above graph and hists 
        {   
            sprintf(hname,"hist_ch%d",j);
            histEnergyArray[j]= new TH1D(hname, hname, nbins, minbin, maxbin_HG);
            graphCal[j] = new TGraphErrors();
            sprintf(hname,"base of ch %d",j);
            hbase[j] = new TH1D(hname,"",1000,-100,200);
            sprintf(hname,"tau of ch %d",j);
            htau[j] = new TH1D(hname,"",1000,0,100000);
            for(int p=0;p<npeaksCal;p++)  
            {
                sprintf(hname,"hist_ch%d_peak%d",j,p);
                histEnergyArrayROI[j][p]= new TH1D(hname, hname, 1000,0.8*peak_HG_Cal[p],1.2*peak_HG_Cal[p]);
            }
            
            trap[j] = new MGWFTrapezoidalFilter(); // Declare trapezoidal filters specific to each channel
            trap[j]->SetRampTime(risetime*us/ns);
            //cout << "rt " << risetime*us/ns << endl;
            trap[j]->SetFlatTime(flattop*us/ns);
            //cout << "ft " << flattop*us/ns << endl;
            trap[j]->SetRestingBaseline(baseline[j]); // This implements the HW map baseline for all cases (external vs HW-map tau)
            if(usetauexternal==true) // Use a tau value specified in the .py script
            {
                trap[j]->SetDecayConstant(extau);
                cout << "extau " << extau << endl;
            }
            else // Use tau value specified in the harware map in the "input" directory
            {
                trap[j]->SetDecayConstant(decaytime[j]);
                cout << decaytime[j] << endl;
            }
        }
        
        
        TChain *t1 = new TChain("MGTree");   // Create chain from data specified runs in .py script
        for(int i=startrun;i<endrun+1;i++)
        {
            sprintf(infile,"OR_run%d",i); // for built data 
            sprintf(infilename,"%s/%s.root",path,infile);
            t1->AddFile(infilename);        
        }
        MGTWaveform* Wave = new MGTWaveform(); // Declare WF object for WF
        MGTWaveform* TrapWave = new MGTWaveform(); // Declare WF object for TrapWF
        MGTEvent* event = 0; //new MGTEvent();
        MGWFLinearFit* linfitbase = new MGWFLinearFit();
        MGWFLinearFit* linfittop = new MGWFLinearFit();
        
        t1->SetBranchAddress("event",&event);
        size_t nentries=t1->GetEntries();
        
        /*
        int chcounter = 1;
        TCanvas cwave;
        cwave.Divide(5,6,0.001,0.001,0);
        TCanvas ctrap;
        ctrap.Divide(5,6,0.001,0.001,0);
        */
        double extrmPt = 0;
        double extrmVl = 0;
        vector<double> smoothWFVector;
        double mean = 0;
        double stdDev = 0;
        double crossing = 0;
        double threshold = 0;
        cout << "Data is not saved... looping over events " << endl;
        for(size_t i=0;i<nentries;i++)
        {
            t1->LoadTree(i);
            t1->GetEntry(i);    
            if(i%50000==0) cout<<"at "<<i<<" of "<<nentries <<endl;

            size_t nwaveforms = event->GetNWaveforms();
            //cout << "WFs in event..: " << nwaveforms << endl;//event->GetNWaveforms()  //GetNDigitizerData() << endl;

            for( size_t i_digit=0; i_digit<event->GetNWaveforms(); i_digit++ )
            {
                //cout << "i "<< i <<" i_digit " << i_digit << endl;
                int channel = event->GetWaveform(i_digit)->GetID(); // GetDigitizerData
                
                //MGTWaveform* wf =  event->GetWaveform(i_digit);
                //cout << "length "<< wf->GetLength() << endl;
                
                int ch = chmap.find(channel)->second; // ch is the index of the channel listed in chmap
                //cout << "channel "<< channel << " ch " << ch << endl; // This shows that there are some instances of ch=0. Can use to confirm channel-ch mapping
                
                if(ch>0)
                {
                    //cout << "looking at detector " << names[ch-1] << " channel "<< channel << " ch " << ch << " baseline " << baseline[ch-1] << endl;
                    ch = ch - 1;  // zero-index the ch index of the channel to match how the index starts with 0 in trap[], hbase[], and htau[]
                
                    Wave = event->GetWaveform(i_digit);
                    if(firstrun==true)
                    {
                        if(lookatwf==true)
                        {
                                TCanvas cwave;
                                cwave.cd();
                                hwave=Wave->GimmeHist(); 
                                hwave->SetLineColor(kBlack);  
                                hwave->DrawCopy(); // the code is currently setup to draw the WF of the first event
                                cwave.Update();
                             
                                
                                trap[ch]->TransformOutOfPlace(*Wave,*TrapWave); // Transform Wave into TrapWave
                                TCanvas ctrap;
                                ctrap.cd();
                                hwave=TrapWave->GimmeHist(); // Get histogram from waveform object
                                cout << "Max Val " << TrapWave->GimmeHist()->GetMaximum() << endl;//->GetBinContent(pickuptime) << endl;
                                cout << "Fixed Val " << TrapWave->GimmeHist()->GetBinContent(pickuptime) << endl;
                                hwave->SetLineColor(kBlack);  
                                hwave->DrawCopy();
                                ctrap.Update();
                                
                                App->Run();
                        }
                        
                        linfitbase->SetFitSamples(200,flatfitsample); // user input
                        linfittop->SetFitSamples(200,topfitsample); // user input
                        linfitbase->PerformFit(*Wave);
                        linfittop->PerformFit(*Wave);
                        
                        double slopebase = linfitbase->GetSlopeTime();
                        double offsetbase = linfitbase->GetYIntercept();
                        double chibase = linfitbase->GetChi2();
                        
                        double slopetop = linfittop->GetSlopeTime();
                        double offsettop = linfittop->GetYIntercept();
                        double chitop = linfittop->GetChi2();
                        
                        hbase[ch]->Fill(offsetbase);
                        htau[ch]->Fill(-(offsettop-offsetbase)/slopetop);
                        //if(i%50000==0) cout << "	tau fill val: " << (-(offsettop-offsetbase)/slopetop) << endl; // test
                    } // END firstrun == true
                    else // these uncalibrated hists do not get filled unless firstrun==false
                    {
                        /////////////////////////////////
                        ///// Estimate t0 of raw WF /////
                        /////////////////////////////////
                        
                        MGWFBaselineRemover* BLR = new MGWFBaselineRemover();
                        BLR->SetStartSample(50);
                        BLR->SetBaselineSamples(250);
                        BLR->TransformInPlace(*Wave);
                        
                        MGWFSavitzkyGolaySmoother* SGS = new MGWFSavitzkyGolaySmoother(24,0,5,"MGWFSavitzkyGolaySmoother"); // (smoothSize,derivativeOrder,polynomialDegree,name)
                        MGTWaveform* smoothWF = new MGTWaveform();
                        SGS->TransformOutOfPlace(*Wave,*smoothWF);
                        MGWFExtremumFinder* EXF = new MGWFExtremumFinder();
                        EXF->SetFindMaximum(true);
                        EXF->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
                        EXF->Transform(smoothWF);
                        extrmPt = EXF->GetTheExtremumPoint();
                        extrmVl = EXF->GetTheExtremumValue();
                        smoothWFVector.resize(0);
                        smoothWFVector = smoothWF->GetVectorData(); // MGWaveform.hh
                        mean = getMean(smoothWFVector);
                        stdDev = getStdDev(smoothWFVector);
                        threshold = mean+3*stdDev;
                        if (threshold<=0)
                        {
                            cout << "WARNING: threshold = " << threshold << " < 0 ... setting threshold to 0" << endl;
                            //cout << "This is b/c of a suspected MGWFTimePointCalculator issue, in which its methods do not like finding 0%-of-max timepoints" << endl;
                            threshold = 1e-20;
                        }
                        MGWFTimePointCalculator* TPC = new MGWFTimePointCalculator();
                        TPC->SetLocalMaximumTime(20000);//(argument in ns) sets upper limit on times to consider
                        TPC->Transform(smoothWF);
                        TPC->AddPoint(threshold/extrmVl);
                        TPC->FindTimePoints(*smoothWF);
                        crossing = TPC->GetFromMaxRiseTime(0);
                        
                        /////////////////////////////////
                        /////////////////////////////////
                        /////////////////////////////////
                        
                        double energy_offline = 0;
                        trap[ch]->TransformOutOfPlace(*smoothWF,*TrapWave);
                        
                        if(fixtrap==false) energy_offline = TrapWave->GimmeHist()->GetMaximum(); // use trapWF maximum
                        else 
                        {
                            energy_offline = TrapWave->GimmeHist()->GetBinContent(pickuptime+crossing); // use trapWF fixed point
                        } 

                        histEnergyArray[ch]->Fill(energy_offline);
                        for (int j = 0; j<npeaksCal; j++)
                        {
                            if(energy_offline>0.8*peak_HG_Cal[j] && energy_offline<1.2*peak_HG_Cal[j])
                            {
                                histEnergyArrayROI[ch][j]->Fill(energy_offline);   
                            }
                        }
                        
                        delete BLR;
                        delete SGS;
                        delete EXF;
                        delete TPC;
                        delete smoothWF;
                    } 
                }
            }
        } // End looping over data
        
        
        if(firstrun==true) // Get baseline and tau results and print them
        {
            //TCanvas* ctest[nchannels]; // test
            for(int i = 0; i<nchannels; i++)
            {
                double base =  hbase[i]->GetBinCenter(hbase[i]->GetMaximumBin());
                double tau =  htau[i]->GetMean();
                cout << names[i] << " " << base << " " << tau << endl;
                
                //ctest[i] = new TCanvas(); // test
                //ctest[i]->cd(); // test
                //htau[i]->Draw(); // test
            } 
            App -> Run(); // End firstrun==true
        }
    
        cout << "uncal'd plots enabled" << endl;
        TCanvas* ctest[nchannels]; // test  for (firstrun=0,lookatwf=0,datasaved=0)
        char uncaldPlotTitle[200];
        for(int i = 14; i<15; i++)//for(int i = 0; i<nchannels; i++)
        {
            ctest[i] = new TCanvas(); // test
            ctest[i]->cd(); // test
            sprintf(uncaldPlotTitle,"%s ch=%d",names[i].c_str(),i);
            histEnergyArray[i]->SetTitle(uncaldPlotTitle);
            histEnergyArray[i]->Draw(); // test
        } 
        App -> Run(); // test
	
        
        TCanvas* cdummy = new TCanvas(); // is this canvas necessary for performing the fit? It doesn't seem to be used. I think calling Fit may automaitcally draw to this canvas too
        cdummy->cd(); // is this canvas necessary for performing the fit? It doesn't seem to be used
        TF1 * fGausFit = new TF1("GausFit", "[0]*Gaus(x,[1],[2]) + [3] + [4]*x ",0,1000);
        for(int k=0; k<nchannels; k++)
        {
            cout << "channel " << k << endl;
            for(int i=0;i<npeaksCal;i++)
            {
                peakpos[i] = histEnergyArrayROI[k][i]->GetBinCenter(histEnergyArrayROI[k][i]->GetMaximumBin()); // gives a guess at ADC value of peak
                cout<< "k " << k << " i " << i << endl;//" peakpos[i] " << peakpos[i] << endl;

                fGausFit->SetRange(0.8*peakpos[i],1.2*peakpos[i]);
                fGausFit->SetLineColor(kBlue);
                fGausFit->SetLineWidth(1);
                fGausFit->SetLineStyle(2);
                fGausFit->SetParameter(0,histEnergyArrayROI[k][i]->GetBinContent(histEnergyArrayROI[k][i]->FindBin(peakpos[i]))); // gaus amp
                fGausFit->SetParameter(1,peakpos[i]); // gaus mean
                fGausFit->SetParameter(2,peakpos[i]/800.); // gaus sigma
                fGausFit->SetParameter(3,0);
                fGausFit->SetParameter(4,0);
                
                histEnergyArrayROI[k][i]->Fit("GausFit","qr+"); // Simple gaus fit to the uncal energy hist zoom-in
                              
                cout << "   energy = " << energy_Cal[i] << " estimated pos " << peakpos[i] << " fit value " << fGausFit->GetParameter(1) << " +/- " << fGausFit->GetParError(1) << endl;
                
                graphCal[k]->SetPoint(i,fGausFit->GetParameter(1),energy_Cal[i]);
                graphCal[k]->SetPointError(i,fGausFit->GetParError(1),0);
            }
        }
        
        //Linear Fit
        TCanvas* ccal = new TCanvas();
        ccal->cd();
        TF1 * fLinFit = new TF1("LinFit", "[0] + [1]*x",0,6000);
        //double slope[nchannels]; // for histCalArray
        //double offset[nchannels]; // for histCalArray
        for(int k=0; k<nchannels; k++)
        {
            fLinFit->SetParameter(0,1);
            fLinFit->SetParameter(1,1);
            fLinFit->SetLineColor(kRed);
            fLinFit->SetLineWidth(1);           
            fLinFit->SetRange(0, maxbin_HG);
            
            graphCal[k]->SetMarkerColor(1+k);
            graphCal[k]->SetMarkerSize(2);
            graphCal[k]->SetMarkerStyle(4);
            graphCal[k]->GetXaxis()->SetLimits(0,maxbin_HG);
            graphCal[k]->GetYaxis()->SetRangeUser(0,3000);
            if(k==0) graphCal[k]->Draw("AP");
            else graphCal[k]->Draw("P same");
            
            graphCal[k]->Fit("LinFit","qr+"); // Fit the calibration map
            slope[k] = fLinFit->GetParameter(1);
            offset[k] = fLinFit->GetParameter(0);
            cout << "channel = " << k << " offset = " << offset[k] << " slope = " << slope[k] << endl; 
        }
        //Apply calibration
        //App->Run();
        
        int flattopint = int(flattop*1000);
        int risetimeint = int(risetime*1000);
        int pickint = int(pickuptime); 
        TFile *f[nchannels];
        for(int k=0; k<nchannels; k++) // Save a calibrated histogram
        {
            int baselineint = int(baseline[k]);
            int tauint = 0;
            if(usetauexternal==true)
            {
                tauint = extau;
            }
            else
            {
                tauint = int(decaytime[k]);
            }
            sprintf(rootname,"./data/histCal%d-%d_%d_%d-%d_%d-%d_%d_%d_%s.root",startrun, endrun, k, risetimeint, flattopint,baselineint,tauint,fixtrap, pickint, ext.c_str()); 
            sprintf(hname,"Channel%d_%s", k, names[k].c_str());
            
            f[k] = new TFile(rootname, "new");
            histCal[k] = new TH1D(hname,hname,nbins,offset[k]+slope[k]*minbin,offset[k]+slope[k]*maxbin_HG); // Declare calibrated histogram with calibrated bounds.
													     // Declaring calibrated bounds, in terms of energy, effectively
													     // calibrates this histogram.
            for(int i = 0; i<nbins; i++) // fill histEnergyArray info into histCal. nbins is the number of bins in histEnergyArray too
            {
                histCal[k]->SetBinContent(i,histEnergyArray[k]->GetBinContent(i));
            }
            histCal[k]->Write(); 
            f[k]->Close();
        }
        /*
        // Cal'd plots
        TCanvas* c = new TCanvas;
        c->Divide(1,4,0.001,0.001,0); // (ncol,nrow,marginx,marginy,color) 
        for(int k = 0; k<nchannels; k++)
        {
            	c->cd(k+1);
                c->SetLogy();
            	//histCal[k] -> GetXaxis() -> SetRangeUser(energy_Cal[p]-10, energy_Cal[p]+10);
            	histCal[k] -> SetLineColor(kBlack);
            	histCal[k] -> Draw("same");
            	//if(k==0) histCal[k] -> Draw(); 
            	//else histCal[k] -> Draw("same");
        }
        */
        App->Run();
           
    } // END datasaved == false
    else // Open the root file saved above
    {

        int flattopint = int(flattop*1000);
        int risetimeint = int(risetime*1000);
        int pickint = int(pickuptime);  
        TFile *f[nchannels];
        for(int k=0; k<nchannels; k++)
        {     
            int baselineint = int(baseline[k]);
            int tauint = 0;
            if(usetauexternal==true)
            {
                tauint = extau;
            }
            else
            {
                tauint = int(decaytime[k]);
            }
            sprintf(rootname,"./data/histCal%d-%d_%d_%d-%d_%d-%d_%d_%d_%s.root",startrun, endrun, k, risetimeint, flattopint,baselineint,tauint,fixtrap,pickint, ext.c_str()); 
            sprintf(hname,"Channel%d_%s", k, names[k].c_str());
            f[k] = new TFile(rootname); // Open the saved root file
            histCal[k] = (TH1D*)f[k]->Get(hname); // Open the saved histCal from within that root file

            for(int n = 0; n<npeaksCal; n++) // histCalArray
            {
                sprintf(hnameCalArray,"Channel%d_%d_%s", k, n, names[k].c_str());
                histCalArray[k][n] = (TH1D*)histCal[k]->Clone(hnameCalArray); 
            }
        }       
         
        // Do fancy fit
        TCanvas* can = new TCanvas();
        can->cd();
        
        TF1* fGausFitFancy[nchannels];
        char fname[100]; 

        for(int k = 0; k<nchannels; k++) // Declare a FancyFit for each channel
        {
            sprintf(fname,"GausFitFancy%d",k);
            fGausFitFancy[k] = new TF1(fname, "[0] + [1]*x + [2]*x*x + [3]*((1.-[6])/sqrt(2.*TMath::Pi())/[5]*exp(-0.5*pow((x-[4])/[5], 2)) + [6]*exp(0.5/[7]/[7] + (x-[4])/[7]/[5])/2./[7]/[5]*TMath::Erfc([7]/sqrt(2.)*(1/[7]/[7] + (x-[4])/[7]/[5])))",0,3000);
            
            fGausFitFancy[k]->SetLineColor(kRed);
            fGausFitFancy[k]->SetParName(0,"const"); //background
            fGausFitFancy[k]->SetParName(1,"1st_coeff"); //background
            fGausFitFancy[k]->SetParName(2,"2nd_coeff"); //background
            
            fGausFitFancy[k]->SetParName(3,"Gaus_Amplitude");   // amplitude of the peak
            fGausFitFancy[k]->SetParName(4,"Gaus_Mean"); // position of the peak
            fGausFitFancy[k]->SetParName(5,"Gaus_sigma");  //width of the peak
            fGausFitFancy[k]->SetParName(6,"frac_expon");   // height of tail
            fGausFitFancy[k]->SetParName(7,"slope_beta");   // length of tail
            
            
            //cout << "channel " << k << endl;
            
            
            for(int i=0;i<npeaks;i++) // Fit npeaks i for each channel k. // had i starting at 1
            {
                //if(k == 0 && i ==1) fGausFitFancy[k]->SetRange(0.95*energy[i],1.0013*energy[i]);
                //if(k == 0 && i ==2) fGausFitFancy[k]->SetRange(0.95*energy[i],1.0008*energy[i]);
                //else fGausFitFancy[k]->SetRange(0.95*energy[i],1.05*energy[i]);
                fGausFitFancy[k]->SetRange(0.95*energy[i],1.05*energy[i]);
                
                fGausFitFancy[k]->SetParameter(0,0); // const
                //fGausFitFancy[k]->SetParLimits(0,0,1000);
                
                fGausFitFancy[k]->SetParameter(1,0); // 1st_coeff
                //fGausFitFancy[k]->SetParLimits(1,-1,1);
                
                fGausFitFancy[k]->SetParameter(2,0); // 2nd_coeff
                //fGausFitFancy[k]->SetParLimits(2,-1,1);
                
                //if(k == 3 && i ==2) fGausFitFancy[k]->SetParameter(3,(histCal[k]->GetBinContent(energy[i]))+10);
                if(k == 1 && i ==2) fGausFitFancy[k]->SetParameter(3,(histCal[k]->GetBinContent(energy[i]))*1.16); // Gaus_Mean
                else fGausFitFancy[k]->SetParameter(3,histCal[k]->GetBinContent(energy[i])); // Gaus_Amplitude
                //fGausFitFancy[k]->SetParameter(3,histCal[k]->GetBinContent(energy[i])); // Gaus_Amplitude

                //if(k == 0 && i ==1) fGausFitFancy[k]->SetParameter(4,energy[i]-1.0); // Gaus_Mean
                //if(k == 0 && i ==2) fGausFitFancy[k]->SetParameter(4,energy[i]+4.0); // Gaus_Mean
                //if(k == 2 && i ==2) fGausFitFancy[k]->SetParameter(4,energy[i]+14.0); // Gaus_Mean
                if(k == 0 && i ==2) fGausFitFancy[k]->SetParameter(4,2616.0); // Gaus_Mean
                else if(k == 1 && i ==2) fGausFitFancy[k]->SetParameter(4,energy[i]+1.0); // Gaus_Mean
                else if(k == 2 && i ==2) fGausFitFancy[k]->SetParameter(4,2615.25); // Gaus_Mean
                else fGausFitFancy[k]->SetParameter(4,energy[i]); // Gaus_Mean
                //fGausFitFancy[k]->SetParameter(4,energy[i]); // Gaus_Mean
                
                fGausFitFancy[k]->SetParameter(5,energy[i]/800.); // Gaus_Sigma

                //if(k == 0 && i ==2) fGausFitFancy[k]->SetParameter(6,0.3); // frac_expon //
                //if(k == 1 && i ==2) fGausFitFancy[k]->SetParameter(6,0.5); // frac_expon // down frac exp up gaus amp up slope beta
                //if(k == 3 && i ==2) fGausFitFancy[k]->SetParameter(6,0.1); // frac_expon // down frac exp down up gaus amp down slope beta
                //else fGausFitFancy[k]->SetParameter(6,0.); // frac_expon 
                fGausFitFancy[k]->SetParameter(6,0.); // frac_expon

                //if(k == 0 && i ==2) fGausFitFancy[k]->SetParameter(7,1.0); // slope_beta
                //if(k == 1 && i ==2) fGausFitFancy[k]->SetParameter(7,2.1); // slope_beta
                //if(k == 3 && i ==2) fGausFitFancy[k]->SetParameter(7,0.5); // slope_beta // down frac exp down up gaus amp down slope beta
                //else fGausFitFancy[k]->SetParameter(7,1.); // slope_beta
                fGausFitFancy[k]->SetParameter(7,1.); // slope_beta
                
                histCalArray[k][i]->Fit(fname,"qr0"); //perform fit  // Fit Channel k // was histCal[k]
                histCalArray[k][i]->Fit(fname,"qr0"); //perform fit
                histCalArray[k][i]->Fit(fname,"qr0"); //perform fit
                histCalArray[k][i]->Fit(fname,"qr0"); //perform fit
                histCalArray[k][i]->Fit(fname,"qr0"); //perform fit
                histCalArray[k][i]->Fit(fname,"qr0"); //perform fit
                histCalArray[k][i]->Fit(fname,"qr0"); //perform fit
                histCalArray[k][i]->Fit(fname,"rq+"); //perform fit
                
                //evaluate fit just for curiosity
                /*
                double chisquare = fGausFitFancy->GetChisquare();
                cout << "chisquare = " << chisquare << endl;
                double ndf = fGausFitFancy->GetNDF();
                cout << "NDF = " << ndf << endl;
                double prob = TMath::Prob(chisquare,ndf);
                cout << "prob = " << prob << endl;
                */
                /*
                for(int p=0;p<8;p++)
                {
                    cout << p << " " << fGausFitFancy[k]->GetParameter(p) << endl;
                }
                */
                
                // calculate FWHM
                /*fGausFitFancy[k]->SetParameter(0,0); // had these uncommented, but seemed unecessary
                fGausFitFancy[k]->SetParameter(1,0);
                fGausFitFancy[k]->SetParameter(2,0);*/
                double fxmax = fGausFitFancy[k]->GetMaximumX(); 
                //cout << "max " << fxmax << endl;
                double fhalfmax = fGausFitFancy[k]->Eval(fxmax)/2.;
                //cout << "half the height " << fhalfmax << endl;
                double fsigma = fGausFitFancy[k]->GetParameter(5); 
                //cout << "sigma from fit " << fsigma << endl;
                double fxlow = fGausFitFancy[k]->GetX(fhalfmax, fxmax-5.*fabs(fsigma), fxmax); //find the x for which f(x) is half the max. Do this for a range left to the center of the peak
                //cout << "x at half max (left)" << fxlow << endl;
                double fxhi = fGausFitFancy[k]->GetX(fhalfmax, fxmax, fxmax+5.*fabs(fsigma)); //find the x for which f(x) is half the max. Do this for a range right to the center of the peak
                //cout << "x at half max (right) " << fxhi << endl;
                
                if(usetauexternal==true)
                {
                    cout << k << " " << names[k] << " " <<  npeakNames[i] << " tau=" << extau << " mean="<<fGausFitFancy[k]->GetParameter(4) << " " << fGausFitFancy[k]->GetParError(4) << " FWHM="<< fxhi - fxlow <<" " << fGausFitFancy[k]->GetParError(5) << endl;
                }
                else cout << k << " " << names[k] << " " << npeakNames[i] << " " << decaytime[k] << " " <<fGausFitFancy[k]->GetParameter(4) << " " << fGausFitFancy[k]->GetParError(4) << " "<< fxhi - fxlow <<" " << fGausFitFancy[k]->GetParError(5) << endl;
                
                //cout << endl;
                
            }
        }
        
        histCalArray[2][2] -> GetXaxis() -> SetRangeUser(energy_Cal[2]-1, energy_Cal[2]+1); // was histCal[k]
        histCalArray[2][2] -> GetXaxis() -> SetRangeUser(energy_Cal[2]-10, energy_Cal[2]+10); // was histCal[k]
        histCalArray[2][2] -> SetLineColor(kBlack);
        histCalArray[2][2]->Draw();
        
        TCanvas* c = new TCanvas;
        c->Divide(npeaksCal,nchannels,0.001,0.001,0);
        //c->cd();
        int divideIndex = 0;
        for(int k = 0; k<nchannels; k++)
        {
            for(int p=0; p<npeaksCal; p++)
            {
                divideIndex++;
                c->cd(divideIndex);
                //histCalArray[k][p] -> TitleSize(0.4);
                //if(k == 0 && p ==1) histCalArray[k][p] -> GetXaxis() -> SetRangeUser(1430.0, 1490.0);
                //else if(k == 0 && p ==2) histCalArray[k][p] -> GetXaxis() -> SetRangeUser(2540.0, 2660.0);
                histCalArray[k][p] -> GetXaxis() -> SetRangeUser(energy_Cal[p]-10, energy_Cal[p]+10);// was histCal[k]
                histCalArray[k][p] -> SetLineColor(kBlack);
                histCalArray[k][p] -> DrawCopy();
            }
        }

        /*
        TCanvas* c2 = new TCanvas;
        c2->Divide(6,4,0.001,0.001,0);
        
        for(int k = 0; k<nchannels; k++)
        {
            c2->cd(k+1);
            histCal[k] -> GetXaxis() -> SetRangeUser(238.5-5, 238.5+5);
            histCal[k] -> SetLineColor(kBlue);
            histCal[k] -> DrawCopy(); 
        }
        */
        
        
        App->Run(); 
    } // END datasaved == true
}

