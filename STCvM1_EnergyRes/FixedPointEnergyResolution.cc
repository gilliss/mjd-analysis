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

using namespace std;

int main (int argc, char* argv[])
{
    cout << "Have " << argc << " arguments:" << endl;
    for (int i = 0; i < argc; ++i) cout << argv[i] << endl;

    if (argc < 11) 
    {
        cout << "too few arguments" << argv[0] << endl;
        return 1;
    }
    
    
    int thechannel = atoi(argv[1]);    // arguments pulled from command issued by python script
    double baseline = atof(argv[2]);
    double tau = atof(argv[3]);
    char path[100];
    sprintf (path,"/global/project/projectdirs/majorana/data/mjd/%s/data/built/%s",argv[4],argv[5]);
    int startrun = atoi(argv[6]);
    int endrun = atoi(argv[7]);
    bool datasaved = atoi(argv[8]);
    bool firstrun = atoi(argv[9]);
    bool lookatwf = atoi(argv[10]);
    double risetime = atof(argv[11]);
    double flattop = atof(argv[12]);
    double ADC = atof(argv[13]);
    bool fixtrap = atoi(argv[14]);
    
    /*
    cout << "These are your settings: " << endl;
    cout << "   system = " << argv[4] << ", partnumber " << argv[5] << endl;
    cout << "   runs = " << startrun << " - " << endrun << endl; 
    cout << "   channel = " << thechannel << endl;
    cout << "   baseline = " << baseline << ", tau = " << tau << "ns" << endl;
    cout << "   risetime = " << risetime << "us, flattop " << flattop << "us" << endl;
    cout << "   peak pos of 2.6 MeV line estimated to be at " << ADC << endl;
    cout << "   using fixed pos of trap filter: " << fixtrap << endl;
    cout << "   first time running = " << firstrun << ", data already saved " << datasaved << ", look at wf " << lookatwf << endl;
    */
    
    
    TApplication *App = new TApplication("App", 0, NULL);    
    char hname[200],rootname[200],infilename[200],infile[200];   
    
    int nchannels=1; // these lines allow for some added automation and added peaks
    int npeaks=3;
    string npeakNames[npeaks];
    TH1D* histCal[nchannels];        
    double peak_HG[npeaks];
    double peak_LG[npeaks];
    double peakpos[npeaks];
    double energy[npeaks];    
    double fitpos[npeaks];
    double fitpos_corr[npeaks];
    double ratio;
    
    string system = argv[4];
    cout << system << endl;
    if(system == "surfmjd")
    {
        double ADC2614keV = ADC; //6855; //user input
        double Energy2614keV = 2614.55;
        ratio = ADC2614keV/Energy2614keV;
        
        energy[0]=583.191;
        peak_HG[0]=ratio*energy[0];
        peak_LG[0]=1/3.*ratio*energy[0];
        
        energy[1]=Energy2614keV-2*511;
        peak_HG[1]=ratio*energy[1];
        peak_LG[1]=1/3.*ratio*energy[1];   
        
        energy[2]=Energy2614keV;
        peak_HG[2]=ADC2614keV; 
        peak_LG[2]=1/3.*ADC2614keV; 
    }
    
    else if(system == "surfst")
    {
        double ADC1460keV = ADC; //user input 
        double Energy1460keV = 1460.859; // 40K
        ratio = ADC1460keV/Energy1460keV;
        
        energy[0]=351.932; // 214Pb
        peak_HG[0]=ratio*energy[0];
        peak_LG[0]=1/3.*ratio*energy[0];
        npeakNames[0]="214Pb";
        
        energy[1]=Energy1460keV; // 40K
        peak_HG[1]=ADC1460keV; 
        peak_LG[1]=1/3.*ADC1460keV;  
        npeakNames[1]="40K";  
        
        energy[2]=2614.55; // 208Tl
        peak_HG[2]=ratio*energy[2]; 
        peak_LG[2]=1/3.*ratio*energy[2];  
        npeakNames[2]="208Tl";
    }
    
    else 
    {
        cout << "There is an issue" << endl;
        exit(-1);
    }
    
    
    /*
    cout << "estimated peak positions " << endl;
    for(int i = 0; i<npeaks; i++)
    {
        cout << "keV: " << energy[i] << " ADC HG: " << peak_HG[i] << endl;
    }
    */
       
     
    if(datasaved==false)
    {    
        TGraphErrors *graphCal[nchannels]; 
        TH1D *histEnergyArray[nchannels]; 
        TH1D *histEnergyArrayROI[nchannels][npeaks]; 
        TH1D* hbase[nchannels];
	TH1D* htau[nchannels];
        TH1D* hwave = new TH1D();
        MGWFTrapezoidalFilter* trap[nchannels];   //This is a trapezoidal filter for each channel
        double ns = 1.e-9;
        double us = 1.e-6;
        double ms = 1.e-3;
        int nbins = 10000;
        int minbin = 0;
        int maxbin_HG = 9000;
        int maxbin_LG = 3000;
        for(int j=0;j<nchannels;j++)  
        {   
            sprintf(hname,"hist_ch%d",j);
            if(j%2==0) histEnergyArray[j]= new TH1D(hname, hname, nbins, minbin, maxbin_HG);
            else histEnergyArray[j]= new TH1D(hname, hname, nbins, minbin, maxbin_LG);
            graphCal[j] = new TGraphErrors();
            sprintf(hname,"base of ch %d",j);
            hbase[j] = new TH1D(hname,"",1000,-100,200);
            sprintf(hname,"tau of ch %d",j);
            htau[j] = new TH1D(hname,"",1000,0,100000);
            for(int p=0;p<npeaks;p++)  
            {
                sprintf(hname,"hist_ch%d_peak%d",j,p);
                if(j%2==0) histEnergyArrayROI[j][p]= new TH1D(hname, hname, 5000,0.8*peak_HG[p],1.2*peak_HG[p]);
                else histEnergyArrayROI[j][p]= new TH1D(hname, hname, 5000,0.8*peak_LG[p],1.2*peak_LG[p]);
            }
            
            trap[j] = new MGWFTrapezoidalFilter();
            trap[j]->SetRampTime(risetime*us/ns); // converts from us to ns // rise time is duration of trap rise
            trap[j]->SetFlatTime(flattop*us/ns);
            trap[j]->SetRestingBaseline(baseline);
            trap[j]->SetDecayConstant(tau);
        }
        
        
        TChain *t1 = new TChain("MGTree");   
        for(int i=startrun;i<endrun+1;i++)
        {
            sprintf(infile,"OR_run%d",i);  
            sprintf(infilename,"%s/%s.root",path,infile);
            t1->AddFile(infilename);        
        }
        MGTWaveform* Wave = new MGTWaveform();
        MGTWaveform* TrapWave = new MGTWaveform();
        MGTEvent* event = new MGTEvent();
        MGWFLinearFit* linfitbase = new MGWFLinearFit();
        MGWFLinearFit* linfittop = new MGWFLinearFit();
        
        t1->SetBranchAddress("event",&event); // one of the built tree's branches holds the "event"
        int nentries=t1->GetEntries();
        if(nentries>5000000) nentries = 5000000;
        
        //nentries = 10;
        for(int i=0;i<nentries;i++) // loop through the entries of the TChain
        {
            t1->GetEntry(i);
            if(i%50000==0) cout<<"at "<<i<<" of "<<nentries <<endl;
	
            for(int i_digit=0; i_digit<event->GetNDigitizerData(); i_digit++) // loop thru the wfs of the current entry
            {
                ///cout << event->GetNDigitizerData() << endl; 	
                int channel = event->GetDigitizerData(i_digit)->GetID();
                cout <<"entry " << i << " wf " << i_digit << " channel " << channel << endl;
                int ch = 0;

                if(channel == thechannel)
                {
                    Wave=event->GetWaveform(i_digit);
                    if(firstrun==true)
                    {
                        
                        if(lookatwf==true)
                        {
                            TCanvas cwave;
                            cwave.cd();
                            hwave=Wave->GimmeHist(); 
                            hwave->SetLineColor(kBlack);
			    char waveTitle[200];
			    sprintf(waveTitle,"Chan %d Event %d WF %d",channel,i,i_digit);
			    hwave->SetTitle(waveTitle);  
                            hwave->DrawCopy();
                            cwave.Update();
                            
                        ////linfitbase->SetFitSamples(200,0); // user input // Tom 9/8/15 this block
                        ////linfittop->SetFitSamples(200,600); //user input
                        ////linfitbase->PerformFit(*Wave);
                        ////linfittop->PerformFit(*Wave);
                       	////double slopebase = linfitbase->GetSlopeTime();
                        ////double offsetbase = linfitbase->GetYIntercept();
                        ////double chibase = linfitbase->GetChi2();
                        ////double slopetop = linfittop->GetSlopeTime();
                        ////double offsettop = linfittop->GetYIntercept();
                        ////double chitop = linfittop->GetChi2();
                        ////cout<<"offsetbase chi "<<offsetbase<<" "<<chibase<<" offsettop chi"<<offsettop<<" "<<chitop<<endl;
			

			////int z = Wave->GimmeHist()->GetMaximumBin();
			////cout<<"WF: Max Bin# "<<z<<" has val "<<Wave->GimmeHist()->GetBinContent(z)<<" which should equal "<<Wave->GimmeHist()->GetMaximum()<<endl; // Tom 9/9/15

                            trap[ch]->TransformOutOfPlace(*Wave,*TrapWave);
                            TCanvas ctrap;
                            ctrap.cd();
                            hwave=TrapWave->GimmeHist(); 
                            hwave->SetLineColor(kBlack);  
                            sprintf(waveTitle,"Chan %d Event %d WF %d",channel,i,i_digit);
                            hwave->SetTitle(waveTitle); 
                            hwave->DrawCopy();
                            ctrap.Update();

                        ////int g = TrapWave->GimmeHist()->GetMaximumBin();
                        ////cout<<"Trap: Max Bin# "<<g<<" has val "<<TrapWave->GimmeHist()->GetBinContent(g)<<" which should equal "<<TrapWave->GimmeHist()->GetMaximum()<<endl; // Tom 9/9/15
       
                            App->Run();
                        }

                        
                        linfitbase->SetFitSamples(200,0); // user input
                        linfittop->SetFitSamples(200,600); //user input
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

		    }
                    
                    
                    
                    double energy_offline = 0;
                    trap[ch]->TransformOutOfPlace(*Wave,*TrapWave);
                    if(fixtrap==false) energy_offline = TrapWave->GimmeHist()->GetMaximum();
                    else 
                    {
                        energy_offline = TrapWave->GimmeHist()->GetBinContent(1400);//700); // ?
                        if(system == "surfst") energy_offline = TrapWave->GimmeHist()->GetBinContent(600); // 850 // ?
                    } 
                    //cout << energy_offline << endl;
                    histEnergyArray[ch]->Fill(energy_offline); // This is the uncalibrated spectrum
                    for (int j = 0; j<npeaks; j++)
                    {
                        if(ch%2==0)
                        {
                            if(energy_offline>0.8*peak_HG[j] && energy_offline<1.2*peak_HG[j])
                            {
                                histEnergyArrayROI[ch][j]->Fill(energy_offline);    // This is an ROI uncalibrated spectrum for initial gaussian fit           
                            }
                        }
                        else
                        {
                            if(energy_offline>0.8*peak_LG[j] && energy_offline<1.2*peak_LG[j])
                            {
                                histEnergyArrayROI[ch][j]->Fill(energy_offline);
                            }
                        }
                    }
                }
            }
        }
        
        if(firstrun==true)   // If firstrun == true, draw uncalibrated spectrum and print baseline and tau
        {
            int flattopint = int(flattop*1000);
            int risetimeint = int(risetime*1000);
            int baselineint = int(baseline);
            int tauint = int(tau);
	    sprintf(rootname,"histUncal%d-%d_%d_%d-%d_%d-%d_%d.root",startrun, endrun, thechannel, risetimeint, flattopint,baselineint,tauint,fixtrap); 
            TFile f(rootname, "new");
       /////////DRAW COMMANDS 
            TCanvas* c = new TCanvas;
            c->cd();
            for(int k = 0; k<nchannels; k++)
            {
                if(k%2==0)
                {
                    histEnergyArray[k] -> SetLineColor(1+k);
                    histEnergyArray[k] -> Draw("same");
		    histEnergyArray[k]->Write(); 
                }
            }
            f.Close();

            for(int i = 0; i<nchannels; i++)
            {
                double base =  hbase[i]->GetBinCenter(hbase[i]->GetMaximumBin());
                double tau =  htau[i]->GetMean();
                cout << base << " " << tau << endl;   
            } 
            App -> Run();
        }
        
        TCanvas* ccal = new TCanvas();
        ccal->cd();
        TF1 * fGausFit = new TF1("GausFit", "[0]*Gaus(x,[1],[2]) + [3] + [4]*x ",0,1000);
        for(int k=0; k<nchannels; k++)
        {
            for(int i=0;i<npeaks;i++)
            {
                peakpos[i] = histEnergyArrayROI[k][i]->GetBinCenter(histEnergyArrayROI[k][i]->GetMaximumBin()); // guess peak by max bin in small window
                fGausFit->SetRange(0.8*peakpos[i],1.2*peakpos[i]);
                fGausFit->SetLineColor(kBlue);
                fGausFit->SetLineWidth(1);
                fGausFit->SetLineStyle(2);
                fGausFit->SetParameter(0,histEnergyArrayROI[k][i]->GetBinContent(histEnergyArrayROI[k][i]->FindBin(peakpos[i])));
                fGausFit->SetParameter(1,peakpos[i]);
                fGausFit->SetParameter(2,peakpos[i]/500.);
                fGausFit->SetParameter(3,0);
                fGausFit->SetParameter(4,0);
                
                histEnergyArrayROI[k][i]->Fit("GausFit","qr+");
                              
                cout << "energy = " << energy[i] << " estimated pos " << peakpos[i] << " fit value " << fGausFit->GetParameter(1) << " +/- " << fGausFit->GetParError(1) << endl;
                
                graphCal[k]->SetPoint(i,fGausFit->GetParameter(1),energy[i]);
                graphCal[k]->SetPointError(i,fGausFit->GetParError(1),0);     // probably want to write these graphs to a root file
            }
        }
       /////////DRAW COMMANDS 
        TCanvas* canArray[npeaks];
        for(int i=0; i<npeaks; i++)
        {
            canArray[i] = new TCanvas();
            canArray[i]->cd();
            for(int k = 0; k<nchannels; k++)
            {
                if(k%2==0)
                {
                    histEnergyArrayROI[k][i] -> SetLineColor(1+k);
                    histEnergyArrayROI[k][i] -> Draw("same"); 
                }
            }    
        }
        
        //Linear Fit
        TF1 * fLinFit = new TF1("LinFit", "[0] + [1]*x",0,6000);
        double slope[nchannels];
        double offset[nchannels];
        for(int k=0; k<nchannels; k++)
        {
            fLinFit->SetParameter(0,1);
            fLinFit->SetParameter(1,1);
            fLinFit->SetLineColor(kRed);
            fLinFit->SetLineWidth(1);           
            if(k%2==0) fLinFit->SetRange(0,maxbin_HG);
            else fLinFit->SetRange(0,maxbin_LG);
            
            graphCal[k]->SetMarkerColor(1+k);
            graphCal[k]->SetMarkerSize(2);
            graphCal[k]->SetMarkerStyle(4);
            graphCal[k]->GetXaxis()->SetLimits(0,maxbin_HG);
            graphCal[k]->GetYaxis()->SetRangeUser(0,3000);
            if(k==0) graphCal[k]->Draw("AP");
            else graphCal[k]->Draw("P same");  // would be good to save this calibration map as a root file
            
            graphCal[k]->Fit("LinFit","qr+");
            slope[k] = fLinFit->GetParameter(1);
            offset[k] = fLinFit->GetParameter(0);
            cout << "offset = " << offset[k] << " slope = " << slope[k] << endl; 
        }
        
        
        //Apply calibration
        TCanvas* ccalhist = new TCanvas();
        ccalhist->cd();
        int flattopint = int(flattop*1000);
        int risetimeint = int(risetime*1000);
        int baselineint = int(baseline);
        int tauint = int(tau);
        
        sprintf(rootname,"histCal%d-%d_%d_%d-%d_%d-%d_%d.root",startrun, endrun, thechannel, risetimeint, flattopint,baselineint,tauint,fixtrap); 
        TFile f(rootname, "new");
        for(int k=0; k<nchannels; k++)
        {
            sprintf(hname,"CalChannel%d", k);
            if(k%2==0) histCal[k] = new TH1D(hname,hname,nbins,offset[k]+slope[k]*minbin,offset[k]+slope[k]*maxbin_HG);
            else histCal[k] = new TH1D(hname,hname,nbins,offset[k]+slope[k]*minbin,offset[k]+slope[k]*maxbin_LG);
						//^^note minbin=0 ADC, and setting bin ranges with y=mx+b y-values b/c it's 
						//EvsADC
						//^^these lines do the calibration simply by rebinning rather than using a "->fill"
						//expression like ADCval*slope+offset=energycal
            for(int i = 0; i<nbins; i++)
            {
                histCal[k]->SetBinContent(i,histEnergyArray[k]->GetBinContent(i));
						//^^histEArray & histCal both have 1000 bins so this'll work
						//Tom 9/10/15
            }
            histCal[k]->Write(); 
            //histCal[k]->SetLineColor(1+k);
            //histCal[k]->Draw("same");
        }
        f.Close();
        App->Run(); // uncommented this, Tom 9/10/15 
           
    } // end of if datasaved == false
    else
    {
        int flattopint = int(flattop*1000);
        int risetimeint = int(risetime*1000);
        int baselineint = int(baseline);
        int tauint = int(tau);
        
        sprintf(rootname,"histCal%d-%d_%d_%d-%d_%d-%d_%d.root",startrun, endrun, thechannel, risetimeint, flattopint,baselineint,tauint,fixtrap); 
        TFile* f = new TFile(rootname);   // setting the file up to be opened
        
        for(int k=0; k<nchannels; k++)
        {     
            sprintf(hname,"CalChannel%d", k);
            histCal[k] = (TH1D*)f->Get(hname); // pulling the histogram from the file
        }       
    }   
    

    //do fancy fit
    TF1 * fGausFitFancy = new TF1("GausFitFancy", "[0] + [1]*x + [2]*x*x + [3]*((1.-[6])/sqrt(2.*TMath::Pi())/[5]*exp(-0.5*pow((x-[4])/[5], 2)) + [6]*exp(0.5/[7]/[7] + (x-[4])/[7]/[5])/2./[7]/[5]*TMath::Erfc([7]/sqrt(2.)*(1/[7]/[7] + (x-[4])/[7]/[5])))",0,3000);
    fGausFitFancy->SetLineColor(kRed);
    fGausFitFancy->SetParName(0,"const"); //background
    fGausFitFancy->SetParName(1,"1st_coeff"); //background
    fGausFitFancy->SetParName(2,"2nd_coeff"); //background
    
    fGausFitFancy->SetParName(3,"Gaus_Amplitude");   // amplitude of the peak
    fGausFitFancy->SetParName(4,"Gaus_Mean"); // position of the peak
    fGausFitFancy->SetParName(5,"Gaus_sigma");  //width of the peak
    fGausFitFancy->SetParName(6,"frac_expon");   // height of tail
    fGausFitFancy->SetParName(7,"slope_beta");   // length of tail

    
    for(int k = 0; k<nchannels; k++)
    {
        cout << "channel " << k << endl;
        for(int i=0;i<npeaks;i++)
        {
            fGausFitFancy->SetRange(0.95*energy[i],1.15*energy[i]);
	    if(i==1) fGausFitFancy->SetRange(0.985*energy[i],1.01*energy[i]);
	    if(i==2) fGausFitFancy->SetRange(0.985*energy[i],1.01*energy[i]); //0.985 for lowerbound
            fGausFitFancy->SetParameter(0,0);
            fGausFitFancy->SetParameter(1,0);
            fGausFitFancy->SetParameter(2,0);
            fGausFitFancy->SetParameter(3,histCal[k]->GetBinContent(energy[i]));
            fGausFitFancy->SetParameter(4,energy[i]);
	    if(i==1) fGausFitFancy->SetParameter(4,energy[i]+1); //0.985 for lowerbound	    
            fGausFitFancy->SetParameter(5,energy[i]/1000.);
            //fGausFitFancy->SetParameter(6,0.5);
            //fGausFitFancy->SetParameter(7,1); 
            fGausFitFancy->SetParameter(6,0.);
            fGausFitFancy->SetParameter(7,1.5);     
            
            histCal[k]->Fit("GausFitFancy","qr+"); //perform fit    // might want to save blow ups of these fits in root files
 
            
            //evaluate fit just for curiosity
            /*
            double chisquare = fGausFitFancy->GetChisquare();
            cout << "chisquare = " << chisquare << endl;
            double ndf = fGausFitFancy->GetNDF();
            cout << "NDF = " << ndf << endl;
            double prob = TMath::Prob(chisquare,ndf);
            cout << "prob = " << prob << endl;
            */
            
            // calculate FWHM
            fGausFitFancy->SetParameter(0,0);
            fGausFitFancy->SetParameter(1,0);
            fGausFitFancy->SetParameter(2,0);
            double fxmax = fGausFitFancy->GetMaximumX(); 
            double fhalfmax = fGausFitFancy->Eval(fxmax)/2.;
            double fsigma = fGausFitFancy->GetParameter(5); 
            double fxlow = fGausFitFancy->GetX(fhalfmax, fxmax-5.*fsigma, fxmax); //find the x for which f(x) is half the max. Do this for a range left to the center of the peak
            double fxhi = fGausFitFancy->GetX(fhalfmax, fxmax, fxmax+5.*fsigma); //find the x for which f(x) is half the max. Do this for a range right to the center of the peak
            
            //if(i==(npeaks-1)) { } // "npeaks-1" gives the 208Tl peak
            
		cout<< npeakNames[i] <<" fit results:"<<endl;
                cout<<" "<<fGausFitFancy->GetParameter(4)<<"+/-"<<fGausFitFancy->GetParError(4)<<" "<<2.35482*fGausFitFancy->GetParameter(5)<<"+/-"<<fGausFitFancy->GetParError(5)<<endl;
                cout<<" "<<fGausFitFancy->GetParameter(4)<<"+/-"<<fGausFitFancy->GetParError(4)<<" "<< fxhi - fxlow <<"+/-~"<<fGausFitFancy->GetParError(5)<<endl;
            
            /*
            cout<<"For peak at "<<fGausFitFancy->GetParameter(4)<<" keV, FWHM (from sigma value of gaussian) = "<<2.35482*fGausFitFancy->GetParameter(5)<<" keV "<< " +/- " << fGausFitFancy->GetParError(5) << " keV" << endl;
            cout<<"For peak at "<<fGausFitFancy->GetParameter(4)<<" keV, FWHM (from xhi - xlow) = "<< fxhi - fxlow <<" keV " << " +/- " << fGausFitFancy->GetParError(5) << " keV" << endl;
            */
            
        }
    }
    
    TCanvas* c = new TCanvas;
    c->cd();

    int flattopint = int(flattop*1000);
    int risetimeint = int(risetime*1000);
    int baselineint = int(baseline);
    int tauint = int(tau);
    sprintf(rootname,"histCalFancyFit%d-%d_%d_%d-%d_%d-%d_%d.root",startrun, endrun, thechannel, risetimeint, flattopint,baselineint,tauint,fixtrap); 
    TFile f(rootname, "new"); /// setup to save the fancyfit histogram

    for(int k = 0; k<nchannels; k++)
    {
        if(k%2==0)
        {
            histCal[k] -> SetLineColor(1+k);
            histCal[k] -> Draw("same");         // might want to save blow ups of this graph's fitted peaks in root files
        }
	histCal[k]->Write(); 
    }
    
    f.Close();

    
    App->Run(); 
    
}

/*
    TCanvas* cSplit = new TCanvas; // This block is for showing each fancy-fitted peak on its own subpad
    cSplit->Divide(npeaks,1);
    for(int k = 0; k<nchannels; k++)
    {
        if(k%2==0)
        {
	    for(int i=0;i<npeaks;i++)
	    {
		cSplit->cd(i+1);
                //histCal[k] -> SetLineColor(1+k);
            	histCal[k]->GetXaxis()->SetRange(0.95*energy[i],1.15*energy[i]);
	    	if(i==2) histCal[k]->GetXaxis()->SetRange(0.985*energy[i],1.01*energy[i]);
                histCal[k] -> DrawClone();	
	    }
        }
    }

*/

