///////////////////////////////
//
// see ResultsFromText.cc for useful plotting reference
// see 2nuBB systematics plotting scripts for useful reference
// Quick & Dirty Method: http://mjwiki.npl.washington.edu/pub/Majorana/AnalysisReports/dirty_and_quick_upper_limit.pdf
//
///////////////////////////////

#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TMultiGraph.h>
#include <TList.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TVectorD.h>
#include <TLegend.h>
#include <TApplication.h>

#include <iostream> // i/o stream, i.e. cout
#include <fstream> // file i/o
#include <stdio.h> // i/o, i.e. sprintf
#include <stdlib.h>
#include <string>

#include <GATDataSet.hh>
#include <MJTChannelMap.hh>
#include "/global/project/projectdirs/majorana/software/sl64/mjsw/mjsw201706Prod/GAT/Apps/DataSetInfo.hh"

using namespace std;

int main (int argc, char* argv[])
{
    TApplication *App = new TApplication("App", 0, NULL);
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Command line arguments and initializations
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 2)
    {
        cout << "Error: Too few arguments" << argv[0] << endl;
        cout << "Usage is: ./___ DS C save" << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // DS#
    int c = atoi(argv[2]); // M#
    int save = atoi(argv[3]); // Save the graph
    if(DS!=0 && DS!=1 && DS!=2 && DS!=3 && DS!=4 && DS!=5 && DS!=6) {cout << "DS pending"<<endl; return 1;}
    if(c!=1 && c!=2) {cout<<"Error: Bad module number"<<endl; return 1;}
    
    // PRINT SETTINGS
    cout<<"------------"<<endl;
    cout<<"Settings:"<<endl;
    cout<<" Plotting DS" << DS << " for M"<<c<< endl;
    cout<<"------------"<<endl;
    
    // INPUT MJTCHANNELMAP FILE
    //TFile *mapFile = new TFile("/global/u2/g/gilliss/2nuBB_Systematics/data/DSChannelMaps.root","READ");
    
    // INPUT DATA FROM TEXT FILE
    char inTextFilePath[100];
    int ds, m, p, d;
    string isotope;
    double line;
    int roiCnt, bgCnt;
    double pkCnt, epkCnt;
    int detStatus;
    double dMkg, dLT;
    
    // INITIALIZATIONS
    int maxMod = 2, maxPos = 7, maxDet = 5;
    double norm_roiCnt = 0., enorm_roiCnt = 0., norm_pkCnt = 0., enorm_pkCnt = 0.;
    double norm_roiCnt_Int = 0., enorm_roiCnt_Int = 0., norm_pkCnt_Int = 0., enorm_pkCnt_Int = 0.;
    double scale = 0;
    int group_i = 0; // counter xAxis location (grouping; For instance a grouping could be a string or a line)
    int nDets = 0; // nDets in grouping
    //map <string, double> xCoordMap; // <label,xCoordForGraph>
    vector <string> xNameVec;
    string xName;
    
    // ISOTOPES TO USE
    map < string , vector<double> > isotopeMap;
    string tempi;
    vector<double> tempe;
    tempi = "228Ac"; tempe.push_back(338.320); tempe.push_back(911.204); tempe.push_back(968.971); tempe.push_back(1588.2);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "108mAg"; tempe.push_back(433.937); tempe.push_back(614.276); tempe.push_back(722.907);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "214Bi"; tempe.push_back(609.312); tempe.push_back(1120.3); tempe.push_back(1764.494); tempe.push_back(2204.1);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "60Co"; tempe.push_back(1173.24); tempe.push_back(1332.5);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "40K"; tempe.push_back(1460.83);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "234mPa"; tempe.push_back(1001.03);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "210Pb"; tempe.push_back(46.539);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "214Pb"; tempe.push_back(351.9321); // progenitor of 214Bi
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "226Ra"; tempe.push_back(186.211);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    tempi = "208Tl"; tempe.push_back(583.191); tempe.push_back(1593.0); tempe.push_back(2105.0); tempe.push_back(2614.53);
    isotopeMap.insert(pair< string , vector<double> >(tempi,tempe));
    tempe.clear();
    //cout<<"isotopeMap.size() = "<<isotopeMap.size()<<endl;
    //cout<<" should be 60Co 1173 1332 "<<isotopeMap.find("60Co")->second[0]<<" "<<isotopeMap.find("60Co")->second[1]<<endl;
    int truncatedLine;
    for(std::map< string , vector<double> >::iterator map_i=isotopeMap.begin(); map_i!=isotopeMap.end(); map_i++)
    {
        for(unsigned int lin_i = 0; lin_i<map_i->second.size(); lin_i++)
        {
            cout<<map_i->first<<" "<<map_i->second[lin_i]<<endl;
            truncatedLine = map_i->second[lin_i];
            xName = map_i->first + " " + to_string(truncatedLine);
            xNameVec.push_back(xName);
        }
    }
    cout<<"------------"<<endl;
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// Plotting and saving
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    // PLOTTING AND SAVING

    TCanvas *canv[2];
    TGraphErrors *graph[2];
    
    canv[0] = new TCanvas(Form("DS%d_C%d_RoiCnt",DS,c),Form("DS%d_C%d_RoiCnt",DS,c));
    canv[1] = new TCanvas(Form("DS%d_C%d_BgSubRoiCnt",DS,c),Form("DS%d_C%d_BgSubRoiCnt",DS,c));
    
    graph[0] = new TGraphErrors();
    graph[0]->SetTitle(Form("DS%d_C%d_RoiCnt",DS,c));
    graph[0]->SetMarkerColor(kBlack);
    graph[0]->SetMarkerStyle(20);
    graph[1] = new TGraphErrors();
    graph[1]->SetTitle(Form("DS%d_C%d_BgSubtractedRoiCnt",DS,c));
    graph[1]->SetMarkerColor(kBlack);
    graph[1]->SetMarkerStyle(20);
    
    group_i = 0;
    for(std::map< string , vector<double> >::iterator map_i=isotopeMap.begin(); map_i!=isotopeMap.end(); map_i++)
    {
        for(unsigned int lin_i = 0; lin_i<map_i->second.size(); lin_i++)
        {
            //cout<<"new line grouping"<<endl;
            sprintf(inTextFilePath,"results/DS%d_M%d_hDetLines.txt",DS,c);
            ifstream inTextFile(inTextFilePath); // scope of loop makes it OK to redefine inTextFile
            while(inTextFile >> m >> p >> d >> isotope >> line >> roiCnt >> bgCnt >> pkCnt >> epkCnt >> detStatus >> dMkg >> dLT)
            {
                //cout<<m<<p<<d<<" "<<isotope<<" "<<line<<" "<<roiCnt<<" "<<bgCnt<<" "<<pkCnt<<" "<<epkCnt<<" "<<detStatus<<" "<<dMkg<<" "<<dLT<<endl;
                if(detStatus==1) // if detector is an isGood detector
                {
                    if(isotope==map_i->first && line==map_i->second[lin_i] /* && !(DS==0 && m==1 && p==1 && d==2 && line==1764.49) && !(DS==0 && m==1 && p==4 && d==1 && line==1764.49) && !(DS==0 && m==1 && p==4 && d==5 && line==1764.49)*/)
                    {
                        //cout<<"   ... "<<m<<p<<d<<" "<<endl;
                        scale = 1/(dMkg*dLT);
                        
                        norm_roiCnt = roiCnt*scale;
                        norm_roiCnt_Int += norm_roiCnt;
                        enorm_roiCnt_Int += sqrt(roiCnt)*scale;
                        
                        norm_pkCnt = pkCnt*scale;
                        norm_pkCnt_Int += norm_pkCnt;
                        enorm_pkCnt_Int += epkCnt*scale; // sqrt() already taken in MultiLine.hh
                        //cout<<"   "<<norm_pkCnt<<" "<<" "<<norm_pkCnt_Int<<" "<<enorm_pkCnt_Int<<endl;
                        
                        nDets++;
                    }
                }
            } // end inTextFile
            if(nDets!=0){
                norm_roiCnt_Int=norm_roiCnt_Int/nDets;
                enorm_roiCnt_Int=enorm_roiCnt_Int/nDets;
                
                norm_pkCnt_Int=norm_pkCnt_Int/nDets;
                enorm_pkCnt_Int=enorm_pkCnt_Int/nDets;
                // Wenqin's Quick & Dirty UL special cases.
                // Can plot the results of these cases if you want
                if(norm_pkCnt_Int > 3*enorm_pkCnt_Int)
                {
                    cout<<"Peak found "<<map_i->first<<" "<<map_i->second[lin_i]<<" "<<norm_pkCnt_Int<<endl;
                }
                if(norm_pkCnt_Int > 0 && norm_pkCnt_Int < 3*enorm_pkCnt_Int)
                {
                    cout<<"90% UL "<<map_i->first<<" "<<map_i->second[lin_i]<<" "<<norm_pkCnt_Int+1.35*enorm_pkCnt_Int<<endl;
                }
                if(norm_pkCnt_Int <= 0)
                {
                    cout<<"90% UL "<<map_i->first<<" "<<map_i->second[lin_i]<<" "<<1.35*enorm_pkCnt_Int<<endl;
                }
            }
            else{ // nDets == 0
                norm_roiCnt_Int=0.;
                enorm_roiCnt_Int=0.;
                
                norm_pkCnt_Int=0.;
                enorm_pkCnt_Int=0.;
            }

            graph[0]->SetPoint(group_i,double(group_i+1),norm_roiCnt_Int);
            graph[0]->SetPointError(group_i,0.0,enorm_roiCnt_Int);
            graph[1]->SetPoint(group_i,double(group_i+1),norm_pkCnt_Int);
            graph[1]->SetPointError(group_i,0.0,enorm_pkCnt_Int);
            
            cout<<map_i->first<<" "<<map_i->second[lin_i]<<" roi "<<norm_roiCnt_Int<<" "<<enorm_roiCnt_Int<<" bgsubroi "<<norm_pkCnt_Int<<" "<<enorm_pkCnt_Int<<" nDets "<<nDets<<endl;
            
            norm_roiCnt_Int=0.;
            enorm_roiCnt_Int=0.;
            norm_pkCnt_Int=0.;
            enorm_pkCnt_Int=0.;
            nDets=0;
            group_i++;
        } // end lin_i; end single grouping
    } // end map_i

    
    // GRAPH
    double bottomMarg = 0.15;
    TAxis* a[2];
    int j = 1;
    int bin_j;

    canv[0]->cd();
    canv[0]->SetBottomMargin(bottomMarg);
    graph[0]->Draw("AP");
    graph[0]->GetYaxis()->SetTitle("Cnt/kg/dy");
    a[0] = graph[0]->GetXaxis();
    j = 1;
    while (j<((a[0]->GetXmax()))) // may need (a->GetXmax())-1
    {
        bin_j = a[0]->FindBin(j);
        a[0]->SetBinLabel(bin_j,xNameVec[j-1].c_str());
        j++;
    }
    canv[0]->Modified();
    canv[0]->Update();
    
    canv[1]->cd();
    canv[1]->SetBottomMargin(bottomMarg);
    graph[1]->Draw("AP");
    graph[1]->GetYaxis()->SetTitle("Cnt/kg/dy");
    a[1] = graph[1]->GetXaxis();
    j = 1;
    while (j<((a[1]->GetXmax()))) // may need (a->GetXmax())-1
    {
        bin_j = a[1]->FindBin(j);
        a[1]->SetBinLabel(bin_j,xNameVec[j-1].c_str());
        j++;
    }
    canv[1]->Modified();
    canv[1]->Update();
    
    if(save==1)
    {
        TFile * f = new TFile("DS_MultiLine_Graphs.root","UPDATE");
        f->cd();
        graph[0]->Write(Form("DS%d_C%d_RoiCnt",DS,c));
        graph[1]->Write(Form("DS%d_C%d_PkCnt",DS,c));
        f->Close();
    }

    cout << "F I N I S H E D" << endl;
    App->Run();
}
//for(int p = 1; p <= maxPos; p++)
//{
//    for(int d = 1; d <= maxDet; d++)
//    {
//        
//    } // d
//} // p
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////
/////////// ...
///////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////