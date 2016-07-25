#define basictracking_cxx

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TGraph.h"
#include "TF2.h"
#include "TArc.h"
#include "TPad.h"
#include "TMinuit.h"
#include "TPDF.h"
#include <string.h>
#include <cstring>
#include <stdlib.h>
#include "TSystem.h"
#include "basictracking.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <fstream>
#include <stdio.h>
#include <TAttFill.h>
#include "TNamed.h"




using namespace std;

TString stra = "0.55GeV_dabc16117205401";                   //file which will contain the plots
TString strb = ".root";                                     //file which will contain the plots
TString strc = "_corr";                                     //file which will contain the plots
TString strd = "_track";                                    //file which will contain the plots
TString filename1 = stra+strc+strb;                         //file which will contain the plots
TString filename2 = stra+strd+strb;                         //file which will contain the plots

// ***************** GLOBAL VARIABLES ********************************
Int_t eventhits3= 0;
Double_t newparam0,newparam1,dnewparam0,dnewparam1;
Double_t chisq = -9999.0,result = -9999.0;
Double_t Riso[150],sigma[150],xpos[150],ypos[150];
Double_t par[6];


// ***************** FUNCTIONS USED ********************************
void minuitfcn (Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag);
Double_t distancefrompoint2line (Double_t& x, Double_t& y,Double_t& slope,Double_t& yIntercept);
Double_t distancefrompoint2line2 (Double_t& x2, Double_t& y2,Double_t& slope2,Double_t& yIntercept2);
Double_t residual (Double_t dis,Double_t radius);
void trackMinimizer (Double_t& b,Double_t &a,Double_t &db,Double_t &da,Double_t chis);
Double_t* line(Double_t& x1, Double_t& x2,Double_t& y1,Double_t& y2);


void basictracking::track()				//START OF MAIN FUNCTION TRACK -->TRACK()!!!!!!!!!!!!!!!
{
    // ***************** VARIABLES, ARRAYS, ARCS, STRINGS ********************************
    Int_t neweventhits3,it,iiter,iter,outiter,truncate10,truncate20,truncate30;
    Double_t prob,chanDis,chanDis2,residual,sumresid,perpYIntercept,perpSlope,hitPosX,hitPosY,alpha,position_y;
    Int_t channel3[150],channelhit3[150],row3[150],column3[150];              // number of channels
    Double_t xposition3[150],yposition3[150],t03[150],t0channel3[150],tot3[150],totchannel3[150],dx[150],newxpos[150],newypos[150],newradius[150];
    Double_t prefitPointX[4],prefitPointY[4],slope[4],yIntercept[4],kSquared[4];
    Double_t polParam[5] = {0.};
    Double_t dslope = 0.001, dyIntercept =0.001;
    Double_t* aa;
    Double_t totoverdx[150][3];
    Double_t totsum,dxsum,temp,totsum10,dxsum10,totsum20,dxsum20,totsum30,dxsum30;
    Double_t array_test1[150];
    Double_t array_test2[150];
    Double_t array_test3[150];

    
    
    
    TH1D* risos = new TH1D("risos","risos",1000,-1.0,9.0); // in mm
    TH1I* multiplicity1 = new TH1I("multiplicity1","multiplicity1",50,0,50);
    TH1I* multiplicity2 = new TH1I("multiplicity2","multiplicity2",50,0,50);
    TH1I* multiplicity3 = new TH1I("multiplicity3","multiplicity3",50,0,50);
    TH1I* multiplicitycheck = new TH1I("multiplicitycheck","multiplicitycheck",50,0,50);
    TH1D* chisqplot = new TH1D("chisqplot","chisqplot",500,0.,50.);
    TH1D* probplot = new TH1D("probplot","probplot",500,0.,1.);
    TH1D* distancefromwire = new TH1D("distancefromwire","distancefromwire",1200,-6.,6.); ///in mm
    TH1D* distancefromisochrone = new TH1D("distancefromisochrone","distancefromisochrone",1000,-1.,1.); ///in mm
    TH1D* meanresiduals = new TH1D("meanresiduals","meanresiduals",1000,0.,2.); ///in mm
    TH2D* residual_vs_r = new TH2D("residual_vs_r","residual_vs_r",1000,-2.,2.,510,0.,5.1); ///in mm
    TH2D* t0_vs_iso = new TH2D("t0_vs_iso","t0_vs_iso",165,0.,165.,1100, -0.0055,0.0055);  ///in meters
    TH2D* t0_vs_iso2 = new TH2D("t0_vs_iso2","t0_vs_iso2",165,0.,165.,550, 0.,0.0055);      ///in meters
    TH1D* resid1 = new TH1D("resid1","resid1",1000,-0.001,0.001); // in meters
    TH1D* resid2 = new TH1D("resid2","resid2",1000,-0.001,0.001); // in meters
    TH1D* resid3 = new TH1D("resid3","resid3",1000,-0.001,0.001); // in meters
    TH1D* resid4 = new TH1D("resid4","resid4",1000,-0.001,0.001); // in meters
    TH1D* resid5 = new TH1D("resid5","resid5",1000,-0.001,0.001); // in meters
    TH1D* resid6 = new TH1D("resid6","resid6",1000,-0.001,0.001); // in meters
    TH1D* resid7 = new TH1D("resid7","resid7",1000,-0.001,0.001); // in meters
    TH1D* resid8 = new TH1D("resid8","resid8",1000,-0.001,0.001); // in meters
    TH1D* resid9 = new TH1D("resid9","resid9",1000,-0.001,0.001); // in meters
    TH1D* resid10 = new TH1D("resid10","resid10",1000,-0.001,0.001); // in meters
    TH1D* resid11 = new TH1D("resid11","resid11",1000,-0.001,0.001); // in meters
    TH1D* resid12 = new TH1D("resid12","resid12",1000,-0.001,0.001); // in meters
    TH1D* resid13 = new TH1D("resid13","resid13",1000,-0.001,0.001); // in meters
    TH1D* resid14 = new TH1D("resid14","resid14",1000,-0.001,0.001); // in meters
    TH1D* resid15 = new TH1D("resid15","resid15",1000,-0.001,0.001); // in meters
    TH1D* resid16 = new TH1D("resid16","resid16",1000,-0.001,0.001); // in meters
    TH1D* resid17 = new TH1D("resid17","resid17",1000,-0.001,0.001); // in meters
    TH1D* resid18 = new TH1D("resid18","resid18",1000,-0.001,0.001); // in meters
    TH1D* resid19 = new TH1D("resid19","resid19",1000,-0.001,0.001); // in meters
    TH1D* resid20 = new TH1D("resid20","resid20",1000,-0.001,0.001); // in meters
    TH1D* iso_error = new TH1D("iso_error","iso_error",20,0.,5.2);
//    TH1F* tot_over_dx = new TH1F("tot_over_dx","tot/dx",100,0.,100.);
//    TH1F* tot_over_dx_10 = new TH1F("tot_over_dx_10","tot/dx_10",100,0.,100.);
//    TH1F* tot_over_dx_20 = new TH1F("tot_over_dx_20","tot/dx_20",100,0.,100.);
//    TH1F* tot_over_dx_30 = new TH1F("tot_over_dx_30","tot/dx_30",100,0.,100.);
    TH1F* tot = new TH1F("tot","tot",1000,0.,1000.);
    TH1F* tot_10 = new TH1F("tot_10","tot_10",1000,0.,1000.);
    TH1F* tot_20 = new TH1F("tot_20","tot_20",1000,0.,1000.);
    TH1F* tot_30 = new TH1F("tot_30","tot_30",1000,0.,1000.);

    
    TH2D* wire_vs_channel = new TH2D("wire_vs_channel","wire_vs_channel",1400,-7.0,7.0,150,0.,150.);  ///in mm
    TH2D* resid_vs_channel = new TH2D("resid_vs_channel","resid_vs_channel",1000,-1.0,1.0,150,0.,150.);  //in mm
    TH1D* proj_y_rt;
    TH1D* new_rt = new TH1D("new_rt","new_rt",155,0.,155.); // in meters
    
    TH1D* resid_per_channel[150] = {0};         /////in mm
    TH1D* wire_per_channel[150] = {0};         /////in mm
    TH2D* rt_per_channel[150] = {0};         /////in mm
    TH2D* rt_per_channel2[150] = {0};         /////in mm
    TH1D* above_resid_per_channel[150] = {0};         /////in mm
    TH1D* below_resid_per_channel[150] = {0};         /////in mm
    TH1D* above_wire_per_channel[150] = {0};         /////in mm
    TH1D* below_wire_per_channel[150] = {0};         /////in mm
    TH1D* resid_per_layer[7] = {0};
    TH1D* wire_per_layer[7] = {0};
    TH2D* new_rt_layer[7] = {0};         /////in mm
    TH2D* new_rt_layer2[7] = {0};         /////in mm
    TH1D* above_resid_per_layer[7] = {0};
    TH1D* below_resid_per_layer[7] = {0};
    TH1D* above_wire_per_layer[7] = {0};
    TH1D* below_wire_per_layer[7] = {0};
    
    
    
    
    TF1* finalfit = new TF1("finalfit","pol4");
    TF1* prefit = new TF1("prefit","pol1");
    TF1* fit = new TF1("fit","pol1");
    prefit -> SetLineColor(kGreen);
    fit -> SetLineColor(kRed);
    TCanvas *c1 = new TCanvas("c1","c1");
    TH2D *plot = new TH2D("event","event",3000,0.0,0.30,1000,0.0,0.075);
    TArc *el1 = new TArc(0.,0.,0.);
    TArc *el2 = new TArc(0.,0.,0.);
    el1->SetLineColor(kBlack);
    el2->SetLineColor(kMagenta);
    el1->SetFillStyle(3000);
    el2->SetFillStyle(3000);
    TString a;
    TString Name1,Name2,Name3,Name4,Name5,Name6,Name7,Name8,Name1a,Name2a,Name3a,Name4a,Name5a,Name6a,Name7a,Name8a;
    gStyle->SetFuncWidth(1);
    TF1* g1 = new TF1("g1","gaus",-0.0055,0.);
    TF1* g2 = new TF1("g2","gaus",0.,0.0055);
    TF1* total = new TF1("total","gaus(0)+gaus(3)",-0.0055,0.0055);
    //gStyle->SetOptStat(0);
    
    for (Int_t i = 0; i< 150; i++) {
        Name1.Form(" Residual of channel %d " ,i);
        resid_per_channel[i] = new TH1D(Name1,Name1,1000,-1.0,1.0);
        Name2.Form(" Distance from wire of channel %d " ,i);
        wire_per_channel[i] = new TH1D(Name2,Name2,1200,-6.0,6.0);
        Name3.Form(" r(t) of channel %d " ,i);
        rt_per_channel[i] = new TH2D(Name3,Name3,155,0.,155.,1100, -5.5,5.5);
        Name4.Form(" abs r(t) of channel %d " ,i);
        rt_per_channel2[i] = new TH2D(Name4,Name4,155,0.,155.,550, 0.,5.5);
        Name5.Form(" Above wire residual of channel %d " ,i);
        above_resid_per_channel[i] = new TH1D(Name5,Name5,1000,-1.0,1.0);
        Name6.Form(" Below wire residual of channel %d " ,i);
        below_resid_per_channel[i] = new TH1D(Name6,Name6,1000,-1.0,1.0);
        Name7.Form(" Above wire distance from wire of channel %d " ,i);
        above_wire_per_channel[i] = new TH1D(Name7,Name7,1200,-6.0,6.0);
        Name8.Form(" Below wire distance from wire of channel %d " ,i);
        below_wire_per_channel[i] = new TH1D(Name8,Name8,1200,-6.0,6.0);
    }
    
    for (Int_t i = 0; i< 7; i++) {
        Name1a.Form(" Residual of layer %d " ,i);
        resid_per_layer[i] = new TH1D(Name1a,Name1a,1000,-1.0,1.0);
        Name2a.Form(" Distance from wire of layer %d " ,i);
        wire_per_layer[i] = new TH1D(Name2a,Name2a,1200,-6.0,6.0);
        Name3a.Form(" new rt curve for layer %d " ,i);
        new_rt_layer[i] = new TH2D(Name3a,Name3a,155,0.,155.,1100, -5.5,5.5);
        Name4a.Form(" new abs rt curve for layer %d " ,i);
        new_rt_layer2[i] = new TH2D(Name4a,Name4a,155,0.,155.,550, 0.,5.5);
        Name5a.Form(" Above wire residual of layer %d " ,i);
        above_resid_per_layer[i] = new TH1D(Name5a,Name5a,1000,-1.0,1.0);
        Name6a.Form(" Below wire residual of layer %d " ,i);
        below_resid_per_layer[i] = new TH1D(Name6a,Name6a,1000,-1.0,1.0);
        Name7a.Form(" Above wire distance from wire of layer %d " ,i);
        above_wire_per_layer[i] = new TH1D(Name7a,Name7a,1200,-6.0,6.0);
        Name8a.Form(" Below wire distance from wire of layer %d " ,i);
        below_wire_per_layer[i] = new TH1D(Name8a,Name8a,1200,-6.0,6.0);
        
    }
    
    
    ifstream polParamData;							     	//declare stream
    polParamData.open("out_Riso_0.55GeV_dabc16117205401.txt");				//open riso txt file
    for(Int_t i = 0; i < 5; i++){
        polParamData >> polParam[i];                       //save in array values from riso txt file
    }
    polParamData.close();									//close stream
    
    
    TFile* treefile = TFile::Open(filename1);		//open file to read data
    TTree* readtree = (TTree*)treefile -> Get("data2");	//get tree from file with the data
    TFile* myfile = new TFile(filename2,"RECREATE");				//create new file to save analysis
    TTree *tree = new TTree("data3","data3");
    tree -> Branch("eventhits3",&eventhits3);
    tree -> Branch("channel3",channel3,"channel3[150]/I");
    tree -> Branch("t03",t03,"t03[150]/D");
    tree -> Branch("tot3",tot3,"tot3[150]/D");
    tree -> Branch("xposition3",xposition3,"xposition3[150]/D");
    tree -> Branch("yposition3",yposition3,"yposition3[150]/D");
    tree -> Branch("row3",row3,"row3[150]/I");
    tree -> Branch("column3",column3,"column3[150]/I");
    tree -> Branch("channelhit3",channelhit3,"channelhit3[150]/I");
    tree -> Branch("t0channel3",t0channel3,"t0channel3[150]/D");
    tree -> Branch("totchannel3",totchannel3,"totchannel3[150]/D");
    tree -> Branch("xpos3",xpos,"xpos[150]/D");
    tree -> Branch("ypos3",ypos,"ypos[150]/D");
    
    Int_t nentries = readtree->GetEntriesFast();				//get entries from tree
    readtree -> SetBranchAddress("eventhits2", &eventhits3);
    readtree -> SetBranchAddress("channel2", &channel3);
    readtree -> SetBranchAddress("t02", &t03);
    readtree -> SetBranchAddress("tot2", &tot3);
    readtree -> SetBranchAddress("xposition2",&xposition3);
    readtree -> SetBranchAddress("yposition2",&yposition3);
    readtree -> SetBranchAddress("row2", &row3);
    readtree -> SetBranchAddress("column2", &column3);
    readtree -> SetBranchAddress("channelhit2", &channelhit3);
    readtree -> SetBranchAddress("t0channel2", &t0channel3);
    readtree -> SetBranchAddress("totchannel2", &totchannel3);
    
    
    cout << "  There are  " << nentries << " entries in this file....... " << endl;
    //c1->Print("event_display.pdf["); // ---------------> CREATION OF THE PDF FILE FOR THE VISUALIZATION OF THE TRACKS !!!!!!!!!!!!!!!!!!!
    for (Int_t lm = 0; lm < nentries; lm++) { 					////// START OF BIG FOR LOOP FOR EVENTS !!!!!!
        readtree -> GetEntry(lm);
        if (lm % 50000 == 0) cout << " ==> " << lm << " events processed........ " << endl;
        if (a != "a"){
            cout << "Do you want to see the next event ? (y , n , a)" << endl;
            cin >> a;
        }
        
        for (Int_t i = 0; i< 150; i++) {     ///// set arrays equal to zero
            Riso[i] = 0.0;
            sigma[i] = 0.0;
            xpos[i] = 0.0;
            ypos[i] = 0.0;
            newxpos[i] = 0.0;
            newypos[i] = 0.0;
            newradius[i] = 0.0;
            dx[i] = 0.0;
            array_test1[i] = 0.0;
            array_test2[i] = 0.0;
            array_test3[i] = 0.0;
        }
        for (Int_t i = 0; i< 4; i++) {            ///// set arrays equal to zero
            prefitPointX[i] = 0.0;
            prefitPointY[i] = 0.0;
            slope[i] = 0.0;
            kSquared[i] = 0.0;
            yIntercept[i] = 0.0;
        }
        
        neweventhits3 = eventhits3;             ///set eventhits3 to neweventhits3
        multiplicity1->Fill(neweventhits3);
        if(a =="y" || a == "a"){    //// start of the if a || y statement !!!!!!!!!!!!!!!
            plot->Draw();
            for(Int_t i = 0; i < eventhits3; i++){        ///start of hits loop per event
                //Riso[i] = (polParam[0]+(polParam[1]*t03[i])+(polParam[2]*pow(t03[i],2))+(polParam[3]*pow(t03[i],3))+(polParam[4]*pow(t03[i],4)));           ///////riso is in meters!!!!!!!!!
                Riso[i] = (-0.000105091+(5.42992e-05*t03[i])+(5.58724e-08*pow(t03[i],2))+(-2.60983e-09*pow(t03[i],3))+(9.13256e-12*pow(t03[i],4)));           ///////riso is in meters!!!!!!!!!
                if(Riso[i] <0.0 && Riso[i] > 0.00035) Riso[i] = 0.000025;
                sigma[i] = (0.000189633 + 2.42154e-05*Riso[i] - 5.96517e-06*Riso[i]*Riso[i]);
                risos->Fill(Riso[i]*1000);
                xpos[i] = xposition3[i];
                ypos[i] = yposition3[i];
                el1->DrawArc(xpos[i],ypos[i],0.005);            ////draw tubes in the map
            }       ///end hits loop per event
            
            it = 0;                                 /// set it to zero
            
            
            for (Int_t j = 0; j<25; j++) {        ///start j loop for columns
                for(Int_t i = 0; i < eventhits3; i++){    ///start hits loop per event
                    if(column3[i] == j && Riso[i] >=0.0 && Riso[i] <= 0.0052){  ///if statement for riso
                        newxpos[it] = xpos[i];
                        newypos[it] = ypos[i];
                        newradius[it] = Riso[i];
                        it++;
                    }       ///end if column and riso
                }			////end loop for hits
            }           ///end j for loop for columns
            
            /////////// ********** FOUR POINTS to draw the prefit ****************************
            prefitPointY[0] = newypos[0] - newradius[0];
            prefitPointY[1] = newypos[0] + newradius[0];
            prefitPointY[2] = newypos[it-1] - newradius[it-1];
            prefitPointY[3] = newypos[it-1] + newradius[it-1];
            prefitPointX[0] = newxpos[0];
            prefitPointX[1] = newxpos[0];
            prefitPointX[2] = newxpos[it-1];
            prefitPointX[3] = newxpos[it-1];
            
            
            for(Int_t i = 0; i < 2; i++){
                if(i == 0){			//first two lines-tracks
                    aa = line(prefitPointX[i],prefitPointX[i+2],prefitPointY[i],prefitPointY[i+2]);
                    kSquared[0] = aa[2];
                    slope[0] = aa[1];
                    yIntercept[0] = aa[0];
                    aa = line(prefitPointX[i],prefitPointX[i+3],prefitPointY[i],prefitPointY[i+3]);
                    kSquared[1] = aa[2];
                    slope[1] = aa[1];
                    yIntercept[1] = aa[0];
                    if (kSquared[1] - kSquared[0] < 0){			//comparing first two lines-tracks and choosing the best according to q^2
                        kSquared[0] = kSquared[1];
                        slope[0] = slope[1];
                        yIntercept[0] = yIntercept[1];
                    }       ////end comparing chi2
                }           ///end if i==1
                if(i == 1){				//third and fourth line-track
                    aa = line(prefitPointX[i],prefitPointX[i+1],prefitPointY[i],prefitPointY[i+1]);
                    kSquared[2] = aa[2];
                    slope[2] = aa[1];
                    yIntercept[2] = aa[0];
                    aa = line(prefitPointX[i],prefitPointX[i+2],prefitPointY[i],prefitPointY[i+2]);
                    kSquared[3] = aa[2];
                    slope[3] = aa[1];
                    yIntercept[3] = aa[0];
                    if (kSquared[3] - kSquared[2] < 0){			//comparing third and fourth line-track and choosing the best according to q^2
                        kSquared[2] = kSquared[3];
                        slope[2] = slope[3];
                        yIntercept[2] = yIntercept[3];
                    }       ////end comparing chi2
                }       ///end if i==1
            }     ///end i for loop
            
            if (kSquared[2] - kSquared[0] < 0){				//choosing the best out of the 4 lines-tracks based on q^2
                kSquared[0] = kSquared[2];
                slope[0] = slope[2];
                yIntercept[0] = yIntercept[2];
            }
            
            
            prefit -> SetParameters(yIntercept[0],slope[0]);       ///set parameters for prefit
            prefit->Draw("same");                   ////draw prefit
            
            outiter = 0;
            for(Int_t i = 0; i < eventhits3; i++){
                if(Riso[i] >=0.0 && Riso[i] <= 0.0052){
                    chanDis = 0.0;
                    chanDis = distancefrompoint2line(xpos[i],ypos[i],slope[0],yIntercept[0]);//distance of track from wire
                    if(TMath::Abs(chanDis) > 0.006 && outiter < 4 && neweventhits3>7){         //if track is more than 7mm far from wire!!!!!!!!!!!!!!!!!!!
                        outiter++;
                        Riso[i] =999.9;
                        neweventhits3--;
                    }       ///end if chanDis
                }       // end if Riso
            }       ///end hits for loop
            
            multiplicity2->Fill(neweventhits3);
            
            
            trackMinimizer (yIntercept[0],slope[0],dyIntercept,dslope,kSquared[0]);		//track construction using minuit function etc
            fit -> SetParameters(newparam0,newparam1);                      ///set parameters for final fit
            fit->Draw("same");              ///draw final track/fit
            
            outiter = 0;
            for(Int_t i = 0; i < eventhits3; i++){
                if(Riso[i] >=0.0 && Riso[i] <= 0.0052){
                    chanDis = 0.0;
                    residual = 0.0;
                    chanDis = distancefrompoint2line(xpos[i],ypos[i],slope[0],yIntercept[0]);//distance of track from wire
                    residual = (chanDis-Riso[i]);
                    if(TMath::Abs(residual) > 0.0006 && outiter < 4 && neweventhits3>7){         //if track is more than 700μm far from isochrone!!!!!!!!!!!!!!!!!!!
                        outiter++;
                        Riso[i] =999.9;
                        neweventhits3--;
                    }       ///end if chanDis
                }
            }
            multiplicity3->Fill(neweventhits3);
            
            if(outiter !=0){
                trackMinimizer (yIntercept[0],slope[0],dyIntercept,dslope,kSquared[0]);		//track construction using minuit function etc
                fit -> SetParameters(newparam0,newparam1);                      ///set parameters for final fit
                fit->Draw("same");              ///draw final track/fit
            }
            
            chisqplot->Fill(result);///fill chi2 plot of final fit
            result = result*(neweventhits3-2);
            prob = 0.0;
            prob = TMath::Prob(result,neweventhits3-2);                         ///calculation of prob value from the chi2
            probplot -> Fill(prob);                                             ////fill prob plot of final fit
            
            sumresid = 0.0;                        /// set variables equal to zero
            iter = 0;
            totsum = 0;
            dxsum = 0;
            totsum10 = 0;
            dxsum10 = 0;
            totsum20 = 0;
            dxsum20 = 0;
            totsum30 = 0;
            dxsum30 = 0;
            temp = 0;
            truncate10 = 0;
            truncate20 = 0;
            truncate30 = 0;
            iiter = 0;
            if(neweventhits3>7){
                // if(prob >0.1){
                for(Int_t i = 0; i < eventhits3; i++){
                    if(Riso[i] >=0.0 && Riso[i] <= 0.0052){
                        el2->DrawArc(xpos[i],ypos[i],Riso[i]);    /////draw ischrones in the map with an if statement
                        iiter++;
                        perpSlope = 0.0;
                        perpYIntercept = 0.0;
                        hitPosX = 0.0;
                        hitPosY = 0.0;
                        alpha = 0.0;
                        chanDis = 0.0;
                        chanDis2 = 0.0;
                        residual = 0.0;
                        position_y = 0.0;
                        chanDis = distancefrompoint2line(xpos[i],ypos[i],newparam1,newparam0);//distance of track from wire
                        chanDis2 = distancefrompoint2line2(xpos[i],ypos[i],newparam1,newparam0);//distance of track from wire
                        residual = (chanDis-Riso[i]);
                        sumresid += TMath::Abs(residual);
                        distancefromwire->Fill(distancefrompoint2line2(xpos[i],ypos[i],newparam1,newparam0));
                        distancefromisochrone->Fill(residual*1000);
                        residual_vs_r->Fill(residual*1000,Riso[i]*1000);
                        t0_vs_iso -> Fill(t03[i],distancefrompoint2line2(xpos[i],ypos[i],newparam1,newparam0));
                        t0_vs_iso2 -> Fill(t03[i],chanDis);
                        position_y = (newparam1*xpos[i] + newparam0);
                        for (Int_t j =0; j<7; j++) {
                            if(row3[i] == j){
                                new_rt_layer[j]->Fill(t03[i],chanDis2*1000);
                                new_rt_layer2[j]->Fill(t03[i],chanDis*1000);
                                wire_per_layer[j]->Fill(chanDis2*1000);
                                resid_per_layer[j]->Fill(residual*1000);
                                if((position_y - ypos[i]) >= 0.0) {
                                    above_resid_per_layer[j]->Fill(residual*1000);
                                    above_wire_per_layer[j]->Fill(chanDis2*1000);
                                }
                                if((position_y - ypos[i]) < 0.0) {
                                    below_resid_per_layer[j]->Fill(residual*1000);
                                    below_wire_per_layer[j]->Fill(chanDis2*1000);
                                }
                            }
                        }
                        wire_vs_channel->Fill(distancefrompoint2line2(xpos[i],ypos[i],newparam1,newparam0)*1000,channel3[i]);
                        resid_vs_channel->Fill(residual*1000,channel3[i]);
                        perpSlope = -1/newparam1;   //making track vertical
                        perpYIntercept = ypos[i] - (perpSlope*xpos[i]); //find where track meets yposition
                        hitPosX = (perpYIntercept - newparam0)/(newparam1-perpSlope);  //x position of the vertical point of the isochrone circle
                        hitPosY = (hitPosX*newparam1)+newparam0;								//y position of the vertical point of the isochrone circle
                        alpha = TMath::ACos(chanDis*1000/5.0);		//pythagoras theorem should replace this line
                        dx[i] = 2*(5.0 * TMath::Sin(alpha));			//distance travelled in the tube
                        if(dx[i]>0.0){
                        totsum += tot3[i];
                        dxsum += dx[i];
                        array_test1[iter] = tot3[i];
                        array_test2[iter] = dx[i];
                        array_test3[iter] = tot3[i]/dx[i];
                       // cout<< " totoverdx is " << array_test3[iter] << " tot is " << array_test1[iter] << " dx is " << array_test2[iter] << " iter is "<< iter << " i is " << i << endl;
                        iter++;
                        }
                        if(Riso[i] <=0.00026) resid1->Fill(residual);
                        if(Riso[i] > 0.00026 && Riso[i] <=0.00052) resid2->Fill(residual);
                        if(Riso[i] > 0.00052 && Riso[i] <=0.00078) resid3->Fill(residual);
                        if(Riso[i] > 0.00078 && Riso[i] <=0.00104) resid4->Fill(residual);
                        if(Riso[i] > 0.00104 && Riso[i] <=0.00130) resid5->Fill(residual);
                        if(Riso[i] > 0.00130 && Riso[i] <=0.00156) resid6->Fill(residual);
                        if(Riso[i] > 0.00156 && Riso[i] <=0.00182) resid7->Fill(residual);
                        if(Riso[i] > 0.00182 && Riso[i] <=0.00208) resid8->Fill(residual);
                        if(Riso[i] > 0.00208 && Riso[i] <=0.00234) resid9->Fill(residual);
                        if(Riso[i] > 0.00234 && Riso[i] <=0.00260) resid10->Fill(residual);
                        if(Riso[i] > 0.00260 && Riso[i] <=0.00286) resid11->Fill(residual);
                        if(Riso[i] > 0.00286 && Riso[i] <=0.00312) resid12->Fill(residual);
                        if(Riso[i] > 0.00312 && Riso[i] <=0.00338) resid13->Fill(residual);
                        if(Riso[i] > 0.00338 && Riso[i] <=0.00364) resid14->Fill(residual);
                        if(Riso[i] > 0.00364 && Riso[i] <=0.00390) resid15->Fill(residual);
                        if(Riso[i] > 0.00390 && Riso[i] <=0.00416) resid16->Fill(residual);
                        if(Riso[i] > 0.00416 && Riso[i] <=0.00442) resid17->Fill(residual);
                        if(Riso[i] > 0.00442 && Riso[i] <=0.00468) resid18->Fill(residual);
                        if(Riso[i] > 0.00468 && Riso[i] <=0.00494) resid19->Fill(residual);
                        if(Riso[i] > 0.00494 && Riso[i] <=0.00520) resid20->Fill(residual);
                        
                        for (Int_t j = 0; j< 150; j++) {
                            if(channel3[i] == j){
                                wire_per_channel[j]->Fill(chanDis2*1000);
                                resid_per_channel[j]->Fill(residual*1000);             ////fill residual per channel
                                rt_per_channel[j]->Fill(t03[i],chanDis2*1000);
                                rt_per_channel2[j]->Fill(t03[i],chanDis*1000);
                                if((position_y - ypos[i]) >= 0.0) {
                                    above_resid_per_channel[j]->Fill(residual*1000);
                                    above_wire_per_channel[j]->Fill(chanDis2*1000);
                                }
                                if((position_y - ypos[i]) < 0.0) {
                                    below_resid_per_channel[j]->Fill(residual*1000);
                                    below_wire_per_channel[j]->Fill(chanDis2*1000);
                                }
                            } // end if channel3[i] == j
                        }   ///end if j for channels
                    }       ///end if riso
                    else continue;
                }           ////end hits for loop
                multiplicitycheck->Fill(iiter);
                meanresiduals->Fill(sumresid*1000/iiter);
                
               // cout << " totsum is " << totsum << " dxsum is " << dxsum << " iter number is " << iter << endl;
               // cout << "------------------------------------------------" << endl;
                Int_t *index = new Int_t[iter];
               // TMath::Sort(iter,array_test3,index,kFALSE); //////sort the totoverdx values
                TMath::Sort(iter,array_test1,index,kFALSE); //////sort the totoverdx values
                
//                for (Int_t i = 0; i<iter; i++) {
//                    cout<< " totoverdx is " << array_test3[index[i]] << " tot is " << array_test1[index[i]] << " dx is " << array_test2[index[i]] << " iter is " << i << " index is " << index[i] << endl;
//                }
//                cout << "------------------------------------------------" << endl;
//
                    truncate10 = 0.9*iter;
                    //cout << " truncate10 is " << truncate10 << " iter number is "<< iter << endl;
                    for(Int_t i = 0; i <truncate10; i++){
                      //  cout<< " totoverdx is " << array_test3[index[i]] << " tot is " << array_test1[index[i]] << " dx is " << array_test2[index[i]] << " iter is " << i << " index is " << index[i] << endl;
                        totsum10 += array_test1[index[i]];
                        dxsum10 += array_test2[index[i]];
                    }
 
              //  cout << "------------------------------------------------" << endl;

                    truncate20 = 0.8*iter;
              //  cout << " truncate20 is " << truncate20 << " iter number is "<< iter << endl;
                    for(Int_t i = 0; i < truncate20; i++){
                   //     cout<< " totoverdx is " << array_test3[index[i]] << " tot is " << array_test1[index[i]] << " dx is " << array_test2[index[i]] << " iter is " << i << " index is " << index[i] << endl;
                        totsum20 += array_test1[index[i]];
                        dxsum20 += array_test2[index[i]];
                    }
                
              //  cout << "------------------------------------------------" << endl;

                    truncate30 = 0.7*iter;
            //    cout << " truncate10 is " << truncate30 << " iter number is "<< iter << endl;
                    for(Int_t i = 0; i < truncate30; i++){
              //          cout<< " totoverdx is " << array_test3[index[i]] << " tot is " << array_test1[index[i]] << " dx is " << array_test2[index[i]] << " iter is " << i << " index is " << index[i] << endl;
                        totsum30 += array_test1[index[i]];
                        dxsum30 += array_test2[index[i]];
                    }
                
                
//                tot_over_dx -> Fill(totsum/dxsum);
//                tot_over_dx_10 -> Fill(totsum10/dxsum10);
//                tot_over_dx_20 -> Fill(totsum20/dxsum20);
//                tot_over_dx_30 -> Fill(totsum30/dxsum30);
                  tot -> Fill(totsum/iter);
                  tot_10 -> Fill(totsum10/truncate10);
                  tot_20 -> Fill(totsum20/truncate20);
                  tot_30 -> Fill(totsum30/truncate30);

            } //end if neweventhits3
            
            
            
            
            //  } ///end if prob
            // c1->Print("event_display.pdf"); /// ---------------> FILL OF THE PDF FILE !!!!!!!!!!!
            tree->Fill();
        } // end of if a || y statement !!!!!!!!!!!!
        else break;
    }			//// END OF BIG FOR LOOP FOR EVENTS !!!!!!!!!!!!!
    
    Double_t newmean = 0.0;
    Double_t newerror = 0.0;
    
    for (Int_t i = 5; i< 155; i+=10) {
        proj_y_rt = t0_vs_iso->ProjectionY("py",i,i);
        proj_y_rt->Fit(g1,"R");
        proj_y_rt->Fit(g2,"R+");
        g1->GetParameters(&par[0]);
        g2->GetParameters(&par[3]);
        proj_y_rt -> Fit("total");
        newmean = (par[4]+TMath::Abs(par[1]))/2;
        newerror = (par[5]+par[2])/2;
        cout << " new mean is " << newmean << " new error is " << newerror << endl;
        new_rt -> SetBinContent(i,newmean);
        new_rt -> SetBinError(i,newerror);
    }
    new_rt->Fit("finalfit");
    fstream out_Riso ("out_Riso_0.55GeV_dabc16117205401.txt", ofstream::out);    ////save parameters from riso calculation
    for(Int_t i = 0; i < 5; i++){
        out_Riso << finalfit -> GetParameter(i) << endl;
    }
    out_Riso.close();
    
    resid1->Fit("gaus");        ///fit plots for isochrone error calculation
    resid2->Fit("gaus");
    resid3->Fit("gaus");
    resid4->Fit("gaus");
    resid5->Fit("gaus");
    resid6->Fit("gaus");
    resid7->Fit("gaus");
    resid8->Fit("gaus");
    resid9->Fit("gaus");
    resid10->Fit("gaus");
    resid11->Fit("gaus");        ///fit plots for isochrone error calculation
    resid12->Fit("gaus");
    resid13->Fit("gaus");
    resid14->Fit("gaus");
    resid15->Fit("gaus");
    resid16->Fit("gaus");
    resid17->Fit("gaus");
    resid18->Fit("gaus");
    resid19->Fit("gaus");
    resid20->Fit("gaus");
    
    iso_error->SetBinContent(1,(resid1->GetFunction("gaus"))->GetParameter(2));  ////get sigma form the fit of the plots for the isochrone error
    iso_error->SetBinContent(2,(resid2->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(3,(resid3->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(4,(resid4->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(5,(resid5->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(6,(resid6->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(7,(resid7->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(8,(resid8->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(9,(resid9->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(10,(resid10->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(11,(resid11->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(12,(resid12->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(13,(resid13->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(14,(resid14->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(15,(resid15->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(16,(resid16->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(17,(resid17->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(18,(resid18->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(19,(resid19->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(20,(resid20->GetFunction("gaus"))->GetParameter(2));
    
    t0_vs_iso2->Fit("finalfit");        ///fit of 4th order polynomial for the new r(t) calculation
    
    
    //c1->Print("event_display.pdf]"); /// ------------------> CLOSE THE PDF FILE !!!!!!!!!!!!
    
    //**************************** WRITING HISTOGRAMS ****************************
    
    risos->Write();
    multiplicity1->Write();
    multiplicity2->Write();
    multiplicity3->Write();
    multiplicitycheck->Write();
    chisqplot->Write();
    probplot->Write();
    distancefromwire->Write();
    distancefromisochrone->Write();
    meanresiduals->Write();
    residual_vs_r->Write();
    t0_vs_iso->Write();
    t0_vs_iso2->Write();
    new_rt->Write();
    wire_vs_channel->Write();
    resid_vs_channel->Write();
    resid1->Write();
    resid2->Write();
    resid3->Write();
    resid4->Write();
    resid5->Write();
    resid6->Write();
    resid7->Write();
    resid8->Write();
    resid9->Write();
    resid10->Write();
    resid11->Write();
    resid12->Write();
    resid13->Write();
    resid14->Write();
    resid15->Write();
    resid16->Write();
    resid17->Write();
    resid18->Write();
    resid19->Write();
    resid20->Write();
    iso_error->Write();
//    tot_over_dx -> Write();
//    tot_over_dx_10 -> Write();
//    tot_over_dx_20 -> Write();
//    tot_over_dx_30 -> Write();
    tot -> Write();
    tot_10 -> Write();
    tot_20 -> Write();
    tot_30 -> Write();
    for (Int_t j = 0; j< 150; j++) {
        wire_per_channel[j]->Write();
        resid_per_channel[j]->Write();
        rt_per_channel[j]->Write();
        rt_per_channel2[j]->Write();
        above_resid_per_channel[j]->Write();
        below_resid_per_channel[j]->Write();
        above_wire_per_channel[j]->Write();
        below_wire_per_channel[j]->Write();
    }
    for (Int_t j = 0; j< 7; j++) {
        new_rt_layer[j]->Write();
        new_rt_layer2[j]->Write();
        wire_per_layer[j]->Write();
        resid_per_layer[j]->Write();
        above_resid_per_layer[j]->Write();
        below_resid_per_layer[j]->Write();
        above_wire_per_layer[j]->Write();
        below_wire_per_layer[j]->Write();
    }
    //************************* WRITING ROOT FILE AND CLOSING ********************
    
    treefile->Close();
    myfile->Write();
    myfile->Close();
    
}           /////END OF MAIN FUNCTION TRACK -->TRACK()!!!!!!!!!!!!!!!!!!


//___________________________________________________function needed for gMinuit______________________________________________________________

void minuitfcn (Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag){
    chisq = 0;
    Double_t dis = 0.0, resid = 0.0;
    Int_t it = 0;
    Double_t aa = par[0];
    Double_t bb = par[1];
    for (Int_t i= 0; i < eventhits3; i++){
        if(Riso[i] >=0.0 && Riso[i] <=0.0052){
            it++;
            dis = distancefrompoint2line (xpos[i],ypos[i],bb,aa);
            resid = residual(dis,Riso[i]);
            //chisq += (resid*resid)/(sigma[i]*sigma[i]);
            chisq += (resid*resid)/(0.000200*0.000200);
        }//riso if
    }//hits loop
    f = chisq/(Double_t(it));
    result = chisq/(Double_t(it-2));
}

//___________________________________________________function for track minimizer______________________________________________________________

void trackMinimizer (Double_t& b,Double_t &a,Double_t &db,Double_t &da,Double_t chis){
    
    TMinuit *gMinuit = new TMinuit(2);
    gMinuit -> SetFCN(minuitfcn);
    gMinuit -> SetPrintLevel(-1);
    
    Double_t arglist[10];
    Int_t ierflag = 0;
    
    arglist[0] = 1;
    gMinuit -> mnexcm("SET ERR", arglist ,1,ierflag);
    
    gMinuit -> mnparm(0,"b",b,db,0,0,ierflag);
    gMinuit -> mnparm(1,"a",a,da,0,0,ierflag);
    
    arglist[0] = 500.;
    arglist[1] = 1.;
    gMinuit -> mnexcm("MIGRAD",arglist,2,ierflag);
    //gMinuit -> mnexcm("MINI",arglist,2,ierflag);
    
    
    gMinuit->GetParameter(0, newparam0, dnewparam0);
    gMinuit->GetParameter(1, newparam1, dnewparam1);
    
    delete gMinuit;
}

// ________________________________________ distance point to line _________________________________________

Double_t distancefrompoint2line (Double_t& x, Double_t& y,Double_t& slope,Double_t& yIntercept){
    return TMath::Abs(y - (slope * x) - yIntercept)/TMath::Sqrt((slope*slope)+1);
}

// ________________________________________ distance point to line 2 _________________________________________

Double_t distancefrompoint2line2 (Double_t& x2, Double_t& y2,Double_t& slope2,Double_t& yIntercept2){
    return (y2 - (slope2 * x2) - yIntercept2)/TMath::Sqrt((slope2*slope2)+1);
}

//__________________________________________residual calculation___________________________________________

Double_t residual (Double_t dis,Double_t radius){
    return TMath::Abs(dis - radius);
}

// __________________________line calculation ________________________________________________________________

Double_t* line(Double_t& x1, Double_t& x2,Double_t& y1,Double_t& y2){
    static Double_t aa[3];
    chisq = 0;
    Double_t dis = 0.0, resid = 0.0;
    Int_t it = 0;
    aa[1]= (y2-y1)/(x2-x1);
    aa[0]= y1 - (aa[1]*x1);
    for (Int_t i= 0; i < eventhits3; i++){
        if(Riso[i] >=0.0 && Riso[i] <= 0.0052){
            it++;
            dis = distancefrompoint2line (xpos[i],ypos[i],aa[1],aa[0]);
            resid = residual(dis,Riso[i]);
            //chisq += (resid*resid)/(sigma[i]*sigma[i]);
            chisq += (resid*resid)/(0.000200*0.000200);
        } //riso if
    } //hits loop
    aa[2] = chisq/(Double_t(it-2));
    return aa;
}
