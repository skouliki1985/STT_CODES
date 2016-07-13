#define myevent_cxx

#include "myevent.h"
#include <TH2.h>
#include <TH1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <stdio.h>
#include <TF1.h>
#include <TMath.h>
#include <TSystem.h>
#include <TChain.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>

using namespace std;


TString stra = "0.55GeV_dabc16117205401";                       //file which will contain the plots
TString strb = ".root";											//file which will contain the plots
TString strc = "_corr";											//file which will contain the plots
TString filename = stra+strb;        							//file which will contain the plots
TString filename1 = stra+strc+strb;        						//file which will contain the plots


void myevent::Loop()						//START OF THE LOOP FUNCTION -->LOOP()
{
    // ******************************************************* variables and arrays ********************************************
    
    Double_t mean[150] = {0.};
    Double_t totCorrection[150] = {0.};
    Double_t t0Correction[150] = {0.};
    Double_t t0[150];
    Double_t t0channel[150];
    Double_t tot[150];
    Double_t totchannel[150];
    Int_t channel[150];
    Int_t channelhit[150];
    Double_t reftime,drifttimes,tottimes,drifttimes1,tottimes1;
    Double_t binMaxContent = 0.0,threshold = 0.0;
    Int_t xbin = 0,binMax = 0;
    Int_t savedata = 1;
    Int_t eventhits,it;
    
    // ******************************************************* histograms ********************************************************
    
    TH1D* reftime0 = new TH1D ("reftime0","reftime0",2000, -2000., 0.);
    TH1D* t0_raw_1 = new TH1D ("t0_raw_1","t0_raw_1",2000, -2000., 0.);
    TH2D* t0_vs_chan_multi = new TH2D("t0_vs_chan_multi","t0 vs channel all hits", 300,-850.,-550.,150,0.,150.);
    TH2D* t0_vs_chan_multi_1 = new TH2D("t0_vs_chan_multi_1","t0 vs channel 1st hits", 300,-850.,-550.,150,0.,150.);
    TH2D* tot_vs_chan_multi = new TH2D("tot_vs_chan_multi","tot vs channel all hits",300,150.,450.,150,0.,150.);
    TH2D* tot_vs_chan_multi_1 = new TH2D("tot_vs_chan_multi_1","tot vs channel 1st hits",300,150.,450.,150,0.,150.);
    TH2D* t0_vs_tot = new TH2D("t0_vs_tot","t0_vs_tot",300,-850.,-550.,300,150.,450.);
    TH2D* t0_vs_tot_1 = new TH2D("t0_vs_tot_1","t0_vs_tot_1",300,-850.,-550.,300,150.,450.);
    TH1D* projX_tot;        /// used only for the calibration method
    TH1D* projX_t0;         /// used only for the calibration method
    
    // ******************************************************* file with histos ********************************************************
    
    TFile *myfile = new TFile(filename,"RECREATE");     // creation of final file with histos
    
    // ******************************************************* general stuff ********************************************************
    
    gStyle->SetOptFit(1111111);
    if(fChain == 0) return;
    Int_t nentries = fChain->GetEntriesFast();
    
    cout << " ************** Total events in file ******************** ======> " << nentries <<"\n";
    
    TTree *tree = new TTree("data1","data1");
    tree -> Branch("eventhits",&eventhits,"eventhits/I");
    tree -> Branch("channel",channel,"channel[150]/I");
    tree -> Branch("t0",t0,"t0[150]/D");
    tree -> Branch("tot",tot,"tot[150]/D");
    tree -> Branch("t0channel",t0channel,"t0channel[150]/D");
    tree -> Branch("totchannel",totchannel,"totchannel[150]/D");
    tree -> Branch("channelhit",channelhit,"channelhit[150]/I");
    
    // ******************************************************* question to save correction files ********************************************************
    
    cout << "Do you want to overwrite the correction files? (1 = yes , 0 = no)" << endl;
    cin >> savedata;
    
    // ******************************************************* main function ********************************************************
    
    for(Int_t jentry=0; jentry < nentries; jentry++) {               // START OF THE MAIN FOR LOOP FOR EVENTS!!!
        Int_t ientry = LoadTree(jentry);
        if(ientry < 0) break;
        fChain->GetEntry(jentry); 				//read each entry of the events
        if(jentry % 50000 == 0)cout <<"  ===========> "<< jentry << " events processed....." << endl;
        
        for(Int_t i=0;i<150;i++){                 // set 0 all arrays
            channel[i] = 0;
            channelhit[i] = 0;
            t0[i] = 0.0;
            tot[i] = 0.0;
            t0channel[i] = 0.0;
            totchannel[i] = 0.0;
        }
        
        eventhits = 0;                          //set equal to 0 the eventhits variable
        it = 0;                                 //set equal to 0 the it variable
        
        for(Int_t i=0;i< totalNTDCChannels;i++){							//for loop for hits per event
            if (TDCChannels_channel[i] == 148){                     // if statement for reference channel 148 times to fill
                reftime = TDCChannels_leadTimes[i][0];              // lead times for channel 148
                reftime0 -> Fill(reftime);                           // fill plot for the times of the reference channel 148
            }   // end if statement
        }   // end for loop for hits per event
        
        
        for(Int_t i=0;i< totalNTDCChannels;i++){							//START OF THE MAIN FOR LOOP FOR HITS PER EVENT!!!
            t0_raw_1->Fill(TDCChannels_leadTimes[i][0]);                 //fill plot of the times without any cuts, really raw data
            for (Int_t j = 0; j < 2 ; j++){
                if(TDCChannels_leadTimes[i][j] >= -850.0 && TDCChannels_leadTimes[i][j]<=-550.0){
                    drifttimes = TDCChannels_leadTimes[i][j];
                    drifttimes1 = TDCChannels_leadTimes[i][0];
                    tottimes = TDCChannels_tots[i][j];
                    tottimes1 = TDCChannels_tots[i][0];
                    t0_vs_chan_multi -> Fill(drifttimes,TDCChannels_channel[i]);
                    t0_vs_chan_multi_1 -> Fill(drifttimes1,TDCChannels_channel[i]);
                    tot_vs_chan_multi -> Fill(tottimes,TDCChannels_channel[i]);
                    tot_vs_chan_multi_1 -> Fill(tottimes1,TDCChannels_channel[i]);
                    t0_vs_tot -> Fill(drifttimes,tottimes);
                    t0_vs_tot_1 -> Fill(drifttimes1,tottimes1);
                    channel[it] = TDCChannels_channel[i];
                    t0[it] = drifttimes1;
                    tot[it] = tottimes1;
                    channelhit[it] = 1;
                    t0channel[it] = drifttimes1;
                    totchannel[it] = tottimes1;
                    it++;
                    eventhits++;
                } ////end first if statement
            }   /////end j loop
        }       ////END OF THE MAIN FOR LOOP FOR HITS PER EVENT!!!!!!!!!!!!
       if (eventhits > 7 && eventhits < 41) tree -> Fill();          //fill tree with data only if we have more than 7 tubes hits per event/track
    }    ///END OF THE MAIN FOR LOOP FOR EVENTS!!!!!!!!!!!!!
    
    
    ///// ********************************* starting correction for times *******************************************
    
    for (Int_t i = 0; i < 150; i++){
        projX_tot = tot_vs_chan_multi_1 ->ProjectionX("_px",i,i);
        projX_tot -> Fit("gaus","q0","q0",130,330);
        if (projX_tot -> GetFunction("gaus")){
            TF1* f1 = projX_tot -> GetFunction("gaus");
            mean[i] = f1 -> GetParameter(1);
            
        }
        else{
            continue;
        }
        
        for (Int_t i = 0; i < 150; i++){
            totCorrection[i] = mean[i] - mean[130];
        }
    
        projX_t0 = t0_vs_chan_multi_1 -> ProjectionX("_px",i,i);
        binMax = projX_t0 -> GetMaximumBin();
        binMaxContent = projX_t0 -> GetBinContent(binMax);
        threshold = binMaxContent * 0.10;
        xbin = projX_t0 -> FindFirstBinAbove(threshold,1);
        t0Correction[i] = (xbin * (-550 + 850)/(300))-850;
        
    }
    
    if (savedata ==1)               ////save correction file for drift times
    {
        fstream out_t0Correction ("out_t0Correction_0.55GeV_dabc16117205401.txt", ofstream::out);
        for (Int_t i = 0; i < 150; i++){
            out_t0Correction << t0Correction[i] << endl;
        }
        out_t0Correction.close();
    }
    
    
    if (savedata ==1)           ////save correction file for tot times times
    {
        fstream out_totCorrection ("out_totCorrection_0.55GeV_dabc16117205401.txt", ofstream::out);
        for (Int_t i = 0; i < 150; i++){
            out_totCorrection << totCorrection[i] << endl;
        }
        out_totCorrection.close();
    }
    
    delete projX_t0;
    delete projX_tot;
    
    //**************************** WRITING HISTOGRAMS ****************************
    
    reftime0->Write();
    t0_raw_1->Write();
    t0_vs_chan_multi -> Write();
    t0_vs_chan_multi_1 -> Write();
    tot_vs_chan_multi -> Write();
    tot_vs_chan_multi_1 -> Write();
    t0_vs_tot -> Write();
    t0_vs_tot_1 -> Write();
    //************************* WRITING ROOT FILE AND CLOSING ********************
    myfile->Write();
    myfile->Close();
}											//END OF THE LOOP FUNCTION -->LOOP()



void myevent::Correction()						//START OF THE CORRECTION FUNCTION -->CORRECTION()
{
    
    // ******************************************************* variables and arrays ********************************************

    Int_t eventhits2,it,it1;
    Double_t xposition2[150],yposition2[150],t02[150],t0channel2[150],tot2[150],totchannel2[150];
    Double_t totCorrection2[150] = {0.};
    Double_t t0Correction2[150] = {0.};
    Int_t row2[150],column2[150],channel2[150],channelhit2[150];
    const Int_t rows = 7,columns = 25;
    Int_t strawMap[rows][columns] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,
        0,99,101,103,105,107,109,111,113,115,117,119,121,123,125,127,129,131,133,135,137,139,141,143,145,
        0,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,
        0,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,
        0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,
        0,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47};
    
    // ******************************************************* histograms ********************************************************

    TH1D* hits_per_channel = new TH1D ("hits_per_channel","hits_per_channel",150, 0., 150.);
    TH2D* t0_vs_chan_corr = new TH2D("t0_vs_chan_corr","t0 vs channel 1st hits corr", 160,-50.,270.,150,0.,150.);
    TH2D* tot_vs_chan_corr = new TH2D("tot_vs_chan_corr","tot vs channel 1st hits corr",300,150.,450.,150,0.,150.);
    TH2D* t0s_vs_tots_corr = new TH2D("t0s_vs_tots_corr","t0s_vs_tots_corr",160,-50.,270.,350,150.,450.);
    TH2D* map1 = new TH2D("map1","map1",300,0.,0.3,300,0.,0.06); // in meters!!!
    TH2D* map2 = new TH2D("map2","map2",30,0.,30.,10,0.,10.); // in meters!!!
    TH1D* t0_allchan_calibration = new TH1D ("t0_allchan_calibration","t0_allchan_calibration",160,-50.,270.);
    TH1D* t0_allchan_calibration_1 = new TH1D ("t0_allchan_calibration_1","t0_allchan_calibration_1",160,-50.,270.);
    TH1D* t0_allchan_calibration_2 = new TH1D ("t0_allchan_calibration_2","t0_allchan_calibration_2",160,-50.,270.);
    TH1D* t0_allchan_calibration_new = new TH1D ("t0_allchan_calibration_new","t0_allchan_calibration_new", 81,0.,162.);
    TH1D* Riso_vs_tdrift = new TH1D ("Riso_vs_tdrift","Isochrone radius vs. drift time",81,0.,162);
    TF1* RisoFit = new TF1("RisoFit","pol4");
    
    
    
    TFile* treefile = TFile::Open(filename);
    TTree* readtree = (TTree*)treefile -> Get("data1");
    
    TFile* myfile = new TFile(filename1,"RECREATE");
    TTree* tree = new TTree("data2","data2");
    tree -> Branch("eventhits2",&eventhits2);
    tree -> Branch("channel2",channel2,"channel2[150]/I");
    tree -> Branch("t02",t02,"t02[150]/D");
    tree -> Branch("tot2",tot2,"tot2[150]/D");
    tree -> Branch("channelhit2",channelhit2,"channelhit2[150]/I");
    tree -> Branch("t0channel2",t0channel2,"t0channel2[150]/D");
    tree -> Branch("totchannel2",totchannel2,"totchannel2[150]/D");
    tree -> Branch("xposition2",xposition2,"xposition2[150]/D");
    tree -> Branch("yposition2",yposition2,"yposition2[150]/D");
    tree -> Branch("row2",row2,"row2[150]/I");
    tree -> Branch("column2",column2,"column2[150]/I");
    
    // stream files
    ifstream totCorrectionData;
    ifstream t0CorrectionData;

    
    // **********************************************************************************************************
    
    Int_t nentries = readtree ->GetEntriesFast();
    cout << " Events in file for correction  ===> " << nentries <<"\n";
    
    
    readtree -> SetBranchAddress("eventhits", &eventhits2);
    readtree -> SetBranchAddress("channel", &channel2);
    readtree -> SetBranchAddress("t0", &t02);
    readtree -> SetBranchAddress("tot", &tot2);
    readtree -> SetBranchAddress("channelhit", &channelhit2);
    readtree -> SetBranchAddress("t0channel", &t0channel2);
    readtree -> SetBranchAddress("totchannel", &totchannel2);
    
    
    totCorrectionData.open("out_totCorrection_0.55GeV_dabc16117205401.txt");     // open input file
    t0CorrectionData.open("out_t0Correction_0.55GeV_dabc16117205401.txt");
    
    for(Int_t i = 0; i < 150; i++){
        totCorrectionData >> totCorrection2[i];
        t0CorrectionData >> t0Correction2[i];
    }
    totCorrectionData.close();
    t0CorrectionData.close();
    
    
    for (Int_t jentry=0; jentry<nentries; jentry++) {         ////start of for loop for events
        readtree -> GetEntry(jentry);
        
        for (Int_t i = 0; i < 150; i++){      /// set arrays equal to 0
            row2[i] = 0;
            column2[i] = 0;
            yposition2[i] = 0.0;
            xposition2[i] = 0.0;
        }
        
        it1=0;
        it=0;
        
        for (Int_t i=0; i < eventhits2; i++){
            hits_per_channel->Fill(channel2[i]);
            t02[i] = t02[i] - t0Correction2[channel2[i]+1];
            tot2[i] = tot2[i];// - totCorrection2[channel2[i]+1];
            t0channel2[i] = t0channel2[i] - t0Correction2[channel2[i]+1];
            totchannel2[i] = totchannel2[i];// - totCorrection2[channel2[i]+1];
            t0_vs_chan_corr->Fill(t02[i],channel2[i]);
            tot_vs_chan_corr->Fill(tot2[i],channel2[i]);
            t0s_vs_tots_corr->Fill(t02[i],tot2[i]);
            
            for (Int_t j = 1; j < rows ; j++){
                for (Int_t k = 1; k < columns ; k++){
                    if(channel2[i] == strawMap[j][k]){
                        yposition2[strawMap[j][k]] = (j*0.00878); ////in meters
                        if(channel2[i] > 0 && channel2[i] < 49){
                            if(channel2[i]%2 !=0) {xposition2[strawMap[j][k]] = 2*k*0.00507;} ////in meters
                            else {xposition2[strawMap[j][k]] = (2*k*0.00507) + 0.00507;} ////in meters
                        }
                        if(channel2[i] > 49 && channel2[i] < 98){
                            if(channel2[i]%2 ==0) {xposition2[strawMap[j][k]] = 2*k*0.00507;} ////in meters
                            else {xposition2[strawMap[j][k]] = (2*k*0.00507) + 0.00507;} ////in meters
                        }
                        if(channel2[i] > 98 && channel2[i] < 147){
                            if(channel2[i]%2 !=0) {xposition2[strawMap[j][k]] = 2*k*0.00507;} ////in meters
                            else {xposition2[strawMap[j][k]] = (2*k*0.00507) + 0.00507;} ////in meters
                        }
                        map1->Fill(xposition2[strawMap[j][k]],yposition2[strawMap[j][k]]); ///plot hitmap
                        map2->Fill(k,j); ///plot hitmap
                        row2[it] = j;
                        column2[it] = k;
                        yposition2[it] = (j*0.00878); ////in meters
                        if(channel2[i] > 0 && channel2[i] < 49){
                            if(channel2[i]%2 !=0) {xposition2[it] = 2*k*0.00507;} ////in meters
                            else {xposition2[it] = (2*k*0.00507) + 0.00507;} ////in meters
                        }
                        if(channel2[i] > 49 && channel2[i] < 98){
                            if(channel2[i]%2 ==0) {xposition2[it] = 2*k*0.00507;} ////in meters
                            else {xposition2[it] = (2*k*0.00507) + 0.00507;} ////in meters
                        }
                        if(channel2[i] > 98 && channel2[i] < 147){
                            if(channel2[i]%2 !=0) {xposition2[it] = 2*k*0.00507;} ////in meters
                            else {xposition2[it] = (2*k*0.00507) + 0.00507;} ////in meters
                        }
                        it++;
                        
                    } //if channel is element of the strawmap
                } // columns for loop
            } //rows for loop
            if(t02[i] >= 0.0 && t02[i] <= 155.0) it1++;
        } ///end for loop for hits
        if(it1 > 7 && it1 < 41) tree->Fill();
    } ////end for loop for events
    treefile->Close();
    
    Double_t N = 0.0, Ni = 0.0, SumNi = 0.0, wuv = 0.0;				////in meters
    Double_t Riso=0.0, Ra = 0.00505, Rwire = 0.00001, Ravalanch = 0.000050;  ////in meters
    
    // ******************************* Riso calculation ***********************************************
    for (Int_t i = 0; i < 150 ; i++){
        t0_allchan_calibration -> Add(t0_vs_chan_corr -> ProjectionX("_px",i,i));
    }
    
    t0_allchan_calibration_1 = (TH1D*)t0_allchan_calibration -> ShowBackground(30,"same");
    t0_allchan_calibration_2 -> Add(t0_allchan_calibration_1,t0_allchan_calibration,-1,1);
    
    for(Int_t i = 0; i < 81; i++ ){
        wuv = t0_allchan_calibration_2 -> GetBinContent(i+26);
        t0_allchan_calibration_new -> SetBinContent(i, wuv);
        N += wuv;
    }
    for(Int_t i = 0; i < 81; i++){
        Ni = t0_allchan_calibration_new -> GetBinContent(i);
        SumNi += Ni;
        Riso = ((SumNi/N)*(Ra - Rwire))+(Ravalanch);
        if(i%5==0){
        Riso_vs_tdrift -> SetBinContent(i,Riso);
        }
    }
    Riso_vs_tdrift -> SetMarkerStyle(3);
    Riso_vs_tdrift -> Fit("RisoFit","q");
    
    fstream out_Riso ("out_Riso_0.55GeV_dabc16117205401.txt", ofstream::out);    ////save parameters from riso calculation
    for(Int_t i = 0; i < 5; i++){
        out_Riso << RisoFit -> GetParameter(i) << endl;
    }
    out_Riso.close();
    
    //**************************** WRITING HISTOGRAMS ****************************

    hits_per_channel->Write();
    map1->Write();
    map2->Write();
    Riso_vs_tdrift->Write();
    t0_vs_chan_corr->Write();
    tot_vs_chan_corr->Write();
    t0s_vs_tots_corr->Write();
    t0_allchan_calibration->Write();
    t0_allchan_calibration_new->Write();
    
    //************************* WRITING ROOT FILE AND CLOSING ********************

    myfile -> Write();
    myfile -> Close();
    
} //END OF THE CORRECTION FUNCTION --> CORRECTION()


