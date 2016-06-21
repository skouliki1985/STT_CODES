
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



void Draw(){
    
//    TFile *thefile1 = TFile::Open("0.55GeV_dabc16117205401_track.root"); // open file
//    TFile *thefile2 = TFile::Open("0.75GeV_dabc16119220215_track.root"); // open file
//    TFile *thefile3 = TFile::Open("1.00GeV_dabc16120173116_track.root"); // open file
//    TFile *thefile4 = TFile::Open("2.95GeV_dabc16116172551_track.root"); // open file
//
//    TCanvas *c = new TCanvas("c","c");
//    c->cd();
//    TH1D * h1 = (TH1D*)thefile1 -> Get("totoverdx");
//    h1->Draw();
//    TH1D * h2 = (TH1D*)thefile2 -> Get("totoverdx");
//    h2->Draw("same");
//    TH1D * h3 = (TH1D*)thefile3 -> Get("totoverdx");
//    h3->Draw("same");
//    TH1D * h4 = (TH1D*)thefile4 -> Get("totoverdx");
//    h4->Draw("same");
    
    TH1D *plot = new TH1D("totoverdx","totoverdx",40,0.,4.0);
    plot-> SetBinContent(6,48.19);
    plot-> SetBinContent(8,42.83);
    plot-> SetBinContent(11,38.17);
    plot-> SetBinContent(30,33.09);

    plot -> Draw();
    

}