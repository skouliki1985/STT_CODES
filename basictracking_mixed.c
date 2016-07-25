//
//  basictracking_mixed.c
//  
//
//  Created by Apostolou Alexandros on 13/07/16.
//
//

#include <stdio.h>
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
int output_ps = 0;   // event display output to postscript file yes=1, no=0
int output_text = 0; // text output of results per event
long int nEvt = 0;   // counter events read in
// event cuts for analysis
int nMinHits = 10;   // >50% effiency (24 hits per layer)
int nMaxHits = 35;   // single track selected

double rStraw = 0.00505;

// arrays for track hits
const Int_t maxHitNumber = 150;
double xPos[maxHitNumber];
double yPos[maxHitNumber];
double Riso[maxHitNumber];
double RisoSigma[maxHitNumber];
double sigma[maxHitNumber];

// cut hits for reasonable isochrones
//double RisoMax =  0.0051;
//double RisoMax =  0.00525;
//double RisoMax =  0.00545;
//double RisoMax =  0.00535;   // PeWi, 12.07. last used, 185µm sigma (160µm)
double RisoMax =  0.0051;
//double RisoMax =  0.00525;
double RisoMin = -0.00035;      // allow small time deviation causing negative radii from polynomial
//double TMaxCut =  175;          // tmax for rt curve, (poly6)
//double TMaxCut =  165;          // tmax for rt curve, (poly6)
//double TMaxCut =  155;          // tmax for rt curve, (poly6)
//double TMaxCut =  145;  // see 2nd legs, average: 140ns +- 15ns  -> +- 100micron spread
//double TMaxCut =  165;  // see 2nd legs, average: 140ns +- 15ns  -> +- 100micron spread
double TMaxCut =  165;  // see 2nd legs, average: 140ns +- 15ns  -> +- 100micron spread

double RisoSetMin =  0.000025;    // set all negative isochrones to 300 micron (assume +-300 micron error), fixed
//double RisoCut1 = 0.0008;       // min rIso (m) cut for 2.nd residual vs channel plot

int eventhits3= 0;
//int eventHitsCut2= 0;

// minuit chi2 fit
double newparam0 = 0.0,newparam1 = 0.0,dnewparam0 = 0.0,dnewparam1 = 0.0;
double chisq = 9999.0,result = 9999.0;
double meanDistance = 9999.0;
//double Riso[maxHitNumber],sigma[maxHitNumber],xPos[maxHitNumber],yPos[maxHitNumber];

// cut outlier hits, wrong/multiple hits
double preFitDistCut = 0.006;          // prefit track distance > 8mm to wire
//double preFitDistCut = 0.006;        // prefit track distance > 8mm to wire
//double OutlierDistCut1 = 0.0010;     // outlier cuts, distance to isochrone > 1000 micron
double OutlierDistCut1 = 0.008;        // outlier cuts, distance to isochrone > 1000 micron
double OutlierDistCut2 = 0.0007;       // outlier cuts, distance to isochrone > 700 micron
double OutlierDistSigmaCut  = 2.5;     // now used, distance for outlier cut (factor x sigma(r))
int nMaxOutliers = 7, nOutlierCut = 6, nOutliers = 0, it1 = 0;
double chanResidual = 0;

// try 2nd position calibration
//#define NSTR 96
// new beam setup (6 layers x 24 straws)
#define NSTR 144
// opposite sign from dist_vs_chan plot
// got 148 micron resolution for 1900V, with all yOffs set to 0
// next it, optimise individual channels
//                        1           5             10             15             20          24

double yOffset[NSTR+1] = {0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};






// ***************** FUNCTIONS USED ********************************
void minuitfcn (Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag);
double distancefrompoint2line (double& x, double& y,double& slope,double& yIntercept);
double distancefrompoint2line2 (double& x2, double& y2,double& slope2,double& yIntercept2);
double residual (double dis,double radius);
void trackMinimizer (double& b,double &a,double &db,double &da,double &chis);
double* line(double& x1, double& x2,double& y1,double& y2);


void basictracking::track()				//START OF MAIN FUNCTION TRACK -->TRACK()!!!!!!!!!!!!!!!
{
    // ***************** VARIABLES, ARRAYS, ARCS, STRINGS ********************************
    int neweventhits3,it,firstiter,seconditer,test;
    double prob,chanDis,residual,sumresid,chanDis2,residual2,sumresid2,perpYIntercept,perpSlope,hitPosX,hitPosY,alpha,totsum,dxsum;
    int channel3[maxHitNumber],channelhit3[maxHitNumber],row3[maxHitNumber],column3[maxHitNumber];              // number of channels
    double xPosition3[maxHitNumber],yPosition3[maxHitNumber],t03[maxHitNumber],t0channel3[maxHitNumber],tot3[maxHitNumber],totchannel3[maxHitNumber],dx[maxHitNumber],newxPos[maxHitNumber],newyPos[maxHitNumber],newradius[maxHitNumber];
    double prefitPointX[4],prefitPointY[4],slope[4],yIntercept[4],kSquared[4];
    double polParam[5] = {0.};
    
    double N = 0.0, Ni = 0.0, SumNi = 0.0, wuv = 0.0;				////in meters
    double Risoch=0.0, Ra = 0.00505, Rwire = 0.00001, Ravalanch = 0.0000;  ////in meters
    
    
    
    // PeWi, 11.07.16 new algo, 3 groups, but above/below wire r-t parameters
    //double Pup[5] = {-6.08987e-05, 5.72153e-05, -2.33973e-08, -2.04856e-09, 8.06522e-12 };
    //double Pdown[5] = {-2.76636e-05, 5.27614e-05, 9.45145e-08, -3.057e-09, 1.06791e-11 };
    
    // PeWi, 12.07.16, iter2, layerposi + riso shift, 165µm (185µm)
    
    double P1up[5]   = {-4.46663e-05,6.03263e-05,-1.19661e-07,-1.01793e-09,4.51796e-12 };
    double P1down[5] = {-2.53343e-05,5.19806e-05,9.17964e-08,-3.05888e-09,1.09956e-11 };
    
    double P2up[5]   = {-4.07661e-05,5.74247e-05,-4.21066e-08,-1.94359e-09, 7.99083e-12};
    double P2down[5] = {1.50478e-05,5.27055e-05,5.4894e-08,-2.32966e-09,7.11056e-12 };
    
    double P3up[5]   = {-6.21554e-05,5.42943e-05,1.9034e-08,-2.08662e-09,7.0285e-12 };
    double P3down[5] = {-1.73207e-05,5.60347e-05,1.07897e-08,-2.50026e-09,9.88379e-12 };
    
    /*
     // PeWi, 11.07.16
     double P1up[5]   = {-4.18741e-05, 5.9459e-05, -8.76039e-08, -1.36968e-09, 5.76252e-12};
     double P1down[5] = {-2.43979e-05, 5.2849e-05,  8.74385e-08, -3.08477e-09, 1.12098e-11};
     
     double P2up[5]   = {-3.5659e-05,  5.70617e-05, -2.30006e-08, -2.17979e-09, 8.89412e-12};
     double P2down[5] = {1.62148e-05,  5.21422e-05,  8.27911e-08, -2.69551e-09, 8.62618e-12};
     
     double P3up[5]   = {-4.51563e-05, 5.51824e-05, 1.02213e-08, -2.0701e-09,  7.1406e-12};
     double P3down[5] = {-2.76859e-05, 5.4289e-05,  4.39799e-08, -2.80029e-09, 1.09395e-11};
     */
    
    
    /*double P1up[5]   = {-7.99205e-05, 5.85984e-05,-5.73561e-08, -1.67123e-09, 6.71924e-12};
     double P1down[5] = {-6.57541e-05, 5.49421e-05, 5.50583e-08, -2.85397e-09, 1.05736e-11};
     
     double P2up[5]   = {-5.61665e-05, 5.71061e-05,-1.3328e-08,  -2.26683e-09, 9.10699e-12};
     double P2down[5] = {-1.83318e-05, 5.25649e-05, 8.78186e-08, -2.83904e-09, 9.42508e-12};
     
     double P3up[5]   = {-6.20864e-05, 5.52805e-05, 1.36255e-08, -2.15837e-09, 7.62948e-12};
     double P3down[5] = {-3.50927e-05, 5.45625e-05, 3.52146e-08, -2.64982e-09, 1.02648e-11};
     */
    
    
    
    /*double P4up[5] =   {-0.00010736,  6.0635e-05, -9.36106e-08, -1.45695e-09, 6.39522e-12};
     double P4down[5] = {-9.17884e-05, 5.47519e-05, 6.57073e-08, -2.9563e-09, 1.09246e-11};
     double P5up[5] =   {-9.58921e-05, 5.81318e-05,-4.14398e-08,- 1.78277e-09, 6.91533e-12 };
     double P5down[5] = {-5.62531e-05, 5.49072e-05, 4.00665e-08, -2.74309e-09, 1.06225e-11};
     double P6up[5] =   {-3.58333e-05, 5.50811e-05, 1.78626e-08, -2.31873e-09, 8.62684e-12};
     double P6down[5] = {-4.3046e-05,  5.39935e-05, 6.76085e-08, -2.90823e-09, 1.06741e-11};
     double P7up[5] =   {-2.75892e-05, 6.34195e-05, -1.57529e-07,-9.88429e-10, 5.07688e-12};
     double P7down[5] = {-6.49887e-05, 4.86653e-05, 2.09173e-07, -4.33507e-09, 1.54244e-11};
     double P8up[5] =   { 5.21076e-05, 5.46187e-05, -4.41175e-08, -1.47933e-09, 5.54743e-12};
     double P8down[5] = {-0.00010458, 5.72759e-05, -7.36817e-09, -2.33645e-09,9.27022e-12};
     */
    
    //double polParam[6] = {0.};
    double dslope = 0.002, dyIntercept = 0.002;
    // PeWi, values above from Alex, my used values below
    //double dslope = 1, dyIntercept = 1;
    double* aa;
    
    // reposition first block by -500 micron
    //   for (int i = 0; i < 60; i++){
    //         yOffset[i] =  -500.;
    //xOffset[i] =  xOffset[i] * 0.000001;
    //   }
    // change from micron to meter unit
    //   for (int i = 0; i < NSTR+1; i++){
    //         yOffset[i] =  yOffset[i] * 0.000001;
    //          //xOffset[i] =  xOffset[i] * 0.000001;
    //    }
    
    
    
    memset(prefitPointX,0,4*sizeof(double));
    memset(prefitPointY,0,4*sizeof(double));
    memset(kSquared,0,4*sizeof(double));
    memset(slope,0,4*sizeof(double));
    memset(yIntercept,0,4*sizeof(double));
    
    memset(Riso,0,maxHitNumber*sizeof(double));
    memset(polParam,0,7*sizeof(double));
    memset(xPos,0,maxHitNumber*sizeof(double));
    memset(yPos,0,maxHitNumber*sizeof(double));
    
    
    gStyle->SetFuncWidth(1);
    
    // my new histos
    TH2F* residual_vs_isochrone_all = new TH2F("residual_vs_isochrones (all ch)","Resid vs. isochr radius (all ch)",540,0.,5.4,400,-2.0,2.0);
    TH2F* residual_vs_isochrone_1 = new TH2F("residual_vs_isochrones (lay1)","Resid vs. isochr radius (lay1)",540,0.,5.4,400,-2.0,2.0);
    TH2F* residual_vs_isochrone_2 = new TH2F("residual_vs_isochrones (lay2)","Resid vs. isochr radius (lay2)",540,0.,5.4,400,-2.0,2.0);
    TH2F* residual_vs_isochrone_3 = new TH2F("residual_vs_isochrones (lay3)","Resid vs. isochr radius (lay3)",540,0.,5.4,400,-2.0,2.0);
    TH2F* residual_vs_isochrone_4 = new TH2F("residual_vs_isochrones (lay4)","Resid vs. isochr radius (lay4)",540,0.,5.4,400,-2.0,2.0);
    TH2F* residual_vs_isochrone_5 = new TH2F("residual_vs_isochrones (lay5)","Resid vs. isochr radius (lay5)",540,0.,5.4,400,-2.0,2.0);
    TH2F* residual_vs_isochrone_6 = new TH2F("residual_vs_isochrones (lay6)","Resid vs. isochr radius (lay6)",540,0.,5.4,400,-2.0,2.0);
    
    TH2F* Wdist_vs_drifttime_all = new TH2F("Wdist_vs_drifttime (all ch)","Track distance vs. drifttime (all ch)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_1 = new TH2F("Wdist_vs_drifttime (lay1)","Track distance vs. drifttime (lay1)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_2 = new TH2F("Wdist_vs_drifttime (lay2)","Track distance vs. drifttime (lay2)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_3 = new TH2F("Wdist_vs_drifttime (lay3)","Track distance vs. drifttime (lay3)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_4 = new TH2F("Wdist_vs_drifttime (lay4)","Track distance vs. drifttime (lay4)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_5 = new TH2F("Wdist_vs_drifttime (lay5)","Track distance vs. drifttime (lay5)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_6 = new TH2F("Wdist_vs_drifttime (lay6)","Track distance vs. drifttime (lay6)",800,-10.,190.,550,-5.5,5.5);
    
    TH2F* Wdist_vs_drifttime_1a = new TH2F("Wdist_vs_drifttimea (lay1)","Track distance vs. drifttime (lay1a)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_2a = new TH2F("Wdist_vs_drifttimea (lay2)","Track distance vs. drifttime (lay2a)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_3a = new TH2F("Wdist_vs_drifttimea (lay3)","Track distance vs. drifttime (lay3a)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_4a = new TH2F("Wdist_vs_drifttimea (lay4)","Track distance vs. drifttime (lay4a)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_5a = new TH2F("Wdist_vs_drifttimea (lay5)","Track distance vs. drifttime (lay5a)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_6a = new TH2F("Wdist_vs_drifttimea (lay6)","Track distance vs. drifttime (lay6a)",800,-10.,190.,550,-5.5,5.5);
    
    TH2F* Wdist_vs_drifttime_1b = new TH2F("Wdist_vs_drifttimeb (lay1)","Track distance vs. drifttime (lay1b)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_2b = new TH2F("Wdist_vs_drifttimeb (lay2)","Track distance vs. drifttime (lay2b)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_3b = new TH2F("Wdist_vs_drifttimeb (lay3)","Track distance vs. drifttime (lay3b)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_4b = new TH2F("Wdist_vs_drifttimeb (lay4)","Track distance vs. drifttime (lay4b)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_5b = new TH2F("Wdist_vs_drifttimeb (lay5)","Track distance vs. drifttime (lay5b)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_6b = new TH2F("Wdist_vs_drifttimeb (lay6)","Track distance vs. drifttime (lay6b)",800,-10.,190.,550,-5.5,5.5);
    
    TH2F* Wdist_vs_drifttime_1c = new TH2F("Wdist_vs_drifttimec (lay1)","Track distance vs. drifttime (lay1c)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_2c = new TH2F("Wdist_vs_drifttimec (lay2)","Track distance vs. drifttime (lay2c)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_3c = new TH2F("Wdist_vs_drifttimec (lay3)","Track distance vs. drifttime (lay3c)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_4c = new TH2F("Wdist_vs_drifttimec (lay4)","Track distance vs. drifttime (lay4c)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_5c = new TH2F("Wdist_vs_drifttimec (lay5)","Track distance vs. drifttime (lay5c)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_6c = new TH2F("Wdist_vs_drifttimec (lay6)","Track distance vs. drifttime (lay6c)",800,-10.,190.,550,-5.5,5.5);
    
    TH2F* Wdist_vs_drifttime_g1 = new TH2F("Wdist_vs_drifttime group1","Track distance vs. drifttime (grp1)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_g2 = new TH2F("Wdist_vs_drifttime group2","Track distance vs. drifttime (grpc)",800,-10.,190.,550,-5.5,5.5);
    TH2F* Wdist_vs_drifttime_g3 = new TH2F("Wdist_vs_drifttime group3","Track distance vs. drifttime (grpc)",800,-10.,190.,550,-5.5,5.5);
    
    
    
    
    TH2F* channel_vs_residuals = new TH2F("channel_vs_residuals","Channel vs Residuals",400,-1.,1.,150,0.,150.);
    TH1F* meanResiduals = new TH1F("meanResiduals","track Residuals in [mm] vs Entries",200,0.,1);
    TH1F* chisquare = new TH1F("chisquare","chi2 vs Entries",200,0.,20);
    TH1F* nhits_outlier = new TH1F ("nhits_outlier","No. outlier hits",20,0.,20.);
    TH1F* nhits_finalfit = new TH1F ("nhits_finalfit","No. of hits in final trackfit",50,0.,50.);
    TH2F* dist_vs_chan = new TH2F ("wdist_vs_chan","wire distance vs channel",150,0.,150.,240,-6.,6.);
    TH2F* riso_vs_chan = new TH2F ("riso_vs_chan","isochrone radius vs channel",150,0.,150.,240,-6.,6.);
    //
    
    TH1D* risos = new TH1D("risos","risos",1000,-1.0,9.0); // in mm
    TH1I* multiplicity = new TH1I("multiplicity","multiplicity",50,0,50);
    TH1I* multiplicity_diff = new TH1I("multiplicity_diff","multiplicity_diff",20,0,20);
    TH1D* chisqplot = new TH1D("chisqplot","chisqplot",500,0.,50.);
    TH1D* distancefromisochrone = new TH1D("distancefromisochrone","distancefromisochrone",1000,-1.,1.); ///in mm
    TH1D* chisqplot2 = new TH1D("chisqplot2","chisqplot2",500,0.,50.);
    TH1D* probplot = new TH1D("probplot","probplot",500,0.,1.);
    TH1D* distancefromisochrone2 = new TH1D("distancefromisochrone2","distancefromisochrone2",1000,-1.,1.); ///in mm
    TH1D* totoverdx = new TH1D("totoverdx","totoverdx",1000,0.0,100.);
    TH2D* residual_vs_r_1 = new TH2D("residual_vs_r_1","residual_vs_r_1",1000,-2.,2.,510,0.,5.1); ///in mm
    TH2D* residual_vs_r_2 = new TH2D("residual_vs_r_2","residual_vs_r_2",1000,-2.,2.,510,0.,5.1); ///in mm
    TH2D* t0_vs_iso = new TH2D("t0_vs_iso","t0_vs_iso",165,0.,165.,1100, -0.0055,0.0055);  ///in m
    TH2D* t0_vs_iso2 = new TH2D("t0_vs_iso2","t0_vs_iso2",165,0.,165.,550, 0.,0.0055);      ///in m
    
    TH2D* t0_vs_iso_pos = new TH2D("t0_vs_iso_pos","t0_vs_iso above wire",165,0.,165.,600, 0.,0.006);  ///in m
    TH2D* t0_vs_iso_neg = new TH2D("t0_vs_iso_neg","t0_vs_iso below wire",165,0.,165.,600, 0.,0.006);      ///in m
    
    TH2D* t0_vs_iso_g1_up = new TH2D("t0_vs_iso_g1_up","t0_vs_iso g1 above wire",170,0.,170.,600, 0.,0.006);  ///in m
    TH2D* t0_vs_iso_g1_do = new TH2D("t0_vs_iso_g1_do","t0_vs_iso g1 below wire",170,0.,170.,600, 0.,0.006);      ///in m
    TH2D* t0_vs_iso_g2_up = new TH2D("t0_vs_iso_g2_up","t0_vs_iso g2 above wire",170,0.,170.,600, 0.,0.006);  ///in m
    TH2D* t0_vs_iso_g2_do = new TH2D("t0_vs_iso_g2_do","t0_vs_iso g2 below wire",170,0.,170.,600, 0.,0.006);      ///in m
    TH2D* t0_vs_iso_g3_up = new TH2D("t0_vs_iso_g3_up","t0_vs_iso g3 above wire",170,0.,170.,600, 0.,0.006);  ///in m
    TH2D* t0_vs_iso_g3_do = new TH2D("t0_vs_iso_g3_do","t0_vs_iso g3 below wire",170,0.,170.,600, 0.,0.006);      ///in m
    
    TH1D* t0_g1_up = new TH1D("t0_g1_up","t0_g1_up",170,0.,170);
    
    TH1D* Riso_vs_tdrift_g1up = new TH1D ("Riso_vs_tdrift_g1up","Isochrone radius vs. drift time grp1 up",81,0.,162);
    TH1D* t0_cali_new = new TH1D ("t0_cali_new","t0_cali_new", 81,0.,162.);
    TF1* RisoFit = new TF1("RisoFit","pol4");
    
    
    /*TH1D* time_g1_up = new TH1D ("time_g1_up","time group1 above wire",200,-10.,190.);
     TH1D* time_g1_do = new TH1D ("time_g1_do","time group1 below wire",200,-10.,190.);
     TH1D* time_g2_up = new TH1D ("time_g2_up","time group2 above wire",200,-10.,190.);
     TH1D* time_g2_do = new TH1D ("time_g2_do","time group2 below wire",200,-10.,190.);
     TH1D* time_g3_up = new TH1D ("time_g3_up","time group3 above wire",200,-10.,190.);
     TH1D* time_g3_do = new TH1D ("time_g3_do","time group3 below wire",200,-10.,190.);
     */
    
    TH1D* resid1 = new TH1D("resid1","resid1",1000,-0.001,0.001); // in m
    TH1D* resid2 = new TH1D("resid2","resid2",1000,-0.001,0.001); // in m
    TH1D* resid3 = new TH1D("resid3","resid3",1000,-0.001,0.001); // in m
    TH1D* resid4 = new TH1D("resid4","resid4",1000,-0.001,0.001); // in m
    TH1D* resid5 = new TH1D("resid5","resid5",1000,-0.001,0.001); // in m
    TH1D* resid6 = new TH1D("resid6","resid6",1000,-0.001,0.001); // in m
    TH1D* resid7 = new TH1D("resid7","resid7",1000,-0.001,0.001); // in m
    TH1D* resid8 = new TH1D("resid8","resid8",1000,-0.001,0.001); // in m
    TH1D* resid9 = new TH1D("resid9","resid9",1000,-0.001,0.001); // in m
    TH1D* resid10 = new TH1D("resid10","resid10",1000,-0.001,0.001); // in m
    TH1D* iso_error = new TH1D("iso_error","iso_error",10,0.,5.);
    
    // r-t isochrone fit for groups and above/below wire
    TF1* rtfitup_g1 = new TF1("rtfitup_g1","pol4");
    TF1* rtfitdo_g1 = new TF1("rtfitdo_g1","pol4");
    TF1* rtfitup_g2 = new TF1("rtfitup_g2","pol4");
    TF1* rtfitdo_g2 = new TF1("rtfitdo_g2","pol4");
    TF1* rtfitup_g3 = new TF1("rtfitup_g3","pol4");
    TF1* rtfitdo_g3 = new TF1("rtfitdo_g3","pol4");
    
    //TF1* finalfit = new TF1("finalfit","pol5");
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
    string a;
    
    //cout << "filenames ********** Isochrone parameters " << filename1u << endl;
    
    //gStyle->SetOptStat(0);
    
    
    ifstream polParamData;							     	//declare stream
    polParamData.open("out_Riso_0.55GeV_dabc16117205401.txt");	//open riso txt file
    cout << "********** Isochrone parameters " << endl;
    //for(int i = 0; i < 6; i++){
    for(int i = 0; i < 5; i++){
        polParamData >> polParam[i];    //save in array values from riso txt file
        cout << "********** P["<< i <<"] "<< polParam[i] << endl;
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
    tree -> Branch("xPosition3",xPosition3,"xPosition3[150]/D");
    tree -> Branch("yPosition3",yPosition3,"yPosition3[150]/D");
    tree -> Branch("row3",row3,"row3[150]/I");
    tree -> Branch("column3",column3,"column3[150]/I");
    tree -> Branch("channelhit3",channelhit3,"channelhit3[150]/I");
    tree -> Branch("t0channel3",t0channel3,"t0channel3[150]/D");
    tree -> Branch("totchannel3",totchannel3,"totchannel3[150]/D");
    tree -> Branch("xPos3",xPos,"xPos[150]/D");
    tree -> Branch("yPos3",yPos,"yPos[150]/D");
    
    int nentries = readtree->GetEntriesFast();				//get entries from tree
    readtree -> SetBranchAddress("eventhits2", &eventhits3);
    readtree -> SetBranchAddress("channel2", &channel3);
    readtree -> SetBranchAddress("t02", &t03);
    readtree -> SetBranchAddress("tot2", &tot3);
    readtree -> SetBranchAddress("xposition2",&xPosition3);
    readtree -> SetBranchAddress("yposition2",&yPosition3);
    readtree -> SetBranchAddress("row2", &row3);
    readtree -> SetBranchAddress("column2", &column3);
    readtree -> SetBranchAddress("channelhit2", &channelhit3);
    readtree -> SetBranchAddress("t0channel2", &t0channel3);
    readtree -> SetBranchAddress("totchannel2", &totchannel3);
    
    
    cout << "  There are  " << nentries << " entries in this file....... " << endl;
    
    //c1->Print("event_display.pdf[","pdf"); // ---------------> START OF THE PDF FILE !!!!!!!!!!!!!!!!!!!
    //  if (output_ps == 1) c1->Print("event_display.ps["); // ---------------> START OF THE PDF FILE !!!!!!!!!!!!!!!!!!!
    
    
    // ******************************************************************************************
    // ****************************** start of big for loop for events !!!!!!  ******************
    for (int lm = 0; lm < nentries; lm++) {
        readtree -> GetEntry(lm);
        nEvt++;
        // first eventcut, min and max number of hits
        // check eventhits 3 (new) or 2 (from read in tree)
        
        if (eventhits3 < nMinHits || eventhits3 > nMaxHits) continue;
        
        if (lm % 50000 == 0) cout << " ==> " << lm << " events processed........ " << endl;
        if (a != "a"){
            cout << "Do you want to see the next event ? (y , n , a)" << endl;
            cin >> a;
        }
        
        
        
        for (int i = 0; i< maxHitNumber; i++) {     ///// set arrays equal to zero
            Riso[i] = 0.0;
            sigma[i] = 0.0;
            xPos[i] = 0.0;
            yPos[i] = 0.0;
            newxPos[i] = 0.0;
            newyPos[i] = 0.0;
            newradius[i] = 0.0;
            dx[i] = 0.0;
        }
        for (int i = 0; i< 4; i++) {            ///// set arrays equal to zero
            prefitPointX[i] = 0.0;
            prefitPointY[i] = 0.0;
            slope[i] = 0.0;
            kSquared[i] = 0.0;
            yIntercept[i] = 0.0;
        }
        
        
        
        neweventhits3 = eventhits3;             ///set eventhits3 to neweventhits3
        
        //        cout << " first eventhits3 is " << neweventhits3 << endl;
        
        if(a =="y" || a == "a"){    //// start of the if a || y statement !!!!!!!!!!!!!!!
            //            cout << " ******* ---------------------> start of event " << lm << " !!!!!!!!!!! " << endl;
            
            plot->Draw();
            
            for(int i = 0; i < eventhits3; i++){
                Riso[i] = (polParam[0]+(polParam[1]*t03[i])+(polParam[2]*pow(t03[i],2))+(polParam[3]*pow(t03[i],3))+(polParam[4]*pow(t03[i],4)));
                ///////riso is in meters!!!!!!!!!
                
                // allow small timing error, but shift to positive radii RisoSetMin
                // Riso unit mm here !
                if (Riso[i] < 0 && Riso[i] > RisoMin) Riso[i] = RisoSetMin;
                //if (Riso[i] < 0 && Riso[i] > RisoMin) Riso[i] = Riso[i]-RisoMin;
                
                //RisoSigma[i] = 250. - (150./5.*Riso[i]*1000.) + (30/(Riso[i]*1000.));    // in micron
                RisoSigma[i] = 250. - (150./5.*Riso[i]*1000.) + (15/(Riso[i]*1000.));    // in micron
                
                // change unit from micron to mm for plotting
                sigma[i] =  RisoSigma[i]*0.001;                    // in mm
                
                risos->Fill(Riso[i]*1000);
                
                xPos[i] = xPosition3[i];
                yPos[i] = yPosition3[i];
                
                // y-re-positioning
                //yPos[i] = yPos[i] + yOffset[channel3[i]]* 0.000001;
                if (channel3[i]>  0 && channel3[i] <  49 && channel3[i]%2 != 0) yPos[i] = yPos[i] -0.00013;
                if (channel3[i]>  0 && channel3[i] <  49 && channel3[i]%2 == 0) yPos[i] = yPos[i] +0.00013;
                if (channel3[i]> 48 && channel3[i] <  97 && channel3[i]%2 != 0) yPos[i] = yPos[i] +0.00013;
                if (channel3[i]> 48 && channel3[i] <  97 && channel3[i]%2 == 0) yPos[i] = yPos[i] -0.00013;
                if (channel3[i]> 97 && channel3[i] < 147 && channel3[i]%2 != 0) yPos[i] = yPos[i] -0.00013;
                if (channel3[i]> 97 && channel3[i] < 147 && channel3[i]%2 == 0) yPos[i] = yPos[i] +0.00013;
                
                
                //		if (channel3[i]> 0 && channel3[i] <=16 ) yPos[i]   = yPos[i] + 0.00020;
                //if (channel3[i] > 0 && channel3[i] < 8 ) yPos[i]   = yPos[i] + 0.00040;
                //if (channel3[i]> 0 && channel3[i] < 9 ) yPos[i]   = yPos[i] + 0.00040;
                //if (channel3[i]> 0 && channel3[i] < 49 && channel3[i]%2 != 0) yPos[i]   = yPos[i]-0.00005;;
                //if (channel3[i]> 0 && channel3[i] < 49 && channel3[i]%2 == 0) yPos[i]   = yPos[i]+0.00020;
                //if (channel3[i]> 48 && channel3[i] < 97 && channel3[i]%2 != 0) yPos[i]  = yPos[i]+0.0001;
                //if (channel3[i]> 48 && channel3[i] < 97 && channel3[i]%2 == 0) yPos[i]  = yPos[i]-0.00010;
                //if (channel3[i]> 97 && channel3[i] < 145 && channel3[i]%2 != 0) yPos[i] = yPos[i]-0.00005;
                //if (channel3[i]> 97 && channel3[i] < 145 && channel3[i]%2 == 0) yPos[i] = yPos[i]+0.00010;
                
                
                el1->DrawArc(xPos[i],yPos[i],rStraw);            ////draw tubes in the map
                if(Riso[i] >=0.0 && Riso[i] <= RisoMax) el2->DrawArc(xPos[i],yPos[i],Riso[i]);    /////draw ischrones in the map
                
                
                //yPos[i] = yPosition3[i] + yOffset[channel3[i]];
                // PeWi, 28.06.16, print
                if (output_text == 1) cout << i <<". ch: " << channel3[i] << " xPos: " << xPos[i]*1000. << " yPos: " << yPos[i]*1000. << " time: " << t03[i] << " riso: " << Riso[i]*1000. << " sigma(µm): " << sigma[i]*1000. <<endl;
                //el1->DrawArc(xPos[i],yPos[i],0.00505);            ////draw tubes in the map
                
                // check if isochrones reasonable
                if (Riso[i] < RisoMin || Riso[i] > RisoMax || t03[i] > TMaxCut) {
                    Riso[i] = 999;      // wrong times
                    neweventhits3--;
                }
                
            }       ///end hits for loop
            
            it = 0;                                 /// set it zero
            
            // *********************** start track pre fit (simple line fit) *********************************
            
            for (int j = 0; j<25; j++) {
                for(int i = 0; i < eventhits3; i++){
                    if(column3[i] == j && Riso[i] >=0.0 && Riso[i] <= RisoMax){
                        newxPos[it] = xPos[i];
                        newyPos[it] = yPos[i];
                        newradius[it] = Riso[i];
                        it++;
                    }       ///end if column and riso
                }			////end loop for hits
            }           ///end j for loop
            
            // calculate coordincate points for pre fit
            prefitPointY[0] = newyPos[0] - newradius[0];
            prefitPointY[1] = newyPos[0] + newradius[0];
            prefitPointY[2] = newyPos[it-1] - newradius[it-1];
            prefitPointY[3] = newyPos[it-1] + newradius[it-1];
            prefitPointX[0] = newxPos[0];
            prefitPointX[1] = newxPos[0];
            prefitPointX[2] = newxPos[it-1];
            prefitPointX[3] = newxPos[it-1];
            
            //            for (int i =0; i<4; i++) {
            //            cout << " prefit point x is " << prefitPointX[i] << " prefitpointY is " << prefitPointY[i] << endl;
            //            }
            //
            
            for(int i = 0; i < 2; i++){
                if(i == 0){			//first two lines-tracks
                    aa = line(prefitPointX[i],prefitPointX[i+2],prefitPointY[i],prefitPointY[i+2]);
                    kSquared[0] = aa[2];
                    slope[0] = aa[1];
                    yIntercept[0] = aa[0];
                    //                    cout << "first chi2 is " << kSquared[0] << endl;
                    aa = line(prefitPointX[i],prefitPointX[i+3],prefitPointY[i],prefitPointY[i+3]);
                    kSquared[1] = aa[2];
                    slope[1] = aa[1];
                    yIntercept[1] = aa[0];
                    //                    cout << "second chi2 is " << kSquared[1] << endl;
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
                    //                    cout << "third chi2 is " << kSquared[2] << endl;
                    aa = line(prefitPointX[i],prefitPointX[i+2],prefitPointY[i],prefitPointY[i+2]);
                    kSquared[3] = aa[2];
                    slope[3] = aa[1];
                    yIntercept[3] = aa[0];
                    //                    cout << "fourth chi2 is " << kSquared[3] << endl;
                    if (kSquared[3] - kSquared[2] < 0){			//comparing third and fourth line-track and choosing the best  according to q^2
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
            
            //   cout << " final prefit chi2 is " << kSquared[0] << endl;
            
            prefit -> SetParameters(yIntercept[0],slope[0]);
            prefit -> Draw("same");
            
            
            //cout << " chi2 of the first line is " << kSquared[0] << endl;
            chisqplot->Fill(kSquared[0]*1e6);
            
            // *********************** track pre fit done (simple line fit) *********************************
            
            
            // ************************************ removing hits far from preFit track (> preFitDist = 8 mm), accidentals, multi tracks
            // ********************
            
            
            for(int i = 0; i < eventhits3; i++){
                if(Riso[i] < RisoMax){
                    chanDis = distancefrompoint2line(xPos[i],yPos[i],slope[0],yIntercept[0]);
                    if(chanDis > preFitDistCut) {
                        Riso[i] = 999;
                        neweventhits3--;
                    }
                }
            }
            
            /*firstiter = 0;
             sumresid = 0.0;
             
             for(int i = 0; i < eventhits3; i++){
             if(Riso[i] >=0.0 && Riso[i] <= RisoMax){
             firstiter++;
             chanDis = 0.0;
             residual = 0.0;
             //chanDis = distancefrompoint2line(xPos[i],yPos[i],slope[0],yIntercept[0]);//distance of track from wire
             chanDis = distancefrompoint2line2(xPos[i],yPos[i],slope[0],yIntercept[0]);//distance of track from wire
             residual =  TMath::Abs(chanDis)-Riso[i];
             distancefromisochrone->Fill(residual*1000);
             residual_vs_r_1->Fill(residual*1000,Riso[i]*1000);
             if(chanDis > preFitDistCut){           //if track is too far in mm!!!!!!!!!!!!!!!!!!!
             Riso[i] = 999.9;
             neweventhits3--;
             }       ///end if chanDis
             
             if (output_text == 1) cout << i <<". ch: " << channel3[i] << " dist wire: " << chanDis*1000 << " rIso (mm): " << Riso[i]*1000. << " residual(mm): " << residual*1000 << endl;
             
             }       // end if Riso
             }       ///end hits for loop
             */
            
            //cout << " final fit chi2 is " << result << endl;
            
            
            //            if (result >= kSquared[0]) {
            //            cout << " event is " << lm << endl;
            //            test++;
            //            }
            
            
            
            // *********** first track fit after prefit
            if (neweventhits3 < nMinHits) {
                if (output_text == 1) cout << "Not enough hits ("<< neweventhits3 << ") --> track fit stopped --> next event" << endl;
                continue;
            }
            
            // TrackMinimizer 1st time
            trackMinimizer (yIntercept[0],slope[0],dyIntercept,dslope,kSquared[0]);		//track construction using minuit function etc
            fit -> SetParameters(newparam0,newparam1);
            fit -> Draw("same");
            
            if (output_text == 1) cout << "1st chi2 fit result: " << chisq*1e6 <<  " mean Dist (mm): " << meanDistance*1000.
                << " nFitHits: " << neweventhits3 << endl;
            
            // 2x outlier removals, 1.2mm then 0.7mm
            // 2 (3) outliers per iteration
            // removing the outlier, maximal 2 outlier hits per track
            // 1st iteration outlier removal <= 2 hits > 1.2 mm  ---------------------------------------------------------------------------------
            nOutliers = 0;
            it1 = 0;
            for(int i = 0; i < eventhits3; i++){
                if(Riso[i] < RisoMax ){
                    chanDis = distancefrompoint2line(xPos[i],yPos[i],newparam1,newparam0);
                    chanResidual = chanDis - Riso[i];
                    if (TMath::Abs(chanResidual) > OutlierDistCut1 && it1 <= nOutlierCut && neweventhits3 > nMinHits){
                        Riso[i] = 999.;
                        it1++;
                        nOutliers++;
                        neweventhits3--;
                    }
                }
            }
            // event cut, reco quality
            if (output_text == 1) cout << "1st chi2 fit 1. outlier check: noutlier = " << it1 << " nFitHits: " << neweventhits3 << endl;
            
            if (neweventhits3 < nMinHits) continue;
            if (nOutliers > nMaxOutliers) continue;
            
            // PeWi, 8.07.16, re-fill Riso with different isochrone parameters above/below wire
            
            for(int i = 0; i < eventhits3; i++){
                chanDis = distancefrompoint2line2(xPos[i],yPos[i],newparam1,newparam0);
                
                // ************************* group with correct r-t 140ns/140ns
                if ( (channel3[i] > 48 && channel3[i] < 65 && channel3[i]%2 !=0) ||
                    (channel3[i] > 64 && channel3[i] < 97 && channel3[i]%2 ==0) ||
                    (channel3[i] > 96 && channel3[i] < 128 ) ||
                    (channel3[i] > 127 && channel3[i] < 145 && channel3[i]%2 ==0 )) {
                    if (chanDis > 0) {
                        Riso[i] = (P1up[0]+(P1up[1]*t03[i])+(P1up[2]*pow(t03[i],2))+(P1up[3]*pow(t03[i],3))+(P1up[4]*pow(t03[i],4)));
                        //if (lm < 100) cout << "event: " << lm << " chandist > 0 " <<  channel3[i] << " " << chanDis << endl;
                    }
                    else {
                        Riso[i] = (P1down[0]+(P1down[1]*t03[i])+(P1down[2]*pow(t03[i],2))+(P1down[3]*pow(t03[i],3))+(P1down[4]*pow(t03[i],4)));
                    }
                }
                // ******************** group with 120ns/160ns r-t below/above wire r-t
                if ( (channel3[i] > 0 && channel3[i] < 48 ) ||
                    (channel3[i] > 48 && channel3[i] < 65 && channel3[i]%2 ==0) ||
                    (channel3[i] > 127 && channel3[i] < 145 && channel3[i]%2 !=0 )) {
                    if (chanDis > 0) {
                        Riso[i] = (P2up[0]+(P2up[1]*t03[i])+(P2up[2]*pow(t03[i],2))+(P2up[3]*pow(t03[i],3))+(P2up[4]*pow(t03[i],4)));
                        //if (lm < 100) cout << "event: " << lm << " chandist > 0 " <<  channel3[i] << " " << chanDis << endl;
                    }
                    else {
                        Riso[i] = (P2down[0]+(P2down[1]*t03[i])+(P2down[2]*pow(t03[i],2))+(P2down[3]*pow(t03[i],3))+(P2down[4]*pow(t03[i],4)));
                        //Riso[i] = Riso[i] + 0.0001;
                    }
                }
                // ******************** group with 150ns/130ns r-t below/above wire r-t
                if (channel3[i] > 64 && channel3[i] < 97 && channel3[i]%2 !=0 ) {
                    if (chanDis > 0) {
                        Riso[i] = (P3up[0]+(P3up[1]*t03[i])+(P3up[2]*pow(t03[i],2))+(P3up[3]*pow(t03[i],3))+(P3up[4]*pow(t03[i],4)));
                        //Riso[i] = Riso[i] + 0.0001;
                    }
                    else {
                        Riso[i] = (P3down[0]+(P3down[1]*t03[i])+(P3down[2]*pow(t03[i],2))+(P3down[3]*pow(t03[i],3))+(P3down[4]*pow(t03[i],4)));
                    }
                }
                
                // PeWi, 11.07.16, correct oberserved shifts in residual distributions
                if (channel3[i] == 96) Riso[i] = Riso[i] - 0.000048;
                if (channel3[i] == 90) Riso[i] = Riso[i] - 0.000026;
                if (channel3[i] == 89) Riso[i] = Riso[i] - 0.000014;
                if (channel3[i] == 86) Riso[i] = Riso[i] - 0.000040;
                if (channel3[i] == 85) Riso[i] = Riso[i] - 0.000014;
                if (channel3[i] == 82) Riso[i] = Riso[i] - 0.000072;
                
                // global residual shift in plot
                Riso[i] = Riso[i] - 0.000025;
                
            }    // end of for, Riso new setting
            
            
            
            
            
            // TrackMinimizer 2nd time
            if(it1 != 0){  // outlier hits found --> refit
                trackMinimizer (yIntercept[0],slope[0],dyIntercept,dslope,kSquared[0]);
                fit -> SetParameters(newparam0,newparam1);
                fit -> SetLineColor(kBlue);
                fit -> Draw("SAME");
                if (output_text == 1) cout << "2nd chi2 re-fit result: " << chisq*1e6 <<  " mean Dist (mm): " << meanDistance*1000.
                    << " nFitHits: " << neweventhits3 << " nOutlierHits: " << nOutliers <<endl;
            }
            
            
            // TrackMinimizer 3rd time, 2nd Outlier removal
            for(int i = 0; i < eventhits3; i++){
                if(Riso[i] < RisoMax ){
                    chanDis = distancefrompoint2line(xPos[i],yPos[i],newparam1,newparam0);
                    chanResidual = chanDis - Riso[i];
                    //if (output_text == 1) cout << "2nd chi2 outlier check: chan " << i << " residual= " <<  chanResidual << " sigma= " << sigma[i] << endl;
                    //if (TMath::Abs(chanResidual) > OutlierDistCut2 && it1 <= nOutlierCut && neweventhits3 > nMinHits){
                    if (TMath::Abs(chanResidual*1000/sigma[i]) > OutlierDistSigmaCut && it1 <= nOutlierCut && neweventhits3 > nMinHits){
                        Riso[i] = 999;
                        it1++;
                        nOutliers++;
                        neweventhits3--;
                    }
                }
            }
            if (output_text == 1) cout << "2nd chi2 fit 2. outlier check: noutlier = " << it1 << " nFitHits: " << neweventhits3 << endl;
            
            
            // event cut, reco quality
            if (neweventhits3 < nMinHits) continue;
            if (nOutliers > nMaxOutliers) continue;
            
            if(it1 != 0){  // outlier hits found --> refit
                trackMinimizer (yIntercept[0],slope[0],dyIntercept,dslope,kSquared[0]);
                fit -> SetParameters(newparam0,newparam1);
                fit -> SetLineColor(kBlue);
                fit -> Draw("SAME");
                if (output_text == 1) cout << "3rd chi2 re-fit result: " << chisq*1e6 <<  " mean Dist (mm): " << meanDistance*1000.
                    << " nFitHits: " << neweventhits3 << " nOutlierHits: " << nOutliers << endl;
            }
            
            //}
            
            
            
            
            // if(result ==0.0) cout << " event is " << lm << endl;
            chisqplot2->Fill(result*1e6);
            prob = 0.0;
            prob = TMath::Prob(result,neweventhits3-2);
            probplot -> Fill(prob);
            
            
            
            
            
            seconditer = 0;
            totsum = 0.0;                           /// set totsum zero
            dxsum = 0.0;
            
            it = 0;
            it1 = 0;
            
            //if(neweventhits3 > nMinHits){
            //  if(result < 20*1e6){
            for(int i = 0; i < eventhits3; i++){
                if(Riso[i] < RisoMax){
                    // cout << " channel is " << channel3[i] << endl;
                    it++;
                    seconditer++;
                    perpSlope = 0.0;
                    perpYIntercept = 0.0;
                    hitPosX = 0.0;
                    hitPosY = 0.0;
                    alpha = 0.0;
                    chanDis2 = 0.0;
                    chanDis = 0.0;
                    residual2 = 0.0;
                    chanDis2 = distancefrompoint2line(xPos[i],yPos[i],newparam1,newparam0);      //distance of track from wire
                    residual2 =  chanDis2-Riso[i];
                    distancefromisochrone2->Fill(residual2*1000);
                    residual_vs_r_2->Fill(residual2*1000,Riso[i]*1000);
                    t0_vs_iso -> Fill(t03[i],distancefrompoint2line2(xPos[i],yPos[i],newparam1,newparam0));
                    //PeWi, 4.07.16 ignore region close to wire for r-t fit
                    if (t03[i]>5 && t03[i] < TMaxCut) t0_vs_iso2 -> Fill(t03[i],chanDis2);
                    
                    
                    perpSlope = -1/newparam1;   //making track vertical
                    perpYIntercept = yPos[i] - (perpSlope*xPos[i]); //find where track meets yPosition
                    hitPosX = (perpYIntercept - newparam0)/(newparam1-perpSlope);  //x position of the vertical point of the isochrone circle
                    hitPosY = (hitPosX*newparam1)+newparam0;								//y position of the vertical point of the isochrone circle
                    alpha = TMath::ACos(chanDis2*1000/5.0);		//pythagoras theorem should replace this line
                    dx[i] = 2*(5.0 * TMath::Sin(alpha));			//distance travelled in the tube
                    totsum += tot3[i];
                    dxsum += dx[i];
                    
                    if(Riso[i] <=0.00052) resid1->Fill(residual2);
                    if(Riso[i] > 0.00052 && Riso[i] <=0.00104) resid2->Fill(residual2);
                    if(Riso[i] > 0.00104 && Riso[i] <=0.00156) resid3->Fill(residual2);
                    if(Riso[i] > 0.00156 && Riso[i] <=0.00208) resid4->Fill(residual2);
                    if(Riso[i] > 0.00208 && Riso[i] <=0.00260) resid5->Fill(residual2);
                    if(Riso[i] > 0.00260 && Riso[i] <=0.00312) resid6->Fill(residual2);
                    if(Riso[i] > 0.00312 && Riso[i] <=0.00364) resid7->Fill(residual2);
                    if(Riso[i] > 0.00364 && Riso[i] <=0.00416) resid8->Fill(residual2);
                    if(Riso[i] > 0.00416 && Riso[i] <=0.00468) resid9->Fill(residual2);
                    if(Riso[i] > 0.00468 && Riso[i] <=RisoMax) resid10->Fill(residual2);
                    
                    //PeWi, 11.07.16, calibrate channel groups and above/below wire
                    chanDis = distancefrompoint2line2(xPos[i],yPos[i],newparam1,newparam0);
                    chanResidual = TMath::Abs(chanDis) - Riso[i];
                    
                    // ************************* group with correct r-t 140ns/140ns
                    if (channel3[i] > 48 && channel3[i] < 65 && channel3[i]%2 !=0) {
                        if (chanDis > 0) {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g1_up -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                        else {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g1_do -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                    }
                    if (channel3[i] > 64 && channel3[i] < 97 && channel3[i]%2 ==0) {
                        if (chanDis > 0) {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g1_up -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                        else {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g1_do -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                    }
                    if (channel3[i] > 96 && channel3[i] < 128 ) {
                        if (chanDis > 0) {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g1_up -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                        else {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g1_do -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                    }
                    if (channel3[i] > 127 && channel3[i] < 145 && channel3[i]%2 ==0 ) {
                        if (chanDis > 0) {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g1_up -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                        else {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g1_do -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                    }
                    
                    // ******************** group with 120ns/160ns r-t below/above wire r-t
                    if (channel3[i] > 0 && channel3[i] < 48 ) {
                        if (chanDis > 0) {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g2_up -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                        else {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g2_do -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                    }
                    if (channel3[i] > 48 && channel3[i] < 65 && channel3[i]%2 ==0) {
                        if (chanDis > 0) {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g2_up -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                        else {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g2_do -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                    }
                    if (channel3[i] > 127 && channel3[i] < 145 && channel3[i]%2 !=0) {
                        if (chanDis > 0) {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g2_up -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                        else {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g2_do -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                    }
                    
                    // ******************** group with 150ns/130ns r-t below/above wire r-t
                    if (channel3[i] > 64 && channel3[i] < 97 && channel3[i]%2 !=0 ) {
                        if (chanDis > 0) {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g3_up -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                        else {
                            if (t03[i]>0 && t03[i] < TMaxCut) t0_vs_iso_g3_do -> Fill(t03[i],TMath::Abs(chanDis));
                        }
                    }
                    
                    
                    residual_vs_isochrone_all -> Fill(Riso[i]*1000,chanResidual*1000);
                    Wdist_vs_drifttime_all -> Fill(t03[i],chanDis*1000);
                    if (channel3[i]> 0 && channel3[i] < 49 && channel3[i]%2 != 0) {
                        residual_vs_isochrone_1 -> Fill(Riso[i]*1000,chanResidual*1000);
                        Wdist_vs_drifttime_1 -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] <= 16) Wdist_vs_drifttime_1a -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 16 && channel3[i] < 33) Wdist_vs_drifttime_1b -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 32 && channel3[i] < 49) Wdist_vs_drifttime_1c -> Fill(t03[i],chanDis*1000);
                    }
                    if (channel3[i]> 0 && channel3[i] < 49 && channel3[i]%2 == 0) {
                        residual_vs_isochrone_2 -> Fill(Riso[i]*1000,chanResidual*1000);
                        Wdist_vs_drifttime_2 -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] <= 16) Wdist_vs_drifttime_2a -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 16 && channel3[i] < 33) Wdist_vs_drifttime_2b -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 32 && channel3[i] < 49) Wdist_vs_drifttime_2c -> Fill(t03[i],chanDis*1000);
                    }
                    
                    if (channel3[i]> 48 && channel3[i] < 97 && channel3[i]%2 != 0) {
                        residual_vs_isochrone_3 -> Fill(Riso[i]*1000,chanResidual*1000);
                        Wdist_vs_drifttime_3 -> Fill(t03[i],chanDis*1000);
                        if (channel3[i]  <= 64) Wdist_vs_drifttime_3a -> Fill(t03[i],chanDis*1000);
                        if (channel3[i]  > 64 && channel3[i] < 81) Wdist_vs_drifttime_3b -> Fill(t03[i],chanDis*1000);
                        if (channel3[i]  > 81 && channel3[i] < 98) Wdist_vs_drifttime_3c -> Fill(t03[i],chanDis*1000);
                    }
                    if (channel3[i]> 48 && channel3[i] < 97 && channel3[i]%2 == 0) {
                        residual_vs_isochrone_4 -> Fill(Riso[i]*1000,chanResidual*1000);
                        Wdist_vs_drifttime_4 -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] <= 64) Wdist_vs_drifttime_4a -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 64 && channel3[i] < 81) Wdist_vs_drifttime_4b -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 81 && channel3[i] < 98) Wdist_vs_drifttime_4c -> Fill(t03[i],chanDis*1000);
                    }
                    if (channel3[i]> 97 && channel3[i] < 145 && channel3[i]%2 != 0) {
                        residual_vs_isochrone_5 -> Fill(Riso[i]*1000,chanResidual*1000);
                        Wdist_vs_drifttime_5 -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] <= 114)  Wdist_vs_drifttime_5a -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 114 && channel3[i] < 130) Wdist_vs_drifttime_5b -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 130 && channel3[i] < 147) Wdist_vs_drifttime_5c -> Fill(t03[i],chanDis*1000);
                    }
                    if (channel3[i]> 97 && channel3[i] < 145 && channel3[i]%2 == 0) {
                        residual_vs_isochrone_6 -> Fill(Riso[i]*1000,chanResidual*1000);
                        Wdist_vs_drifttime_6 -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] <= 114) Wdist_vs_drifttime_6a -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 114 && channel3[i] < 130) Wdist_vs_drifttime_6b -> Fill(t03[i],chanDis*1000);
                        if (channel3[i] > 130 && channel3[i] < 147) Wdist_vs_drifttime_6c -> Fill(t03[i],chanDis*1000);
                    }
                    channel_vs_residuals -> Fill(chanResidual*1000,channel3[i]);
                    dist_vs_chan -> Fill (channel3[i],chanDis*1000);
                    if (chanDis < 0) {
                        riso_vs_chan -> Fill (channel3[i],-Riso[i]*1000);
                    }
                    else {
                        riso_vs_chan -> Fill (channel3[i],Riso[i]*1000);
                    }
                    
                    if (output_text == 1) cout << it <<". ch: "<< channel3[i] <<" dist wire: "<< chanDis*1000 <<" rIso (mm): "<< Riso[i]*1000.<<" residual(mm): "<<chanResidual*1000 << endl;
                    
                    //         cout << " NEW channel is " << channel3[i] << " NEW distance from wire is " << chanDis*1000 << " NEW residual is " << (chanDis - Riso[i])*1000 << endl;
                }       //         end if riso
            }               //         end hits for loop
            
            
            totoverdx->Fill(totsum/dxsum);
            meanResiduals -> Fill(meanDistance*1000.);
            chisquare -> Fill(chisq*1e6);
            nhits_finalfit -> Fill(neweventhits3);
            nhits_outlier -> Fill(nOutliers);
            //}       ///end if result
            if (output_text == 1) cout << "chi2 outlier final fit result: " << chisq*1e6 <<  " mean Dist (mm): " << meanDistance*1000.
                << " nFitHits: " << neweventhits3 << " nOutlier " <<  nOutliers << endl;
            
            
            
            //c1->Print("event_display.pdf"); /// ---------------> FILL OF THE PDF FILE !!!!!!!!!!!
            // if (output_ps == 1 ) c1->Print("event_display.ps"); /// ---------------> FILL OF THE PDF FILE !!!!!!!!!!!
            multiplicity->Fill(neweventhits3);
            multiplicity_diff->Fill(eventhits3-neweventhits3);
            // if(seconditer <nMinHits) cout << " event is " << lm << endl;
            // my histos, fill fit result
            
            tree->Fill();
            
            // end of if a || y statement !!!!!!!!!!!!
        }
        else break;
        //        cout << " -----------> end of event " << lm << " !!!!!!!! " << endl;
    }			////end of big for loop for events !!!!!!!!!!!!!
    
    resid1->Fit("gaus");
    resid2->Fit("gaus");
    resid3->Fit("gaus");
    resid4->Fit("gaus");
    resid5->Fit("gaus");
    resid6->Fit("gaus");
    resid7->Fit("gaus");
    resid8->Fit("gaus");
    resid9->Fit("gaus");
    resid10->Fit("gaus");
    
    iso_error->SetBinContent(1,(resid1->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(2,(resid2->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(3,(resid3->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(4,(resid4->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(5,(resid5->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(6,(resid6->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(7,(resid7->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(8,(resid8->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(9,(resid9->GetFunction("gaus"))->GetParameter(2));
    iso_error->SetBinContent(10,(resid10->GetFunction("gaus"))->GetParameter(2));
    
    
    //    cout << "sigma is " << (resid1->GetFunction("gaus"))->GetParameter(2) << endl;
    cout << "************ Isochrone r-t fit ***************" << endl;
    t0_vs_iso2->Fit("finalfit");
    cout << "************ Isochrone r-t fit (above wire) ***************" << endl;
    t0_vs_iso_pos->Fit("finalfit");
    cout << "************ Isochrone r-t fit (below wire) ***************" << endl;
    t0_vs_iso_neg->Fit("finalfit");
    
    
    t0_vs_iso_g1_up ->Fit("rtfitup_g1");
    t0_vs_iso_g1_do ->Fit("rtfitup_g1");
    t0_vs_iso_g2_up ->Fit("rtfitup_g2");
    t0_vs_iso_g2_do ->Fit("rtfitup_g2");
    t0_vs_iso_g3_up ->Fit("rtfitup_g3");
    t0_vs_iso_g3_do ->Fit("rtfitup_g3");
    
    
    //***************** PeWi, 13.07.16, re-fit r-t from time spectra
    t0_g1_up -> Add(t0_vs_iso_g1_up->ProjectionX("_px"));
    /*t0_vs_iso_g1_do->ProjectionX();
     t0_vs_iso_g2_up->ProjectionX();
     t0_vs_iso_g2_do->ProjectionX();
     t0_vs_iso_g3_up->ProjectionX();
     t0_vs_iso_g3_do->ProjectionX();
     */
    //t0_vs_iso_g1_up_prox -> Add(t0_vs_iso_g1_up->ProjectionX("_px))
    //for(int j = 0; j < 81; j++ ){
    N = 0;
    for(int j = 0; j < 171; j++ ){
        //wuv =  t0_vs_iso_g1_up->ProjectionX("_px") -> GetBinContent(j+10);
        wuv =  t0_vs_iso_g1_up->ProjectionX("_px") -> GetBinContent(j);
        //wuv =  t0_g1_up -> GetBinContent(j+10);
        t0_cali_new -> SetBinContent(j, wuv);
        N += wuv;
        cout <<"Bin "<<j<<" entries "<<wuv<<" sumentries "<<N<<endl;
    }
    //for(int j = 0; j < 81; j++){
    for(int j = 0; j < 171; j++){
        Ni = t0_cali_new -> GetBinContent(j);
        SumNi += Ni;
        Risoch = ((SumNi/N)*(Ra - Rwire))+(Ravalanch);
        if(j%5==0){
            Riso_vs_tdrift_g1up -> SetBinContent(j,Risoch);
        }
        cout <<"Bin "<<j<<" entries "<<Ni<<" Riso "<<Risoch<<endl;
    }
    Riso_vs_tdrift_g1up -> SetMarkerStyle(3);
    Riso_vs_tdrift_g1up -> Fit("RisoFit","q");
    
    // ******************************* end re-fit r-t **********************
    
    /*time_g1_up -> Add(t0_vs_iso_g1_up -> ProjectionX("_px",1,600));
     time_g1_do -> Add(t0_vs_iso_g1_do -> ProjectionX("_px",1,600));
     time_g2_up -> Add(t0_vs_iso_g2_up -> ProjectionX("_px",1,600));
     time_g2_do -> Add(t0_vs_iso_g2_do -> ProjectionX("_px",1,600));
     time_g3_up -> Add(t0_vs_iso_g3_up -> ProjectionX("_px",1,600));
     time_g3_do -> Add(t0_vs_iso_g3_do -> ProjectionX("_px",1,600));
     */
    
    
    
    //    cout << "*************************************" << test << endl;
    
    //c1->Print("event_display.pdf]"); /// ------------------> CLOSE THE PDF FILE !!!!!!!!!!!!
    // if (output_ps ==1) c1->Print("event_display.ps]"); /// ------------------> CLOSE THE PDF FILE !!!!!!!!!!!!
    
    risos->Write();
    multiplicity->Write();
    multiplicity_diff->Write();
    chisqplot->Write();
    distancefromisochrone->Write();
    chisqplot2->Write();
    probplot->Write();
    distancefromisochrone2->Write();
    totoverdx->Write();
    residual_vs_r_1->Write();
    residual_vs_r_2->Write();
    t0_vs_iso->Write();
    t0_vs_iso2->Write();
    
    t0_vs_iso_g1_up ->Write();
    t0_vs_iso_g1_do ->Write();
    t0_vs_iso_g2_up ->Write();
    t0_vs_iso_g2_do ->Write();
    t0_vs_iso_g3_up ->Write();
    t0_vs_iso_g3_do ->Write();
    
    t0_cali_new ->Write();
    Riso_vs_tdrift_g1up -> Write();
    t0_g1_up ->Write();
    /*time_g1_up -> Write();
     time_g1_do -> Write();
     time_g2_up -> Write();
     time_g2_do -> Write();
     time_g3_up -> Write();
     time_g3_do -> Write();
     */
    
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
    iso_error->Write();
    
    // my new histos
    
    t0_vs_iso_pos->Write();
    t0_vs_iso_neg->Write();
    
    residual_vs_isochrone_all -> Write();
    residual_vs_isochrone_1 -> Write();
    residual_vs_isochrone_2 -> Write();
    residual_vs_isochrone_3 -> Write();
    residual_vs_isochrone_4 -> Write();
    residual_vs_isochrone_5 -> Write();
    residual_vs_isochrone_6 -> Write();
    
    Wdist_vs_drifttime_all -> Write();
    Wdist_vs_drifttime_1 -> Write();
    Wdist_vs_drifttime_2 -> Write();
    Wdist_vs_drifttime_3 -> Write();
    Wdist_vs_drifttime_4 -> Write();
    Wdist_vs_drifttime_5 -> Write();
    Wdist_vs_drifttime_6 -> Write();
    
    Wdist_vs_drifttime_1a -> Write();
    Wdist_vs_drifttime_2a -> Write();
    Wdist_vs_drifttime_3a -> Write();
    Wdist_vs_drifttime_4a -> Write();
    Wdist_vs_drifttime_5a -> Write();
    Wdist_vs_drifttime_6a -> Write();
    
    Wdist_vs_drifttime_1b -> Write();
    Wdist_vs_drifttime_2b -> Write();
    Wdist_vs_drifttime_3b -> Write();
    Wdist_vs_drifttime_4b -> Write();
    Wdist_vs_drifttime_5b -> Write();
    Wdist_vs_drifttime_6b -> Write();
    
    Wdist_vs_drifttime_1c -> Write();
    Wdist_vs_drifttime_2c -> Write();
    Wdist_vs_drifttime_3c -> Write();
    Wdist_vs_drifttime_4c -> Write();
    Wdist_vs_drifttime_5c -> Write();
    Wdist_vs_drifttime_6c -> Write();
    
    channel_vs_residuals -> Write();
    dist_vs_chan -> Write();
    riso_vs_chan -> Write();
    meanResiduals -> Write();
    chisquare -> Write();
    nhits_finalfit -> Write();
    nhits_outlier -> Write();
    //
    
    treefile->Close();
    myfile->Write();
    myfile->Close();
    
}           /////END OF MAIN FUNCTION TRACK -->TRACK()!!!!!!!!!!!!!!!!!!


//___________________________________________________function needed for gMinuit______________________________________________________________

//void minuitfcn (Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag){
void minuitfcn (Int_t& npar, Double_t* gin, Double_t& result, Double_t* par, Int_t iflag){
    chisq = 0.0;
    double dis = 0.0, resid = 0.0, mresid = 0.0;
    int it = 0;
    double aa = par[0];
    double bb = par[1];
    for (int i= 0; i < eventhits3; i++){
        if(Riso[i] >=0.0 && Riso[i] <= RisoMax){
            it++;
            dis = distancefrompoint2line (xPos[i],yPos[i],bb,aa);
            resid = residual(dis,Riso[i]);
            // cout << " residual from 2nd line is " << resid << endl;
            mresid += resid*resid;
            chisq += (resid*resid)/(sigma[i]*sigma[i]);
        }    //riso if
    }        //hits loop
    //cout << " chi2 from 2nd line is " << chisq << endl;
    //f = chisq/(double(it-2));
    result = chisq/(double(it));
    chisq = chisq/(double(it)-2);
    meanDistance = TMath::Sqrt(mresid/(double(it)));
}

//___________________________________________________function for track minimizer______________________________________________________________

void trackMinimizer (double& b,double &a,double &db,double &da,double &chis){
    
    //TMinuit *gMinuit = new TMinuit(5);    was before from Alex
    TMinuit *gMinuit = new TMinuit(2);    // 2 parameters fit
    gMinuit -> SetFCN(minuitfcn);
    gMinuit -> SetPrintLevel(-1);
    
    Double_t arglist[10];
    Int_t ierflag = 0;
    
    
    //arglist[0] = 1;  // intermediate FCNs strategy
    //arglist[0] = 2;   // extensive FCNs strategy, try to improve minimum
    //gMinuit -> mnexcm("SET STR",arglist,1,ierflag);
    
    
    arglist[0] = 1;   // chi2 fit
    gMinuit -> mnexcm("SET ERR", arglist ,1,ierflag);
    
    gMinuit -> mnparm(0,"b",b,db,0,0,ierflag);
    gMinuit -> mnparm(1,"a",a,da,0,0,ierflag);
    
    arglist[0] = 500.;
    arglist[1] = 1.;
    gMinuit -> mnexcm("MIGRAD",arglist,2,ierflag);
    gMinuit -> mnexcm("MINI",arglist,2,ierflag);
    //gMinuit -> mnexcm("MIG",arglist,2,ierflag);
    //gMinuit -> mnexcm("IMP",arglist,1,ierflag);
    gMinuit->GetParameter(0, newparam0, dnewparam0);
    gMinuit->GetParameter(1, newparam1, dnewparam1);
    
    //    if (ierflag > 0) cout << "Problem in Fit" << endl;
    
    delete gMinuit;
}

// ________________________________________ distance point to line _________________________________________

double distancefrompoint2line (double& x, double& y,double& slope,double& yIntercept){
    return TMath::Abs(y - (slope * x) - yIntercept)/TMath::Sqrt((slope*slope)+1);
}

// ________________________________________ distance point to line _________________________________________

double distancefrompoint2line2 (double& x2, double& y2,double& slope2,double& yIntercept2){
    return -(y2 - (slope2 * x2) - yIntercept2)/TMath::Sqrt((slope2*slope2)+1);
}

//__________________________________________residual calculation___________________________________________

double residual (double dis,double radius){
    return (dis - radius);
    //return TMath::Abs(dis - radius);
}

// __________________________line calculation ________________________________________________________________

double* line(double& x1, double& x2,double& y1,double& y2){
    static double aa[3];
    chisq = 0;
    double dis = 0.0, resid = 0.0, mresid = 0.0;
    int it = 0;
    aa[1]= (y2-y1)/(x2-x1);
    aa[0]= y1 - (aa[1]*x1);
    for (int i= 0; i < eventhits3; i++){
        if(Riso[i] >=0.0 && Riso[i] <= RisoMax){
            dis = distancefrompoint2line (xPos[i],yPos[i],aa[1],aa[0]);
            resid = residual(dis,Riso[i]);
            //cout << " residual from 1st line is " << resid << endl;
            it++;
            mresid += (resid*resid);
            chisq += (resid*resid)/(sigma[i]*sigma[i]);
            //chisq += (resid*resid)/(0.000150*0.000150);
        } //riso if
    } //hits loop
    //cout << " chi2 from 1st line is " << chisq << endl;
    aa[2] = chisq/(double(it-2));
    return aa;
}





