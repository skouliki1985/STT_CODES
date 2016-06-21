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
int eventhits3= 0;
double newparam0 = 0.0,newparam1 = 0.0,dnewparam0 = 0.0,dnewparam1 = 0.0;
double chisq = 9999.0,result = 9999.0;
double Riso[150],sigma[150],xpos[150],ypos[150];


// ***************** FUNCTIONS USED ********************************
void minuitfcn (Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag);
double distancefrompoint2line (double& x, double& y,double& slope,double& yIntercept);
double distancefrompoint2line2 (double& x2, double& y2,double& slope2,double& yIntercept2);
double residual (double dis,double radius);
void trackMinimizer (double& b,double &a,double &db,double &da);
double* line(double& x1, double& x2,double& y1,double& y2);


void basictracking::track()				//START OF MAIN FUNCTION TRACK -->TRACK()!!!!!!!!!!!!!!!
{
    // ***************** VARIABLES, ARRAYS, ARCS, STRINGS ********************************
    int neweventhits3,it,firstiter,seconditer,test;
    double prob,chanDis,residual,sumchanDis,sumresid,chanDis2,residual2,sumchanDis2,sumresid2,perpYIntercept,perpSlope,hitPosX,hitPosY,alpha,totsum,dxsum;
    int channel3[150],channelhit3[150],row3[150],column3[150];              // number of channels
    double xposition3[150],yposition3[150],t03[150],t0channel3[150],tot3[150],totchannel3[150],dx[150],newxpos[150],newypos[150],newradius[150];
    double prefitPointX[4],prefitPointY[4],slope[4],yIntercept[4],kSquared[4];
    double polParam[5] = {0.};
    double dslope = 0.002, dyIntercept = 0.002;
    double* aa;
    
    gStyle->SetFuncWidth(1);
    
    
    
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
    TH1D* resid_per_channel[150] = {0};         /////in mm
    
    
    
    
    
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
    TString Name;
    
    for (int i = 0; i< 150; i++) {     ///// set arrays equal to zero
        Name.Form(" Residual of channel %d " ,i);
        resid_per_channel[i] = new TH1D(Name,Name,1000,-1.0,1.0);
    }
    
    
    //gStyle->SetOptStat(0);
    
    
    ifstream polParamData;							     	//declare stream
    polParamData.open("out_Riso_0.55GeV_dabc16117205401.txt");				//open riso txt file
    for(int i = 0; i < 5; i++){
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
    
    int nentries = readtree->GetEntriesFast();				//get entries from tree
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
    //    c1->Print("event_display.pdf["); // ---------------> START OF THE PDF FILE !!!!!!!!!!!!!!!!!!!
    for (int lm = 0; lm < nentries; lm++) { 					////// start of big for loop for events !!!!!!
        readtree -> GetEntry(lm);
        if (lm % 50000 == 0) cout << " ==> " << lm << " events processed........ " << endl;
        if (a != "a"){
            cout << "Do you want to see the next event ? (y , n , a)" << endl;
            cin >> a;
        }
        
        
        
        
        
        for (int i = 0; i< 150; i++) {     ///// set arrays equal to zero
            Riso[i] = 0.0;
            sigma[i] = 0.0;
            xpos[i] = 0.0;
            ypos[i] = 0.0;
            newxpos[i] = 0.0;
            newypos[i] = 0.0;
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
                t03[i] = t03[i] - 7.0;
                Riso[i] = (polParam[0]+(polParam[1]*t03[i])+(polParam[2]*pow(t03[i],2))+(polParam[3]*pow(t03[i],3))+(polParam[4]*pow(t03[i],4)));///////riso is in meters!!!!!!!!!
                //Riso[i] = (-8.19526e-06+(5.52853e-05*t03[i])+(7.56188e-08*pow(t03[i],2))+(-3.31433e-09*pow(t03[i],3))+(1.28182e-11*pow(t03[i],4)));///////riso is in meters!!!!!!!!!
                sigma[i] = (2.15420e-04 - 4.40523e-05*Riso[i] + 7.10446e-06*Riso[i]*Riso[i]);
                risos->Fill(Riso[i]*1000);
                xpos[i] = xposition3[i];
                ypos[i] = yposition3[i];
                //                cout << " channel is " << channel3[i] << " xpos is " << xpos[i] << " ypos is " << ypos[i] << " t0 is " << t03[i] << " riso is " << Riso[i] << endl;
                el1->DrawArc(xpos[i],ypos[i],0.005);            ////draw tubes in the map
                if(Riso[i] >=0.0 && Riso[i] <= 0.0052) el2->DrawArc(xpos[i],ypos[i],Riso[i]);    /////draw ischrones in the map
            }       ///end hits for loop
            
            it = 0;                                 /// set it zero
            
            
            for (int j = 0; j<25; j++) {
                for(int i = 0; i < eventhits3; i++){
                    if(column3[i] == j && Riso[i] >=0.0 && Riso[i] <= 0.0052){
                        newxpos[it] = xpos[i];
                        newypos[it] = ypos[i];
                        newradius[it] = Riso[i];
                        it++;
                    }       ///end if column and riso
                }			////end loop for hits
            }           ///end j for loop
            
            prefitPointY[0] = newypos[0] - newradius[0];
            prefitPointY[1] = newypos[0] + newradius[0];
            prefitPointY[2] = newypos[it-1] - newradius[it-1];
            prefitPointY[3] = newypos[it-1] + newradius[it-1];
            prefitPointX[0] = newxpos[0];
            prefitPointX[1] = newxpos[0];
            prefitPointX[2] = newxpos[it-1];
            prefitPointX[3] = newxpos[it-1];
            
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
            prefit->Draw("same");
            
            
            //cout << " chi2 of the first line is " << kSquared[0] << endl;
            chisqplot->Fill(kSquared[0]);
            
            firstiter = 0;
            sumchanDis = 0.0;
            sumresid = 0.0;
            
            
            for(int i = 0; i < eventhits3; i++){
                if(Riso[i] >=0.0 && Riso[i] <= 0.0052){
                    firstiter++;
                    chanDis = 0.0;
                    residual = 0.0;
                    chanDis = distancefrompoint2line(xpos[i],ypos[i],slope[0],yIntercept[0]);//distance of track from wire
                    sumchanDis += chanDis;
                    residual = (chanDis-Riso[i]);
                    sumresid += TMath::Abs(residual);
                    distancefromisochrone->Fill(residual*1000);
                    residual_vs_r_1->Fill(residual*1000,Riso[i]*1000);
                    if(TMath::Abs(residual) > 0.0007){         //if track is too far in mm!!!!!!!!!!!!!!!!!!!
                        Riso[i] = -999.9;
                        neweventhits3--;
                    }       ///end if chanDis
                    //                    cout << " channel is " << channel3[i] << " distance from wire is " << chanDis*1000 << " residual is " << (chanDis - Riso[i])*1000 << endl;
                }       // end if Riso
            }       ///end hits for loop
            
            
            //            cout << " final fit chi2 is " << result << endl;
            
            
            //            if (result >= kSquared[0]) {
            //            cout << " event is " << lm << endl;
            //            test++;
            //            }
            
            if (neweventhits3>7) {
                trackMinimizer (yIntercept[0],slope[0],dyIntercept,dslope);		//track construction using minuit function etc
                fit -> SetParameters(newparam0,newparam1);
                fit->Draw("same");
            }
            
            //cout << " chi2 of the second line is " << result << endl;
            
            // if(result ==0.0) cout << " event is " << lm << endl;
            
//            for(int i = 0; i < eventhits3; i++){
//                if(Riso[i] >=0.0 && Riso[i] <= 0.0052){
//                    chanDis2 = 0.0;
//                    residual2 = 0.0;
//                    chanDis2 = distancefrompoint2line(xpos[i],ypos[i],newparam1,newparam0);//distance of track from wire
//                    residual2 = (chanDis2-Riso[i]);
//                    if(TMath::Abs(residual2) > 0.0007){         //if track is too far in mm!!!!!!!!!!!!!!!!!!!
//                        Riso[i] = -999.9;
//                        neweventhits3--;
//                    }       ///end if chanDis
//                }
//            }
            
            chisqplot2->Fill(result);
            prob = 0.0;
            prob = TMath::Prob(result,neweventhits3-2);
            probplot -> Fill(prob);
            

            
            seconditer = 0;
            sumchanDis2 = 0.0;
            sumresid2 = 0.0;
            totsum = 0.0;                           /// set totsum zero
            dxsum = 0.0;
            
            
            if(neweventhits3 > 7){
                if(result < 2){
                    for(int i = 0; i < eventhits3; i++){
                        if(Riso[i] >=0.0 && Riso[i] <= 0.0052){
                            // cout << " channel is " << channel3[i] << endl;
                            seconditer++;
                            perpSlope = 0.0;
                            perpYIntercept = 0.0;
                            hitPosX = 0.0;
                            hitPosY = 0.0;
                            alpha = 0.0;
                            chanDis2 = 0.0;
                            residual2 = 0.0;
                            chanDis2 = distancefrompoint2line(xpos[i],ypos[i],newparam1,newparam0);//distance of track from wire
                            sumchanDis2 += chanDis2;
                            residual2 = (chanDis2-Riso[i]);
                            sumresid2 += TMath::Abs(residual2);
                            distancefromisochrone2->Fill(residual2*1000);
                            residual_vs_r_2->Fill(residual2*1000,Riso[i]*1000);
                            t0_vs_iso -> Fill(t03[i],distancefrompoint2line2(xpos[i],ypos[i],newparam1,newparam0));
                            t0_vs_iso2 -> Fill(t03[i],chanDis2);
                            perpSlope = -1/newparam1;   //making track vertical
                            perpYIntercept = ypos[i] - (perpSlope*xpos[i]); //find where track meets yposition
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
                            if(Riso[i] > 0.00468 && Riso[i] <=0.00520) resid10->Fill(residual2);
                            for (int j = 0; j< 150; j++) {
                                if(channel3[i] == j){
                                    resid_per_channel[j]->Fill(residual2*1000);
                                }
                            }
                            
                            
                            
                            //                    cout << " NEW channel is " << channel3[i] << " NEW distance from wire is " << chanDis*1000 << " NEW residual is " << (chanDis - Riso[i])*1000 << endl;
                        }       ///end if riso
                    }           ////end hits for loop
                    
                    totoverdx->Fill(totsum/dxsum);
                }       ///end if result
                
                //         c1->Print("event_display.pdf"); /// ---------------> FILL OF THE PDF FILE !!!!!!!!!!!
                multiplicity->Fill(neweventhits3);
                multiplicity_diff->Fill(eventhits3-neweventhits3);
                // if(seconditer <7) cout << " event is " << lm << endl;
            }    ///end if neweventhits3
            
            
            tree->Fill();
            
            
        } // end of if a || y statement !!!!!!!!!!!!
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
    t0_vs_iso2->Fit("finalfit");
    
    //    cout << "*************************************" << test << endl;
    
    //   c1->Print("event_display.pdf]"); /// ------------------> CLOSE THE PDF FILE !!!!!!!!!!!!
    
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
    for (int j = 0; j< 150; j++) {
        resid_per_channel[j]->Write();
    }
    
    treefile->Close();
    myfile->Write();
    myfile->Close();
    
}           /////END OF MAIN FUNCTION TRACK -->TRACK()!!!!!!!!!!!!!!!!!!


//___________________________________________________function needed for gMinuit______________________________________________________________

void minuitfcn (Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag){
    chisq = 0;
    double dis = 0.0, resid = 0.0;
    int it = 0;
    double aa = par[0];
    double bb = par[1];
    for (int i= 0; i < eventhits3; i++){
        if(Riso[i] >=0.0 && Riso[i] <= 0.0052){
            it++;
            dis = distancefrompoint2line (xpos[i],ypos[i],bb,aa);
            resid = residual(dis,Riso[i]);
            // cout << " residual from 2nd line is " << resid << endl;
            //chisq += (resid*resid)/(sigma[i]*sigma[i]);
            chisq += (resid*resid)/(0.000150*0.000150);
        }//riso if
    }//hits loop
    //cout << " chi2 from 2nd line is " << chisq << endl;
    f = chisq/(double(it-2));
    result = f;
    //cout << "***************************" << result << endl;
}

//___________________________________________________function for track minimizer______________________________________________________________

void trackMinimizer (double& b,double &a,double &db,double &da){
    
    TMinuit *gMinuit = new TMinuit(5);
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
    return (y2 - (slope2 * x2) - yIntercept2)/TMath::Sqrt((slope2*slope2)+1);
}

//__________________________________________residual calculation___________________________________________

double residual (double dis,double radius){
    return TMath::Abs(dis - radius);
}

// __________________________line calculation ________________________________________________________________

double* line(double& x1, double& x2,double& y1,double& y2){
    static double aa[3];
    chisq = 0;
    double dis = 0.0, resid = 0.0;
    int it = 0;
    aa[1]= (y2-y1)/(x2-x1);
    aa[0]= y1 - (aa[1]*x1);
    for (int i= 0; i < eventhits3; i++){
        if(Riso[i] >=0.0 && Riso[i] <= 0.0052){
            dis = distancefrompoint2line (xpos[i],ypos[i],aa[1],aa[0]);
            resid = residual(dis,Riso[i]);
            //cout << " residual from 1st line is " << resid << endl;
            it++;
            //chisq += (resid*resid)/(sigma[i]*sigma[i]);
            chisq += (resid*resid)/(0.000150*0.000150);
        } //riso if
    } //hits loop
    //cout << " chi2 from 1st line is " << chisq << endl;
    aa[2] = chisq/(double(it-2));
    return aa;
}





