
#ifndef myevent_h
#define myevent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "EventIII.h"
#include "TDCHit.h"
#include "ADCHit.h"
#include <TObject.h>

const Int_t kMaxTDCChannels = 203;

class myevent {
    public :
    TTree          *fChain;
    Int_t           fCurrent;
    UInt_t          fUniqueID;
    UInt_t          fBits;
    Int_t           totalNTDCChannels;
    Int_t           TDCChannels_;
    UInt_t          TDCChannels_fUniqueID[kMaxTDCChannels];
    UInt_t          TDCChannels_fBits[kMaxTDCChannels];
    Int_t           TDCChannels_channel[kMaxTDCChannels];
    Double_t        TDCChannels_leadTime1[kMaxTDCChannels];
    Double_t        TDCChannels_trailTime1[kMaxTDCChannels];
    Double_t        TDCChannels_tot1[kMaxTDCChannels];
    Double_t        TDCChannels_referenceDiff1[kMaxTDCChannels];
    Double_t        TDCChannels_leadTimes[kMaxTDCChannels][100];
    Double_t        TDCChannels_trailTimes[kMaxTDCChannels][100];
    Double_t        TDCChannels_tots[kMaxTDCChannels][100];
    Double_t        TDCChannels_referenceDiffs[kMaxTDCChannels][100];
    Int_t           TDCChannels_hitsNum[kMaxTDCChannels];
    
    
    TBranch        *b_eventIII_fUniqueID;
    TBranch        *b_eventIII_fBits;
    TBranch        *b_eventIII_totalNTDCChannels;
    TBranch        *b_eventIII_TDCChannels_;
    TBranch        *b_TDCChannels_fUniqueID;
    TBranch        *b_TDCChannels_fBits;
    TBranch        *b_TDCChannels_channel;
    TBranch        *b_TDCChannels_leadTime1;
    TBranch        *b_TDCChannels_trailTime1;
    TBranch        *b_TDCChannels_tot1;
    TBranch        *b_TDCChannels_referenceDiff1;
    TBranch        *b_TDCChannels_leadTimes;
    TBranch        *b_TDCChannels_trailTimes;
    TBranch        *b_TDCChannels_tots;
    TBranch        *b_TDCChannels_referenceDiffs;
    TBranch        *b_TDCChannels_hitsNum;
    
    myevent(TTree *tree=0);
    virtual ~myevent();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual void 	  Correction();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myevent_cxx
myevent::myevent(TTree *tree) : fChain(0)
{
    if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dabc16117205401.hld_times_hits.root");       ///path of the root file to analyze
        if (!f || !f->IsOpen()) {
            f = new TFile("dabc16117205401.hld_times_hits.root");          ///path of the root file to analyze
        }
        f->GetObject("T",tree);
        
    }
    Init(tree);
}

myevent::~myevent()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t myevent::GetEntry(Long64_t entry)
{
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t myevent::LoadTree(Long64_t entry)
{
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void myevent::Init(TTree *tree)
{
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
    
    fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_eventIII_fUniqueID);
    fChain->SetBranchAddress("fBits", &fBits, &b_eventIII_fBits);
    fChain->SetBranchAddress("totalNTDCChannels", &totalNTDCChannels, &b_eventIII_totalNTDCChannels);
    fChain->SetBranchAddress("TDCChannels", &TDCChannels_, &b_eventIII_TDCChannels_);
    fChain->SetBranchAddress("TDCChannels.fUniqueID", TDCChannels_fUniqueID, &b_TDCChannels_fUniqueID);
    fChain->SetBranchAddress("TDCChannels.fBits", TDCChannels_fBits, &b_TDCChannels_fBits);
    fChain->SetBranchAddress("TDCChannels.channel", TDCChannels_channel, &b_TDCChannels_channel);
    fChain->SetBranchAddress("TDCChannels.leadTime1", TDCChannels_leadTime1, &b_TDCChannels_leadTime1);
    fChain->SetBranchAddress("TDCChannels.trailTime1", TDCChannels_trailTime1, &b_TDCChannels_trailTime1);
    fChain->SetBranchAddress("TDCChannels.tot1", TDCChannels_tot1, &b_TDCChannels_tot1);
    fChain->SetBranchAddress("TDCChannels.referenceDiff1", TDCChannels_referenceDiff1, &b_TDCChannels_referenceDiff1);
    fChain->SetBranchAddress("TDCChannels.leadTimes[100]", TDCChannels_leadTimes, &b_TDCChannels_leadTimes);
    fChain->SetBranchAddress("TDCChannels.trailTimes[100]", TDCChannels_trailTimes, &b_TDCChannels_trailTimes);
    fChain->SetBranchAddress("TDCChannels.tots[100]", TDCChannels_tots, &b_TDCChannels_tots);
    fChain->SetBranchAddress("TDCChannels.referenceDiffs[100]", TDCChannels_referenceDiffs, &b_TDCChannels_referenceDiffs);
    fChain->SetBranchAddress("TDCChannels.hitsNum", TDCChannels_hitsNum, &b_TDCChannels_hitsNum);
    Notify();
}

Bool_t myevent::Notify()
{
    return kTRUE;
}

void myevent::Show(Long64_t entry)
{
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t myevent::Cut(Long64_t entry)
{
    return 1;
}
#endif
