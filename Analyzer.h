/*
 Copyright: 
 2010 daid KAHL, OTA shinsuke, HASHIMOTO takashi 
 2011, 2012, 2013 daid KAHL

 This file is part of crabat

 crabat is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 crabat is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with crabat.  If not, see <http://www.gnu.org/licenses/>.
   
 In order to publish results obtained from crabat or any
 future derivative software, the authors must make reference
 to the paper describing it (in preparation): D. Kahl et al.
 ``Algorithms for Active Target Data Analysis'' 2013.  This is an 
 additional term to the license as stipulated in section 7b.
 Publication includes, but is not limited to: books, peer-
 reviewed journals, conference proceedings, annual reports, etc.
 
*/

/**
@file Analyzer.h
@author daid
@version 1.0
@brief Defines the Analyzer class and members
@date 09 Nov 2011 23:19:47  
*/

//The base code was derived from ROOT as follows:
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 22 01:17:19 2010 by ROOT version 5.26/00
// from TTree rawdata/rawdata
// found on file: hoge.root
//////////////////////////////////////////////////////////

#ifndef Analyzer_h
#define Analyzer_h

#include <iostream>
#include <TClonesArray.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include <TCutG.h>
#include <TGraph.h>
#include <vector>
using std::vector;

const Int_t kMaxfadc = 437;

/**
 * @class TArtFadcHit
 * @brief Array object for flash ADC (TPC) data
 * @details Contains all the variables needed for TPC data analysis\n
 * We call TClonesArray to implement it.
 */
class TArtFadcHit : public TObject {
public:
   static const int kMaxSample = 4000;
   int fID;
   int fHit;
   int fNSample;
   int fTime;
   Int_t *fSample; //[fNSample]
   Int_t *fClock; //[fNSample]
   TArtFadcHit () {
     fSample = new Int_t[kMaxSample];
     fClock  = new Int_t[kMaxSample];
   }
   virtual ~TArtFadcHit() {
     if (fSample) delete[] fSample;
     if (fClock) delete[]  fClock;
     fSample = fClock = 0;
   }
   virtual void Clear(Option_t*) {
     if (fSample) delete[]  fSample;
     if (fClock) delete[]  fClock;
     fSample = fClock = 0;
   }
   ClassDef(TArtFadcHit,2);
};

ClassImp(TArtFadcHit);

class TDetector {
//class TDetector : public TObject {
public :
  TDetector(){}
  virtual ~TDetector(){}
   ClassDef(TDetector,2);
   virtual vector<vector<Short_t> > SetTpcMap();
   virtual vector<vector<Double_t> > SetTpcBgain();
   virtual vector<Double_t> SetGeoXB();
   virtual vector<Double_t> SetGeoXC();
};

ClassImp(TDetector);
/** 
 * @class Analyzer
 * @brief Class for analyzing the data.
 * @details Creates the TChain and branch addresses for the passed tree.\n
 * Creates all the basic data variables.\n
 * ADC and TDC data are just done with simple variables.\n
 * Flash ADC data is handled with running TClonesArray on TArtFadcHit.\n
 */
class Analyzer {
private:
   TFile *histofile;

public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Define any constants -- note in class structure we cannot use static const style
   // as it is forbidden by ISO C++
   //#define        rf 59.70; // real RF value
   // Declaration of leaf types
   Double_t        fSiE[18][9];   //
   Double_t        fSiEcal[18][9];   //
   Double_t        fSiEmax[18];
   Double_t        fSiEStripmax[12];
   Short_t        fSiEStripmax_ch[12];
   Double_t        fSiT[18];   //
   Double_t        fSiTcal[18];   //
   Double_t	   trig_dly_tpcB[48];
   Double_t	   trig_dly_rf[2];
   Double_t        trig_dly_ppac[2][4];
   Bool_t 	   SiIsHit[18]; // wants to know if an SSD is hit based on fSiT>-100
   Bool_t 	   SsdOR; // want to know if the SSD-OR trigger was invoked based on 1a,2a,3a,5a,6b
   UShort_t	   ssd_mult; // SSD Multiplicity
   UShort_t        ssdhit[5][2]; // which ssd is triggered
   UShort_t        ssd_evts[6];
   Double_t        fPPAC[2][5];
   Double_t        fPPACcal[2][5];
   Double_t        fRF[2];
   Double_t        fRFcal[2];
   Double_t        tof[2]; 
   Double_t        Tof; 
   Bool_t          PpacIsHit[2];
  // TPC pad mapping 
   //Short_t tpc_ch[2][144];
   Bool_t          PadIsHit[144];
   Double_t	   fSampleRebin[1334]; // kMaxSample/3 elements
   UShort_t	   rebinCounter;
   Double_t	   rebinTemp;
   Bool_t 	   truepeak;
   Double_t        fSampleMax[144][20][20];
   Double_t        fSampleInt[144][20][2];
   Double_t        fClockMax[144][20][20];
   Double_t        fTimeMax[144][20][20];
   UShort_t           MultiHit[144][20];
   //UShort_t           fPulseNo[144];
   Short_t           fPulseNo[144];
   Double_t        fHitMax[144];
   Double_t        baseline[144][20];
   Double_t        baseline_dev[144][20];
   Double_t        trough;
   UShort_t           pad[144]; // place holder for the pad in physics loops
   // Beam pad
   UShort_t        TracksB[48];
   Bool_t	   PadIsHitB[48][20];
   Double_t        QratioB[48][20];
   Double_t        XpadB[48][20];
   Double_t        XpadBraw[48][20];
   Double_t        YpadB[48][20];
   Double_t        ZpadB[48][20];
   Double_t        dEpadB[48][20];
   Double_t        dEpadBcalib[48][20];
   Double_t        TpadB[48][20];
   UShort_t        padsHitB[20];
   // Left pad
   UShort_t        TracksL[8];
   Double_t        XpadL[8][20];
   Double_t        YpadL[8][20];
   Double_t        ZpadL[8][20];
   Double_t        dEpadL[8][20];
   Double_t        TpadL[8][20];
   Double_t        dEpadLTotal[20];
   UShort_t        padsHitL[20];
   // Right pad
   UShort_t        TracksR[8];
   Double_t        XpadR[8][20];
   Double_t        YpadR[8][20];
   Double_t        ZpadR[8][20];
   Double_t        dEpadR[8][20];
   Double_t        TpadR[8][20];
   Double_t        dEpadRTotal[20];
   UShort_t        padsHitR[20];
   // Center pad
   UShort_t        TracksC[8];
   Bool_t	   PadIsHitC[8][20];
   Double_t        QratioC[8][20];
   Double_t        XpadC[8][20];
   Double_t        YpadC[8][20];
   Double_t        ZpadC[8][20];
   Double_t        dEpadC[8][20];
   Double_t        TpadC[8][20];
   Double_t        dEpadCTotal[20];
   UShort_t        padsHitC[20];
   // ppac calib?
   Double_t PpacX[2], PpacY[2];
   Double_t TargetX, TargetY;
   TF1 *fTargetX ;
   TF1 *fTargetY ;
   
   TClonesArray *arr;
   
   // List of branches
   TBranch        *b_fSiE;   //!
   TBranch        *b_fSiT;   //!
   TBranch        *b_fPPAC;   //!
   TBranch        *b_fRF;   //!
   TBranch        *b_fadc_;   //!

   Analyzer(TChain *tree=0);
   virtual ~Analyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   void Init(TChain *tree);
   //virtual Short_t **SetTpcMap();
   virtual void SetRF();
   virtual UShort_t SetPpacSepZ();
   virtual Bool_t WindowCut(Double_t x=99., Double_t y=99.);
   //virtual Bool_t WindowCut(Double_t &x, Double_t &y);
   virtual void     Loop(Int_t run=0,
     Bool_t flag_raw=0, Bool_t flag_detail=0, 
     Bool_t flag_ssd=0, Bool_t flag_strip=0, 
     Bool_t flag_ppac=0, 
     Bool_t flag_tpc=0);
   virtual UShort_t PpacXYCalib(Bool_t flag_detail=0, Bool_t flag_ppac=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef Analyzer_cxx

Analyzer::Analyzer(TChain *tree)
{
// if parameter tree is not specified (or zero), exit 
   if (tree == 0) {
     return;
   }
   Init(tree);
}

/**
 * @brief Default destructor for Analyzer
 * @details deletes the Chain of the current file and the TClonesArray arr TArtFadcHit.
 */
Analyzer::~Analyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   delete arr;
}

Int_t Analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analyzer::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   
   fChain->SetBranchAddress("ssdE", fSiE, &b_fSiE);
   fChain->SetBranchAddress("ssdT", fSiT, &b_fSiT);

   fChain->SetBranchAddress("ppac", fPPAC, &b_fPPAC);
   fChain->SetBranchAddress("rf", fRF, &b_fRF);

   arr = new TClonesArray("TArtFadcHit");
   arr->Clear();
   fChain->SetBranchAddress("fadc", &arr);
   fChain->GetBranch("fadc")->SetAutoDelete(kFALSE);

   Notify();
}

Bool_t Analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef Analyzer_cxx
