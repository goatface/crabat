/*
 Copyright: 
 2010 daid KAHL, OTA shinsuke, HASHIMOTO takashi, YAMAGUCHI hidetoshi 
 2011, 2012 daid KAHL

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
*/

/**
@file Analyzer.cxx
@author daid
@version 1.0
@brief Analyzer for the run data.
@details Takes the Chain and loops over all entries.\n
Fill histograms, do physics.
@date 07 Nov 2011 18:57:28 
*/

//////////////////////////////////////
//    analyze the data here!
/////////////////////////////////////


// If run from ROOT, exit with message
#if defined (__CINT__)
{
cout << "This is C++ code, not a ROOT macro!" << endl;
cout << "To continue, run:" << endl;
cout << "make" << endl;
cout << "./run" << endl;
gApplication->Terminate();
}
#endif

#define Analyzer_cxx

#include "Analyzer.h" // all other headers are here
#include "Calibration.cxx"
#include "Analyzer_config.cxx"

Bool_t gate_30s=false,gate_29p=false,gate_17f=false;
Bool_t proton_beam=false;
Bool_t ribeam_ssd=false;
Bool_t goodRF_30s=false, goodPpacX_30s=false;
Bool_t goodRF_29p=false, goodPpacX_29p=false;
Bool_t goodRF_17f=false, goodPpacX_17f=false;

// for individual pulse analysis of selected cases
int pulsetrack=0;
int pulsetrack1a=0;
int pulsetrack2a=0;
int pulsetrack3a=0;
bool specialtrack=false;
// peak finding rebin size
const int rebin=4;
	    
	    double centerprojsum=0.;
	    int ytestcounter=0;
	    double ypointscounter=0;

double XpadBcalib[2];

const float padsize = 4.; // pad size mm
//const float padsize = (196./48.); // pad size mm

/**
 * @brief Loops over the TChain and entries
 * @details Oh my god there are so many details here.
 */
void Analyzer::Loop(Int_t run, 
     Bool_t flag_raw, Bool_t flag_detail, 
     Bool_t flag_ssd, Bool_t flag_strip, 
     Bool_t flag_ppac, 
     Bool_t flag_tpc)
{
  cout << "Analyzer::Loop processing run number " << run << endl;

  if (fChain == 0) return;
  //Fast function returns kBigNumber when it is TChain
  //Long64_t nentries = fChain->GetEntriesFast(); // I think this caused segfault -- daid
  Long64_t nentries = fChain->GetEntries();
  cout << "Entries: " << nentries << endl;

  Long64_t nbytes = 0, nb = 0;

  //Fancy color palette
  gStyle->SetPalette(1,0);

  //Calibration params

  Calibration *ssd_strip_calib=new Calibration("ssd_strip_calib.dat",97);//12x8
  Calibration *ssd_padE_calib=new Calibration("ssd_padE_calib.dat",19);
  Calibration *ssd_padT_calib=new Calibration("ssd_padT_calib.dat",19);
  Calibration *ppac_calib=new Calibration("ppac_calib.dat",10);
  Calibration *rf_calib=new Calibration("rf_calib.dat",2);
  Calibration *rf_ds_calib=new Calibration("rf_delay_calib.dat",2);
  Calibration *beam_x_left_calib=new Calibration("beam_x_left_calib.dat",48);
  Calibration *beam_x_right_calib=new Calibration("beam_x_right_calib.dat",48);
  Calibration *hg_center_calib=new Calibration("hg_center_calib.dat",8);
 
  TDetector MSTPC;
  vector<vector<Short_t> > tpc_ch=MSTPC.TDetector::SetTpcMap();
  vector<vector<Double_t> > padBgain=MSTPC.TDetector::SetTpcBgain();

  if (flag_detail){
    if (flag_strip){
      ssd_strip_calib->ScanFile();
      ssd_strip_calib->ShowEntries();
    }
    if (flag_ssd){
      ssd_padE_calib->ScanFile();
      ssd_padE_calib->ShowEntries();
      ssd_padT_calib->ScanFile();
      ssd_padT_calib->ShowEntries();
    }
    if (flag_ppac){
      ppac_calib->ScanFile();
      ppac_calib->ShowEntries();

      // PPAC downscale coincidence delay (nearly 805 ns, but here in ch)
      //there should be no need for these data? 07 Dec 2011 01:01:38 
      trig_dly_ppac[0][0]=8246.;
      trig_dly_ppac[0][1]=8244.;
      trig_dly_ppac[0][2]=8244.;
      trig_dly_ppac[0][3]=8242.;
      trig_dly_ppac[1][0]=8247.;
      trig_dly_ppac[1][1]=8246.;
      trig_dly_ppac[1][2]=8250.;
      trig_dly_ppac[1][3]=8248.;

      rf_calib->ScanFile();
      rf_calib->ShowEntries();
      rf_ds_calib->ScanFile();
      rf_ds_calib->ShowEntries();
      
      SetRF();

    }
    if (flag_tpc){
      
      beam_x_right_calib->ScanFile();
      beam_x_right_calib->ShowEntries();
      beam_x_left_calib->ScanFile();
      beam_x_left_calib->ShowEntries();
      
      hg_center_calib->ScanFile();
      hg_center_calib->ShowEntries();

    }
  } // end if: flag_detail
  //int ssdTime[18] = {0};
  for (int i=0;i<6;i++) ssd_evts[i]=0;
  static const UShort_t PpacSepZ = SetPpacSepZ();
  // target projection function
  // PPACa set as 0 position
  fTargetX = new TF1("fTargetX","[0]+((x/[2])*([1]-[0]))",0,1000); 
  // PPACb set as 0 -- idential results if we use the target position which subtracts [2] for calling
  // 	or we could make perfectly identical results if we summed that into the function itself
  //fTargetX = new TF1("fTargetX","[0]+((x+[2])*(([1]-[0])/[2]))",0,1000); 
  fTargetX->SetParNames("PPACaX","PPACbX","PPAC Separation");
  fTargetY = new TF1("fTargetY","[0]+((x/[2])*([1]-[0]))",0,1000);
  fTargetY->SetParNames("PPACaY","PPACbY","PPAC Separation");
   
   //Loop on entries
   for (Long64_t jentry=0; jentry<nentries-1;jentry++) {//-1...needed when the last event is not complete.
   // cout << jentry << endl;
      if (jentry%100==0 || jentry==nentries-2) {
         cout << " " << jentry << "/" << nentries-1 << "\r"<< flush;
      }
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

     //       arr->Clear(); 


      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //      cout << nb << "Bytes, Entry:"<< ientry << ","<<jentry << endl;
      //      fChain->Show(jentry);


      // if (Cut(ientry) < 0) continue;
	
      //Can begin filling any histrograms you defined
      
      //trigger conditions - even though using SSD data, any event may care about this
      //need to go through run.h and this code for needless use of flag_ssd and so on.
      for (UShort_t i=0; i<18; i++) {
        SiIsHit[i] = false;
        if ( fSiT[i] > -90 ) {
          SiIsHit[i] = true;
        } // end if : fSiT[i] > 0
      } // end for i
      // Set up the SSD-OR trigger condition
      SsdOR = false; // SSD-OR is always false until true
      ssd_mult = 0; // SSD Multiplicity is 0 until we know otherwise
      for ( UShort_t i=0;i<5;i++){ //which SSDs triggered?
        // SSD numbers
        if (i==0) ssdhit[i][0]=0;
        if (i==1) ssdhit[i][0]=1;
        if (i==2) ssdhit[i][0]=2;
        if (i==3) ssdhit[i][0]=4;
        if (i==4) ssdhit[i][0]=11;
        // By default, none are hit
        ssdhit[i][1]=0; 
      }
      // Experiment trigger condition before run 1008: 1a,2a,3a,6b
      if ( run < 1008 && ( SiIsHit[0] || SiIsHit[1] || SiIsHit[2] || SiIsHit[11] )) SsdOR = true; // SSD-OR trigger started the event
      // Experiment trigger condition from run 1008 and later: 1a,2a,3a,5a,6b
      if ( run > 1007 && ( SiIsHit[0] || SiIsHit[1] || SiIsHit[2] || SiIsHit[4] || SiIsHit[11] )) SsdOR = true; // SSD-OR trigger started the event
      // was 4a in the trigger?!  check the dE-E telescopes!
      // see online log Vol II p. 90 for record, but we don't know *when* it happened (so far did not check)

      if ( SsdOR ) { // Get the SSD trigger multiplicity
        if ( SiIsHit[0] ) {
          ssd_mult++;
          ssdhit[0][1]=1;
        }
        if ( SiIsHit[1] ) {
          ssd_mult++;
          ssdhit[1][1]=1;
        }
        if ( SiIsHit[2] ) {
          ssd_mult++;
          ssdhit[2][1]=1;
        }
        if ( SiIsHit[4] && run > 1007 ) {
          ssd_mult++;
          ssdhit[3][1]=1;
        }
        if ( SiIsHit[11] ) {
          ssd_mult++;
          ssdhit[4][1]=1;
        }
      }

      //PPAC
      if (flag_ppac) { // process PPAC?
	if (flag_raw) { // raw PPAC data?
	  for (UShort_t i=0;i<2;i++) if (fRF[i] > 0. && !SsdOR) hRf_ds[i]->Fill(fRF[i]);
	  for (UShort_t i=0;i<2;i++) if (fRF[i] > 0. && SsdOR) hRf_ssd[i]->Fill(fRF[i]);
          
          //trigger timing
	  if (1){
	  //if (!SsdOR){
	  //raw PPAC data
          if (fPPAC[0][0]>0.) hPpac0TX1->Fill(fPPAC[0][0]);
          if (fPPAC[0][1]>0.) hPpac0TX2->Fill(fPPAC[0][1]);
          if (fPPAC[0][2]>0.) hPpac0TY1->Fill(fPPAC[0][2]);
          if (fPPAC[0][3]>0.) hPpac0TY2->Fill(fPPAC[0][3]);
          if (fPPAC[1][0]>0.) hPpac1TX1->Fill(fPPAC[1][0]);
          if (fPPAC[1][1]>0.) hPpac1TX2->Fill(fPPAC[1][1]);
          if (fPPAC[1][2]>0.) hPpac1TY1->Fill(fPPAC[1][2]);
          if (fPPAC[1][3]>0.) hPpac1TY2->Fill(fPPAC[1][3]);
          }
	} // end if: flag_raw
	
	if (flag_detail) { // detailed PPAC data?
          //Calibration 
	  for (UShort_t i=0;i<2;i++) {
	    if (fRF[i] > 0) {
	      if (SsdOR) {
	        fRFcal[i]=rf_calib->Calib(i,fRF[i]);
	      }else{  
	        fRFcal[i]=rf_ds_calib->Calib(i,fRF[i]);
	        // the delay calibration can lead to negative RF
	        // in that case, wrap-around the true RF.
	        if (fRFcal[i]<0.) fRFcal[i]+=rf;
	      }
	      hRfcal[i]->Fill(fRFcal[i]);
	    }
	  }
          for(UShort_t i=0;i<2;i++){
	  PpacIsHit[i] = false;
	  if (fPPAC[i][0]>0. && fPPAC[i][1]>0. && fPPAC[i][2]>0. && fPPAC[i][3]>0.) PpacIsHit[i]=true;
            for(UShort_t j=0;j<4;j++){
              fPPACcal[i][j]=0.;
	      //check here for PPACb X centering issue 18 Jul 2011 17:37:49 
	      if ( fPPAC[i][j] > 0. )  {
	        if (!SsdOR) fPPAC[i][j]-=trig_dly_ppac[i][j];
		fPPACcal[i][j]=ppac_calib->Calib(i*5+j,fPPAC[i][j]);
	      }
            }    
          }
          tof[0]=tof[1]=-100.;
	  if (PpacXYCalib(flag_detail,flag_ppac)) cout << "PPAC XY Calibration failed!" << endl;
          if (tof[0]>0. && tof[1]>0. && PpacIsHit[0] && PpacIsHit[1] ) Tof=tof[1]-tof[0]; 
          //PpacX[1]=PpacX[1]-2; 
	  // fill basic PPAC histograms
	  if (PpacIsHit[0] && PpacIsHit[1]){
	    hRF0Tof->Fill(fRFcal[0],Tof);
            hRF1Tof->Fill(fRFcal[1],Tof);
	  }
          if (PpacIsHit[0]) hPpac0XY->Fill(PpacX[0],PpacY[0]);
          if (PpacIsHit[1]) hPpac1XY->Fill(PpacX[1],PpacY[1]);
          if (fRF[0] > 0.)
          {
	    //if (1){
	    if (PpacIsHit[0] ){
              //if ((fSiTcal[0]-tof[0])>290 && (fSiTcal[0]-tof[0])<295.5) { // 30S careful, can crash production runs, used for energyloss
                hPpac0XRF0->Fill(PpacX[0],fRFcal[0]);
	      //}
            }
	    //if (1){
	    if (PpacIsHit[1] ){
	      hPpac1XRF0->Fill(PpacX[1],fRFcal[0]);
	    }
          }
          if (fRF[1] > 0.)
          {
	    //if (1 ){
	    if (PpacIsHit[0] ){
              hPpac0XRF1->Fill(PpacX[0],fRFcal[1]);
	    }
	    //if (1 ){ // goatface 03 Sep 2012 19:21:14 
	    if (PpacIsHit[1] ){
              hPpac1XRF1->Fill(PpacX[1],fRFcal[1]);
	    }
          }
  fTargetX->SetParameters(PpacX[0],PpacX[1],PpacSepZ);
  fTargetY->SetParameters(PpacY[0],PpacY[1],PpacSepZ);
  TargetX=-100.;
  TargetY=-100.;
  if (PpacIsHit[0] && PpacIsHit[1]){
              // need a check to never allow fTarget < 0 ... this is extrapolation not interpolation!
	      //TargetX=fTargetX->Eval(452.+388.4);
	      //TargetY=fTargetY->Eval(452.+388.4); // PPACa center to window
	      TargetX=fTargetX->Eval(452);
	      TargetY=fTargetY->Eval(452); // PPACa center to window
    
  }
	  // Turn off all the gates for a new event
          gate_30s=false;
	  gate_29p=false;
          gate_17f=false;
          goodRF_30s=false;
	  goodPpacX_30s=false;
          goodRF_29p=false;
	  goodPpacX_29p=false;
          goodRF_17f=false;
	  goodPpacX_17f=false;
	  // Set the gates for 30S/29P RI beams production run
	  // 30S done again w/ new PPAC calib 20 Sep 2011 17:22:52 
	  // 30S redone with RFcal 16 Jan 2012 17:32:44 
	  if (PpacIsHit[0] && PpacIsHit[1]){
	    if (  ((fRFcal[0]>=20) && (fRFcal[0]<=30)) || ((fRFcal[1]>=53)&&(fRFcal[1]<=60)) || (fRFcal[1]<=3) ) 
	      {goodRF_30s=true;} //30S gate
            if ((PpacX[0]>=-9 && PpacX[0]<=8) && (PpacX[1]>=-10 && PpacX[1]<=12)) {goodPpacX_30s=true;} //30S gate
            if (goodRF_30s && goodPpacX_30s) {gate_30s=true;}
            // PPAC X was redone 20 Sep 2011 17:23:07 
	    // 29P redone with RFcal 16 Jan 2012 17:35:10 
	    if (  ((fRFcal[0]>=6) && (fRFcal[0]<=15)) || ((fRFcal[1]>=39)&&(fRFcal[1]<=48)) ) 
	       {goodRF_29p=true;} //29P gate
            if ((PpacX[0]>=-4 && PpacX[0]<=9) && (PpacX[1]>=-5 && PpacX[1]<=14)) {goodPpacX_29p=true;} //29P gate
            if (goodRF_29p && goodPpacX_29p) {gate_29p=true;}
            // mystery ion...it should be some light ion at 0deg SSD in coincidence with 30S or 29P...check run 1027
	    // mystery ion RFcal 16 Jan 2012 17:42:46
	    // PPAC X redone 16 Jan 2012 17:45:43 
	    if (((fRFcal[0]>=34) && (fRFcal[0]<=43)) || ((fRFcal[1]>=6)&&(fRFcal[1]<=15))) {goodRF_17f=true;} 
            if ((PpacX[0]>=-9 && PpacX[0]<=8) && (PpacX[1]>=-11 && PpacX[1]<=2)) {goodPpacX_17f=true;} 
            if (goodRF_17f && goodPpacX_17f) {gate_17f=true;}
          }
/*	  
	  // 30S energy loss case (WF and F3 slits changed slightly) 24 Jan 2012 13:29:18  
      
	  if (PpacIsHit[0] && PpacIsHit[1]){
	    if (  ((fRFcal[0]>=19) && (fRFcal[0]<=29)) || ((fRFcal[1]>=52)&&(fRFcal[1]<=60)) || (fRFcal[1]<=2) ) 
	       {goodRF_30s=true;} //30S gate
            if ((PpacX[0]>=-9 && PpacX[0]<=13) && (PpacX[1]>=-7 && PpacX[1]<=17)) {goodPpacX_30s=true;} //30S gate
            if (goodRF_30s && goodPpacX_30s) {gate_30s=true;}
            //if (goodRF_30s && goodPpacX_30s && (WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)))) {gate_30s=true;}
            // PPAC X was redone 20 Sep 2011 17:23:07 
	    // 29P redone with RFcal 16 Jan 2012 17:35:10 
	    if (  ((fRFcal[0]>=3) && (fRFcal[0]<=14)) || ((fRFcal[1]>=36)&&(fRFcal[1]<=47)) ) 
	       {goodRF_29p=true;} //29P gate
            if ((PpacX[0]>=-5 && PpacX[0]<=14) && (PpacX[1]>=-5 && PpacX[1]<=19)) {goodPpacX_29p=true;} //29P gate
            if (goodRF_29p && goodPpacX_29p) {gate_29p=true;}
            //if (goodRF_29p && goodPpacX_29p && (WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)))) {gate_29p=true;}
	    //if (gate_30s) for (int ii=200; ii<800; ii=ii+10){
	   //   hPpacProj_3D->Fill(fTargetX->Eval(ii),fTargetY->Eval(ii),ii);
	  //  }
	  }
*/
          // mystery ion...it should be some light ion at 0deg SSD in coincidence with 30S or 29P...check run 1027
	  // mystery ion RFcal 16 Jan 2012 17:42:46
	  // PPAC X redone 16 Jan 2012 17:45:43 
//	  if (((fRFcal[0]>=34) && (fRFcal[0]<=43)) || ((fRFcal[1]>=6)&&(fRFcal[1]<=15))) {goodRF_17f=true;} 
//          if ((PpacX[0]>=-9 && PpacX[0]<=8) && (PpacX[1]>=-11 && PpacX[1]<=2)) {goodPpacX_17f=true;} 
//          if (goodRF_17f && goodPpacX_17f) {gate_17f=true;}
          
          if (gate_30s) {
	    hPpacToF_30s->Fill(Tof);
	      hPpac0XRF0_30s->Fill(PpacX[0],fRFcal[0]);
	      if (!SsdOR) hPpac0XRF0_30s_ds->Fill(PpacX[0],fRFcal[0]);
              hPpac1XRF0_30s->Fill(PpacX[1],fRFcal[0]);
	      hPpac0XRF1_30s->Fill(PpacX[0],fRFcal[1]);
              hPpac1XRF1_30s->Fill(PpacX[1],fRFcal[1]);
	      //hTargetXY_30s->Fill(TargetX,TargetY); // PPACa center to window
	      if (WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)))  hTargetXY_30s->Fill(TargetX,TargetY); // PPACa center to window
	      //if (WindowCut(TargetX,TargetY))  hTargetXYcut_30s->Fill(TargetX,TargetY); // PPACa center to window
	  }
          if (gate_29p) {
	    hPpacToF_29p->Fill(Tof);
            //if ((fSiTcal[0]-tof[0])>295.5 && (fSiTcal[0]-tof[0])<303) { // 29P careful, can crash production runs, used for energyloss
            //if ((fSiTcal[0]-tof[0])>290 && (fSiTcal[0]-tof[0])<295.5) { // 30S careful, can crash production runs, used for energyloss
	      hPpac0XRF0_29p->Fill(PpacX[0],fRFcal[0]);
              hPpac1XRF0_29p->Fill(PpacX[1],fRFcal[0]);
	      hPpac0XRF1_29p->Fill(PpacX[0],fRFcal[1]);
              hPpac1XRF1_29p->Fill(PpacX[1],fRFcal[1]);
	    //}
	      //hTargetXY_29p->Fill(TargetX,TargetY); // PPACa center to window
	      if (WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)))  hTargetXY_29p->Fill(TargetX,TargetY); // PPACa center to window
	      //if (WindowCut(TargetX,TargetY))  hTargetXY_29p->Fill(TargetX,TargetY); // PPACa center to window
	      //hTargetXY_29p->Fill(fTargetX->Eval(654.4),fTargetY->Eval(654.4));
	      //hTargetXY_29p->Fill(fTargetX->Eval(452),fTargetY->Eval(452));
	  }
	  if (gate_17f){
	      hPpac0XRF0_17f->Fill(PpacX[0],fRFcal[0]);
              hPpac0XRF1_17f->Fill(PpacX[0],fRFcal[1]);
	  }
	  if (gate_30s){
	        hPpac0XRF1cut->Fill(PpacX[0],fRFcal[1]);
	        hPpac1XRF1cut->Fill(PpacX[1],fRFcal[1]);
	  }
	  //if (PpacIsHit[0] && PpacIsHit[1]) hTargetXY->Fill(fTargetX->Eval(654.4),fTargetY->Eval(654.4));
	  //if (PpacIsHit[0] && PpacIsHit[1]) {
	  //if (PpacIsHit[0] && PpacIsHit[1] && gate_30s) {
	  if(0){
	  //if (PpacIsHit[0] && PpacIsHit[1] && ribeam_ssd && gate_30s) {
	    hTargetXY->Fill(TargetX,TargetY);
	    //if (WindowCut(TargetX,TargetY))  hTargetXYcut->Fill(TargetX,TargetY); // PPACa center to window

	  }

        } // end if: flag_detail
      } // end if: flag_ppac
      //SSD Pad
      if (flag_ssd) { // process SSD data?
        //SiIsHit was defined here 28 Nov 2011 12:13:51 
	if (flag_raw) { // process raw SSD data?
          for (UShort_t i=0; i<18; i++) {
            if (SiIsHit[i]) {
              hSiPadHit->Fill(i);
              hSiPadE->Fill(i,fSiE[i][0]);
              //if (SsdOR && (i==0||i==1||i==2||i==4||i==11)) // checking consistency with scalars for SSD-OR
	      //if (fSiE[i][0]>200)
	      hSiT->Fill(i,fSiT[i]);
              //cout << "fSiT[" << i << "]:" << fSiT[i] << endl;
	      padE_ch[i]->Fill(fSiE[i][0]);
	      //if (fSiE[0][4]>200 || fSiE[0][5]>200 || fSiE[0][3]>200) padE_ch[i]->Fill(fSiE[i][0]);
	      padT_ch[i]->Fill(fSiT[i]);
            } // end if: SiIsHit
          } // end for: i
        } // end if: flag_raw
        
        if (flag_detail){
          proton_beam=false;
          //if  (( (fSiE[1][1] > 890.) && (fSiE[1][1]<950.) ) || ( (fSiE[1][2] > 920.) && (fSiE[1][2]<1000.)  ) || ((fSiE[1][3] > 920.) && (fSiE[1][3] < 1000.))){// proton beam (I assume?)
          if  (( (fSiE[1][2] > 920.) && (fSiE[1][2]<1000.) ) ){// proton beam (I assume?)
          //if  (( fSiE[1][1] > 890. && fSiE[1][1] < 950. ) || ( fSiE[1][2] > 920. && fSiE[1][2] < 1000. ) || ( fSiE[1][3] > 920. && fSiE < 980.)){// proton beam (I assume?)
            proton_beam=true;
          }
          //all data in fSiE[18][9],fSiT[18]
          for (UShort_t i=0; i<18; i++) {
            fSiEcal[i][0]=-100.;
            fSiTcal[i]=-100.;
	    if (SiIsHit[i] == true ){
              //fSiTcal[i]=0.09765625*fSiT[i]; //calibration from ch to ns
            //if (SiIsHit[i] == true && proton_beam == true){
              //if (i==1)
              fSiTcal[i]=ssd_padT_calib->Calib(i,fSiT[i]);
              hSiTcal->Fill(i,fSiTcal[i]);
	      padT_cal_ch[i]->Fill(fSiTcal[i]);
              //hSiT->Fill(1,fSiT[1]);
              //Fill individual ssd data 
              if (fSiE[i][0] > 0 ){ // don't fill the histogram with junk
                fSiEcal[i][0]=ssd_padE_calib->CalibMeV(i,fSiE[i][0]);
                //pad_ch[i]->Fill(fSiE[i][0]);
                // BUG gate_30s is false up to here!!
                if (gate_30s) { // this will only do something useful if ppac_detail == 1
                   pad_ch_gated[i]->Fill(fSiE[i][0]);
                }
                hSiPadEcal->Fill(i,fSiEcal[i][0]);	
              } // end if: fSiE[i][0] > 0
            } // end if: SiIsHit

          	// this is all garbage: clean it out
                // anyway it needs an ssd_strip_detail check, otherwise can seg fault
		/*
		if (i<6) {//dE-E... which is which telescope?
                  if ((fRF[0] > 240 && fRF[0] < 330)||(fRF[0] > 470 && fRF[0] < 540 )){
                         dE_EBeamGateRF0[i]->Fill(fSiE[i][0],fSiE[i+6][0]);
                         dE_EBeamGateRFAll->Fill(fSiE[i][0],fSiE[i+6][0]);
                      }
                  if ((fRF[1] > 240 && fRF[1] < 340)||(fRF[1]>480)||(fRF[1]<130)){
                         dE_EBeamGateRF1[i]->Fill(fSiE[i][0],fSiE[i+6][0]);
                         dE_EBeamGateRFAll->Fill(fSiE[i][0],fSiE[i+6][0]);
                      }
                }*/
                 
          } // end for: i<18
/*
            for (Int_t i=0; i<6; i++){
                	if (fSiEmax[i] > 0. && fSiEmax[i+6] > 0. && fSiEcal[i+12][0] == 0.){
            	dE_Estrip[i]->Fill(fSiEmax[i],(fSiEmax[i]+fSiEmax[i+6]));}
                	else if(fSiEmax[i] > 0. && fSiEmax[i+6] > 0. && fSiEcal[i+12][0] > 0.){
            	dE_Estrip[i]->Fill(fSiEmax[i],(fSiEmax[i]+fSiEmax[i+6]+fSiEcal[i+12][0]));	}
                  if (fSiEcal[i][0] > 0. && fSiEcal[i+6][0] > 0. && fSiEcal[i+12][0]==0.){
          	dE_E[i]->Fill(fSiEcal[i][0],(fSiEcal[i][0]+fSiEcal[i+6][0]));}
                  else if (fSiEcal[i][0] > 0. && fSiEcal[i+6][0] > 0. && fSiEcal[i+12][0] > 0.){
          	dE_E[i]->Fill(fSiEcal[i][0],(fSiEcal[i][0]+fSiEcal[i+6][0]+fSiEcal[i+12][0]));}
            
            	}
  */    
        // SsdOR was set true here 28 Nov 2011 12:14:07 
	}//end if: flag_detail
      } // end if: ssd
     
     if (flag_strip){
        // SSD Strip
        ribeam_ssd=false;
        if (flag_raw) { // process raw SSD strip data?
          for(UShort_t i=0;i<12;i++){
            for (UShort_t j=1; j<9; j++) {
              // here goatface 28 Aug 2012 11:43:43 
	      if (fSiE[i][j]>4000) fSiE[i][j]=-90.;
	      if (i==0 && (j==7 || j==8)){
	        if (fSiE[i][j]< 3000 && fSiE[i][j]>2000) {
		  ribeam_ssd=true;
		}
	      }
	      hSiStripE->Fill(i*8+j,fSiE[i][j]);
	      //if (i==0 && (j==7 || j==8) && fSiE[0][7]>500 && fSiE[0][8]>500) continue;
	      strip_ch[i*8+j]->Fill(fSiE[i][j]); 
            } // end for: j
	  } // end for: i
	  //if (fSiE[0][8]>500 && fSiE[0][7]>500) cout << "." << endl;
	} // end if: ssd_strip_raw
	if (flag_detail){
          for(UShort_t i=0;i<12;i++){
            fSiEmax[i]=0.;
            for (UShort_t j=1; j<9; j++) {
              //if (SsdOR && gate_30s)
	      fSiEcal[i][j]=ssd_strip_calib->CalibMeV(i*8+j,fSiE[i][j]);
              if (fSiE[i][j]>fSiEmax[i] && fSiE[i][j]<4000.) {//find maximum strip
              //if (fSiEcal[i][j]>fSiEmax[i] && fSiE[i][j]<4000.) {//find maximum strip
                //fSiEmax[i]=fSiEcal[i][j];
                fSiEmax[i]=fSiE[i][j];
		fSiEStripmax[i]=fSiE[i][j]; // goatface 29 Aug 2012 16:44:27 
		fSiEStripmax_ch[i]=j; // goatface 29 Aug 2012 16:44:27 
              }
              //Fill each strip data 
              hSiStripEcal->Fill(i*8+j,fSiEcal[i][j]);
              //if ( fSiE[i][j] > 0 && fSiE[i][j] < 4000 && proton_beam){ // just for high and low pedestals (high part is too tight)
              //if ( fSiE[i][j] > 0 && fSiE[i][j] < 4097){ // just for high and low pedestals (high part is too tight)
              //if ( (SiIsHit[i]) && (fSiE[i][j] > -100) ){ // just for high and low pedestals (high part is too tight)
              //if ( (fSiE[i][j] > -100 && !SsdOR) ){
	      //if ( (fSiE[i][j] > -100) && fSiE[i][j]<4000 && gate_29p && SsdOR ){
              if ( (fSiE[i][j] > -100 ) ){
                strip_cal_ch[i*8+j]->Fill(fSiEcal[i][j]);
              }
            } //end for: j
            hSiEmax->Fill(i,fSiEmax[i]);
          } //end for: i
        } // end if: flag_detail 
      } // end if: flag_strip

      if (flag_detail && flag_ssd && flag_strip){
         
/*
        // telescopes at zero degrees have cabling switched somehow...
        //if (fSiEmax[0] > 0. && fSiEmax[7] > 0. && SiIsHit[0]) // 1X + 2Y 
        if (fSiEmax[0] > 0. && fSiEmax[7] > 0.) // 1X + 2Y 
          dE_Estrip[0]->Fill(fSiEmax[0],(fSiEmax[0]+fSiEmax[7]+fSiEcal[12][0]));
        //if (fSiEmax[1] > 0. && fSiEmax[6] > 0. && SiIsHit[1]) // 2X + 1Y + 1c
        if (fSiEmax[1] > 0. && fSiEmax[6] > 0.) // 2X + 1Y + 1c
          dE_Estrip[1]->Fill(fSiEmax[1],(fSiEmax[1]+fSiEmax[6]+fSiEcal[13][0]));
        //if (fSiE[2][0] > 0. && fSiE[8][0] > 0. && SiIsHit[2])
        if (fSiE[2][0] > 0. && fSiE[8][0] > 0.)
        //if (fSiEcal[2][0] > 0. && fSiEcal[8][0] > 0.) // causes strange problems...
          dE_Estrip[2]->Fill(fSiEcal[2][0],(fSiEcal[2][0]+fSiEcal[8][0]));
        
        //if (fSiEmax[3] > 0. && fSiEmax[9] > 0. && SiIsHit[3] )
        if (fSiEmax[3] > 0. && fSiEmax[9] > 0. )
          dE_Estrip[3]->Fill(fSiEmax[3],(fSiEmax[3]+fSiEmax[9]));
  */    


        if (fSiEmax[0] > 0. && fSiEmax[7] > 0. ) // 1X + 2Y 
          dE_Estrip[0]->Fill(fSiEmax[0],(fSiEmax[0]+fSiEmax[7]+fSiEcal[12][0]));
        if (fSiEmax[1] > 0. && fSiEmax[6] > 0. ) // 2X + 1Y + 1c
          dE_Estrip[1]->Fill(fSiEmax[1],(fSiEmax[1]+fSiEmax[6]+fSiEcal[13][0]));
        if (fSiEcal[2][0] > 0. && fSiEcal[8][0] > 0.) 
          dE_Estrip[2]->Fill(fSiEcal[2][0],(fSiEcal[2][0]+fSiEcal[8][0]+fSiEcal[14][0]));
        if (fSiEmax[3] > 0. && fSiEmax[9] > 0.) 
          dE_Estrip[3]->Fill(fSiEmax[3],(fSiEmax[3]+fSiEmax[9]));
        for (int i=0;i<8;i++){
          if (SiIsHit[0] && fSiE[0][i+1]>180 && fSiE[0][0]>-90)
	  padE_strip[i]->Fill(fSiE[0][0],fSiE[0][i+1]);
	}
      } // end if: flag_detail & flag_ssd & flag_strip

      // PPAC and SSD
      if (flag_detail && flag_ppac){
	if (SsdOR){ // physics event
	  if (fRF[0] > 0. ) {
	    if (PpacIsHit[0] ){
	      hPpac0XRF0ssd->Fill(PpacX[0],fRFcal[0]);
	    }
	    if (PpacIsHit[1] ){
              hPpac1XRF0ssd->Fill(PpacX[1],fRFcal[0]);
              hPpac1XYcut->Fill(PpacX[1],PpacY[1]);
	    }
	  }
	  if (fRF[1] > 0. ) {
	    if (PpacIsHit[0] ){
	      hPpac0XRF1ssd->Fill(PpacX[0],fRFcal[1]);
            }
	    if (PpacIsHit[1] ){
	      hPpac1XRF1ssd->Fill(PpacX[1],fRFcal[1]);
	    }
	  }
	}
	else{ // downscale event
	  if (fRF[0] > 0.){
	    if (PpacIsHit[0] ){
	      hPpac0XRF0ds->Fill(PpacX[0],fRFcal[0]);
	    }
	    if (PpacIsHit[1] ){
              hPpac1XRF0ds->Fill(PpacX[1],fRFcal[0]);
	    }
	    if (gate_30s){
	    // can fill here
	    }
	  }
	  if (fRF[1] > 0.){
	    if (PpacIsHit[0] ){
	      hPpac0XRF1ds->Fill(PpacX[0],fRFcal[1]);
	    }
	    if (PpacIsHit[1] ){
              hPpac1XRF1ds->Fill(PpacX[1],fRFcal[1]);
	    }
	    if (gate_30s){
	      //hPpac1XRF1cut->Fill(PpacX[1],fRF[1]);
	    }
	  }
	}

      } // end if : flag_detail & flag_ppac

      if (flag_detail && flag_ssd && flag_ppac) { // PPAC and SSD pad detailed analysis?
	for(UShort_t i=0;i<18;i++){
	    if (SiIsHit[i]){
	      if (gate_30s)  {
	        //hPpac0TpadT_30s_ch[i]->Fill(fSiTcal[i]-tof[0]);
	        hPpac1TpadT_30s_ch[i]->Fill(fSiTcal[i]-tof[1]);
              }
	      if (gate_29p) {
	        //hPpac0TpadT_29p_ch[i]->Fill(fSiTcal[i]-tof[0]);
	        hPpac1TpadT_29p_ch[i]->Fill(fSiTcal[i]-tof[1]);
              }
	    }
	}
      } // end if: flag_detail & flag_ssd & flag_ppac
      if (flag_detail && flag_ssd && flag_strip && flag_ppac) { // PPAC and SSD strip detailed analysis?
        for (UShort_t i=0;i<1;i++){
            //if ((fSiTcal[0]-tof[0])>296 && (fSiTcal[0]-tof[0])<300.5 && fSiEStripmax[i]>2000) { // carful, can crash production runs, used for energyloss
            //if (gate_29p) {
            //if (gate_29p && fSiEStripmax[i]>2000) {
//	    if (PpacIsHit[0] && PpacIsHit[1]){
	    if (WindowCut (fTargetX->Eval(452.),fTargetY->Eval(452.))&& PpacIsHit[0] && PpacIsHit[1]){
	        if (fSiE[0][7]>1000) hPpac0TpadT_29p_ch[i]->Fill(tof[0]);
	        //if (fSiE[0][7]>1000) hPpac0TpadT_29p_ch[i]->Fill(fSiTcal[i]-tof[0]);
	       strip_ch_gated29p[i*8+fSiEStripmax_ch[i]]->Fill(fSiEStripmax[i]);
	       if (PpacIsHit[0] && PpacIsHit[1] && fSiE[0][7]>1000 && gate_29p && fSiE[0][8]<1000) (hTargetXY_29p->Fill(TargetX,TargetY)); // PPACa center to window
	       //if (gate_29p && fSiE[0][7]>1000 && fSiEStripmax_ch[i]==7) (hTargetXY_29p->Fill(TargetX,TargetY)); // PPACa center to window
	      //if (WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)))  hTargetXYcut->Fill(TargetX,TargetY); // PPACa center to window
	      //if (gate_29p && fSiE[0][7]>1000 && fSiEStripmax_ch[i]==7 && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)))  hTargetXYcut->Fill(TargetX,TargetY); // PPACa center to window
	      if (gate_29p && fSiE[0][8]>1000 && fSiE[0][7]<1000)  hTargetXYcut->Fill(TargetX,TargetY); // PPACa center to window
            //}
            //if (gate_30s) {
            //if (gate_30s && fSiEStripmax[i]>2000) {
            //if ((fSiTcal[0]-tof[0])>292 && (fSiTcal[0]-tof[0])<296 && fSiEStripmax[i]>2000 && gate_30s) {
	       strip_ch_gated30s[i*8+fSiEStripmax_ch[i]]->Fill(fSiEStripmax[i]);
	        if (fSiE[0][8]>1000)hPpac0TpadT_30s_ch[i]->Fill(tof[0]);
	        //if (fSiE[0][8]>1000)hPpac0TpadT_30s_ch[i]->Fill(fSiTcal[i]-tof[0]);
	      
	        if (fSiE[0][8]>1000 && fSiE[0][7]<1000 && PpacIsHit[0] && PpacIsHit[1] && gate_30s) hTargetXY_30s->Fill(TargetX,TargetY); // PPACa center to window
	        //if (gate_30s && fSiE[0][8]>1000 && fSiEStripmax_ch[i]==8) hTargetXY_30s->Fill(TargetX,TargetY); // PPACa center to window
	        //if (gate_30s && fSiE[0][8]>1000 && fSiEStripmax_ch[i]==8 && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.))) hTargetXYcut_30s->Fill(TargetX,TargetY); // PPACa center to window
	        if (gate_30s && fSiE[0][8]<1000 && fSiE[0][7]>1000 && fSiEStripmax_ch[i]==8) hTargetXYcut_30s->Fill(TargetX,TargetY); // PPACa center to window
	    }
         /* 
	  for (UShort_t j=1;j<9;j++){
            if (gate_29p) {
	       strip_ch_gated29p[i*8+j]->Fill(fSiE[i][j]);
            }
            if (gate_30s) {
               if (j==8 && fSiE[i][j]>2000 && fSiE[i][7] > 2000) cout << "Multihit!"  << endl;
               strip_ch_gated30s[i*8+j]->Fill(fSiE[i][j]);
            }
            if (gate_17f) {
               strip_ch_gated17f[i*8+j]->Fill(fSiE[i][j]);
            }
	  } // end for: j
	  */
	} // end for: i
      } // end if: flag_detail & flag_ssd & flag_strip & flag_ppac

      //GEM-MSTPC
      if (flag_tpc){  // analyze TPC data?
        
	TArtFadcHit *fadc;
        Int_t nadclines = arr->GetEntriesFast();
	//cout << "nadclines: " << nadclines << endl;
        // Init the FADC Max arrays
        if (flag_detail) { // process TPC detailed analysis? 
	  for(UShort_t j=0;j<144;j++){
	    PadIsHit[j]=false;
            fPulseNo[j]=0;  
            for(UShort_t k=0;k<20;k++){ // fPulseNo
	      baseline[j][k] = 0.;
	      baseline_dev[j][k] = 0.;
              MultiHit[j][k] = 1; // Multiple hit number 
              fSampleInt[j][k][0]=0.;
              fSampleInt[j][k][1]=0.;
              for(UShort_t l=0;l<20;l++){ // MultiHit
              fSampleMax[j][k][l] = 0.; 
              fClockMax[j][k][l] = 0.;
	      fTimeMax[j][k][l] = 0.;
              }    
            }    
          }
	  for (UShort_t j=0;j<1334;j++){
            fSampleRebin[j]=0.;
	  }
	} // end if: flag_detail

	for (Int_t i=0; i<nadclines; i++) { // loop over the number of ADC lines
	  fadc=(TArtFadcHit *) arr->At(i);
          UShort_t id=fadc->fID; // pad ID (range 0 to 143)
          if (id>143 || id<0) cout << "Error in Pad ID: " << id << endl;
	  fPulseNo[id]=fadc->fHit;  // Pulse number if there are events, never 0
          //if (fPulseNo[id]>1) cout << "fPulseNo[" << id << "]: " << fPulseNo[id] << endl;
	  if (flag_raw) {// process TPC raw data?
	   // if (SiIsHit[2] && ssd_mult==1){
	    for (Int_t j=0; j<fadc->fNSample; j++) {
              hGemMstpcAll->Fill(fadc->fClock[j], fadc->fSample[j]);
              hGemMstpcCh->Fill(id,fadc->fSample[j]);
              //if (SiIsHit[0] && fSiE[0][0]>170)
	      // playing with the baseline and what goes on with fHit for the HG Gems 14 Mar 2013 18:20:07 
	      //if (fadc->fSample[j]>550) hTpc_dE_ch[id]->Fill(fadc->fClock[j]+fadc->fTime, fadc->fSample[j]);
	      hTpc_dE_ch[id]->Fill(fadc->fClock[j], fadc->fSample[j]);
              hfSample[id]->Fill(fadc->fSample[j]);
              if (id < 96){//Main (low gain) beam pads
                   hBeamPadCh->Fill(id,fadc->fSample[j]); 
              }
	      else {//high gain pads
		 hGemMstpcSide->Fill(fadc->fClock[j], fadc->fSample[j]);
	      } // end if: id<96
	    } // end for: j fNSample
	 //   } // end if ssdishit
	  } // end if: flag_raw
	 
	 /*
	  // this is for viewing 500 sample tracks one-by-one for QA 
	  //if (id==122 && pulsetrack<500){
	  if (id==37 && pulsetrack<500){ //beam check
	  //if (id==115 && pulsetrack<500){
	    if ((fadc->fSample[0]!=0.) && (fadc->fSample[1]!=0.) && (fadc->fSample[2]!=0.) && (fadc->fSample[3]!=0.) && (fadc->fSample[4]!=0.)){
	      //if (flag_detail && flag_ppac) { // high gain left check
	      //if (flag_detail && SiIsHit[2]) { // high gain left check
	      if (flag_detail && flag_ppac && gate_30s) { // beam check
	        for (Int_t j=0; j<fadc->fNSample; j++) {
                  hTpcPulse[pulsetrack]->Fill(fadc->fClock[j], fadc->fSample[j]);
		  specialtrack=true;
	        }
	      pulsetrack++;
	      }
	    }
	  else continue;
	  }
	  */

          if (flag_detail) { // process TPC detailed analysis?
	    if (fPulseNo[id]!=0) PadIsHit[id] = true;
	    else continue; // junk; truncate
	    
            //baseline check 
	    for(UShort_t j=0; j<5;j++){ // get the first five pulses
              if (fadc->fSample[j]>0.) baseline[id][fPulseNo[id]] = baseline[id][fPulseNo[id]] + fadc->fSample[j]; 
            } // end for j: baseline determination
            
	    // make sure there are 5 pulses
	    if ((fadc->fSample[0]!=0.) && (fadc->fSample[1]!=0.) && (fadc->fSample[2]!=0.) && (fadc->fSample[3]!=0.) && (fadc->fSample[4]!=0.)){
	      baseline[id][fPulseNo[id]]=baseline[id][fPulseNo[id]]/5.;
	      // find maximum deviation of the first five samples from the baseline
	      Double_t maxdeviation=0.;
	      for (UShort_t j=0;j<5;j++){
	        if (TMath::Abs(fadc->fSample[j]-baseline[id][fPulseNo[id]])>maxdeviation) maxdeviation=TMath::Abs(fadc->fSample[j]-baseline[id][fPulseNo[id]]);
	      }
              baseline_dev[id][fPulseNo[id]]=(TMath::Abs(maxdeviation/baseline[id][fPulseNo[id]]));
	    }
	    else continue; // junk; truncate
	    if((fadc->fSample[0]) < 400 || fadc->fSample[0]>600) continue; // the pulse is junk or partly missed or buried
            
	    // daid's moving average rebinning
            if (fadc->fNSample > 30 ){  // make sure the pulse width is reasonable
              rebinCounter=0;
	      rebinTemp=0.;
	      for (UShort_t j=0; j<fadc->fNSample; j++) { // loop all the FADC samples
	        rebinTemp+=fadc->fSample[j];
		rebinCounter++;
		if (rebinCounter==rebin){
		  fSampleRebin[j/rebin]=(rebinTemp/rebin);
		  rebinCounter=0;
	          rebinTemp=0.;
		}
		if (rebinCounter>rebin || rebinCounter<0) cout << "fadc rebinCounter is broken!" << endl;
	      }
            /*
	    // original peak finder
            trough=1e6; //init trough
              for (UShort_t j=3; j<fadc->fNSample-3; j++) { // loop all the FADC samples
                //hGemMstpcAll->Fill(fadc->fClock[j], fadc->fSample[j]);
                if (fadc->fSample[j]<trough && fadc->fSample[j]>0.) trough=fadc->fSample[j];
                // this call is strange here! moved it up until I remember why it was put specifically here
                // daid 13 Apr 2012 18:09:30 
                // if((fadc->fSample[j]) < 400) continue; // 500 should be the baseline ?
                  if(fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] > fadc->fSample[j]) continue; // is the present fSample larger than fSampleMax?
                    if( !(  // be sure it is a true peak, by looking +3 and -3 samples
                        (fadc->fSample[j-3] < fadc->fSample[j-2]) &&
                        (fadc->fSample[j-2] < fadc->fSample[j-1]) &&
                        (fadc->fSample[j-1] < fadc->fSample[j])   &&
                        (fadc->fSample[j]   > fadc->fSample[j+1]) &&
                        (fadc->fSample[j+1] > fadc->fSample[j+2]) &&
                        (fadc->fSample[j+2] > fadc->fSample[j+3])
                        )) continue; // no peak -- skip
                      // true peak: fill the FADC Max arrays
                      if (TMath::Abs(trough-baseline[id][fPulseNo[id]])>(baseline_dev[id][fPulseNo[id]]*baseline[id][fPulseNo[id]])){
                        // trough and baseline don't agree!
                        hTrough[id]->Fill(trough-baseline[id][fPulseNo[id]]);
                      }
                      fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] = fadc->fSample[j]-baseline[id][fPulseNo[id]];
                      fClockMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] = fadc->fClock[j];
                      fTimeMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] = fadc->fTime; // we should modify fTime somehow for each peak?
                      trough=1e6; // reset trough
                      if (specialtrack) {
                        cout << "pulsetrack: " << pulsetrack << endl;
                        cout << "fSampleMax[" << id << "][" << fPulseNo[id] << "][" << MultiHit[id][fPulseNo[id]] << "]: " << fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] << endl;
                      }
                      MultiHit[id][fPulseNo[id]]++; // valid peak: increment the hit number
              } // end for j: FADC samples
              specialtrack=false;
 		*/


	    
	    // daid's moving average peak finder
              truepeak=false; // reset the true peak event
	      for (UShort_t j=2; j<((fadc->fNSample/rebin)-2);j++){ // loop all the FADC samples

	      //for (UShort_t j=2; j<((fadc->fNSample/rebin)-(rebin-1));j++){ // loop all the FADC samples
		if (rebinCounter>=rebin){ rebinCounter=0; } // probably we don't need this here... 09 Apr 2013 20:52:55 
	        if(fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] > fSampleRebin[j]) { continue; } // is the present fSample larger than fSampleMax?
                if( !(  // be sure it is a true peak, by looking +2 and -2 samples
                  (fSampleRebin[j-2] < fSampleRebin[j-1]) && 
                  (fSampleRebin[j-1] < fSampleRebin[j]  ) && 
		  (fSampleRebin[j] > fSampleRebin[j+1]  ) &&
		  (fSampleRebin[j+1] > fSampleRebin[j+2]) &&
		  (fSampleRebin[j]-baseline[id][fPulseNo[id]]>10)
	          )) { continue; } // no peak -- skip
	        // peak!
                // find the highest local peak within the rebinned region
		rebinCounter=0;
		for (int i=0;i<rebin;i++){
		  if( fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] < fadc->fSample[j*rebin+i]-baseline[id][fPulseNo[id]])
	          {
                    rebinCounter=i; // get the right rebin adjustment
		    fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] = fadc->fSample[j*rebin+i]-baseline[id][fPulseNo[id]];
		  }
		}
		// now we have the exact peak location
		fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] = fadc->fSample[j*rebin+rebinCounter]-baseline[id][fPulseNo[id]];
                fClockMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] = fadc->fClock[j*rebin+rebinCounter];
                fTimeMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] = fadc->fTime; // we should modify fTime somehow for each peak?
//	        trough=1e6; // reset trough
		 // output data on the 500 sample tracks one-by-one
		/*
		//if (specialtrack && id==122 ) {
		//if (specialtrack && id==115 ) {
		//if (0) {
		if (specialtrack) {
		  cout << "pulsetrack: " << pulsetrack << endl; 
		  cout << "fSampleMax[" << id << "][" << fPulseNo[id] << "][" << MultiHit[id][fPulseNo[id]] << "]: " << fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] << endl;
		  //cout << "fSampleRebin[" << j << "]: " << fSampleRebin[j] << endl;
		}
		*/

                MultiHit[id][fPulseNo[id]]++; // valid peak: increment the multiple hit number
		rebinCounter++; // don't need it?! 09 Apr 2013 20:54:02 
		truepeak=true;
              } // end for j: FADC samples
	      specialtrack=false;
            


	      // peak integrator
              // for now, we throw away events which have multiple hits both here and in later analysis
	      // that is because at present we lack the coding capability to deconvolute overlapping peaks
              if (truepeak && MultiHit[id][fPulseNo[id]]-1==1){
                for (UShort_t j=1; j<fadc->fNSample; j++) { // loop all the FADC samples
                if (fadc->fSample[j]<baseline[id][fPulseNo[id]] ) continue;
		  fSampleInt[id][fPulseNo[id]][1]=fSampleInt[id][fPulseNo[id]][1]+fadc->fSample[j]-baseline[id][fPulseNo[id]];
	        }
	      }
            



	    } //end if: fNSample > 30
          
	  /*
	  if (fPulseNo[id]>1) {
	    cout << "fSampleMax[" << id << "][1][1]: " << fSampleMax[id][1][1]<< endl;
	    cout << "fSampleMax[" << id << "][2][1]: " << fSampleMax[id][2][1]<< endl;
	  }
	  */

          } // end if: flag_detail
	} // end for i: nadclines
   
        if (flag_detail) { // process TPC detailed analysis?
	  //*  *****************************
          //   START GEM-MSTPC PHYSICS HERE!
          //*  *****************************
          // init arrays
	  XpadBcalib[0]=0.;
	  XpadBcalib[1]=0.;
	  for (UShort_t i=0;i<20;i++){
	    padsHitB[i]=0;
	    padsHitL[i]=0;
	    padsHitC[i]=0;
	  }
	  for (UShort_t i=0;i<48;i++){
            TracksB[i]=0; 
            for (UShort_t j=0;j<20;j++) {
              XpadB[i][j]=-1000.;
              YpadB[i][j]=0.;
              ZpadB[i][j]=0.;
              dEpadB[i][j]=0.;
              TpadB[i][j]=0.;
            }
          }
          for (UShort_t i=0;i<8;i++){
            TracksC[i]=0; 
            TracksL[i]=0; 
            TracksR[i]=0; 
            for (UShort_t j=0;j<20;j++) {
              XpadC[i][j]=0.;
              XpadR[i][j]=0.;
              XpadL[i][j]=0.;
              YpadC[i][j]=0.;
              YpadR[i][j]=0.;
              YpadL[i][j]=0.;
              ZpadC[i][j]=0.;
              ZpadR[i][j]=0.;
              ZpadL[i][j]=0.;
              dEpadC[i][j]=0.;
              dEpadR[i][j]=0.;
              dEpadL[i][j]=0.;
              TpadC[i][j]=0.;
              TpadR[i][j]=0.;
              TpadL[i][j]=0.;
            }
          }
          // BEAM
          for(UShort_t i=0;i<96;i++){ //beam_itop: // all right Beam pads
            for(UShort_t j=0;j<96;j++){ //beam_jtop: // all left Beam pads
              if((tpc_ch[0][i]==-1)||(tpc_ch[1][j]==-1) || // tpc_ch will get us the pad number (0 to 47)
              (tpc_ch[0][i] != tpc_ch[1][j]) || 
              (PadIsHit[i] == false ) || (PadIsHit[j] == false )) continue; // only run for the same pads; make sure both pads are hit
              pad[i] = tpc_ch[0][i]; // place holder for the loop
              for(UShort_t k=1;k<=fPulseNo[i];k++){ //beam_ktop: // loop on left side pulse number
	        //for(UShort_t l=1;l<2;l++){ // loop on left side hit number
	        for(UShort_t l=1;l<=MultiHit[i][k];l++){ // loop on left side hit number
                  for(UShort_t m=1;m<=fPulseNo[j];m++){ //beam_mtop: // loop on right side pulse number
                    //for(UShort_t n=2;n<2; n++){ // loop on right side hit number
                    for(UShort_t n=1;n<=MultiHit[j][m]; n++){ // loop on right side hit number
                      if (TracksB[pad[i]] >= 20 ) continue; // don't screw up memory; TracksB[i] cannot access more than 19
                      if(!(fSampleMax[i][k][l] > 0.1) || !(fSampleMax[j][m][n] > 0.1)  || // both samples are nonzero?
                         !(fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	              //if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
                      TpadB[pad[i]][TracksB[pad[i]]] = (fTimeMax[i][k][l]+fTimeMax[j][m][n])/2.;
	              // pad 6 seems to collect a lot of the data from pad 7, and thus shows roughly 2x number of real tracks
		      /* // this was from an array error 21 Mar 2013 20:09:53 
		      if (pad[i]==6 && TracksB[6]>0) 
	                for (UShort_t pad6extra=0;pad6extra<TracksB[6];pad6extra++){
	                  if (TMath::Abs(TpadB[6][TracksB[6]]-TpadB[6][pad6extra])<=1){
	                    TpadB[6][TracksB[6]]=0.;
			    if ((m+1)<=fPulseNo[j]){ m++; goto beam_mtop;}
	                    else if ((k+1)<=fPulseNo[i]){ k++; goto beam_ktop;}
	                    else if ((j+1)<96) { j++; goto beam_jtop;}
	                    else if ((i+1)<96) { i++; goto beam_itop;}
	                    else {goto beam_end;}
			  }
	                }*/
			// cm don't use
	              /*XpadB[pad[i]][TracksB[pad[i]]] = 5.5*((fSampleMax[i][k][l]-fSampleMax[j][m][n])/
	                (fSampleMax[i][k][l]+fSampleMax[j][m][n]));*/
	              // mm
		      // if runs are gain 1 and if channel is 1, or if channel is 6...need to do something different 
		      // calibrated X
		      /*
		      if (run < 1005 && (pad[i]==1 || pad[i]==6)) {
                        if (pad[i]==1){
		        XpadBcalib[0]=fSampleInt[i][k][l]*1.02296;
		        XpadBcalib[1]=fSampleInt[j][m][n]*0.97701;
			}
                        if (pad[i]==6){
		        XpadBcalib[0]=fSampleInt[i][k][l]*1.0956;
		        XpadBcalib[1]=fSampleInt[j][m][n]*0.90440;
			}
		      }
                      else {
		        XpadBcalib[0]=beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l]);
		        XpadBcalib[1]=beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n]);
		      }
		      */
		       //calibrated
		      /*
		      //if (fSampleInt[i][k][l]>100 && fSampleInt[j][m][n]>100)
		      if (XpadBcalib[0]!=0)
		      XpadB[pad[i]][TracksB[pad[i]]] = 55.*((XpadBcalib[0]-XpadBcalib[1])/
	                (fSampleInt[i][k][l]+fSampleInt[j][m][n]));
	              */
		      
			 //raw data
		      //integrated
		      XpadB[pad[i]][TracksB[pad[i]]] = 55.*((fSampleInt[i][k][l]-fSampleInt[j][m][n])/
	                (fSampleInt[i][k][l]+fSampleInt[j][m][n])); 
		      //max 
		      //XpadB[pad[i]][TracksB[pad[i]]] = 55.*((fSampleMax[i][k][l]-fSampleMax[j][m][n])/
	              //  (fSampleMax[i][k][l]+fSampleMax[j][m][n])); 
		   
		     /* // better above
		      XpadB[pad[i]][TracksB[pad[i]]] = 55.*((beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l])-beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n]))/
	                (fSampleInt[i][k][l]+fSampleInt[j][m][n])); */
		      YpadB[pad[i]][TracksB[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	                + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	              if (!SsdOR) YpadB[pad[i]][TracksB[pad[i]]]=YpadB[pad[i]][TracksB[pad[i]]]-39.;
	              ZpadB[pad[i]][TracksB[pad[i]]] = padsize*(pad[i]-23)-0.21;
                      //dEpadB[pad[i]][TracksB[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]); // peak 
                      
		      //if (fSampleInt[i][k][l]+fSampleInt[j][m][n]>10) // integrated
		      // calibrated
		      //dEpadB[pad[i]][TracksB[pad[i]]] = (XpadBcalib[0]+XpadBcalib[1]); // integrated
		      //raw
		      //dEpadB[pad[i]][TracksB[pad[i]]] = (fSampleInt[i][k][l]+fSampleInt[j][m][n]); // integrated
		      dEpadB[pad[i]][TracksB[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]); // max
                      TracksB[pad[i]]++; // increase the number of tracks for that pad
		      padsHitB[TracksB[pad[i]]]++;
                    } // end for n: right side hit number
	          } // end m
	        } // end l
	      } // end k
            } // end for j: loop over left Beam pads 
          } // end for i: loop over left Beam pads 
	  //beam_end:

	  // Center
          for(UShort_t i=120;i<128;i++) // all right Center pads
            for(UShort_t j=96;j<104;j++){ // all left Center pads
              if((tpc_ch[0][i]==-1)||(tpc_ch[1][j]==-1) || // tpc_ch will get us the pad number (0 to 7)
              (tpc_ch[0][i] != tpc_ch[1][j]) || 
              (PadIsHit[i] == false ) || (PadIsHit[j] == false )) continue; // only run for the same pads; make sure both pads are hit
                pad[i] = tpc_ch[0][i]; // place holder for the loop
                //if (j==99) cout << pad[i] << endl;
		//if (j==103) pad[i]=0;
                //if (j==101) pad[i]=2;
		for(UShort_t k=1;k<=fPulseNo[i];k++) // loop on left side pulse number
	  	//for(UShort_t l=1;l<2;l++) // loop on left side hit number
	  	for(UShort_t l=1;l<=MultiHit[i][k];l++) // loop on left side hit number
                    for(UShort_t m=1;m<= fPulseNo[j];m++) // loop on right side pulse number
                      //for(UShort_t n=1;n<2; n++){ // loop on right side hit number
                      for(UShort_t n=1;n<=MultiHit[j][m]; n++){ // loop on right side hit number
                        if (TracksC[pad[i]] >= 20 ) continue; // don't screw up memory; TracksC[i] cannot access more than 19
	  	        if(!(fSampleMax[i][k][l] > 0.1) || !(fSampleMax[j][m][n] > 0.1)  || // both samples are nonzero?
                           !(fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	  		// zero suppression threshold is often too low, making the clock values wrong!  14 Mar 2013 18:24:14 
			//if(TMath::Abs((fClockMax[i][k][l])-(fClockMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
			if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
			// this seems to work for the original ridf2root....but it's throwing events away?
			//if (fSampleMax[i][k][l] < 10 || fSampleMax[j][m][n]<10.) continue;
			//if (fSampleMax[i][k][l]+fSampleMax[j][m][n]<70.) continue;
	  		if (pad[i]==0) continue; // for X checking, this pad data is way to crazy 11 Feb 2013 22:35:59 
			// X by peak
			//if (j==101) cout << "time left: " << fTimeMax[j][m][n]+fClockMax[j][m][n]  << " time right: " << fTimeMax[i][k][l]+fClockMax[i][k][l] << endl;
			//if (j==97 || j==100) cout << "time " << j << " : "  << fTimeMax[j][m][n] << endl;
			XpadC[pad[i]][TracksC[pad[i]]] = 55.*((fSampleMax[i][k][l]-fSampleMax[j][m][n])/
	  		  (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
			//if (j==101) cout << l << " " << n << " time left: " << fTimeMax[j][m][n]+fClockMax[j][m][n]  << " time right: " << fTimeMax[i][k][l]+fClockMax[i][k][l] << 
			//"\t" << fTimeMax[j][m][n] << " " << fClockMax[j][m][n] << " " << fTimeMax[i][k][l] << " "  << fClockMax[i][k][l] << "\t" << XpadC[pad[i]][TracksC[pad[i]]] << " \t" << fSampleMax[j][m][n] << "\t" << fSampleMax[i][k][l] << endl;
		        //if (pad[i]==5) XpadC[5][TracksC[5]]-=5.; // checking how pad 5 alters projection 19 Feb 2013 14:50:20 
		        //if (pad[i]==0) XpadC[0][TracksC[0]]+=85.; // checking how pad 0 inclusion alters projection 19 Feb 2013 14:50:20 
			// X by integration
			//XpadC[pad[i]][TracksC[pad[i]]] = 55.*((fSampleInt[i][k][l]-fSampleInt[j][m][n])/
	                //  (fSampleInt[i][k][l]+fSampleInt[j][m][n])); 
                        YpadC[pad[i]][TracksC[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	  		  + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	  		ZpadC[pad[i]][TracksC[pad[i]]] = 9.575 + padsize*(pad[i]); 
			//dEpadC[pad[i]][TracksC[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
			// quick and dirty calibration for CNS interview Feb 2013
			//if (pad[i]==7 || pad[i]==8) continue;
			//double center_gain[8]={3.49,1.13,1,1,1,1.04,0,0};
			//dEpadC[pad[i]][TracksC[pad[i]]] = center_gain[pad[i]]*(fSampleMax[i][k][l]+fSampleMax[j][m][n]);
			//calib
			if (fSampleMax[i][k][l]+fSampleMax[j][m][n]>40. ) //goatface 28 Feb 2013 20:20:31 emis
			//dEpadC[pad[i]][TracksC[pad[i]]] = hg_center_calib->Calib(pad[i],(fSampleMax[i][k][l]+fSampleMax[j][m][n]));
			//raw
			dEpadC[pad[i]][TracksC[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
                        TracksC[pad[i]]++; // increase the number of tracks for that pad
			padsHitC[TracksC[pad[i]]]++;
                      } // end for n: right side hit number
            } // end for j: loop over left Center pads 
          
	  // Left -- note: we are always in tpc_ch[1][x] because this is the left side!
          for(UShort_t i=112;i<120;i++) // all upstream Left pads
            for(UShort_t j=104;j<112;j++){ // all downstream Left pads
              if((tpc_ch[1][i]==-1)||(tpc_ch[1][j]==-1) || // tpc_ch will get us the pad number (0 to 7)
              (tpc_ch[1][i] != tpc_ch[1][j]) || 
              (PadIsHit[i] == false ) || (PadIsHit[j] == false )) continue; // only run for the same pads; make sure both pads are hit
                pad[i] = tpc_ch[1][i]; // place holder for the loop
                //if (i==118) pad[i]=1;
		for(UShort_t k=1;k<=fPulseNo[i];k++) // loop on left side pulse number
	  	//for(UShort_t l=1;l<2;l++) // loop on left side hit number ONCE
	  	for(UShort_t l=1;l<=MultiHit[i][k];l++) // loop on left side hit number
                    for(UShort_t m=1;m<= fPulseNo[j];m++) // loop on right side pulse number
                      //for(UShort_t n=1;n<2; n++){ // loop on right side hit number ONCE
                      for(UShort_t n=1;n<=MultiHit[j][m]; n++){ // loop on right side hit number
                        if (TracksL[pad[i]] >= 20 ) continue; // don't screw up memory; TracksL[i] cannot access more than 19
	  	        if(!(fSampleMax[i][k][l] > 0.1) || !(fSampleMax[j][m][n] > 0.1)  || // both samples are nonzero?
                           (!fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	  		if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
			if (fSampleMax[i][k][l]+fSampleMax[j][m][n]<80.) continue;
	  		// what's the units of Z?
			ZpadL[pad[i]][TracksL[pad[i]]] = 117.5*((fSampleMax[j][m][n]-fSampleMax[i][k][l])/
	  		  (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
                        YpadL[pad[i]][TracksL[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	  		  + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	  		XpadL[pad[i]][TracksL[pad[i]]] = -7.621 - padsize*(pad[i]); 
                        dEpadL[pad[i]][TracksL[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
                        TracksL[pad[i]]++; // increase the number of tracks for that pad
			padsHitL[TracksL[pad[i]]]++;
                      } // end for n: right side hit number
            } // end for j: loop over downstream Left pads 
          
	  // Right
          for(UShort_t i=136;i<144;i++) // all upstream Right pads
            for(UShort_t j=128;j<136;j++){ // all downstream Right pads
              if((tpc_ch[0][i]==-1)||(tpc_ch[0][j]==-1) || // tpc_ch will get us the pad number (0 to 7)
              (tpc_ch[0][i] != tpc_ch[0][j]) || 
              (PadIsHit[i] == false ) || (PadIsHit[j] == false )) continue; // only run for the same pads; make sure both pads are hit
                pad[i] = tpc_ch[0][i]; // place holder for the loop
                for(UShort_t k=1;k<=fPulseNo[i];k++) // loop on left side pulse number
	  	for(UShort_t l=1;l<=MultiHit[i][k];l++) // loop on left side hit number
                    for(UShort_t m=1;m<= fPulseNo[j];m++) // loop on right side pulse number
                      for(UShort_t n=1;n<=MultiHit[j][m]; n++){ // loop on right side hit number
                        if (TracksR[pad[i]] >= 20 ) continue; // don't screw up memory; TracksR[i] cannot access more than 19
	  	        if(!(fSampleMax[i][k][l] > 0.1) || !(fSampleMax[j][m][n] > 0.1)  || // both samples are nonzero?
                           !(fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	  	        if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
	  	        ZpadR[pad[i]][TracksR[pad[i]]] = 117.5*((fSampleMax[i][k][l]-fSampleMax[j][m][n])/
	  	          (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
                        YpadR[pad[i]][TracksR[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	  	          + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	  	        XpadR[pad[i]][TracksR[pad[i]]] = 7.621 + padsize*(pad[i]); 
                        dEpadR[pad[i]][TracksR[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
                        TracksR[pad[i]]++; // increase the number of tracks for that pad
			padsHitR[TracksR[pad[i]]]++;
                    } // end for n: right side hit number
            } // end for j: loop over downstream Right pads 
            

	  //beam
			 //PpacX[1]+=0.5;
          for(UShort_t i=0;i<48;i++){ // loop over all pads (now by pad number)
                         //TargetX=fTargetX->Eval((452.+82.+(i*padsize)+2.1));
			 //hTargetXZ->Fill(TargetX,i);
	    for(UShort_t j=0;j<1;j++){ // loop for first track only
	      //for(UShort_t j=0;j<TracksB[i];j++){ // loop for each track
	      //if (i==4 && TracksB[4]>1) cout << "TpadB[4][" << j << "]: " << TpadB[4][j] << endl;
	      //if (i==6) cout << "TpadB[6][" << j << "]: " << TpadB[6][j] << endl;
	      //cout << "TracksB[" << i << "]: " << TracksB[i] << endl;
	      if (!SsdOR) hBraggB->Fill(i,dEpadB[i][j]);
	      //implement me in Analyzer_config! 
	      double shift;
	      double paddistance;
	      const double theta=.01980872946019476473; //measured one -- seems correct 17 Apr 2013 17:32:51 atan((10.5*0.5)/265)
	      //const double theta=.02083032003621708494; //atan((0.5*7.5)/180) 14 May 2013 15:01:24 
	      //Int_t rotpoint=133; // in mm, from beginning of pad, estimated from geometry needing ~150mm to shift by 6 in Y, -17 mm to front of pad from foot
	      //Int_t rotpoint=341; // = (6.8/tan(.01980872946019476473))-(0.5*4) for first pad
	      //Int_t rotpoint=330; // this comes from the geometry 14 May 2013 15:07:46
	      Int_t rotpoint = 286 ; //pure geometry
	      //Int_t rotpoint=281; // = (5.56/tan(.01980872946019476473))-(0.5*4) for first pad ... 5.56 comes from unshifted residual Y in'cept
	      //Int_t rotpoint=266; // in mm, from beginning of pad, estimated from geometry needing ~150mm to shift by 6 in Y, -17 mm to front of pad from foot
	      paddistance=((i+0.5)*padsize)+rotpoint;
	      shift=paddistance*TMath::Tan(theta);
              //shift=0.;
	     TargetX=TargetY=-999.;
	     if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && gate_30s && XpadB[i][j]>-999. && dEpadB[i][j]>50){
	     //if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && gate_29p){
                         //checking TargetX works as expected
			 //TargetX=PpacX[0]+(((452+83+(i*padsize)+(padsize/2))/PpacSepZ)*(PpacX[1]-PpacX[0]));
			 // testing retarded stuff: 23 Apr 2013 16:28:43 
			 //TargetX=fTargetX->Eval((452.+83.+(i*padsize)+(padsize/2)));
			 // correct one I think:
			 TargetX=fTargetX->Eval((452.+83.+(i*padsize)+(padsize/2)));
                         TargetY=fTargetY->Eval((452.+83.+(i*padsize)+padsize/2)); // PPACa center to window, inactive region, plus pad
			 //to test averages (pad 6 has a lot of data for no obvious reason) 
			 //if (i<24  && i!=6) hPadXB1geo->Fill(i,(TargetX-(XpadB[i][j]+shift+5.4)));
			 //if (i<24) hPadXB1geo->Fill(i,TargetX);
			 //hPadXBgeo->Fill(i,TargetX); // simple target projection by PPAC for reference
			 hPadXBgeo_30s->Fill(i,((XpadB[i][j]+shift)-TargetX)); // proper residual
			 //hPadXBgeo->Fill(i,((XpadB[i][j])-TargetX)); // proper residual w/o shifting
			 //hPadXBgeo->Fill(i,(TargetX-(XpadB[i][j]+shift)));
			 hPadXB_30s->Fill(i,(XpadB[i][j])-TargetX);
			 //hPadXB->Fill(i,(XpadB[i][j]+shift));
			 //if (i<24) hPadXB1geo->Fill(i,(TargetX-(XpadB[i][j]+shift+5.4)));
	                 hPadYB->Fill(i,YpadB[i][j]);
			 hPadYBgeo->Fill(i,TargetY);
			 
             }
	     TargetX=TargetY=-999.;
	     if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && gate_29p && XpadB[i][j]>-999. && dEpadB[i][j]>50){
			 TargetX=fTargetX->Eval((452.+83.+(i*padsize)+(padsize/2)));
                         TargetY=fTargetY->Eval((452.+83.+(i*padsize)+padsize/2)); // PPACa center to window, inactive region, plus pad
			 hPadXBgeo_29p->Fill(i,TargetX); // proper residual
			 hPadXB_29p->Fill(i,(XpadB[i][j])+shift);
	     
	     }
	     TargetX=TargetY=-999.;
	     if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && (gate_29p||gate_30s)  && XpadB[i][j]>-999. && dEpadB[i][j]>50){
			 TargetX=fTargetX->Eval((452.+83.+(i*padsize)+(padsize/2)));
                         TargetY=fTargetY->Eval((452.+83.+(i*padsize)+padsize/2)); // PPACa center to window, inactive region, plus pad
			 hPadXBgeo->Fill(i,((XpadB[i][j]+shift)-TargetX)); // proper residual
			 hPadXB->Fill(i,(XpadB[i][j])-TargetX);
	     
	     }
	     //hPadXB->Fill(i,XpadB[i][j]);
	     if (!SsdOR) {

	     }
	     if (SsdOR) hPadB_3D->Fill(XpadB[i][j],YpadB[i][j],i);
	       //if (ssd_detail && !SsdOR) { //downscale beam condition
	       //  hBraggB_ds_ch[i]->Fill(dEpadB[i][j]);
	       //}
	       //if (ssd_detail && SsdOR) { // ssd-or condition
	       //  hBraggB_ssd_ch[i]->Fill(dEpadB[i][j]);
	       //}
             if (flag_ppac && gate_30s) {
	       
	       //goatface 29 Oct 2012 17:52:36 
	       
	       if (dEpadB[i][j]>10){
	       //if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && gate_30s && XpadB[i][j]>-999. && dEpadB[i][j]>50){
	         hBraggB_30s->Fill(i,dEpadB[i][j]);
	         hBraggB_3D->Fill(i,j,dEpadB[i][j]);
	         //hPadXB_30s->Fill(i,XpadB[i][j]);
	         if (!SsdOR) { //downscale beam condition
	           hBraggB_30s_ds_ch[i]->Fill(dEpadB[i][j]);
	           hBraggB_30s_ds->Fill(i,dEpadB[i][j]);
	           hBraggB_ds_ch[i]->Fill(dEpadB[i][j]);
	           hBraggB_3D_ds->Fill(i,j,dEpadB[i][j]);
	         } // end if : !SsdOR
	         // PUT PHYSICS HERE
	         if ( SsdOR) { // ssd-or condition
	           hBraggB_30s_ssd_ch[i]->Fill(dEpadB[i][j]);
	           hBraggB_30s_ssd->Fill(i,dEpadB[i][j]);
	           hPadYB_30s->Fill(i,YpadB[i][j]);
	           hBraggB_ssd_ch[i]->Fill(dEpadB[i][j]);
	           hBraggB_3D_ssd->Fill(i,j,dEpadB[i][j]);
	           
	           UShort_t k=i-1; // channel to check against
	           if (k==31) k--; // these channels are always bad
	           if (run < 1006){
	             if (k==40) k--;
	           if (i==24) k=21;
	           }
	           if (run > 1005) { 
	             if ( i==24) k=20;
	             if ( k==47) k--; // inclusion makes no sense because it's the last pad.
	           }
	           if (i>15) { // below this is junk, also avoids seg faults for i<0
	             Float_t deltaE=(dEpadB[i][j]-dEpadB[k][j])/i; 
	             //if (1) { // view all
	             if (deltaE>6) { // condition of reaction point
	               hBraggB_rp->Fill(i,deltaE);
	               
	             //work from here!
	             //fSiIsHit-TpadC[i][j]
	             
	             //cout << i << endl;
	             }
	           }
	           // Separate the Bragg peaks by the track number. 
	           hBraggB_30s_jch[j]->Fill(i,dEpadB[i][j]);
	         } // end if : SsdOR
	         if (!SsdOR) { // d/s condition
	         } // end if : !SsdOR
	       } // end if dEpadB>10
             } // end if : gate_30s
	     if (flag_ppac && gate_29p) {
	       if (dEpadB[i][j]>10){
	     //if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.))){
	     //can check window cut
	         hBraggB_29p->Fill(i,dEpadB[i][j]);
	         if (!SsdOR) { //downscale beam condition
	           hBraggB_29p_ds->Fill(i,dEpadB[i][j]);
	           hBraggB_29p_ds_ch[i]->Fill(dEpadB[i][j]);
                 }	
	         if ( SsdOR) { // ssd-or condition
	           hBraggB_29p_ssd->Fill(i,dEpadB[i][j]);
	           hBraggB_29p_ssd_ch[i]->Fill(dEpadB[i][j]);
	         }
	       } // end if: dEpadB > 10
	     } // end if: gate_29p
             if (gate_17f){
	       if (dEpadB[i][j]>10){
	         hBraggB_17f->Fill(i,dEpadB[i][j]);
                 if (SsdOR) {
	           hBraggB_17f_ssd_ch[i]->Fill(dEpadB[i][j]);
	         }
	         else{ // downscale
	           hBraggB_17f_ds_ch[i]->Fill(dEpadB[i][j]);

	         }
	       } // end if: dEpadB > 10
	     }
            } // end for j: track loop
	  } // end for i : pad number loop
          
	 
	  //right
	  for (int i=0;i<20;i++)
	    dEpadRTotal[i]=0; // 
	  for(UShort_t i=0;i<8;i++){ // loop over all pads (now by pad number)
            for(UShort_t j=0;j<1;j++){ // just one track
            //for(UShort_t j=0;j<TracksL[i];j++){ // loop for each track
	      if ( dEpadR[i][j]>1.){// 
	        dEpadRTotal[j] = dEpadRTotal[j] + dEpadR[i][j];
	        hBraggR->Fill(i,dEpadR[i][j]);
	        if (SiIsHit[4] && dEpadR[i][j]>1.){
	          hPadYR->Fill(i,YpadR[i][j]);
	          hPadZR[0]->Fill(i,ZpadR[i][j]);
	          //if (SiIsHit[3] && ssd_mult==1){
	          /*
		  if (SiIsHit[2] && pulsetrack3a<100 && pulsetrack3a>0 && padsHitL[1]>1 && gate_30s){
	            hPadZL[pulsetrack3a]->Fill(i,ZpadL[i][j]);
                  } 
		  */
                }
	      }
	    }	
	  }
	  //left
	  for (int i=0;i<20;i++)
	    dEpadLTotal[i]=0; // 
	  for(UShort_t i=0;i<8;i++){ // loop over all pads (now by pad number)
            //for(UShort_t j=0;j<1;j++){ // just one track
            for(UShort_t j=0;j<TracksL[i];j++){ // loop for each track
	      if ( dEpadL[i][j]>1.){// 
	        dEpadLTotal[j] = dEpadLTotal[j] + dEpadL[i][j];
	        hBraggL->Fill(i,dEpadL[i][j]);
	        //if (SiIsHit[2] && dEpadL[i][j]>1.){
	        if (SiIsHit[2] && dEpadL[i][j]>1.&& fSiE[2][1]>200 && fSiEStripmax_ch[2]==1){
	          hPadYL->Fill(i,YpadL[i][j]);
	          hPadZL[0]->Fill(i,ZpadL[i][j]);
	          //if (SiIsHit[3] && ssd_mult==1){
	          if (SiIsHit[2] && pulsetrack3a<100 && pulsetrack3a>0 && padsHitL[1]>1 && gate_30s){
	            hPadZL[pulsetrack3a]->Fill(i,ZpadL[i][j]);
                  }   
                }
	      }
	    }	
	  }
	  if (SiIsHit[2] && padsHitL[1]>1 && gate_30s) pulsetrack3a++;
	  if (SiIsHit[2] && gate_30s) ssd_evts[2]++;
	 //center
	  for (int i=0;i<20;i++)
	    dEpadCTotal[i]=0; // 
	        // test rotation of high gain center from beam pad rotation parameters 28 Jan 2013 17:26:29 
		double shift;
		double paddistance;
		//const double theta=.03961897402027128139; //measured one	
		const double theta=.01980872946019476473; //measured one -- seems correct 17 Apr 2013 17:32:51 
		// old one 29 Jan 2013 23:35:34 
		//Int_t rotpoint=36; // imaginary pads (should be the same as using 12 for beam gem, with adding all 48 beam gem pads)
			 //[Int_t rotpoint=340.5; // in mm, 133+192+15.5
			 //Int_t rotpoint=548.5; // in mm, 341+192+15.5
			 Int_t rotpoint=493.5; // in mm, 286+192+15.5
			 //from beginning of pad, estimated from geometry needing ~150mm to shift by 6 in Y, -17 mm to front of pad from foot
	  
	  for(UShort_t i=0;i<8;i++){ // loop over all pads (now by pad number)
		//testing rotation 28 Jan 2013 17:27:52 
		//paddistance=((i+rotpoint+0.5)*padsize+15); // rotation origin is before high gain gem, so positive, and 15 is distance from beam gem to central gem
		// maybe easier to read
		paddistance=((i+0.5)*padsize)+rotpoint;
		shift=paddistance*TMath::Tan(theta);
	    //for(UShort_t j=0;j<2;j++){ // first two tracks debugging pad 0 13 Mar 2013 22:49:18 
	    for(UShort_t j=0;j<1;j++){ // first track
                XpadC[i][j]+=shift; 
            //for(UShort_t j=0;j<TracksC[i];j++){ // loop for each track
	      if ( dEpadC[i][j]>1.){// 
	      //if ( dEpadC[i][j]>1. && proton_beam ){// proton beam (I assume?)
	      //if ( dEpadC[i][j]>1. && ((fSiE[1][1]>110.) || ( fSiE[1][2]>110.) || ( fSiE[1][3]>120.))){// proton beam (I assume?)
	      //if ( dEpadC[i][j]>1. && ((fSiE[1][1]<1100. && fSiE[1][1]>900.) || (fSiE[1][2]<1100. && fSiE[1][2]>900.) || (fSiE[1][3]<1100. && fSiE[1][3]>900.))){// proton beam (I assume?)
              //if (dEpadC[i][j]>1. && (SiIsHit[0] == true || SiIsHit[1] == true)) { // 
	        dEpadCTotal[j] = dEpadCTotal[j] + dEpadC[i][j];
	        
		if ( padsHitC[1]>1  ) hBraggC->Fill(i,dEpadC[i][j]);
	        //if (gate_17f && SiIsHit[1] && padsHitC[1]>1  ) hBraggC->Fill(i,dEpadC[i][j]);
	        if (padsHitC[1]>1){
	        //if (SiIsHit[0] && ssd_mult==1){
	          //if (1){
		  //if (SiIsHit[0]) {
		    //XpadC[i][j]+=shift;
		    //cout << "shift :" << i << " = " << shift << endl;
		  if (fSiE[0][3]>200) {
		    if (XpadC[i][j]>-40 && XpadC[i][j]<68)
		    hPadXC1asum3->Fill(i,XpadC[i][j]);
		  }
		  if (fSiE[0][4]>200) {
		    if (XpadC[i][j]>-40 && XpadC[i][j]<68)
		    hPadXC1asum4->Fill(i,XpadC[i][j]);
		  }
		  if (fSiE[0][5]>200) {
		    if (XpadC[i][j]>-40 && XpadC[i][j]<68)
		    hPadXC1asum5->Fill(i,XpadC[i][j]);
		  }
		    //hPadXC1asum->Fill(i,XpadC[i][j]+shift);
	            if (SiIsHit[0])
		    hPadYC1a->Fill(i,YpadC[i][j]);
		  
	          if (fSiE[1][4]>200) {
	          //if (SiIsHit[1]) {
		    hPadXC2asum->Fill(i,XpadC[i][j]);
	            hPadYC2a->Fill(i,YpadC[i][j]);
		  }
	          /*
	          for(int k=1;k<9;k++){ 
		    if (fSiE[0][k]>300 && fSiE[0][k]<4000)
		    hPadXC[k]->Fill(i,XpadC[i][j]);
                  }*/
                }
		if (SiIsHit[0] && pulsetrack1a<100 && padsHitC[1]>1 && //gate_30s &&
	         fSiE[0][8]>0 ) {
		//(fSiE[0][1]>0 || fSiE[0][2]>0 || fSiE[0][3]>0 || fSiE[0][4]>0 || 
		 //fSiE[0][5]>0 || fSiE[0][6]>0 || fSiE[0][7]>0 || fSiE[0][8]>0 )){
		//if (SiIsHit[0] && pulsetrack2<100 && pulsetrack2>0 && padsHitC[1]>1){
		  hPadXC1a[pulsetrack1a]->Fill(i,XpadC[i][j]);
		}
		if (SiIsHit[1] && pulsetrack2a<100 && padsHitC[1]>1 && gate_17f ){//&&
	         //fSiE[1][3]>0 ) {
		 //(fSiE[1][1]>0 || fSiE[1][2]>0 || fSiE[1][3]>0 || fSiE[1][4]>0 || 
		 //fSiE[1][5]>0 || fSiE[1][6]>0 || fSiE[1][7]>0 || fSiE[1][8]>0 )){
		  hPadXC2a[pulsetrack2a]->Fill(i,XpadC[i][j]);
                }   
	      }
	    }
	  }
	  // ssd projection test 06 Feb 2013 22:43:54 
	  //if (dEpadC[3][0]>1 && dEpadC[4][0]>1 && fSiE[0][5]>200){
            //double ytest = (((XpadC[4][0]-XpadC[3][0])/4.)*(67.-14.))+XpadC[3][0];
	  if (dEpadC[1][0]>10 && dEpadC[5][0]>10 && fSiE[0][3]>200){
            //double ytest = (((XpadC[5][0]-XpadC[1][0])/16.)*(-290+6.))+XpadC[1][0];
            //double ytest = (((XpadC[4][0]-XpadC[3][0])/4.)*(99.-14.))+XpadC[3][0];
            //double ytest = (((XpadC[5][0]-XpadC[1][0])/16.)*(99.-6.))+XpadC[1][0];
	    //cout << "Projection\t" << XpadC[5][0] << "\t" << XpadC[1][0] << "\t" << ytest << endl;
	    //centerprojsum+=ytest;
	    //ytestcounter++;
	  }
          if (padsHitC[1]>1 && fSiE[0][5]>200  && fSiEStripmax_ch[0]==5){
	    // some while or for loop which removes points? 19 Feb 2013 13:43:52 
	    int points = padsHitC[1];
	    TGraph *Graph;
	    Graph = new TGraph(points);
            int jj=0;
	    for (int ii=0;ii<8;ii++){
	      if(TracksC[ii]>=1){
	        Graph->SetPoint(jj,(2+(ii*4)),XpadC[ii][0]);
		//cout << "point: " << jj << "\t" << (2+(ii*4)) << "\t" << XpadC[ii][0] << endl;
	        jj++;
	      }
	    }
	    TF1 *linear;
	    linear = new TF1("linear","pol1");
	    Graph->Fit("linear","Q");
	    double linpar[2];
	    linear->GetParameters(&linpar[0]);
	    if (1){
	    //if (linear->GetChisquare()<10){
	      //double ytest = ((-270.25)*linpar[1])+linpar[0]; // alpha source plane see Vol 5 pp. 4-5
	      double ytest = ((98.15)*linpar[1])+linpar[0]; // SSD plane see Vol 5 pp.4-5
	      //double ytest = ((-50)*linpar[1])+linpar[0]; // trying to find the X=0 intercept upstream
	      //double ytest = linpar[1]*99+linpar[0];
	      centerprojsum+=ytest;
	      ytestcounter++;
	      ypointscounter+=points;
	    }
	    delete linear;
	    delete Graph;
	  }


	  if (SiIsHit[0] &&  padsHitC[1]>1 && //gate_30s && 
	  fSiE[0][8]>0 ) {
         //  (fSiE[0][1]>0 || fSiE[0][2]>0 || fSiE[0][3]>0 || fSiE[0][4]>0 || 
         //  fSiE[0][5]>0 || fSiE[0][6]>0 || fSiE[0][7]>0 || fSiE[0][8]>0 )){
	    pulsetrack1a++;
	  }
	  if (SiIsHit[1] &&  padsHitC[1]>1 && gate_17f ){//&& 
//	  fSiE[1][3]>0 ) {
          // (fSiE[1][1]>0 || fSiE[1][2]>0 || fSiE[1][3]>0 || fSiE[1][4]>0 || 
          // fSiE[1][5]>0 || fSiE[1][6]>0 || fSiE[1][7]>0 || fSiE[1][8]>0 )){
	    pulsetrack2a++;
	  }
	  if (SiIsHit[0] && gate_30s ) ssd_evts[0]++;
	  if (SiIsHit[1] && gate_30s ) ssd_evts[1]++;
          
	  /*
	  for (UShort_t i=0;i<3;i++){
	    if (dEpadCTotal > 1. && fSiE[1][i+1]>10. && fSiE[1][i+1]<4000.){
	      hTpc_ch_dE[i+9]->Fill(fSiEcal[1][i+1],dEpadCTotal);	
	    }
	  }
          */


	 //multiplicity
	 //if (flag_ppac){
	 if(1){
	   if (padsHitB[1]>0 && PpacIsHit[0] && PpacIsHit[1]) hMultB->Fill(padsHitB[1]);
	   if (padsHitL[1]>0 &&  SiIsHit[2] && dEpadLTotal[0]>50.) hMultL_3a->Fill(padsHitL[1]);
	   //if (padsHitL[1]>0 && PpacIsHit[0] && PpacIsHit[1] && SiIsHit[2] && dEpadLTotal[0]>50.) hMultL_3a->Fill(padsHitL[1]);
	   if (flag_ssd && dEpadRTotal[0]>1 && SiIsHit[0]){ 
             hBraggR_E[0]->Fill(dEpadRTotal[0],fSiE[0][0]);
	     for (int i=1;i<9;i++){
               if (padsHitR[1]==i) {
	         hBraggR_E[i]->Fill(dEpadRTotal[0],fSiE[0][0]);
	         }
	     }
	   }
	   //if (flag_ssd && dEpadCTotal[0]>1 && fSiE[0][4]>120){ 
	   //if (flag_ssd && dEpadCTotal[0]>1 && fSiEmax[0]>50){ 
	   if (flag_ssd && dEpadCTotal[0]>0.5 && SiIsHit[0]){ 
             //hBraggC_E[0]->Fill(dEpadCTotal[0],fSiE[0][4]);
             //hBraggC_E[0]->Fill(dEpadCTotal[0],fSiEmax[0]); // alpha-calib-post
             //hBraggC_E[0]->Fill(dEpadCTotal[0],fSiE[0][0]);  // production
	     for (int i=1;i<9;i++){
	       if (dEpadC[i-1][0] > 80){ 
	         //hBraggC_E[i]->Fill(dEpadC[i-1][0],fSiEmax[0]);
	         //if (fSiE[0][4]>=200 && i==7) hBraggC_Edown[i]->Fill(dEpadC[i-1][0],fSiE[0][4]);
	         if (fSiE[0][4]>=290 && fSiE[0][4]<=350) hBraggC_Edown[i]->Fill(dEpadC[i-1][0],fSiE[0][4]);
	         if (fSiE[0][4]>=360 && fSiE[0][4]<=430) hBraggC_Eup[i]->Fill(dEpadC[i-1][0],fSiE[0][4]);
	         //if (fSiE[0][4]>200 && i==7) cout << dEpadC[6][0] << endl;
		 //hBraggC_E[i]->Fill(dEpadC[i-1][0],fSiEmax[0]);
	         //hBraggC_E[i]->Fill(dEpadCTotal[0],fSiEmax[0]); // alpha calib post
	         //hBraggC_E[i]->Fill(dEpadCTotal[0],fSiE[0][0]); //production
               }
	       if (padsHitC[1]==i) {
	         // CNS Interview Feb 2013 dirty way
		 hBraggC_E[0]->Fill((dEpadCTotal[0]*(6/i)),fSiE[0][0]);
	       }
	     }
	   }
	   if (flag_ssd && dEpadLTotal[0]>1 && SiIsHit[2]){ 
	   //if (flag_ssd && dEpadCTotal[0]>1 && SiIsHit[0]){ 
             hBraggL_E[0]->Fill(dEpadLTotal[0],fSiE[2][0]);
	     for (int i=1;i<9;i++){
               if (padsHitL[1]==i) {
	         hBraggL_E[i]->Fill(dEpadLTotal[0],fSiE[2][0]);
	         }
	     }
	   }
	   if (gate_30s && SsdOR){
	     if (padsHitB[1]>0 ) hMultB_30s->Fill(padsHitB[1]);
	     
	   }
	   if (1){
	     if (dEpadCTotal[0]>10.){
	     //if (!(dEpadC[5][1]>5) && !(padsHitC[1]==1)){
	       //for (int i=0;i<8;i++){
	       
	         //if (padsHitC[1]>0 &&  SiIsHit[0] && gate_30s) hMultC_30s_1a->Fill(padsHitC[1]);
	         if (padsHitC[1]>0 && SiIsHit[0] && fSiE[0][4]>190) hMultC_30s_1a->Fill(padsHitC[1]);
	         if (padsHitC[1]>0 && SiIsHit[1] && gate_30s) hMultC_30s_2a->Fill(padsHitC[1]);
	       //}
	     }
	     if (dEpadLTotal[0]>50.){
	       if (padsHitL[1]>0 && SiIsHit[2] && gate_30s) hMultL_30s_3a->Fill(padsHitL[1]);
	     }
	   }
	 }

// end multiplicity

        } // end if: flag_detail

      } // end if: flag_tpc

  // arr->Clear(); 
   fTargetX->ReleaseParameter(0);
   fTargetX->ReleaseParameter(1);
   fTargetX->ReleaseParameter(2);
   fTargetY->ReleaseParameter(0);
   fTargetY->ReleaseParameter(1);
   fTargetY->ReleaseParameter(2);
   
   } // close loop on jentry
   cout << "ssd_evts[0]: " << ssd_evts[0] << endl;
   cout << "pulsetrack1a: " << pulsetrack1a << endl;
   cout << "ssd_evts[1]: " << ssd_evts[1] << endl;
   cout << "pulsetrack2a: " << pulsetrack2a << endl;
   cout << "ssd_evts[2]: " << ssd_evts[2] << endl;
   cout << "pulsetrack3a: " << pulsetrack3a << endl;
   cout << "projected center: " << centerprojsum/ytestcounter << endl;
   cout << "No events in projection: " << ytestcounter << endl;
   cout << "Average no. pads: " << ypointscounter/ytestcounter << endl;
   delete fTargetX;
   delete fTargetY;
   delete ssd_strip_calib;
   delete ssd_padE_calib;
   delete ssd_padT_calib;
   delete ppac_calib;
   delete rf_calib;
   delete rf_ds_calib;
   delete beam_x_right_calib;
   delete beam_x_left_calib;
} // close Analyzer::Loop
