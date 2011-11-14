/*
 Copyright: 
 2010 daid KAHL, OTA shinsuke, HASHIMOTO takashi, YAMAGUCHI hidetoshi 
 2011 daid KAHL

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

Bool_t gate_30s=false,gate_29p=false,gate_31s=false;
Bool_t proton_beam=false;
Bool_t goodRF_30s=false, goodPpacX_30s=false;
Bool_t goodRF_29p=false, goodPpacX_29p=false;
Bool_t goodRF_31s=false, goodPpacX_31s=false;

// TPC pad mapping 
Short_t tpc_ch[2][144]={

{ -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,01:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,02:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,03:DISABLED
  40,41,42,43,44,45,46,47,  // SET 0,04:Beam right downstream
  32,33,34,35,36,37,38,39,  // SET 0,05:Beam right downstream
  24,25,26,27,28,29,30,31,  // SET 0,06:Beam right downstream
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,07:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,08:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,09:DISABLED
  16,17,18,19,20,21,22,23,  // SET 0,10:Beam right upstream
  8,9,10,11,12,13,14,15,    // SET 0,11:Beam right upstream
  0,1,2,3,4,5,6,7,          // SET 0,12:Beam right upstream
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,13:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,14:DISABLED
  0,1,2,3,4,5,6,7,          // SET 0,15:Left upstream
  7,6,5,4,3,2,1,0,          // SET 0,16:Center right      should it be?: // 0,1,2,3,4,5,6,7, // Center right
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,17:DISABLED
  0,1,2,3,4,5,6,7},         // SET 0,18:Right upstream

 { 31,30,29,28,27,26,25,24, // SET 1,01:Beam left downstream
  39,38,37,36,35,34,33,32,  // SET 1,02:Beam left downstream
  47,46,45,44,43,42,41,40,  // SET 1,03:Beam left downstream
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,04:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,05:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,06:DISABLED 
  7,6,5,4,3,2,1,0,          // SET 1,07:Beam left upstream
  15,14,13,12,11,10,9,8,    // SET 1,08:Beam left upstream
  23,22,21,20,19,18,17,16,  // SET 1,09:Beam left upstream
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,10:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,11:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,12:DISABLED
  0,1,2,3,4,5,6,7,          // SET 1,13:Center left
  0,1,2,3,4,5,6,7,          // SET 1,14:Left downstream
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,15:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,16:DISABLED
  0,1,2,3,4,5,6,7           // SET 1,17:Right downstream
  -1,-1,-1,-1,-1,-1,-1,-1}  // SET 1,18:DISABLED
};

//internal calibration for beam pad
Double_t padBgain[2][48]={
{ // runs 1001-1005
  1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,// 0-7
  1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 8-15
  0.712629, .69137310,0.706643,0.732711,.74086605,0.655,1.,1., // 16-23
  0.708874,0.706074,0.70031,0.701813,0.696598713,0.70119423,0.6988148,1., //24-31
  0.69839515,0.693352,0.69736996,0.696989,0.697088,0.696793,0.69828,0.696932,  //32-39
  1.,0.65594417,0.701976,0.701372,0.706313,0.70744137,.70876,0.68},  //40-47
  //46 is weird here even though the same method is used...
  //1.,0.65594417,0.701976,0.701372,0.706313,0.70744137,.77039460,0.68},  //40-47

{ //other runs
  1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 0-7
  1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 8-15
  1.  ,1.10,1.  ,1.  ,0.95,1.  ,1.  ,1.  , // 16-23
  1.  ,1.  ,1.  ,1.  ,0.97,0.99,0.94,1.  , // 24-31
  0.97,1.  ,0.98,1.  ,1.  ,1.  ,1.  ,1.  , // 32-39
  1.03,1.03,1.  ,1.  ,1.  ,0.97,0.90,1.  // 40-47
}};


/*
Double_t padBgain[48]={
1.,1.,1.,1.,1.,1.,1.,1.,
1.,1.,1.,1.,1.,1.,1.,1.,
1.,1.,1.,1.,1.,1.,1.,1.,
1.,1.,1.,1.,1.,1.,1.,1.,
1.,1.,1.,1.,1.,1.,1.,1.,
1.,1.,1.,1.,1.,1.,1.,1.};
*/

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

   Calibration *ssd_strip_calib=new Calibration("ssd_strip_calib.dat",96);//12x8
   Calibration *ssd_pad_calib=new Calibration("ssd_pad_calib.dat",19);
   Calibration *ppac_calib=new Calibration("ppac_calib.dat",10);
   
   if (flag_detail){
     if (flag_strip){
       ssd_strip_calib->ScanFile();
       ssd_strip_calib->ShowEntries();
     }
     if (flag_ssd){
       ssd_pad_calib->ScanFile();
       ssd_pad_calib->ShowEntries();
     }
     if (flag_ppac){
       ppac_calib->ScanFile();
       ppac_calib->ShowEntries();
       //PPACa -- check these values! 24 Jan 2011 11:56:55 
       PpacOffset[0][0]=0.92; //18 Jul 2011 16:56:44 
       PpacOffset[0][1]=1.58; //18 Jul 2011 16:56:57 
       PpacOffsetLine[0][0]=0.0; //18 Jul 2011 16:58:08 
       PpacOffsetLine[0][1]=-0.22; //18 Jul 2011 16:58:10 
       PpacPositionGain[0][0]=0.6200; // 18 Jul 2011 16:54:57 
       PpacPositionGain[0][1]=0.6210; // 18 Jul 2011 16:54:59 
       PpacOffsetGeometry[0][0]=-0.60;
       PpacOffsetGeometry[0][1]=0.0; // checked 08 Jul 2011 16:00:27 
       //PpacOffsetGeometry[1][1]=1.83;
       //PPACb -- check these values! 24 Jan 2011 11:56:52 
       PpacOffset[1][0]=0.17; //18 Jul 2011 16:57:56 
       PpacOffset[1][1]=0.11; //18 Jul 2011 16:57:54 
       PpacOffsetLine[1][0]=-4.3; //18 Jul 2011 16:58:03 
       PpacOffsetLine[1][1]=-1.0; //18 Jul 2011 16:58:04 
       PpacPositionGain[1][0]=0.6205; //18 Jul 2011 16:55:04 
       PpacPositionGain[1][1]=0.6125; //18 Jul 2011 16:55:06 
       PpacOffsetGeometry[1][0]=0.0;
       PpacOffsetGeometry[1][1]=0.0; // checked 08 Jul 2011 16:00:22 
       //PpacOffsetGeometry[1][1]=-3.5;
     }
   } // end if: flag_detail
   //int ssdTime[18] = {0};

   //Loop on entries
   for (Long64_t jentry=0; jentry<nentries-1;jentry++) {//-1...needed when the last event is not complete.
   // cout << jentry << endl;
      if (jentry%100==0 || jentry==nentries-2) {
         cout << jentry << "/" << nentries-1 << "\r"<< flush;
      }
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      //      arr->Clear(); 


      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //      cout << nb << "Bytes, Entry:"<< ientry << ","<<jentry << endl;
      //      fChain->Show(jentry);


      // if (Cut(ientry) < 0) continue;
	
      //Can begin filling any histrograms you defined
      
      //SSD Pad
      if (flag_ssd) { // process SSD data?
        for (UShort_t i=0; i<18; i++) {
          SiIsHit[i] = false;
          if ( fSiT[i] > -90 ) {
            SiIsHit[i] = true;
          } // end if : fSiT[i] > 0
        } // end for i
        if (flag_raw) { // process raw SSD data?
          for (UShort_t i=0; i<18; i++) {
            if (SiIsHit[i]) {
              hSiPadHit->Fill(i);	
              hSiPadE->Fill(i,fSiE[i][0]);
              hSiT->Fill(i,fSiT[i]);
              pad_ch[i]->Fill(fSiEcal[i][0]);
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
            if (SiIsHit[i] == true ){
              fSiTcal[i]=0.09765625*fSiT[i]; //calibration from ch to ns
            //if (SiIsHit[i] == true && proton_beam == true){
              //if (i==1)
              hSiTcal->Fill(i,fSiT[i]);
              //hSiT->Fill(1,fSiT[1]);
              //Fill individual ssd data 
              if (fSiE[i][0] > 0 ){ // don't fill the histogram with junk
                fSiEcal[i][0]=ssd_pad_calib->CalibMeV(i,fSiE[i][0]);
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
        }//end if: flag_detail
      } // end if: ssd
     
     if (flag_strip){
        // SSD Strip
        if (flag_raw) { // process raw SSD strip data?
          for(UShort_t i=0;i<12;i++){
            for (UShort_t j=1; j<9; j++) {
              hSiStripE->Fill(i*8+j,fSiE[i][j]);
              strip_ch[i*8+j]->Fill(fSiE[i][j]); 
            } // end for: j
	  } // end for: i
	} // end if: ssd_strip_raw
	if (flag_detail){
          for(UShort_t i=0;i<12;i++){
            fSiEmax[i]=0;
            for (UShort_t j=1; j<9; j++) {
              fSiEcal[i][j]=ssd_strip_calib->CalibMeV(i*8+j,fSiE[i][j]);
              if (fSiEcal[i][j]>fSiEmax[i] && fSiE[i][j]<4000.) {//find maximum strip
                fSiEmax[i]=fSiEcal[i][j];
              }
              //Fill each strip data 
              hSiStripEcal->Fill(i*8+j,fSiEcal[i][j]);
              //if ( fSiE[i][j] > 0 && fSiE[i][j] < 4000 && proton_beam){ // just for high and low pedestals (high part is too tight)
              //if ( fSiE[i][j] > 0 && fSiE[i][j] < 4097){ // just for high and low pedestals (high part is too tight)
              //if ( (SiIsHit[i]) && (fSiE[i][j] > -100) ){ // just for high and low pedestals (high part is too tight)
              //if ( (fSiE[i][j] > -100 && !SsdOR) ){
	      if ( (fSiE[i][j] > -100) ){
                strip_cal_ch[i*8+j]->Fill(fSiEcal[i][j]);
              }
              /*
              if (gate_29p) {
                 strip_ch_gated29p[i*8+j]->Fill(fSiE[i][j]);
              }
              if (gate_30s) {
                 strip_ch_gated30s[i*8+j]->Fill(fSiE[i][j]);
              }*/
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
  
      } // end if: flag_detail & flag_ssd & flag_strip

      //PPAC
      if (flag_ppac) { // process PPAC?
	if (flag_raw) { // raw PPAC data?
	  for (UShort_t i=0;i<2;i++) if (fRF[i] > 0) hRf[i]->Fill(fRF[i]);
          
          //raw PPAC data
          if (fPPAC[0][0]>0.) hPpac0TX1->Fill(fPPAC[0][0]);
          if (fPPAC[0][1]>0.) hPpac0TX2->Fill(fPPAC[0][1]);
          if (fPPAC[0][2]>0.) hPpac0TY1->Fill(fPPAC[0][2]);
          if (fPPAC[0][3]>0.) hPpac0TY2->Fill(fPPAC[0][3]);
          if (fPPAC[1][0]>0.) hPpac1TX1->Fill(fPPAC[1][0]);
          if (fPPAC[1][1]>0.) hPpac1TX2->Fill(fPPAC[1][1]);
          if (fPPAC[1][2]>0.) hPpac1TY1->Fill(fPPAC[1][2]);
          if (fPPAC[1][3]>0.) hPpac1TY2->Fill(fPPAC[1][3]);
        
	} // end if: flag_raw
	
	if (flag_detail) { // detailed PPAC data?
          //Calibration 
          
          for(UShort_t i=0;i<2;i++){
	  PpacIsHit[i] = false;
	  if (fPPAC[i][0]>0. && fPPAC[i][1]>0. && fPPAC[i][2]>0. && fPPAC[i][3]>0.) PpacIsHit[i]=true;
            for(UShort_t j=0;j<4;j++){
              fPPACcal[i][j]=0.;
	      //check here for PPACb X centering issue 18 Jul 2011 17:37:49 
	      if ( fPPAC[i][j] > 0. )  {
	        fPPACcal[i][j]=ppac_calib->Calib(i*5+j,fPPAC[i][j]);
	      }
            }    
          }

          // PPAC calibration for loop
          for (UShort_t i=0;i<2;i++){
	    //average time
	    tof[i] = 0.25*(fPPACcal[i][0]+fPPACcal[i][1]+fPPACcal[i][2]+fPPACcal[i][3]);
	    dxy[0]=(fPPACcal[i][0]-fPPACcal[i][1]);
            dxy[1]=(fPPACcal[i][2]-fPPACcal[i][3]);
            for (UShort_t ii = 0; ii < 2; ii++) {
             // time calibration
              dxy[ii] -=PpacOffset[i][ii];
              dxy[ii] -= PpacOffsetLine[i][ii];
             // position calibration
              dxy[ii] *= PpacPositionGain[i][ii];
              dxy[ii] -= PpacOffsetGeometry[i][ii];
            }
              PpacX[i]=dxy[0];
              if (i==0 ){PpacY[i]=dxy[1];} // PPACaY is normal
              if (i==1 ){PpacY[i]=(dxy[1]*-1.0);} // PPACbY cables were inverted
          } //end PPAC calibration for loop
          if (tof[0]>0. && tof[1]>0. && PpacIsHit[0] && PpacIsHit[1] ) Tof=tof[1]-tof[0]; 
          
	  // fill basic PPAC histograms
          hRF0Tof->Fill(fRF[0],Tof);
          hRF1Tof->Fill(fRF[1],Tof);
          hPpac0XY->Fill(PpacX[0],PpacY[0]);
          hPpac1XY->Fill(PpacX[1],PpacY[1]);
          if (fRF[0] > 0)
          {
            hPpac0XRF0->Fill(PpacX[0],fRF[0]);
            hPpac1XRF0->Fill(PpacX[1],fRF[0]);
          }
          if (fRF[1] > 0)
          {
            hPpac0XRF1->Fill(PpacX[0],fRF[1]);
            hPpac1XRF1->Fill(PpacX[1],fRF[1]);
          }
          
	  // Turn off all the gates for a new event
          gate_30s=false;
	  gate_29p=false;
          gate_31s=false;
          goodRF_30s=false;
	  goodPpacX_30s=false;
          goodRF_29p=false;
	  goodPpacX_29p=false;
          goodRF_31s=false;
	  goodPpacX_31s=false;
          // Set the gates for 30S/29P RI beams
	  // 30S done again w/ new PPAC calib 20 Sep 2011 17:22:52  
          if (  ((fRF[0]>=250) && (fRF[0]<=330)) || ((fRF[0]>=460)&&(fRF[0]<=540)) ||
	    ((fRF[1]>=255)&&(fRF[1]<=335)) || ((fRF[1]>=480)&&(fRF[1]<=530)) || ((fRF[1]>=95)&&(fRF[1]<=120))
	    ) {goodRF_30s=true;} //30S gate
          if ((PpacX[0]>=-9 && PpacX[0]<=8) && (PpacX[1]>=-10 && PpacX[1]<=10)) {goodPpacX_30s=true;} //30S gate
          if (goodRF_30s && goodPpacX_30s) {gate_30s=true;}
          // PPAC X was redone 20 Sep 2011 17:23:07 
	  if (  ((fRF[0]>=160) && (fRF[0]<=210)) || ((fRF[0]>=375)&&(fRF[0]<=425)) ||
	    ((fRF[1]>=160)&&(fRF[1]<=220)) || ((fRF[1]>=380)&&(fRF[1]<=440))
	     ) {goodRF_29p=true;} //29P gate
          if ((PpacX[0]>=-4 && PpacX[0]<=9) && (PpacX[1]>=-5 && PpacX[1]<=14)) {goodPpacX_29p=true;} //29P gate
          if (goodRF_29p && goodPpacX_29p) {gate_29p=true;}
          // don't think this is 31S...but maybe something...
	  if (((fRF[0]>=355) && (fRF[0]<=430)) || ((fRF[1]>=140)&&(fRF[1]<=220))) {goodRF_31s=true;} //30S gate
          if ((PpacX[0]>=-6 && PpacX[0]<=2) && (PpacX[1]>=-5 && PpacX[1]<=3)) {goodPpacX_31s=true;} //30S gate
          if (goodRF_31s && goodPpacX_31s) {gate_31s=true;}
          
          if (gate_30s) {
	    hPpacToF_30s->Fill(Tof);
	      hPpac0XRF0_30s->Fill(PpacX[0],fRF[0]);
              hPpac1XRF0_30s->Fill(PpacX[1],fRF[0]);
	      hPpac0XRF1_30s->Fill(PpacX[0],fRF[1]);
              hPpac1XRF1_30s->Fill(PpacX[1],fRF[1]);
	  }
          if (gate_29p) {
	    hPpacToF_29p->Fill(Tof);
	      hPpac0XRF0_29p->Fill(PpacX[0],fRF[0]);
              hPpac1XRF0_29p->Fill(PpacX[1],fRF[0]);
	      hPpac0XRF1_29p->Fill(PpacX[0],fRF[1]);
              hPpac1XRF1_29p->Fill(PpacX[1],fRF[1]);
	  }
	  if (gate_30s){
	        hPpac0XRF1cut->Fill(PpacX[0],fRF[1]);
	        hPpac1XRF1cut->Fill(PpacX[1],fRF[1]);
	  }
        } // end if: flag_detail
      } // end if: flag_ppac

      // PPAC and SSD
      if (flag_detail && flag_ssd && flag_ppac){
	if (SsdOR){ // physics event
	  if (fRF[0] > 0. ) {
	    hPpac0XRF0ssd->Fill(PpacX[0],fRF[0]);
            hPpac1XRF0ssd->Fill(PpacX[1],fRF[0]);
            hPpac1XYcut->Fill(PpacX[1],PpacY[1]);
	  }
	  if (fRF[1] > 0. ) {
	    hPpac0XRF1ssd->Fill(PpacX[0],fRF[1]);
            hPpac1XRF1ssd->Fill(PpacX[1],fRF[1]);
	  }
	}
	else{ // downscale event
	  if (fRF[0] > 0.){
	    hPpac0XRF0ds->Fill(PpacX[0],fRF[0]);
            hPpac1XRF0ds->Fill(PpacX[1],fRF[0]);
	    if (gate_30s){
	    // can fill here
	    }
	  }
	  if (fRF[1] > 0.){
	    hPpac0XRF1ds->Fill(PpacX[0],fRF[1]);
            hPpac1XRF1ds->Fill(PpacX[1],fRF[1]);
	    if (gate_30s){
	      //hPpac1XRF1cut->Fill(PpacX[1],fRF[1]);
	    }
	  }
	}

      } // end if : flag_detail & flag_ssd & flag_ppac

      if (flag_detail && flag_strip && flag_ppac) { // PPAC and SSD Strip detailed analysis?
	for(UShort_t i=0;i<18;i++){
	    if (SiIsHit[i]){
	      if (gate_30s)  hPpac0TpadT_30s_ch[i]->Fill(fSiT[i]-tof[0]);
	      if (gate_29p) hPpac0TpadT_29p_ch[i]->Fill(fSiT[i]-tof[0]);
	    }
	}
      } // end if: flag_detail & flag_strip & flag_ppac

      //GEM-MSTPC
      if (flag_tpc){  // analyze TPC data?
        
	TArtFadcHit *fadc;
        Int_t nadclines = arr->GetEntriesFast();
        // Init the FADC Max arrays
        if (flag_detail) { // process TPC detailed analysis? 
	  for(UShort_t j=0;j<144;j++){
	    PadIsHit[j]=false;
            fNHit[j]=0;  
            for(UShort_t k=0;k<20;k++){ // fNHit
	      baseline[j][k] = 0.;
	      baseline_dev[j][k] = 0.;
              HitNo[j][k] = 1; // Multiple hit number 
              for(UShort_t l=0;l<20;l++){ // HitNo
              fSampleMax[j][k][l] = 0.; 
              fClockMax[j][k][l] = 0.;
	      fTimeMax[j][k][l] = 0.;
              }    
            }    
          }
	} // end if: flag_detail

	for (Int_t i=0; i<nadclines; i++) { // loop over the number of ADC lines
	  fadc=(TArtFadcHit *) arr->At(i);
          UShort_t id=fadc->fID; // pad ID (range 0 to 143)
          fNHit[id]=fadc->fHit;  // Pulse number: if hit never 0
          
	  if (flag_raw) {// process TPC raw data?
	    for (Int_t j=0; j<fadc->fNSample; j++) {
              hGemMstpcAll->Fill(fadc->fClock[j], fadc->fSample[j]);
              hGemMstpcCh->Fill(id,fadc->fSample[j]);
              hTpc_dE_ch[id]->Fill(fadc->fClock[j], fadc->fSample[j]);
              hfSample[id]->Fill(fadc->fSample[j]);
              if (id < 96){//Main (low gain) beam pads
                   hBeamPadCh->Fill(id,fadc->fSample[j]); 
              }
	      else {//high gain pads
                 hGemMstpcSide->Fill(fadc->fClock[j], fadc->fSample[j]);
	      } // end if: id<96
	    } // end for: j fNSample
	  } // end if: flag_raw

          if (flag_detail) { // process TPC detailed analysis?
	    if (fNHit[id]!=0) PadIsHit[id] = true;
	    else continue; // junk; truncate
	    
            //baseline check 
	    for(UShort_t j=0; j<5;j++){ // get the first five pulses
              if (fadc->fSample[j]>0.) baseline[id][fNHit[id]] = baseline[id][fNHit[id]] + fadc->fSample[j]; 
            } // end for j: baseline determination
            
	    // make sure there are 5 pulses
	    if ((fadc->fSample[0]!=0.) && (fadc->fSample[1]!=0.) && (fadc->fSample[2]!=0.) && (fadc->fSample[3]!=0.) && (fadc->fSample[4]!=0.)){
	      baseline[id][fNHit[id]]=baseline[id][fNHit[id]]/5.;
	      // find maximum deviation of the first five samples from the baseline
	      Double_t maxdeviation=0.;
	      for (UShort_t j=0;j<5;j++){
	        if (TMath::Abs(fadc->fSample[j]-baseline[id][fNHit[id]])>maxdeviation) maxdeviation=TMath::Abs(fadc->fSample[j]-baseline[id][fNHit[id]]);
	      }
              baseline_dev[id][fNHit[id]]=(TMath::Abs(maxdeviation/baseline[id][fNHit[id]]));
	    }
	    else continue; // junk; truncate
	
            // peak finder
	    trough=1e6; //init trough
            if (fadc->fNSample > 30 )  // make sure the pulse width is reasonable
              for (UShort_t j=3; j<fadc->fNSample-3; j++) { // loop all the FADC samples
                //hGemMstpcAll->Fill(fadc->fClock[j], fadc->fSample[j]);
	        if (fadc->fSample[j]<trough && fadc->fSample[j]>0.) trough=fadc->fSample[j];
	        if((fadc->fSample[j]) < 400) continue; // 500 should be the baseline ?
	          if(fSampleMax[id][fNHit[id]][HitNo[id][fNHit[id]]] > fadc->fSample[j]) continue; // is the present fSample larger than fSampleMax?
                    if( !(  // be sure it is a true peak, by looking +3 and -3 samples
		        (fadc->fSample[j-3] < fadc->fSample[j-2]) && 
                        (fadc->fSample[j-2] < fadc->fSample[j-1]) &&
                        (fadc->fSample[j-1] < fadc->fSample[j])   &&
                        (fadc->fSample[j]   > fadc->fSample[j+1]) &&
                        (fadc->fSample[j+1] > fadc->fSample[j+2]) &&
                        (fadc->fSample[j+2] > fadc->fSample[j+3])
			)) continue; // no peak -- skip
                      // true peak: fill the FADC Max arrays
                      if (TMath::Abs(trough-baseline[id][fNHit[id]])>(baseline_dev[id][fNHit[id]]*baseline[id][fNHit[id]])){ 
	                // trough and baseline don't agree!
	                hTrough[id]->Fill(trough-baseline[id][fNHit[id]]);
	              }
                      fSampleMax[id][fNHit[id]][HitNo[id][fNHit[id]]] = fadc->fSample[j]-baseline[id][fNHit[id]];
                      fClockMax[id][fNHit[id]][HitNo[id][fNHit[id]]] = fadc->fClock[j];
                      fTimeMax[id][fNHit[id]][HitNo[id][fNHit[id]]] = fadc->fTime; // we should modify fTime somehow for each peak?
                      HitNo[id][fNHit[id]]++; // valid peak: increment the hit number
	              trough=1e6; // reset trough
              } // end for j: FADC samples
          } // end if: flag_detail

	} // end for i: nadclines
   
        if (flag_detail) { // process TPC detailed analysis?
	  //*  *****************************
          //   START GEM-MSTPC PHYSICS HERE!
          //*  *****************************
          // init arrays
	  for (UShort_t i=0;i<48;i++){
            TracksB[i]=0; 
            for (UShort_t j=0;j<20;j++) {
              XpadB[i][j]=0.;
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
          
	  // Beam
          for(UShort_t i=0;i<96;i++) // all right Beam pads
            for(UShort_t j=0;j<96;j++){ // all left Beam pads
              if((tpc_ch[0][i]==-1)||(tpc_ch[1][j]==-1) || // tpc_ch will get us the pad number (0 to 47)
              (tpc_ch[0][i] != tpc_ch[1][j]) || 
              (PadIsHit[i] == false ) || (PadIsHit[j] == false )) continue; // only run for the same pads; make sure both pads are hit
                pad[i] = tpc_ch[0][i]; // place holder for the loop
                for(UShort_t k=1;k<=fNHit[i];k++) // loop on left side pulse number
	  	  for(UShort_t l=1;l<=HitNo[i][k];l++) // loop on left side hit number
                    for(UShort_t m=1;m<= fNHit[j];m++) // loop on right side pulse number
                      for(UShort_t n=1;n<=HitNo[j][m]; n++){ // loop on right side hit number
                        if (TracksB[pad[i]] >= 20 ) continue; // don't screw up memory; TracksB[i] cannot access more than 19
                        if(!(fSampleMax[i][k][l] > 0.1) || !(fSampleMax[j][m][n] > 0.1)  || // both samples are nonzero?
                           !(fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	  	        if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
	  	        XpadB[pad[i]][TracksB[pad[i]]] = 5.5*((fSampleMax[i][k][l]-fSampleMax[j][m][n])/
	  	          (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
                        YpadB[pad[i]][TracksB[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	  	          + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	  	        ZpadB[pad[i]][TracksB[pad[i]]] = 0.42*(pad[i]-23)-0.21;
                        dEpadB[pad[i]][TracksB[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
                        TpadB[pad[i]][TracksB[pad[i]]] = (fTimeMax[i][k][l]+fTimeMax[j][m][n])/2.;
                        TracksB[pad[i]]++; // increase the number of tracks for that pad
                      } // end for n: right side hit number
            } // end for j: loop over left Beam pads 

	  // Center
          for(UShort_t i=120;i<128;i++) // all right Center pads
            for(UShort_t j=96;j<104;j++){ // all left Center pads
              if((tpc_ch[0][i]==-1)||(tpc_ch[1][j]==-1) || // tpc_ch will get us the pad number (0 to 7)
              (tpc_ch[0][i] != tpc_ch[1][j]) || 
              (PadIsHit[i] == false ) || (PadIsHit[j] == false )) continue; // only run for the same pads; make sure both pads are hit
                pad[i] = tpc_ch[0][i]; // place holder for the loop
                for(UShort_t k=1;k<=fNHit[i];k++) // loop on left side pulse number
	  	for(UShort_t l=1;l<=HitNo[i][k];l++) // loop on left side hit number
                    for(UShort_t m=1;m<= fNHit[j];m++) // loop on right side pulse number
                      for(UShort_t n=1;n<=HitNo[j][m]; n++){ // loop on right side hit number
                        if (TracksC[pad[i]] >= 20 ) continue; // don't screw up memory; TracksC[i] cannot access more than 19
	  	        if(!(fSampleMax[i][k][l] > 0.1) || !(fSampleMax[j][m][n] > 0.1)  || // both samples are nonzero?
                           !(fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	  		if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
	  		XpadC[pad[i]][TracksC[pad[i]]] = 5.0*((fSampleMax[i][k][l]-fSampleMax[j][m][n])/
	  		  (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
                        YpadC[pad[i]][TracksC[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	  		  + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	  		ZpadC[pad[i]][TracksC[pad[i]]] = 9.575 + 0.42*(pad[i]); 
                        dEpadC[pad[i]][TracksC[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
                        TracksC[pad[i]]++; // increase the number of tracks for that pad
                      } // end for n: right side hit number
            } // end for j: loop over left Center pads 
          
	  // Left
          for(UShort_t i=112;i<120;i++) // all upstream Left pads
            for(UShort_t j=104;j<112;j++){ // all downstream Left pads
              if((tpc_ch[0][i]==-1)||(tpc_ch[1][j]==-1) || // tpc_ch will get us the pad number (0 to 7)
              (tpc_ch[0][i] != tpc_ch[1][j]) || 
              (PadIsHit[i] == false ) || (PadIsHit[j] == false )) continue; // only run for the same pads; make sure both pads are hit
                pad[i] = tpc_ch[0][i]; // place holder for the loop
                for(UShort_t k=1;k<=fNHit[i];k++) // loop on left side pulse number
	  	for(UShort_t l=1;l<=HitNo[i][k];l++) // loop on left side hit number
                    for(UShort_t m=1;m<= fNHit[j];m++) // loop on right side pulse number
                      for(UShort_t n=1;n<=HitNo[j][m]; n++){ // loop on right side hit number
                        if (TracksL[pad[i]] >= 20 ) continue; // don't screw up memory; TracksL[i] cannot access more than 19
	  	        if(!(fSampleMax[i][k][l] > 0.1) || !(fSampleMax[j][m][n] > 0.1)  || // both samples are nonzero?
                           (!fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	  		if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
	  		ZpadL[pad[i]][TracksL[pad[i]]] = 11.75*((fSampleMax[i][k][l]-fSampleMax[j][m][n])/
	  		  (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
                        YpadL[pad[i]][TracksL[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	  		  + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	  		XpadL[pad[i]][TracksL[pad[i]]] = -7.621 - 0.42*(pad[i]); 
                        dEpadL[pad[i]][TracksL[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
                        TracksL[pad[i]]++; // increase the number of tracks for that pad
                      } // end for n: right side hit number
            } // end for j: loop over downstream Left pads 
          
	  // Right
          for(UShort_t i=136;i<144;i++) // all upstream Right pads
            for(UShort_t j=128;j<136;j++){ // all downstream Right pads
              if((tpc_ch[0][i]==-1)||(tpc_ch[1][j]==-1) || // tpc_ch will get us the pad number (0 to 7)
              (tpc_ch[0][i] != tpc_ch[1][j]) || 
              (PadIsHit[i] == false ) || (PadIsHit[j] == false )) continue; // only run for the same pads; make sure both pads are hit
                pad[i] = tpc_ch[0][i]; // place holder for the loop
                for(UShort_t k=1;k<=fNHit[i];k++) // loop on left side pulse number
	  	for(UShort_t l=1;l<=HitNo[i][k];l++) // loop on left side hit number
                    for(UShort_t m=1;m<= fNHit[j];m++) // loop on right side pulse number
                      for(UShort_t n=1;n<=HitNo[j][m]; n++){ // loop on right side hit number
                        if (TracksR[pad[i]] >= 20 ) continue; // don't screw up memory; TracksR[i] cannot access more than 19
	  	        if(!(fSampleMax[i][k][l] > 0.1) || !(fSampleMax[j][m][n] > 0.1)  || // both samples are nonzero?
                           !(fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	  	        if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
	  	        ZpadR[pad[i]][TracksR[pad[i]]] = 11.75*((fSampleMax[i][k][l]-fSampleMax[j][m][n])/
	  	          (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
                        YpadR[pad[i]][TracksR[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	  	          + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	  	        XpadR[pad[i]][TracksR[pad[i]]] = 7.621 + 0.42*(pad[i]); 
                        dEpadR[pad[i]][TracksR[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
                        TracksR[pad[i]]++; // increase the number of tracks for that pad
                    } // end for n: right side hit number
            } // end for j: loop over downstream Right pads 
          
	  //beam
          for(UShort_t i=0;i<48;i++){ // loop over all pads (now by pad number)
	    for(UShort_t j=0;j<TracksB[i];j++){ // loop for each track
	     if ( run < 1006 ) dEpadB[i][j] = dEpadB[i][j] * padBgain[0][i]; // artificial internal calibration
	     if ( run > 1005 ) dEpadB[i][j] = dEpadB[i][j] * padBgain[1][i]; // artificial internal calibration
             hBraggB->Fill(i,dEpadB[i][j]);
	     hPadXB->Fill(i,XpadB[i][j]);
	     hPadYB->Fill(i,YpadB[i][j]);
	     if (SsdOR) hPadB_3D->Fill(XpadB[i][j],YpadB[i][j],i);
	       //if (ssd_detail && !SsdOR) { //downscale beam condition
	       //  hBraggB_ds_ch[i]->Fill(dEpadB[i][j]);
	       //}
	       //if (ssd_detail && SsdOR) { // ssd-or condition
	       //  hBraggB_ssd_ch[i]->Fill(dEpadB[i][j]);
	       //}
             if (gate_30s) {
	       
	       hBraggB_30s->Fill(i,dEpadB[i][j]);
	       hBraggB_ch[i]->Fill(dEpadB[i][j]);
	       hBraggB_3D->Fill(i,j,dEpadB[i][j]);
	       hPadXB_30s->Fill(i,XpadB[i][j]);
	       if (flag_ssd && !SsdOR) { //downscale beam condition
	         hBraggB_ds_ch[i]->Fill(dEpadB[i][j]);
	         hBraggB_3D_ds->Fill(i,j,dEpadB[i][j]);
	       } // end if : flag_ssd && !SsdOR
	       // PUT PHYSICS HERE
	       if (flag_ssd && SsdOR) { // ssd-or condition
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
	       } // end if : flag_ssd && SsdOR
	       if (flag_ssd && !SsdOR) { // d/s condition
	       } // end if : flag_ssd && !SsdOR
             } // end if : gate_30s
	     if (gate_29p) {
	       hBraggB_29p->Fill(i,dEpadB[i][j]);
	     } // end if: gate_29p
             if (gate_31s) hBraggB_31s->Fill(i,dEpadB[i][j]);
            } // end for j: track loop
	  } // end for i : pad number loop
          
	  //center
	  Double_t dEpadCTotal=0; // declare this somewhere else and just reset it here 04 Apr 2011 16:49:04 
	  
	  for(UShort_t i=0;i<8;i++){ // loop over all pads (now by pad number)
            //for(UShort_t j=0;j<1;j++){ // first track
            for(UShort_t j=0;j<TracksC[i];j++){ // loop for each track
	      if ( dEpadC[i][j]>1.){// 
	      //if ( dEpadC[i][j]>1. && proton_beam ){// proton beam (I assume?)
	      //if ( dEpadC[i][j]>1. && ((fSiE[1][1]>110.) || ( fSiE[1][2]>110.) || ( fSiE[1][3]>120.))){// proton beam (I assume?)
	      //if ( dEpadC[i][j]>1. && ((fSiE[1][1]<1100. && fSiE[1][1]>900.) || (fSiE[1][2]<1100. && fSiE[1][2]>900.) || (fSiE[1][3]<1100. && fSiE[1][3]>900.))){// proton beam (I assume?)
              //if (dEpadC[i][j]>1. && (SiIsHit[0] == true || SiIsHit[1] == true)) { // 
	        dEpadCTotal = dEpadCTotal + dEpadC[i][j];
	        hBraggC->Fill(i,dEpadC[i][j]);
	        hPadXC->Fill(i,XpadC[i][j]);
	      }
	    }
	  }
          /*
	  for (UShort_t i=0;i<3;i++){
	    if (dEpadCTotal > 1. && fSiE[1][i+1]>10. && fSiE[1][i+1]<4000.){
	      hTpc_ch_dE[i+9]->Fill(fSiEcal[1][i+1],dEpadCTotal);	
	    }
	  }
          */
        } // end if: flag_detail

      } // end if: flag_tpc

   } // close loop on jentry

   delete ssd_strip_calib;
   delete ssd_pad_calib;
   delete ppac_calib;
} // close Analyzer::Loop
