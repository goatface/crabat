/*
 Copyright: 
 2010, 2011, 2012, 2013 daid KAHL

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
@file run.h
@author daid
@version 1.0
@brief Global flag initalization, histogram memory control and initalization.
@details Basically allocates, initalizes and deletes memory for all histograms.
@date 09 Nov 2011 22:50:31 
*/
//Memory allocation / freeing for crabat
//This appeared to speed up processing time;
//however, no benchmarks were taken
//(previously some of this was done inline in Analyzer.cxx)
//Anyway, it's a little more modular feeling this way

#include <TH2.h>
#include <TH3.h>


/**
@brief Prints the invocation usage to standard error.
@details The function is called whenever we get strange runtime options.\n
Prints to standard error the options accepted by the binary run.\n
*/
void Usage() {
  fprintf(stderr, "\nrun processes the data analysis\n");
  fprintf(stderr, "Must analyze at least one data type (raw and/or detailed)\n");
  fprintf(stderr, "Must analyze at least one detector type (ssd, strip, ppac, tpc)\n");
  fprintf(stderr, "Copyright 2011 daid KAHL (and others, see AUTHORS)\n");
  fprintf(stderr, "Usage: run [-R RunNumber] [-B ] [-H header] [-I path] [-r ] [-d ] [-s ] [-x ] [-p ] [-t ]\n");
  fprintf(stderr, "Description of options:\n");
  fprintf(stderr, "  -R <int>\t: Call run on RunNumber <int> (can be entered manually at runtime also)\n");
  fprintf(stderr, "  -B\t\t: Call run in non-interactive Batch mode.  Requires -R flag.\n");
  fprintf(stderr, "  -H <string>\t: Use run header <string> instead of run header in run.conf\n");
  fprintf(stderr, "  -I <string>\t: Use ROOT Tree input directory path <string> instead of reading from run.conf\n");
  fprintf(stderr, "  -r\t\t: Raw data histograms\n");
  fprintf(stderr, "  -d\t\t: Detailed data histograms\n");
  fprintf(stderr, "  -s\t\t: Analyze SSD pad data\n");
  fprintf(stderr, "  -x\t\t: Analyze SSD strip data\n");
  fprintf(stderr, "  -p\t\t: Analyze PPAC data\n");
  fprintf(stderr, "  -t\t\t: Analyze TPC data\n");
};


// read any arguments passed to run at runtime
// Never change these here!
Bool_t flag_Run = false; // flag_Run controls automated mode (no user input)
Bool_t flag_Batch = false; // controls Batch mode - must be used with -R flag as well
Bool_t flag_Header = false; // pass a header rather than read it from run.conf
Bool_t flag_Inpath = false; // pass the ROOT Tree Input Path (directory) rather than read from run.conf
Bool_t flag_raw=false, flag_detail=false;
Bool_t flag_ssd=false, flag_strip=false, flag_ppac=false, flag_tpc=false; 

//Histograms
TH1F *hRf_ds[2]; 
TH1F *hRf_ssd[2]; 
TH1F *hRfcal[2]; 

TH1F *hSiStripHit ;
TH1F *hSiPadHit ;
TH1F *hfSample[144]; 
TH1F *hTrough[144]; 
TH2F *hBeamPadCh ;
TH2F *hGemMstpcCh ;

TH2F *hSiStripE ;
TH2F *hSiStripEcal ;
TH2F *hSiPadE ;
TH2F *hSiEmax ;
TH2F *hSiPadEcal ;
//Each ssd channel
TH1F *strip_ch[97];
TH1F *strip_cal_ch[97];
TH1F *strip_ch_gated30s[17];
TH1F *strip_ch_gated29p[17];
TH1F *strip_ch_gated17f[17];
TH1F *padE_ch[18];
TH1F *padT_ch[18];
TH1F *padT_cal_ch[18];
TH1F *pad_ch_gated[18];

//DeltaE-E
TH2F *dE_E[6];
TH2F *padE_strip[8];
TH2F *dE_Estrip[6];
TH2F *dE_EBeamGateRF0[6];
TH2F *dE_EBeamGateRF1[6];
TH2F *dE_EBeamGateRFAll ;
TH2F *hSiT ;
TH2F *hSiTcal ;

TH1F *hPpac0TpadT_30s_ch[18];
TH1F *hPpac1TpadT_30s_ch[18];
TH1F *hPpac0TpadT_29p_ch[18];
TH1F *hPpac1TpadT_29p_ch[18];
TH1F *hPpacToF_30s ;
TH1F *hPpacToF_29p ;

TH2F *hRF0Tof ;
TH2F *hRF1Tof ;

TH1F *hPpac0TX1 ;
TH1F *hPpac0TX2 ;
TH1F *hPpac0TY1 ;
TH1F *hPpac0TY2 ;
TH1F *hPpac1TX1 ;
TH1F *hPpac1TX2 ;
TH1F *hPpac1TY1 ;
TH1F *hPpac1TY2 ;

TH2F *hPpac0XY ;
TH2F *hPpac1XY ;
TH2F *hPpac1XYcut ;
TH2F *hTargetXY ;
TH2F *hTargetXZ ;
TH2F *hTargetXYcut ;
TH2F *hTargetXY_30s ;
TH2F *hTargetXYcut_30s ;
TH2F *hTargetXY_29p ;
TH2F *hPpac0XRF0 ;
TH2F *hPpac0XRF0_17f ;
TH2F *hPpac0XRF1_17f ;
TH2F *hPpac0XRF0_29p ;
TH2F *hPpac0XRF1_29p ;
TH2F *hPpac0XRF0_30s ;
TH2F *hPpac0XRF0_30s_ds ;
TH2F *hPpac0XRF1_30s ;
TH2F *hPpac0XRF0ds ;
TH2F *hPpac0XRF1ds ;
TH2F *hPpac0XRF0ssd ;
TH2F *hPpac0XRF1ssd ;
TH2F *hPpac1XRF0 ;
TH2F *hPpac1XRF0_29p ;
TH2F *hPpac1XRF1_29p ;
TH2F *hPpac1XRF0_30s ;
TH2F *hPpac1XRF1_30s ;
TH2F *hPpac1XRF0ds ;
TH2F *hPpac1XRF1ds ;
TH2F *hPpac1XRF0ssd ;
TH2F *hPpac1XRF1ssd ;
TH2F *hPpac0XRF1 ;
TH2F *hPpac0XRF1cut ;
TH2F *hPpac1XRF1cut ;
TH2F *hPpac1XRF1 ;

//TH2F *hGemMstpcAll ;
TH2F *hGemMstpcAll ;
TH2F *hGemMstpcSide ;
TH2F *hPadC ;
TH2F *hTpcPulse[500];

TH2F *hTpc_dE_ch[144];

// low gain pads

TH1F *hBraggB_30s_ssd_ch[48];
TH1F *hBraggB_30s_ds_ch[48];
TH1F *hBraggB_29p_ssd_ch[48];
TH1F *hBraggB_29p_ds_ch[48];
TH1F *hBraggB_17f_ssd_ch[48];
TH1F *hBraggB_17f_ds_ch[48];
TH1F *hBraggB_ds_ch[48];
TH1F *hBraggB_ssd_ch[48];

TH2F *hBraggB_30s_jch[15]; 
TH2F *hBraggB_rp ;
TH2F *hBraggB_30s ;
TH2F *hBraggB_30s_ssd ;
TH2F *hBraggB_30s_ds ;
TH2F *hBraggB_29p ;
TH2F *hBraggB_29p_ssd ;
TH2F *hBraggB_29p_ds ;
TH2F *hBraggB_17f ;
TH2F *hBraggB; 

// pad multiplicity
TH1F *hMultB;
TH1F *hMultB_30s;
TH1F *hMultC_30s_1a;
TH1F *hMultC_30s_2a;
TH1F *hMultL_30s_3a;
TH1F *hMultL_3a;

TH2F *hPadXB ;
TH2F *hPadXB_30s ;
TH2F *hPadXB_29p ;
TH2F *hPadXBgeo ;
TH2F *hPadXBgeo_30s ;
TH2F *hPadXBgeo_29p ;

TH2F *hPadYB ;
TH2F *hPadYBgeo ;
TH2F *hPadYB_30s ;

TH3F *hPpacProj_3D ;

TH3F *hPadB_3D ;

TH3F *hBraggB_3D ;
TH3F *hBraggB_3D_ds ;
TH3F *hBraggB_3D_ssd ;

// high gain pads
TH2F *hBraggL ;
TH2F *hBraggL_E[9] ;
TH2F *hBraggR_E[9] ;
TH2F *hBraggR ;
TH2F *hBraggC ;
TH2F *hBraggC_Eup[9];
TH2F *hBraggC_Edown[9];
TH2F *hBraggC_E[1];
//TH2F *hBraggC ;

//TGraph *gPadXC;
TH2F *hPadXC1asum3;
TH2F *hPadXC1asum4;
TH2F *hPadXC1asum5;
TH2F *hPadXC2asum;
TH2F *hPadXC1a[100] ;
TH2F *hPadXC2a[100] ;
TH2F *hPadZL[100] ;
TH2F *hPadZR[100] ;
TH2F *hPadYC1a;
TH2F *hPadYC2a;
TH2F *hPadYL;
TH2F *hPadYR;

/**
@brief Initalize the histograms based on which data are to be processed
@details This function must be called before any histograms are Fill'd
All memory allocation for histograms is completed previously in the run header\n
Here we define the histograms and set all relevant histogramming options.
*/
void HistInit(){ // histogram initalization
   
  char name[100];
  int xbin,ybin;
  double xmin,xmax,ymin,ymax;
  if (flag_ssd){ 
    
    if (flag_raw){ 
       //hSiPadHit goes here
      hSiPadHit = new TH1F("hSiPadHit","hSiPadHit",18,-0.5,17.5);
      hSiT = new TH2F("hSiT","hSiT",18,-0.5,17.5, 30000,0.,30000.);
      hSiT->SetOption("COL");
      hSiPadE = new TH2F("hSiPadE","hSiPadE",18,-0.5,17.5, 410,0.,4100.);
      hSiPadE->SetOption("COL");   
      for (UShort_t i=0;i<18;i++) {
        sprintf(name,"padE_ch%d",i);
        padE_ch[i]= new TH1F(name, name, 4000,0,4000);
        sprintf(name,"padT_ch%d",i);
        padT_ch[i]= new TH1F(name, name, 300,0,30000);
      } 
    } // end if: flag_raw
    
    if (flag_detail){
      hSiTcal = new TH2F("hSiTcal","hSiTcal",18,-0.5,17.5, 30000,0.,3000.);
      hSiTcal->SetOption("COL");
      hSiTcal->GetXaxis()->SetTitle("Detector number");
      hSiTcal->GetXaxis()->CenterTitle(true);
      hSiTcal->GetYaxis()->SetTitle("Time (ns)");
      hSiTcal->GetYaxis()->CenterTitle(true);
      hSiPadEcal = new TH2F("hSiPadEcal","hSiPadEcal",18,-0.5,17.5, 410,0.,4100.);
      for (UShort_t i=0;i<18;i++) {
        sprintf(name,"padT_cal_ch%d",i);
        padT_cal_ch[i]= new TH1F(name, name, 2000,0,2000);
      } 
    } // end if: flag_detail
    
  } // end if: flag_ssd
  
  if (flag_strip) {
  
    if (flag_raw){
      hSiStripHit = new TH1F("hSiStripHit","hSiStripHit",96,-0.5,95.5);
      hSiStripE = new TH2F("hSiStripE","hSiStripE",96,.5,96.5, 410,0.,4100.);
      hSiStripE->SetOption("COL");
      for (UShort_t i=1;i<97;i++) {
        sprintf(name,"strip_ch%d",i);
        strip_ch[i]= new TH1F(name, name, 4096,0,4096);
      }
    } // end if: flag_raw
  
    if (flag_detail){
      hSiStripEcal = new TH2F("hSiStripEcal","hSiStripEcal",96,-0.5,95.5, 400,0.,20.);
      hSiEmax = new TH2F("hSiEmax","hSiEmax",12,-0.5,11.5, 410,0.,4100.);
      hSiEmax->SetOption("COL");
      for (UShort_t i=1;i<97;i++) {
        sprintf(name,"strip_cal_ch%d",i);
        strip_cal_ch[i]= new TH1F(name, name, 300,0,30);
      }
    } // end if: flag_detail

  } // end if: flag_strip
  
  if (flag_detail && flag_ssd && flag_strip){
    for (UShort_t i=0;i<6;i++) {
      sprintf(name,"dE_E%d",i);
      dE_E[i]= new TH2F(name, name, 300,0.,30., 1000,0.,100.);
      dE_E[i]->SetOption("COL");
      sprintf(name,"dE_Estrip%d",i);
      dE_Estrip[i]= new TH2F(name, name, 300,0.,30., 1000,0.,100.);
      dE_Estrip[i]->SetOption("COL");
    }
    for (UShort_t i=0;i<8;i++) {
      sprintf(name,"padE_strip%d",i);
      padE_strip[i]= new TH2F(name, name, 500,100,600, 500,100,600);
      padE_strip[i]->SetOption("COL");
    }

  } // end if: flag_detail & flag_ssd & flag_strip

  if (flag_ppac){
  
    if(flag_raw){
      for (UShort_t i=0;i<2;i++){
        sprintf(name,"hRf_ds%d",i);
        hRf_ds[i] = new TH1F(name,name,1000,0,1000);
        sprintf(name,"hRf_ssd%d",i);
        hRf_ssd[i] = new TH1F(name,name,1000,0,1000);
      }
      hPpac0TX1 = new TH1F("hPpac0TX1","PPACa X1",23000,0.,23000.);
      hPpac0TX1->GetXaxis()->SetTitle("Raw timing (ch)");
      hPpac0TX1->GetXaxis()->CenterTitle(true);
      hPpac0TX1->GetYaxis()->SetTitle("Counts");
      hPpac0TX1->GetYaxis()->CenterTitle(true);
      hPpac0TX2 = new TH1F("hPpac0TX2","PPACa X2",23000,0.,23000.);
      hPpac0TX2->GetXaxis()->SetTitle("Raw timing (ch)");
      hPpac0TX2->GetXaxis()->CenterTitle(true);
      hPpac0TX2->GetYaxis()->SetTitle("Counts");
      hPpac0TX2->GetYaxis()->CenterTitle(true);
      hPpac0TY1 = new TH1F("hPpac0TY1","PPACa Y1",23000,0.,23000.);
      hPpac0TY1->GetXaxis()->SetTitle("Raw timing (ch)");
      hPpac0TY1->GetXaxis()->CenterTitle(true);
      hPpac0TY1->GetYaxis()->SetTitle("Counts");
      hPpac0TY1->GetYaxis()->CenterTitle(true);
      hPpac0TY2 = new TH1F("hPpac0TY2","PPACa Y2",23000,0.,23000.);
      hPpac0TY2->GetXaxis()->SetTitle("Raw timing (ch)");
      hPpac0TY2->GetXaxis()->CenterTitle(true);
      hPpac0TY2->GetYaxis()->SetTitle("Counts");
      hPpac0TY2->GetYaxis()->CenterTitle(true);
      hPpac1TX1 = new TH1F("hPpac1TX1","PPACb X1",23000,0.,23000.);
      hPpac1TX1->GetXaxis()->SetTitle("Raw timing (ch)");
      hPpac1TX1->GetXaxis()->CenterTitle(true);
      hPpac1TX1->GetYaxis()->SetTitle("Counts");
      hPpac1TX1->GetYaxis()->CenterTitle(true);
      hPpac1TX2 = new TH1F("hPpac1TX2","PPACb X2",23000,0.,23000.);
      hPpac1TX2->GetXaxis()->SetTitle("Raw timing (ch)");
      hPpac1TX2->GetXaxis()->CenterTitle(true);
      hPpac1TX2->GetYaxis()->SetTitle("Counts");
      hPpac1TX2->GetYaxis()->CenterTitle(true);
      hPpac1TY1 = new TH1F("hPpac1TY1","PPACb Y1",23000,0.,23000.);
      hPpac1TY1->GetXaxis()->SetTitle("Raw timing (ch)");
      hPpac1TY1->GetXaxis()->CenterTitle(true);
      hPpac1TY1->GetYaxis()->SetTitle("Counts");
      hPpac1TY1->GetYaxis()->CenterTitle(true);
      hPpac1TY2 = new TH1F("hPpac1TY2","PPACb Y2",23000,0.,23000.);
      hPpac1TY2->GetXaxis()->SetTitle("Raw timing (ch)");
      hPpac1TY2->GetXaxis()->CenterTitle(true);
      hPpac1TY2->GetYaxis()->SetTitle("Counts");
      hPpac1TY2->GetYaxis()->CenterTitle(true);
    } // end if: flag_raw
    
    if (flag_detail){
      for (UShort_t i=0;i<2;i++){
        sprintf(name,"hRfcal%d",i);
        hRfcal[i] = new TH1F(name,name,300,0.,60.);
        hRfcal[i]->GetXaxis()->SetTitle("RF (ns)");
	hRfcal[i]->GetXaxis()->CenterTitle(true);
        hRfcal[i]->GetYaxis()->SetTitle("counts");
	hRfcal[i]->GetYaxis()->CenterTitle(true);
      }
      hRF0Tof = new TH2F("hRF0Tof","RF0 vs. PPAC Tof",300,0.,60.,200,0.,20.);
      hRF0Tof->SetOption("COL");
      hRF0Tof->GetXaxis()->SetTitle("RF (ns)");
      hRF0Tof->GetXaxis()->CenterTitle(true);
      hRF0Tof->GetYaxis()->SetTitle("Time of Flight (ns)");
      hRF0Tof->GetYaxis()->CenterTitle(true);
      
      hRF1Tof = new TH2F("hRF1Tof","RF1 vs. PPAC Tof",300,0.,60.,200,0.,20.);
      hRF1Tof->SetOption("COL");
      hRF1Tof->GetXaxis()->SetTitle("RF (ns)");
      hRF1Tof->GetXaxis()->CenterTitle(true);
      hRF1Tof->GetYaxis()->SetTitle("Time of Flight (ns)");
      hRF1Tof->GetYaxis()->CenterTitle(true);
      
      hPpac0XY = new TH2F("hPpac0XY","PPACa X vs. Y",160,-40.,40.,160,-40.,40.);
      hPpac0XY->SetOption("COL");
      hPpac0XY->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XY->GetXaxis()->CenterTitle(true);
      hPpac0XY->GetYaxis()->SetTitle("Y Position (mm)");
      hPpac0XY->GetYaxis()->CenterTitle(true);
      
      hPpac1XY = new TH2F("hPpac1XY","PPACb X vs. Y",160,-40.,40.,160,-40.,40.);
      hPpac1XY->SetOption("COL");   
      hPpac1XY->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XY->GetXaxis()->CenterTitle(true);
      hPpac1XY->GetYaxis()->SetTitle("Y Position (mm)");
      hPpac1XY->GetYaxis()->CenterTitle(true);
      
      hPpac1XYcut = new TH2F("hPpac1XYcut","PPACb X vs. Y",160,-40.,40.,160,-40.,40.);
      hPpac1XYcut->SetOption("COL");   
      hPpac1XYcut->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XYcut->GetXaxis()->CenterTitle(true);
      hPpac1XYcut->GetYaxis()->SetTitle("Y Position (mm)");
      hPpac1XYcut->GetYaxis()->CenterTitle(true);
      
      hTargetXY = new TH2F("hTargetXY","Target X vs. Y",80,-40.,40.,80,-40.,40.);
      hTargetXY->SetOption("COL");   
      hTargetXY->GetXaxis()->SetTitle("X Position (mm)");
      hTargetXY->GetXaxis()->CenterTitle(true);
      hTargetXY->GetYaxis()->SetTitle("Y Position (mm)");
      hTargetXY->GetYaxis()->CenterTitle(true);
      
      hTargetXZ = new TH2F("hTargetXZ","Target X vs. Z",160,-40.,40.,48,-0.5,47.5);
      hTargetXZ->SetOption("COL");   
      hTargetXZ->GetXaxis()->SetTitle("X Position (mm)");
      hTargetXZ->GetXaxis()->CenterTitle(true);
      hTargetXZ->GetYaxis()->SetTitle("Pad Number");
      hTargetXZ->GetYaxis()->CenterTitle(true);
      
      hTargetXYcut = new TH2F("hTargetXYcut","Target X vs. Y window cut",320,-40.,40.,320,-40.,40.);
      hTargetXYcut->SetOption("COL");   
      hTargetXYcut->GetXaxis()->SetTitle("X Position (mm)");
      hTargetXYcut->GetXaxis()->CenterTitle(true);
      hTargetXYcut->GetYaxis()->SetTitle("Y Position (mm)");
      hTargetXYcut->GetYaxis()->CenterTitle(true);
      
      hTargetXY_30s = new TH2F("hTargetXY_30s","^{30}S Target X vs. Y",320,-40.,40.,320,-40.,40.);
      hTargetXY_30s->SetOption("COL");   
      hTargetXY_30s->GetXaxis()->SetTitle("X Position (mm)");
      hTargetXY_30s->GetXaxis()->CenterTitle(true);
      hTargetXY_30s->GetYaxis()->SetTitle("Y Position (mm)");
      hTargetXY_30s->GetYaxis()->CenterTitle(true);
      
      hTargetXYcut_30s = new TH2F("hTargetXYcut_30s","^{30}S Target X vs. Y window cut",320,-40.,40.,320,-40.,40.);
      hTargetXYcut_30s->SetOption("COL");   
      hTargetXYcut_30s->GetXaxis()->SetTitle("X Position (mm)");
      hTargetXYcut_30s->GetXaxis()->CenterTitle(true);
      hTargetXYcut_30s->GetYaxis()->SetTitle("Y Position (mm)");
      hTargetXYcut_30s->GetYaxis()->CenterTitle(true);
      
      hTargetXY_29p = new TH2F("hTargetXY_29p","^{29}P Target X vs. Y",320,-40.,40.,320,-40.,40.);
      hTargetXY_29p->SetOption("COL");   
      hTargetXY_29p->GetXaxis()->SetTitle("X Position (mm)");
      hTargetXY_29p->GetXaxis()->CenterTitle(true);
      hTargetXY_29p->GetYaxis()->SetTitle("Y Position (mm)");
      hTargetXY_29p->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0 = new TH2F("hPpac0XRF0","PPACa X vs. RF0",320,-40.,40.,300,0.,60.);
      hPpac0XRF0->SetOption("COL");
      hPpac0XRF0->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0->GetXaxis()->CenterTitle(true);
      hPpac0XRF0->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF0->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1 = new TH2F("hPpac0XRF1","PPACa X vs. RF1",320,-40.,40.,300,0.,60.);
      hPpac0XRF1->SetOption("COL");
      hPpac0XRF1->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1->GetXaxis()->CenterTitle(true);
      hPpac0XRF1->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF1->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0 = new TH2F("hPpac1XRF0","PPACb X vs. RF0",320,-40.,40.,300,0.,60.);
      hPpac1XRF0->SetOption("COL");
      hPpac1XRF0->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0->GetXaxis()->CenterTitle(true);
      hPpac1XRF0->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF0->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0_29p = new TH2F("hPpac0XRF0_29p","PPACa X vs. RF0 29P",160,-40.,40.,300,0.,60.);
      hPpac0XRF0_29p->SetOption("COL");
      hPpac0XRF0_29p->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0_29p->GetXaxis()->CenterTitle(true);
      hPpac0XRF0_29p->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF0_29p->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1_29p = new TH2F("hPpac0XRF1_29p","PPACa X vs. RF1 29P",160,-40.,40.,300,0.,60.);
      hPpac0XRF1_29p->SetOption("COL");
      hPpac0XRF1_29p->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1_29p->GetXaxis()->CenterTitle(true);
      hPpac0XRF1_29p->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF1_29p->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0_17f = new TH2F("hPpac0XRF0_17f","PPACa X vs. RF0 17F",160,-40.,40.,300,0.,60.);
      hPpac0XRF0_17f->SetOption("COL");
      hPpac0XRF0_17f->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0_17f->GetXaxis()->CenterTitle(true);
      hPpac0XRF0_17f->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF0_17f->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1_17f = new TH2F("hPpac0XRF1_17f","PPACa X vs. RF1 17F",160,-40.,40.,300,0.,60.);
      hPpac0XRF1_17f->SetOption("COL");
      hPpac0XRF1_17f->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1_17f->GetXaxis()->CenterTitle(true);
      hPpac0XRF1_17f->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF1_17f->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0_30s = new TH2F("hPpac0XRF0_30s","PPACa X vs. RF0 30S",160,-40.,40.,300,0.,60.);
      hPpac0XRF0_30s->SetOption("COL");
      hPpac0XRF0_30s->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0_30s->GetXaxis()->CenterTitle(true);
      hPpac0XRF0_30s->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF0_30s->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0_30s_ds = new TH2F("hPpac0XRF0_30s_ds","PPACa X vs. RF0 30S DS",160,-40.,40.,300,0.,60.);
      hPpac0XRF0_30s_ds->SetOption("COL");
      hPpac0XRF0_30s_ds->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0_30s_ds->GetXaxis()->CenterTitle(true);
      hPpac0XRF0_30s_ds->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF0_30s_ds->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1_30s = new TH2F("hPpac0XRF1_30s","PPACa X vs. RF1 30S",160,-40.,40.,300,0.,60.);
      hPpac0XRF1_30s->SetOption("COL");
      hPpac0XRF1_30s->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1_30s->GetXaxis()->CenterTitle(true);
      hPpac0XRF1_30s->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF1_30s->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0_29p = new TH2F("hPpac1XRF0_29p","PPACb X vs. RF0 29P",160,-40.,40.,300,0.,60.);
      hPpac1XRF0_29p->SetOption("COL");
      hPpac1XRF0_29p->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0_29p->GetXaxis()->CenterTitle(true);
      hPpac1XRF0_29p->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF0_29p->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1_29p = new TH2F("hPpac1XRF1_29p","PPACb X vs. RF1 29P",160,-40.,40.,300,0.,60.);
      hPpac1XRF1_29p->SetOption("COL");
      hPpac1XRF1_29p->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1_29p->GetXaxis()->CenterTitle(true);
      hPpac1XRF1_29p->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF1_29p->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0_30s = new TH2F("hPpac1XRF0_30s","PPACb X vs. RF0 30S",160,-40.,40.,300,0.,60.);
      hPpac1XRF0_30s->SetOption("COL");
      hPpac1XRF0_30s->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0_30s->GetXaxis()->CenterTitle(true);
      hPpac1XRF0_30s->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF0_30s->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1_30s = new TH2F("hPpac1XRF1_30s","PPACb X vs. RF1 30S",160,-40.,40.,300,0.,60.);
      hPpac1XRF1_30s->SetOption("COL");
      hPpac1XRF1_30s->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1_30s->GetXaxis()->CenterTitle(true);
      hPpac1XRF1_30s->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF1_30s->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1cut = new TH2F("hPpac0XRF1cut","PPACa X vs. RF1",160,-40.,40.,300,0.,60.);
      hPpac0XRF1cut->SetOption("COL");
      hPpac0XRF1cut->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1cut->GetXaxis()->CenterTitle(true);
      hPpac0XRF1cut->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF1cut->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1cut = new TH2F("hPpac1XRF1cut","PPACb X vs. RF1",160,-40.,40.,300,0.,60.);
      hPpac1XRF1cut->SetOption("COL");
      hPpac1XRF1cut->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1cut->GetXaxis()->CenterTitle(true);
      hPpac1XRF1cut->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF1cut->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1 = new TH2F("hPpac1XRF1","PPACb X vs. RF1",320,-40.,40.,300,0.,60.);
      hPpac1XRF1->SetOption("COL");
      hPpac1XRF1->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1->GetXaxis()->CenterTitle(true);
      hPpac1XRF1->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF1->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0ds = new TH2F("hPpac0XRF0ds","PPACa X vs. RF0 downscale",160,-40.,40.,300,0.,60.);
      hPpac0XRF0ds->SetOption("COL");
      hPpac0XRF0ds->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0ds->GetXaxis()->CenterTitle(true);
      hPpac0XRF0ds->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF0ds->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1ds = new TH2F("hPpac0XRF1ds","PPACa X vs. RF1 downscale",160,-40.,40.,300,0.,60.);
      hPpac0XRF1ds->SetOption("COL");
      hPpac0XRF1ds->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1ds->GetXaxis()->CenterTitle(true);
      hPpac0XRF1ds->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF1ds->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0ssd = new TH2F("hPpac0XRF0ssd","PPACa X vs. RF0 SSD-OR",160,-40.,40.,300,0.,60.);
      hPpac0XRF0ssd->SetOption("COL");
      hPpac0XRF0ssd->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0ssd->GetXaxis()->CenterTitle(true);
      hPpac0XRF0ssd->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF0ssd->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1ssd = new TH2F("hPpac0XRF1ssd","PPACa X vs. RF1 SSD-OR",160,-40.,40.,300,0.,60.);
      hPpac0XRF1ssd->SetOption("COL");
      hPpac0XRF1ssd->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1ssd->GetXaxis()->CenterTitle(true);
      hPpac0XRF1ssd->GetYaxis()->SetTitle("RF (ns)");
      hPpac0XRF1ssd->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0ds = new TH2F("hPpac1XRF0ds","PPACb X vs. RF0 downscale",160,-40.,40.,300,0.,60.);
      hPpac1XRF0ds->SetOption("COL");
      hPpac1XRF0ds->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0ds->GetXaxis()->CenterTitle(true);
      hPpac1XRF0ds->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF0ds->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1ds = new TH2F("hPpac1XRF1ds","PPACb X vs. RF1 downscale",160,-40.,40.,300,0.,60.);
      hPpac1XRF1ds->SetOption("COL");
      hPpac1XRF1ds->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1ds->GetXaxis()->CenterTitle(true);
      hPpac1XRF1ds->GetYaxis()->SetTitle("RF (ns");
      hPpac1XRF1ds->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0ssd = new TH2F("hPpac1XRF0ssd","PPACb X vs. RF0 SSD-OR",160,-40.,40.,300,0.,60.);
      hPpac1XRF0ssd->SetOption("COL");
      hPpac1XRF0ssd->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0ssd->GetXaxis()->CenterTitle(true);
      hPpac1XRF0ssd->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF0ssd->GetYaxis()->CenterTitle(true);

      hPpac1XRF1ssd = new TH2F("hPpac1XRF1ssd","PPACb X vs. RF1 SSD-OR",160,-40.,40.,300,0.,60.);
      hPpac1XRF1ssd->SetOption("COL");
      hPpac1XRF1ssd->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1ssd->GetXaxis()->CenterTitle(true);
      hPpac1XRF1ssd->GetYaxis()->SetTitle("RF (ns)");
      hPpac1XRF1ssd->GetYaxis()->CenterTitle(true);

      hPpacToF_30s = new TH1F("hPpacToF_30s","hPpacToF_30s",200,0.,20.);
      hPpacToF_29p = new TH1F("hPpacToF_29p","hPpacToF_29p",200,0.,20.);
      
      hPpacProj_3D = new TH3F("hPpacProj_3D","hPpacProj_3D",10,10,30.,10,-10.,10.,60,200,800);
      hPpacProj_3D->SetOption("GLCOL");
    } // end if: flag_detail
  
  } // end if: flag_ppac
  
  if (flag_detail && flag_ssd && flag_ppac){
    for (UShort_t i=0;i<18;i++) {
      sprintf(name,"hPpac0TpadT_30s_ch%d",i);
      hPpac0TpadT_30s_ch[i]= new TH1F(name, name, 15000,0.,1500.);
      sprintf(name,"hPpac1TpadT_30s_ch%d",i);
      hPpac1TpadT_30s_ch[i]= new TH1F(name, name, 15000,0.,1500.);
      sprintf(name,"hPpac0TpadT_29p_ch%d",i);
      hPpac0TpadT_29p_ch[i]= new TH1F(name, name, 15000,0,1500.);
      sprintf(name,"hPpac1TpadT_29p_ch%d",i);
      hPpac1TpadT_29p_ch[i]= new TH1F(name, name, 15000,0,1500.);
      sprintf(name,"pad_ch_gated%d",i);
      pad_ch_gated[i]= new TH1F(name, name, 800,0,4100);
    }
  } // end if: flag_detail & flag_ssd & flag_ppac
  
  if (flag_ssd && flag_tpc ){
    //wtf why is this existing and empty?
  } // end if flag_ssd & flag_tpc
  
  if (flag_detail && flag_ssd && flag_strip && flag_ppac){
    for (UShort_t i=1;i<17;i++) {
      sprintf(name,"strip_ch_gated30s%d",i);     
      strip_ch_gated30s[i]= new TH1F(name, name, 4096,0,4096);
      sprintf(name,"strip_ch_gated29p%d",i);     
      strip_ch_gated29p[i]= new TH1F(name, name, 4096,0,4096);
      sprintf(name,"strip_ch_gated17f%d",i);     
      strip_ch_gated17f[i]= new TH1F(name, name, 4096,0,4096);
    }
  } // end if: flag_detail & flag_strip & flag_ppac
  
  if (flag_detail && flag_ssd && flag_strip && flag_ppac){
    dE_EBeamGateRFAll = new TH2F("dE_EBeamGateRFAll", "dE_EBeamGateRFAll", 205,0,4100, 205,0,4100);
    dE_EBeamGateRFAll->SetOption("COL");
    for (UShort_t i=0;i<6;i++) {
      sprintf(name,"dE_EBeamGateRF0_%d",i);
      dE_EBeamGateRF0[i]= new TH2F(name, name, 300,0,30, 300,0,30);
      dE_EBeamGateRF0[i]->SetOption("COL");
      sprintf(name,"dE_EBeamGateRF1_%d",i);
      dE_EBeamGateRF1[i]= new TH2F(name, name, 300,0,30, 300,0,30);
      dE_EBeamGateRF1[i]->SetOption("COL");
    }
  } // end if: flag_detail & flag_ssd & flag_strip & flag_ppac
  
  if (flag_tpc){
  
    if (flag_raw){
      hGemMstpcAll = new TH2F("hGemMstpcAll","GEM-MSTPC fadc all channels/events",200,0.,500.,300,0.,2000.);
      hGemMstpcAll->SetOption("COLZ");
      hGemMstpcSide = new TH2F("hGemMstpcSide","GEM-MSTPC fadc high gain channels/events",200,0.,500.,300,0.,2000.);
      hGemMstpcSide->SetOption("COLZ");
      hBeamPadCh = new TH2F("hBeamPadCh","hBeamPadCh",96,0.,96.,2000,0.,2000.);
      hGemMstpcCh = new TH2F("hGemMstpcCh","hGemMstpcCh",144,0.,144.,2000,0.,2000.);
      for (UShort_t i=0;i<144;i++) {
        sprintf(name,"hfSample%d",i);
        hfSample[i] = new TH1F(name,name,1000,0,1500);
        sprintf(name,"hTpc_dE_ch%d",i);
        hTpc_dE_ch[i]= new TH2F(name, name, 300,0,300, 200,0,2000);
        hTpc_dE_ch[i]->SetOption("COLZ");
      }
    } // end if: flag_raw
  
    if (flag_detail){ 
    
      xbin=8;xmin=0.5;xmax=8.5;
      hMultL_3a = new TH1F("hMultL_3a","Left Pad Multiplicity Gated SSD-OR",xbin,xmin,xmax);
      hMultL_3a->GetXaxis()->SetTitle("No. Pads Hit");
      hMultL_3a->GetXaxis()->CenterTitle(true);
      
      xbin=48;xmin=-0.5;xmax=47.5;ybin=3000;ymin=0.;ymax=30000.;
      hBraggB = new TH2F("hBraggB","Beam Bragg Curves",xbin,xmin,xmax,ybin,ymin,ymax);
      hBraggB->SetOption("COLZ");
      hBraggB->GetXaxis()->SetTitle("Pad No.");
      hBraggB->GetXaxis()->CenterTitle(true);
      hBraggB->GetYaxis()->SetTitle("#DeltaE (arb.)");
      hBraggB->GetYaxis()->CenterTitle(true);
       
      xbin=8;xmin=-0.5;xmax=7.5;ybin=300;ymin=0.;ymax=3000.;
      hBraggL = new TH2F("hBraggL","hBraggL",xbin,xmin,xmax,ybin,ymin,ymax);
      hBraggL->SetOption("COLZ");
      hBraggR = new TH2F("hBraggR","hBraggR",xbin,xmin,xmax,ybin,ymin,ymax);
      hBraggR->SetOption("COLZ");
      hBraggC = new TH2F("hBraggC","hBraggC",xbin,xmin,xmax,ybin,ymin,ymax);
      hBraggC->SetOption("COLZ");
      
      xbin=8;xmin=-0.5;xmax=7.5;ybin=114;ymin=-40.;ymax=74.;
      hPadXC1asum3 = new TH2F("hPadXC1asum3","hPadXC1asum3",xbin,xmin,xmax,ybin,ymin,ymax);
      hPadXC1asum3->SetOption("COLZ");
      hPadXC1asum4 = new TH2F("hPadXC1asum4","hPadXC1asum4",xbin,xmin,xmax,ybin,ymin,ymax);
      hPadXC1asum4->SetOption("COLZ");
      hPadXC1asum5 = new TH2F("hPadXC1asum5","hPadXC1asum5",xbin,xmin,xmax,ybin,ymin,ymax);
      hPadXC1asum5->SetOption("COLZ");
      hPadXC2asum = new TH2F("hPadXC2asum","hPadXC2asum",xbin,xmin,xmax,ybin,ymin,ymax);
      hPadXC2asum->SetOption("COLZ");
      for (UShort_t i=0;i<100;i++) {
        sprintf(name,"hPadXC1a%d",i);
        hPadXC1a[i] = new TH2F(name,name,xbin,xmin,xmax,ybin,ymin,ymax);
        hPadXC1a[i]->SetOption("COLZ");
	sprintf(name,"hPadXC2a%d",i);
        hPadXC2a[i] = new TH2F(name,name,xbin,xmin,xmax,ybin,ymin,ymax);
        hPadXC2a[i]->SetOption("COLZ");
      }
      for (UShort_t i=0;i<100;i++) {
        sprintf(name,"hPadZL%d",i);
        hPadZL[i] = new TH2F(name,name,8,-0.5,7.5,234,-117,117.);
        hPadZL[i]->SetOption("COLZ");
        sprintf(name,"hPadZR%d",i);
        hPadZR[i] = new TH2F(name,name,8,-0.5,7.5,234,-117,117.);
        hPadZR[i]->SetOption("COLZ");
      }
      hPadYC1a = new TH2F("hPadYC1a","hPadYC1a",8,-0.5,7.5,100,0.,300.);
      hPadYC1a->SetOption("COLZ");
      hPadYC2a = new TH2F("hPadYC2a","hPadYC2a",8,-0.5,7.5,100,0.,300.);
      hPadYC2a->SetOption("COLZ");
      hPadYL = new TH2F("hPadYL","hPadYL",8,-0.5,7.5,100,0.,300.);
      hPadYL->SetOption("COLZ");
      hPadYR = new TH2F("hPadYR","hPadYR",8,-0.5,7.5,100,0.,300.);
      hPadYR->SetOption("COLZ");

      hPadXB = new TH2F("hPadXB","Beam X Residual",48,-0.5,47.5,500,-50.,50.);
      hPadXB->SetOption("COLZ");
      hPadXB->GetXaxis()->SetTitle("Pad No.");
      hPadXB->GetXaxis()->CenterTitle(true);
      hPadXB->GetYaxis()->SetTitle("#DeltaX position (mm)");
      hPadXB->GetYaxis()->CenterTitle(true);
      
      hPadYB = new TH2F("hPadYB","Beam Y",48,-0.5,47.5,300,0.,300.);
      hPadYB->SetOption("COLZ");
      hPadYB->GetXaxis()->SetTitle("Pad No.");
      hPadYB->GetXaxis()->CenterTitle(true);
      hPadYB->GetYaxis()->SetTitle("Y position (arb.)");
      hPadYB->GetYaxis()->CenterTitle(true);
      
      hPadB_3D = new TH3F("hPadB_3D","hPadB_3D",3,-2.,1.,10,70.,115.,48,-0.5,47.5);
      hPadB_3D->SetOption("GLCOL");

      for (UShort_t i=0;i<144;i++) {
        sprintf(name,"hTrough%d",i);
        hTrough[i] = new TH1F(name,name,1500,-500.,1000.);
      } 

    } // end if: flag_detail
  
  } // end if: flag_tpc

  if (flag_ssd && flag_tpc && flag_detail){
      xbin=100;xmin=0;xmax=10000;ybin=400;ymin=0.;ymax=1600.;
        sprintf(name,"hBraggC_E0");
        hBraggC_E[0] = new TH2F(name,name,xbin,xmin,xmax,ybin,ymin,ymax);
        hBraggC_E[0] ->SetOption("COLZ");
      for (UShort_t i=1;i<9;i++) {
        sprintf(name,"hBraggC_Edown%d",i);
        hBraggC_Edown[i] = new TH2F(name,name,xbin,xmin,xmax,ybin,ymin,ymax);
        hBraggC_Edown[i] ->SetOption("COLZ");
        sprintf(name,"hBraggC_Eup%d",i);
        hBraggC_Eup[i] = new TH2F(name,name,xbin,xmin,xmax,ybin,ymin,ymax);
        hBraggC_Eup[i] ->SetOption("COLZ");
      }
      for (UShort_t i=0;i<9;i++) {
        sprintf(name,"hBraggL_E%d",i);
        hBraggL_E[i] = new TH2F(name,name,xbin,xmin,xmax,ybin,ymin,ymax);
        hBraggL_E[i] ->SetOption("COLZ");
      }
      for (UShort_t i=0;i<9;i++) {
        sprintf(name,"hBraggR_E%d",i);
        hBraggR_E[i] = new TH2F(name,name,xbin,xmin,xmax,ybin,ymin,ymax);
        hBraggR_E[i] ->SetOption("COLZ");
      }
  
  }
  if (flag_ppac && flag_tpc && flag_detail){
    
    xbin=48;xmin=0.5;xmax=48.5;
    hMultB = new TH1F("hMultB","Beam Pad Multiplicity Gated on Beam",xbin,xmin,xmax);
    hMultB->GetXaxis()->SetTitle("No. Pads Hit");
    hMultB->GetXaxis()->CenterTitle(true);
    
    hMultB_30s = new TH1F("hMultB_30s","Beam Pad Multiplicity Gated on ^{30}S",xbin,xmin,xmax);
    hMultB_30s->GetXaxis()->SetTitle("No. Pads Hit");
    hMultB_30s->GetXaxis()->CenterTitle(true);

    xbin=8;xmin=0.5;xmax=8.5;
    
    hMultC_30s_1a = new TH1F("hMultC_30s_1a","Center Pad Multiplicity Gated on ^{30}S and 1a",xbin,xmin,xmax);
    hMultC_30s_1a->GetXaxis()->SetTitle("No. Pads Hit");
    hMultC_30s_1a->GetXaxis()->CenterTitle(true);

    hMultC_30s_2a = new TH1F("hMultC_30s_2a","Center Pad Multiplicity Gated on ^{30}S and 2a",xbin,xmin,xmax);
    hMultC_30s_2a->GetXaxis()->SetTitle("No. Pads Hit");
    hMultC_30s_2a->GetXaxis()->CenterTitle(true);

    hMultL_30s_3a = new TH1F("hMultL_30s_3a","Left Pad Multiplicity Gated on ^{30}S and 3a",xbin,xmin,xmax);
    hMultL_30s_3a->GetXaxis()->SetTitle("No. Pads Hit");
    hMultL_30s_3a->GetXaxis()->CenterTitle(true);

    //xbin=48;xmin=-0.5;xmax=47.5;ybin=200;ymin=0.;ymax=2000.;
    xbin=48;xmin=-0.5;xmax=47.5;ybin=3000;ymin=0.;ymax=30000.;
    hBraggB_30s = new TH2F("hBraggB_30s","Beam Bragg Curve Gated on ^{30}S",xbin,xmin,xmax,ybin,ymin,ymax);
    hBraggB_30s->SetOption("COLZ");
    hBraggB_30s->GetXaxis()->SetTitle("Pad No.");
    hBraggB_30s->GetXaxis()->CenterTitle(true);
    hBraggB_30s->GetYaxis()->SetTitle("#DeltaE (arb.)");
    hBraggB_30s->GetYaxis()->CenterTitle(true);
    
    hBraggB_29p = new TH2F("hBraggB_29p","hBraggB_29p",xbin,xmin,xmax,ybin,ymin,ymax);
    hBraggB_29p->SetOption("COLZ");
    hBraggB_17f = new TH2F("hBraggB_17f","hBraggB_17f",xbin,xmin,xmax,ybin,ymin,ymax);
    hBraggB_17f->SetOption("COLZ");

    hPadYBgeo = new TH2F("hPadYBgeo","Beam Y",48,-0.5,47.5,80,-20.,20.);
    hPadYBgeo->SetOption("COLZ");
    hPadYBgeo->GetXaxis()->SetTitle("Pad No.");
    hPadYBgeo->GetXaxis()->CenterTitle(true);
    hPadYBgeo->GetYaxis()->SetTitle("Y position projected (mm)");
    hPadYBgeo->GetYaxis()->CenterTitle(true);
      
    hPadXBgeo = new TH2F("hPadXBgeo","Beam X Residual (corrected)",48,-0.5,47.5,500,-50.,50.);
    hPadXBgeo->SetOption("COLZ");
    hPadXBgeo->GetXaxis()->SetTitle("Pad No.");
    hPadXBgeo->GetXaxis()->CenterTitle(true);
    hPadXBgeo->GetYaxis()->SetTitle("#DeltaX position (mm)");
    hPadXBgeo->GetYaxis()->CenterTitle(true);
      
    hPadXBgeo_30s = new TH2F("hPadXBgeo_30s","^{30}S Beam X Residual (corrected)",48,-0.5,47.5,500,-50.,50.);
    hPadXBgeo_30s->SetOption("COLZ");
    hPadXBgeo_30s->GetXaxis()->SetTitle("Pad No.");
    hPadXBgeo_30s->GetXaxis()->CenterTitle(true);
    hPadXBgeo_30s->GetYaxis()->SetTitle("#DeltaX position (mm)");
    hPadXBgeo_30s->GetYaxis()->CenterTitle(true);
      
    hPadXBgeo_29p = new TH2F("hPadXBgeo_29p","^{29}P Beam X Residual (corrected)",48,-0.5,47.5,500,-50.,50.);
    hPadXBgeo_29p->SetOption("COLZ");
    hPadXBgeo_29p->GetXaxis()->SetTitle("Pad No.");
    hPadXBgeo_29p->GetXaxis()->CenterTitle(true);
    hPadXBgeo_29p->GetYaxis()->SetTitle("#DeltaX position (mm)");
    hPadXBgeo_29p->GetYaxis()->CenterTitle(true);
      
    hPadXB_30s = new TH2F("hPadXB_30s","^{30}S Beam X Residual",48,-0.5,47.5,500,-50.,50.);
    hPadXB_30s->SetOption("COLZ");
    hPadXB_30s->GetXaxis()->SetTitle("Pad No.");
    hPadXB_30s->GetXaxis()->CenterTitle(true);
    hPadXB_30s->GetYaxis()->SetTitle("#DeltaX position (mm)");
    hPadXB_30s->GetYaxis()->CenterTitle(true);
    
    hPadXB_29p = new TH2F("hPadXB_29p","^{29}P Beam X Residual",48,-0.5,47.5,500,-50.,50.);
    hPadXB_29p->SetOption("COLZ");
    hPadXB_29p->GetXaxis()->SetTitle("Pad No.");
    hPadXB_29p->GetXaxis()->CenterTitle(true);
    hPadXB_29p->GetYaxis()->SetTitle("#DeltaX position (mm)");
    hPadXB_29p->GetYaxis()->CenterTitle(true);
    
    hPadYB_30s = new TH2F("hPadYB_30s","Beam Y Gated on ^{30}S",48,-0.5,47.5,300,0.,300.);
    hPadYB_30s->SetOption("COLZ");
    hPadYB_30s->GetXaxis()->SetTitle("Pad No.");
    hPadYB_30s->GetXaxis()->CenterTitle(true);
    hPadYB_30s->GetYaxis()->SetTitle("Y position (cm)");
    hPadYB_30s->GetYaxis()->CenterTitle(true);
    
    for (UShort_t j=0;j<15;j++){
      sprintf(name,"hBraggB_30s_jch%d",j);
      hBraggB_30s_jch[j] = new TH2F(name,name,xbin,xmin,xmax,ybin,ymin,ymax);
      hBraggB_30s_jch[j]->SetOption("COLZ");
    }
    
    for (UShort_t j=0;j<500;j++){
      sprintf(name,"hTpcPulse%d",j);
      hTpcPulse[j] = new TH2F(name,name,33,0.,100.,600,400.,1000.);
      //hTpcPulse[j] = new TH2F(name,name,66,0.,200.,600,400.,2000.);
      hTpcPulse[j]->GetXaxis()->SetTitle("fClock (20 ns / bin)");
      hTpcPulse[j]->GetXaxis()->CenterTitle(true);
      hTpcPulse[j]->GetYaxis()->SetTitle("fSample (arb)");
      hTpcPulse[j]->GetYaxis()->CenterTitle(true);
    }
    
    hBraggB_3D = new TH3F("hBraggB_3D","hBraggB_3D",48,-0.5,47.5,4,-0.5,3.5,400,0.,2000.);
    hBraggB_3D->SetOption("GLCOL");
    hBraggB_3D_ds = new TH3F("hBraggB_3D_ds","hBraggB_3D_ds",48,-0.5,47.5,4,-0.5,3.5,400,0.,2000.);
    hBraggB_3D_ds->SetOption("GLCOL");
    hBraggB_3D_ssd = new TH3F("hBraggB_3D_ssd","hBraggB_3D_ssd",48,-0.5,47.5,4,-0.5,3.5,400,0.,2000.);
    hBraggB_3D_ssd->SetOption("GLCOL");
      
    hBraggB_30s_ds = new TH2F("hBraggB_30s_ds","Beam Bragg Curve Gated on ^{30}S DS",xbin,xmin,xmax,ybin,ymin,ymax);
    hBraggB_30s_ds->SetOption("COLZ");
    hBraggB_30s_ds->GetXaxis()->SetTitle("Pad No.");
    hBraggB_30s_ds->GetXaxis()->CenterTitle(true);
    hBraggB_30s_ds->GetYaxis()->SetTitle("#DeltaE (arb.)");
    hBraggB_30s_ds->GetYaxis()->CenterTitle(true);
    
    hBraggB_30s_ssd = new TH2F("hBraggB_30s_ssd","Beam Bragg Curve Gated on ^{30}S SSD-OR",xbin,xmin,xmax,ybin,ymin,ymax);
    hBraggB_30s_ssd->SetOption("COLZ");
    hBraggB_30s_ssd->GetXaxis()->SetTitle("Pad No.");
    hBraggB_30s_ssd->GetXaxis()->CenterTitle(true);
    hBraggB_30s_ssd->GetYaxis()->SetTitle("#DeltaE (arb.)");
    hBraggB_30s_ssd->GetYaxis()->CenterTitle(true);
    
    hBraggB_29p_ds = new TH2F("hBraggB_29p_ds","Beam Bragg Curve Gated on ^{29}P DS",xbin,xmin,xmax,ybin,ymin,ymax);
    hBraggB_29p_ds->SetOption("COLZ");
    hBraggB_29p_ds->GetXaxis()->SetTitle("Pad No.");
    hBraggB_29p_ds->GetXaxis()->CenterTitle(true);
    hBraggB_29p_ds->GetYaxis()->SetTitle("#DeltaE (arb.)");
    hBraggB_29p_ds->GetYaxis()->CenterTitle(true);
    
    hBraggB_29p_ssd = new TH2F("hBraggB_29p_ssd","Beam Bragg Curve Gated on ^{29}P SSD-OR",xbin,xmin,xmax,ybin,ymin,ymax);
    hBraggB_29p_ssd->SetOption("COLZ");
    hBraggB_29p_ssd->GetXaxis()->SetTitle("Pad No.");
    hBraggB_29p_ssd->GetXaxis()->CenterTitle(true);
    hBraggB_29p_ssd->GetYaxis()->SetTitle("#DeltaE (arb.)");
    hBraggB_29p_ssd->GetYaxis()->CenterTitle(true);
    
    hBraggB_rp = new TH2F("hBraggB_rp","Beam Scattering Location for ^{30}S",xbin,xmin,xmax,ybin,ymin,ymax);
    hBraggB_rp->SetOption("COLZ");
    hBraggB_rp->GetXaxis()->SetTitle("Pad No.");
    hBraggB_rp->GetXaxis()->CenterTitle(true);
    hBraggB_rp->GetYaxis()->SetTitle("(#DeltaE_{i}-#DeltaE_{i-1})/Pad i");
    hBraggB_rp->GetYaxis()->CenterTitle(true);
      
    for (UShort_t i=0;i<48;i++) { // Beam Bragg pads
      sprintf(name,"hBraggB_ds_ch%d",i);
      hBraggB_ds_ch[i] = new TH1F(name,name,2000,0,2000);
      hBraggB_ds_ch[i]->GetXaxis()->SetTitle("#DeltaE (arb.)");
      hBraggB_ds_ch[i]->GetXaxis()->CenterTitle(true);
      hBraggB_ds_ch[i]->GetYaxis()->SetTitle("Counts");
      hBraggB_ds_ch[i]->GetYaxis()->CenterTitle(true);
      sprintf(name,"hBraggB_ssd_ch%d",i);
      hBraggB_ssd_ch[i] = new TH1F(name,name,2000,0,2000);
      hBraggB_ssd_ch[i]->GetXaxis()->SetTitle("#DeltaE (arb.)");
      hBraggB_ssd_ch[i]->GetXaxis()->CenterTitle(true);
      hBraggB_ssd_ch[i]->GetYaxis()->SetTitle("Counts");
      hBraggB_ssd_ch[i]->GetYaxis()->CenterTitle(true);
    }
    for (UShort_t i=0;i<48;i++) {
      sprintf(name,"hBraggB_30s_ssd_ch%d",i);
      hBraggB_30s_ssd_ch[i] = new TH1F(name,name,ybin,ymin,ymax);
      sprintf(name,"hBraggB_30s_ds_ch%d",i);
      hBraggB_30s_ds_ch[i] = new TH1F(name,name,ybin,ymin,ymax);
      sprintf(name,"hBraggB_29p_ssd_ch%d",i);
      hBraggB_29p_ssd_ch[i] = new TH1F(name,name,ybin,ymin,ymax);
      sprintf(name,"hBraggB_29p_ds_ch%d",i);
      hBraggB_29p_ds_ch[i] = new TH1F(name,name,ybin,ymin,ymax);
      sprintf(name,"hBraggB_17f_ssd_ch%d",i);
      hBraggB_17f_ssd_ch[i] = new TH1F(name,name,ybin,ymin,ymax);
      sprintf(name,"hBraggB_17f_ds_ch%d",i);
      hBraggB_17f_ds_ch[i] = new TH1F(name,name,ybin,ymin,ymax);
    }

  } // end if: flag_tpc & flag_ppac & flag_detail

} // close HistInit

/**
@brief Write all the histograms.
@details This function should be called after all histograms have been Fill'd in Analyzer.\n
It will write all histograms that were initalized based on the data options.\n
*/
void HistWrite() {
  // WRITE THE HISTOGRAMS

  if (flag_ssd) { // process SSD?
    if (flag_raw) { // write the SSD pad raw histograms
      hSiT->Write(); 
      hSiPadE->Write(); 
      hSiPadHit->Write();  
      for (UShort_t i=0;i<18;i++) {
         padE_ch[i]->Write();
         padT_ch[i]->Write();
      } // end for : i
    } // end if: flag_raw
    
    if (flag_detail) { // write the SSD pad detailed histograms
      hSiTcal->Write(); 
      hSiPadEcal->Write();
      for (UShort_t i=0;i<18;i++) {
         padT_cal_ch[i]->Write();
      }
    } // end if: flag_detail
  } //end if: flag_ssd

  if(flag_strip){
    if (flag_raw){ // write the SSD strip raw histograms
      hSiStripE->Write(); 
      hSiStripHit->Write(); 
      for (UShort_t i=1;i<97;i++) {
         strip_ch[i]->Write();
       }
    } // end if: flag_raw
    
    if (flag_detail){ // write the SSD strip detailed histograms
      hSiEmax->Write(); 
      for (UShort_t i=1;i<97;i++) {
         strip_cal_ch[i]->Write();
      }
    } // end if : flag_detail
  } //end if: flag_strip

  if (flag_detail && flag_ssd && flag_strip){
      for (UShort_t i=0;i<6;i++) {
         dE_E[i]->Write();
         dE_Estrip[i]->Write();
      }
      for (UShort_t i=0;i<8;i++) {
         padE_strip[i]->Write();
      }
  } // end if: flag_detail & flag_ssd & flag_strip

  if (flag_ppac){ // process PPAC?
    if (flag_raw) { // write the PPAC raw histograms
       for (UShort_t i=0;i<2;i++) hRf_ssd[i]->Write();
       for (UShort_t i=0;i<2;i++) hRf_ds[i]->Write();
       hPpac0TX1->Write();
       hPpac0TX2->Write();
       hPpac0TY1->Write();
       hPpac0TY2->Write();
       hPpac1TX1->Write();
       hPpac1TX2->Write();
       hPpac1TY1->Write();
       hPpac1TY2->Write();
    } // end if: flag_raw

    if (flag_detail) { // write the PPAC detailed histograms
       for (UShort_t i=0;i<2;i++) hRfcal[i]->Write();
       hRF0Tof->Write();
       hRF1Tof->Write();
       hPpac0XY->Write();
       hPpac1XY->Write();  
       hPpac1XYcut->Write();  
       hTargetXY->Write();  
       hTargetXZ->Write();  
       hTargetXYcut->Write();  
       hTargetXY_30s->Write();  
       hTargetXYcut_30s->Write();  
       hTargetXY_29p->Write();  
       hPpac0XRF0->Write(); 
       hPpac1XRF0->Write(); 
       hPpac0XRF1->Write(); 
       hPpac0XRF1cut->Write(); 
       hPpac1XRF1cut->Write(); 
       hPpac1XRF1->Write();
       hPpac0XRF0_17f->Write(); 
       hPpac0XRF1_17f->Write(); 
       hPpac0XRF0_29p->Write(); 
       hPpac0XRF1_29p->Write(); 
       hPpac0XRF0_30s->Write(); 
       hPpac0XRF0_30s_ds->Write(); 
       hPpac0XRF1_30s->Write(); 
       hPpac1XRF0_29p->Write(); 
       hPpac1XRF1_29p->Write(); 
       hPpac1XRF0_30s->Write(); 
       hPpac1XRF1_30s->Write(); 
       hPpac0XRF0ds->Write(); 
       hPpac0XRF1ds->Write(); 
       hPpac0XRF0ssd->Write(); 
       hPpac0XRF1ssd->Write(); 
       hPpac1XRF0ds->Write(); 
       hPpac1XRF1ds->Write(); 
       hPpac1XRF0ssd->Write(); 
       hPpac1XRF1ssd->Write(); 
       hPpacToF_30s->Write();
       hPpacToF_29p->Write();
       hPpacProj_3D->Write();
    } // end if: flag_detail
  } // end if: flag_ppac

  if (flag_ssd && flag_ppac){ // process PPAC & SSD?
    if (flag_detail) { // write PPAC / SSD detailed histograms
       for (UShort_t i=0;i<18;i++){
         hPpac0TpadT_30s_ch[i]->Write();
         hPpac1TpadT_30s_ch[i]->Write();
         hPpac0TpadT_29p_ch[i]->Write();
         hPpac1TpadT_29p_ch[i]->Write();
       }
    } // end if: flag_detail
  } // end if: flag_ssd & flag_ppac
  if (flag_ssd && flag_tpc && flag_detail){
    for (int i=0;i<9;i++) {
      if (i==0) hBraggC_E[i]->Write();
      else {
      hBraggC_Eup[i]->Write();
      hBraggC_Edown[i]->Write();
      }
      hBraggL_E[i]->Write();
      hBraggR_E[i]->Write();
    }
  }
  if (flag_detail && flag_ssd && flag_strip && flag_ppac){ // process PPAC & SSD strip details?
    for (int i=1;i<17;i++){
      strip_ch_gated30s[i]->Write();
      strip_ch_gated29p[i]->Write();
      strip_ch_gated17f[i]->Write();
    } // enf for: i
  } // end if: flag_detail && flag_strip && flag_ppac
  if (flag_tpc) { // process TPC?
    if (flag_raw) { // write TPC raw histograms
       hGemMstpcAll->Write();
       hGemMstpcSide->Write();
       hGemMstpcCh->Write();
       hBeamPadCh->Write();
       for (UShort_t i=0;i<144;i++) {
          hfSample[i]->Write();
          hTpc_dE_ch[i]->Write();
       }
    } // end if: flag_raw
    if (flag_detail) { // write TPC detailed histograms
      hMultL_3a->Write();
      //hPadC->Write();
      hBraggL->Write();
      hBraggR->Write();
      hBraggC->Write();
      hPadXC1asum3->Write();
      hPadXC1asum4->Write();
      hPadXC1asum5->Write();
      hPadXC2asum->Write();
      for (UShort_t i=0;i<100;i++) {
        hPadXC1a[i]->Write();
        hPadXC2a[i]->Write();
      }
      for (UShort_t i=0;i<100;i++) {
        hPadZL[i]->Write();
        hPadZR[i]->Write();
      }
      hPadYC1a->Write();
      hPadYC2a->Write();
      hPadYL->Write();
      hPadYR->Write();
      if (flag_ppac){
        hBraggB->Write();
        hPadXB->Write();
        hPadYB->Write();
        hPadB_3D->Write();
        for (UShort_t i=0;i<144;i++) {
           hTrough[i]->Write();
           if (i<48){
             hBraggB_30s_ssd_ch[i]->Write();
             hBraggB_30s_ds_ch[i]->Write();
             hBraggB_29p_ssd_ch[i]->Write();
             hBraggB_29p_ds_ch[i]->Write();
             hBraggB_17f_ssd_ch[i]->Write();
             hBraggB_17f_ds_ch[i]->Write();
           }
        }
      } // end if: flag_ppac
    } // end if: flag_detail
  } // end if: flag_tpc
  
  if (flag_ppac && flag_tpc && flag_detail){
    hPadXB_30s->Write();
    hPadXB_29p->Write();
    hPadYBgeo->Write();
    hPadXBgeo->Write();
    hPadXBgeo_30s->Write();
    hPadXBgeo_29p->Write();
    hPadYB_30s->Write();
    hBraggB_30s->Write();
    hMultB->Write();
    hMultB_30s->Write();
    hMultC_30s_1a->Write();
    hMultC_30s_2a->Write();
    hMultL_30s_3a->Write();
    //hBraggB_30s2->Write();
    hBraggB_29p->Write();
    hBraggB_17f->Write();
    for (UShort_t i=0;i<15;i++) hBraggB_30s_jch[i]->Write();
    hBraggB_3D->Write();

    hBraggB_30s_ds->Write();
    hBraggB_30s_ssd->Write();
    hBraggB_29p_ds->Write();
    hBraggB_29p_ssd->Write();
    hBraggB_rp->Write();
    hBraggB_3D_ds->Write();
    hBraggB_3D_ssd->Write();
    for (UShort_t i=0;i<48;i++){
      hBraggB_ds_ch[i]->Write();
      hBraggB_ssd_ch[i]->Write();
    }
    for (UShort_t i=0;i<500;i++){
      hTpcPulse[i]->Write();
    }
  } // end if: flag_ppac && flag_tpc && flag_detail

} // close HistWrite


/**
@brief delete all the histograms
@details This function should be called after all histograms have been Write'n.\n
In principle it is not totally necessary to do this kind of garbage collection.\n
However, for memory checking with valgrind, this is helpful.
*/
void HistClean(){ // function to free histogram memory
  delete hSiStripHit ;
  delete hSiPadHit ;
  delete hBeamPadCh ;
  delete hGemMstpcCh ;
  delete hSiStripE ;
  delete hSiStripEcal ;
  delete hSiPadE ;
  delete hSiEmax ;
  delete hSiPadEcal ;
  delete dE_EBeamGateRFAll ;
  delete hSiT ;
  delete hSiTcal ;
  delete hRF0Tof;
  delete hRF1Tof;
  delete hPpac0TX1;
  delete hPpac0TX2;
  delete hPpac0TY1;
  delete hPpac0TY2;
  delete hPpac1TX1;
  delete hPpac1TX2;
  delete hPpac1TY1;
  delete hPpac1TY2;
  delete hPpac0XY ;
  delete hPpac1XY ;
  delete hPpac1XYcut ;
  delete hTargetXY ;
  delete hTargetXZ ;
  delete hTargetXYcut ;
  delete hTargetXY_29p ;
  delete hTargetXY_30s ;
  delete hTargetXYcut_30s ;
  delete hPpac0XRF0 ;
  delete hPpac0XRF0_17f ;
  delete hPpac0XRF1_17f ;
  delete hPpac0XRF0_29p ;
  delete hPpac0XRF1_29p ;
  delete hPpac0XRF0_30s ;
  delete hPpac0XRF0_30s_ds ;
  delete hPpac0XRF1_30s ;
  delete hPpac0XRF0ds ;
  delete hPpac0XRF1ds ;
  delete hPpac0XRF0ssd ;
  delete hPpac0XRF1ssd ;
  delete hPpac1XRF0 ;
  delete hPpac1XRF0_29p ;
  delete hPpac1XRF1_29p ;
  delete hPpac1XRF0_30s ;
  delete hPpac1XRF1_30s ;
  delete hPpac1XRF0ds ;
  delete hPpac1XRF1ds ;
  delete hPpac1XRF0ssd ;
  delete hPpac1XRF1ssd ;
  delete hPpac0XRF1 ;
  delete hPpac0XRF1cut ;
  delete hPpac1XRF1cut ;
  delete hPpac1XRF1 ;
  delete hPpacProj_3D ;
  delete hGemMstpcAll ;
  delete hGemMstpcSide ;
  delete hBraggB_30s ;
  delete hMultB;
  delete hMultB_30s;
  delete hMultC_30s_1a;
  delete hMultC_30s_2a;
  delete hMultL_30s_3a;
  delete hMultL_3a;
  delete hBraggB_30s_ds ;
  delete hBraggB_30s_ssd ;
  delete hPadXB ;
  delete hPadXB_30s ;
  delete hPadXB_29p ;
  delete hPadXBgeo;
  delete hPadXBgeo_30s;
  delete hPadXBgeo_29p;
  delete hPadYB ;
  delete hPadYB_30s ;
  delete hBraggB_3D ;
  delete hBraggB_3D_ds ;
  delete hBraggB_3D_ssd ;
  delete hPadB_3D;
  delete hBraggB_29p ;
  delete hBraggB_29p_ssd ;
  delete hBraggB_29p_ds ;
  delete hBraggB_17f ;
  delete hBraggB ;
  delete hBraggL ;
  delete hBraggR ;
  delete hBraggC ;
    for (int i=0;i<9;i++){ 
      if (i==0) delete hBraggC_E[i] ;
      else {
      delete hBraggC_Eup[i] ;
      delete hBraggC_Edown[i] ;
      }
      delete hBraggL_E[i] ;
      delete hBraggR_E[i] ;
    }
  delete hBraggB_rp;
  for (int i=0;i<2;i++) delete hRf_ds[i]; 
  for (int i=0;i<2;i++) delete hRf_ssd[i]; 
  for (int i=0;i<2;i++) delete hRfcal[i]; 
  for (int j=0;j<15;j++) delete hBraggB_30s_jch[j];
  for (int j=0;j<500;j++) delete hTpcPulse[j];
  for (int i=0;i<97;i++){
    delete strip_ch[i];
    delete strip_cal_ch[i];
  }
  for (int i=0;i<17;i++){
    delete strip_ch_gated30s[i];
    delete strip_ch_gated29p[i];
    delete strip_ch_gated17f[i];
  }
  delete hPpacToF_30s;
  delete hPpacToF_29p;
  for (int i=0;i<18;i++){
    delete padE_ch[i];
    delete padT_ch[i];
    delete padT_cal_ch[i];
    delete hPpac0TpadT_30s_ch[i];
    delete hPpac0TpadT_29p_ch[i];
    delete pad_ch_gated[i];
  }
  for (int i=0;i<6;i++){
    delete dE_E[i];
    delete dE_Estrip[i];
    delete dE_EBeamGateRF0[i];
    delete dE_EBeamGateRF1[i];
  }
  for (int i=0;i<8;i++){
    delete padE_strip[i];
  }
  for (int i=0;i<144;i++){
    delete hfSample[i]; 
    delete hTrough[i];
    delete hTpc_dE_ch[i];
  }
  for (int i=0;i<48;i++){
    delete hBraggB_30s_ssd_ch[i];
    delete hBraggB_30s_ds_ch[i];
    delete hBraggB_29p_ssd_ch[i];
    delete hBraggB_29p_ds_ch[i];
    delete hBraggB_17f_ssd_ch[i];
    delete hBraggB_17f_ds_ch[i];
    delete hBraggB_ds_ch[i];
    delete hBraggB_ssd_ch[i];
  }
} // close HistClean    
