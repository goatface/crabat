/*
 Copyright: 
 2010, 2011 daid KAHL

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
TH1F *hRf[2]; 

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
TH1F *strip_ch_gated30s[97];
TH1F *strip_ch_gated29p[97];
TH1F *pad_ch[18];
TH1F *pad_ch_gated[18];

//DeltaE-E
TH2F *dE_E[6];
TH2F *dE_Estrip[6];
TH2F *dE_EBeamGateRF0[6];
TH2F *dE_EBeamGateRF1[6];
TH2F *dE_EBeamGateRFAll ;
TH2F *hSiT ;
TH2F *hSiTcal ;

TH1F *hPpac0TpadT_30s_ch[18];
TH1F *hPpac0TpadT_29p_ch[18];
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
TH2F *hPpac0XRF0 ;
TH2F *hPpac0XRF0_29p ;
TH2F *hPpac0XRF1_29p ;
TH2F *hPpac0XRF0_30s ;
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

TH2F *hTpc_dE_ch[144];

// low gain pads

TH1F *hBraggB_ch[48];
TH1F *hBraggB_ds_ch[48];
TH1F *hBraggB_ssd_ch[48];

TH2F *hBraggB_30s_jch[15]; 
TH2F *hBraggB_rp ;
TH2F *hBraggB_30s ;
TH2F *hBraggB_30s_ssd ;
TH2F *hBraggB_30s_ds ;
TH2F *hBraggB_29p ;
TH2F *hBraggB_31s ;
TH2F *hBraggB; 

TH2F *hPadXB ;
TH2F *hPadXB_30s ;

TH2F *hPadYB ;
TH2F *hPadYB_30s ;

TH3F *hPadB_3D ;

TH3F *hBraggB_3D ;
TH3F *hBraggB_3D_ds ;
TH3F *hBraggB_3D_ssd ;

// high gain pads
TH2F *hBraggL ;
TH2F *hBraggR ;
TH2F *hBraggC ;
//TH2F *hBraggC ;

TH2F *hPadXC ;

/**
@brief Initalize the histograms based on which data are to be processed
@details This function must be called before any histograms are Fill'd
All memory allocation for histograms is completed previously in the run header\n
Here we define the histograms and set all relevant histogramming options.
*/
void HistInit(){ // histogram initalization
   
  char name[100];
  
  if (flag_ssd){ 
    
    if (flag_raw){ 
       //hSiPadHit goes here
      hSiPadHit = new TH1F("hSiPadHit","hSiPadHit",18,-0.5,17.5);
      hSiT = new TH2F("hSiT","hSiT",18,-0.5,17.5, 30000,0.,30000.);
      hSiT->SetOption("COL");
      hSiPadE = new TH2F("hSiPadE","hSiPadE",18,-0.5,17.5, 410,0.,4100.);
      hSiPadE->SetOption("COL");   
      for (UShort_t i=0;i<18;i++) {
        sprintf(name,"pad_ch%d",i);
        pad_ch[i]= new TH1F(name, name, 3950,50,4000);
      } 
    } // end if: flag_raw
    
    if (flag_detail){
      hSiTcal = new TH2F("hSiTcal","hSiTcal",18,-0.5,17.5, 30000,0.,30000.);
      hSiTcal->SetOption("COL");
      hSiPadEcal = new TH2F("hSiPadEcal","hSiPadEcal",18,-0.5,17.5, 410,0.,4100.);
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

  } // end if: flag_detail & flag_ssd & flag_strip

  if (flag_ppac){
  
    if(flag_raw){
      for (UShort_t i=0;i<2;i++){
        sprintf(name,"hRf%d",i);
        hRf[i] = new TH1F(name,name,1000,0,1000);
      }
      hPpac0TX1 = new TH1F("hPpac0TX1","hPpac0TX1",22000,0.,22000.);
      hPpac0TX2 = new TH1F("hPpac0TX2","hPpac0TX2",22000,0.,22000.);
      hPpac0TY1 = new TH1F("hPpac0TY1","hPpac0TY1",22000,0.,22000.);
      hPpac0TY2 = new TH1F("hPpac0TY2","hPpac0TY2",22000,0.,22000.);
      hPpac1TX1 = new TH1F("hPpac1TX1","hPpac1TX1",22000,0.,22000.);
      hPpac1TX2 = new TH1F("hPpac1TX2","hPpac1TX2",22000,0.,22000.);
      hPpac1TY1 = new TH1F("hPpac1TY1","hPpac1TY1",22000,0.,22000.);
      hPpac1TY2 = new TH1F("hPpac1TY2","hPpac1TY2",22000,0.,22000.);
    } // end if: flag_raw
    
    if (flag_detail){
      hRF0Tof = new TH2F("hRF0Tof","RF0 vs. PPAC Tof",600,0.,600.,200,0.,20.);
      hRF0Tof->SetOption("COL");
      hRF0Tof->GetXaxis()->SetTitle("RF (arb)");
      hRF0Tof->GetXaxis()->CenterTitle(true);
      hRF0Tof->GetYaxis()->SetTitle("Time of Flight (ns)");
      hRF0Tof->GetYaxis()->CenterTitle(true);
      
      hRF1Tof = new TH2F("hRF1Tof","RF1 vs. PPAC Tof",600,0.,600.,200,0.,20.);
      hRF1Tof->SetOption("COL");
      hRF1Tof->GetXaxis()->SetTitle("RF (arb)");
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
      
      hPpac0XRF0 = new TH2F("hPpac0XRF0","PPACa X vs. RF0",160,-40.,40.,600,0.,600.);
      hPpac0XRF0->SetOption("COL");
      hPpac0XRF0->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0->GetXaxis()->CenterTitle(true);
      hPpac0XRF0->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF0->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1 = new TH2F("hPpac0XRF1","PPACa X vs. RF1",160,-40.,40.,600,0.,600.);
      hPpac0XRF1->SetOption("COL");
      hPpac0XRF1->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1->GetXaxis()->CenterTitle(true);
      hPpac0XRF1->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF1->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0 = new TH2F("hPpac1XRF0","PPACb X vs. RF0",160,-40.,40.,600,0.,600.);
      hPpac1XRF0->SetOption("COL");
      hPpac1XRF0->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0->GetXaxis()->CenterTitle(true);
      hPpac1XRF0->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF0->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0_29p = new TH2F("hPpac0XRF0_29p","PPACa X vs. RF0 29P",160,-40.,40.,600,0.,600.);
      hPpac0XRF0_29p->SetOption("COL");
      hPpac0XRF0_29p->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0_29p->GetXaxis()->CenterTitle(true);
      hPpac0XRF0_29p->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF0_29p->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1_29p = new TH2F("hPpac0XRF1_29p","PPACa X vs. RF1 29P",160,-40.,40.,600,0.,600.);
      hPpac0XRF1_29p->SetOption("COL");
      hPpac0XRF1_29p->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1_29p->GetXaxis()->CenterTitle(true);
      hPpac0XRF1_29p->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF1_29p->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0_30s = new TH2F("hPpac0XRF0_30s","PPACa X vs. RF0 30S",160,-40.,40.,600,0.,600.);
      hPpac0XRF0_30s->SetOption("COL");
      hPpac0XRF0_30s->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0_30s->GetXaxis()->CenterTitle(true);
      hPpac0XRF0_30s->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF0_30s->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1_30s = new TH2F("hPpac0XRF1_30s","PPACa X vs. RF1 30S",160,-40.,40.,600,0.,600.);
      hPpac0XRF1_30s->SetOption("COL");
      hPpac0XRF1_30s->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1_30s->GetXaxis()->CenterTitle(true);
      hPpac0XRF1_30s->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF1_30s->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0_29p = new TH2F("hPpac1XRF0_29p","PPACb X vs. RF0 29P",160,-40.,40.,600,0.,600.);
      hPpac1XRF0_29p->SetOption("COL");
      hPpac1XRF0_29p->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0_29p->GetXaxis()->CenterTitle(true);
      hPpac1XRF0_29p->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF0_29p->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1_29p = new TH2F("hPpac1XRF1_29p","PPACb X vs. RF1 29P",160,-40.,40.,600,0.,600.);
      hPpac1XRF1_29p->SetOption("COL");
      hPpac1XRF1_29p->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1_29p->GetXaxis()->CenterTitle(true);
      hPpac1XRF1_29p->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF1_29p->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0_30s = new TH2F("hPpac1XRF0_30s","PPACb X vs. RF0 30S",160,-40.,40.,600,0.,600.);
      hPpac1XRF0_30s->SetOption("COL");
      hPpac1XRF0_30s->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0_30s->GetXaxis()->CenterTitle(true);
      hPpac1XRF0_30s->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF0_30s->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1_30s = new TH2F("hPpac1XRF1_30s","PPACb X vs. RF1 30S",160,-40.,40.,600,0.,600.);
      hPpac1XRF1_30s->SetOption("COL");
      hPpac1XRF1_30s->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1_30s->GetXaxis()->CenterTitle(true);
      hPpac1XRF1_30s->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF1_30s->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1cut = new TH2F("hPpac0XRF1cut","PPACa X vs. RF1",160,-40.,40.,600,0.,600.);
      hPpac0XRF1cut->SetOption("COL");
      hPpac0XRF1cut->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1cut->GetXaxis()->CenterTitle(true);
      hPpac0XRF1cut->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF1cut->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1cut = new TH2F("hPpac1XRF1cut","PPACb X vs. RF1",160,-40.,40.,600,0.,600.);
      hPpac1XRF1cut->SetOption("COL");
      hPpac1XRF1cut->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1cut->GetXaxis()->CenterTitle(true);
      hPpac1XRF1cut->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF1cut->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1 = new TH2F("hPpac1XRF1","PPACb X vs. RF1",160,-40.,40.,600,0.,600.);
      hPpac1XRF1->SetOption("COL");
      hPpac1XRF1->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1->GetXaxis()->CenterTitle(true);
      hPpac1XRF1->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF1->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0ds = new TH2F("hPpac0XRF0ds","PPACa X vs. RF0 downscale",160,-40.,40.,600,0.,600.);
      hPpac0XRF0ds->SetOption("COL");
      hPpac0XRF0ds->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0ds->GetXaxis()->CenterTitle(true);
      hPpac0XRF0ds->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF0ds->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1ds = new TH2F("hPpac0XRF1ds","PPACa X vs. RF1 downscale",160,-40.,40.,600,0.,600.);
      hPpac0XRF1ds->SetOption("COL");
      hPpac0XRF1ds->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1ds->GetXaxis()->CenterTitle(true);
      hPpac0XRF1ds->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF1ds->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF0ssd = new TH2F("hPpac0XRF0ssd","PPACa X vs. RF0 SSD-OR",160,-40.,40.,600,0.,600.);
      hPpac0XRF0ssd->SetOption("COL");
      hPpac0XRF0ssd->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF0ssd->GetXaxis()->CenterTitle(true);
      hPpac0XRF0ssd->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF0ssd->GetYaxis()->CenterTitle(true);
      
      hPpac0XRF1ssd = new TH2F("hPpac0XRF1ssd","PPACa X vs. RF1 SSD-OR",160,-40.,40.,600,0.,600.);
      hPpac0XRF1ssd->SetOption("COL");
      hPpac0XRF1ssd->GetXaxis()->SetTitle("X Position (mm)");
      hPpac0XRF1ssd->GetXaxis()->CenterTitle(true);
      hPpac0XRF1ssd->GetYaxis()->SetTitle("RF (arb)");
      hPpac0XRF1ssd->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0ds = new TH2F("hPpac1XRF0ds","PPACb X vs. RF0 downscale",160,-40.,40.,600,0.,600.);
      hPpac1XRF0ds->SetOption("COL");
      hPpac1XRF0ds->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0ds->GetXaxis()->CenterTitle(true);
      hPpac1XRF0ds->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF0ds->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF1ds = new TH2F("hPpac1XRF1ds","PPACb X vs. RF1 downscale",160,-40.,40.,600,0.,600.);
      hPpac1XRF1ds->SetOption("COL");
      hPpac1XRF1ds->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1ds->GetXaxis()->CenterTitle(true);
      hPpac1XRF1ds->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF1ds->GetYaxis()->CenterTitle(true);
      
      hPpac1XRF0ssd = new TH2F("hPpac1XRF0ssd","PPACb X vs. RF0 SSD-OR",160,-40.,40.,600,0.,600.);
      hPpac1XRF0ssd->SetOption("COL");
      hPpac1XRF0ssd->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF0ssd->GetXaxis()->CenterTitle(true);
      hPpac1XRF0ssd->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF0ssd->GetYaxis()->CenterTitle(true);

      hPpac1XRF1ssd = new TH2F("hPpac1XRF1ssd","PPACb X vs. RF1 SSD-OR",160,-40.,40.,600,0.,600.);
      hPpac1XRF1ssd->SetOption("COL");
      hPpac1XRF1ssd->GetXaxis()->SetTitle("X Position (mm)");
      hPpac1XRF1ssd->GetXaxis()->CenterTitle(true);
      hPpac1XRF1ssd->GetYaxis()->SetTitle("RF (arb)");
      hPpac1XRF1ssd->GetYaxis()->CenterTitle(true);

      hPpacToF_30s = new TH1F("hPpacToF_30s","hPpacToF_30s",200,0.,20.);
      hPpacToF_29p = new TH1F("hPpacToF_29p","hPpacToF_29p",200,0.,20.);
    } // end if: flag_detail
  
  } // end if: flag_ppac
  
  if (flag_detail && flag_ssd && flag_ppac){
    for (UShort_t i=0;i<18;i++) {
      sprintf(name,"hPpac0TpadT_30s_ch%d",i);
      hPpac0TpadT_30s_ch[i]= new TH1F(name, name, 15000,0.,1500.);
      sprintf(name,"hPpac0TpadT_29p_ch%d",i);
      hPpac0TpadT_29p_ch[i]= new TH1F(name, name, 15000,0,1500.);
      sprintf(name,"pad_ch_gated%d",i);
      pad_ch_gated[i]= new TH1F(name, name, 800,0,4100);
    }
  } // end if: flag_detail & flag_ssd & flag_ppac
  
  if (flag_ssd && flag_tpc ){
    //wtf why is this existing and empty?
  } // end if flag_ssd & flag_tpc
  
  if (flag_detail && flag_strip && flag_ppac){
    for (UShort_t i=1;i<97;i++) {
      sprintf(name,"strip_ch_gated30s%d",i);     
      strip_ch_gated30s[i]= new TH1F(name, name, 4096,0,4096);
      sprintf(name,"strip_ch_gated29p%d",i);     
      strip_ch_gated29p[i]= new TH1F(name, name, 4096,0,4096);
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
        hTpc_dE_ch[i]= new TH2F(name, name, 200,0,200, 200,0,2000);
        hTpc_dE_ch[i]->SetOption("COLZ");
      }
    } // end if: flag_raw
  
    if (flag_detail){ 
      hBraggB = new TH2F("hBraggB","Beam Bragg Curves",48,-0.5,47.5,200,0.,2000.);
      hBraggB->SetOption("COLZ");
      hBraggB->GetXaxis()->SetTitle("Pad No.");
      hBraggB->GetXaxis()->CenterTitle(true);
      hBraggB->GetYaxis()->SetTitle("#DeltaE (arb.)");
      hBraggB->GetYaxis()->CenterTitle(true);
       
      hBraggL = new TH2F("hBraggL","hBraggL",8,-0.5,7.5,200,0.,2000.);
      hBraggL->SetOption("COLZ");
      hBraggR = new TH2F("hBraggR","hBraggR",8,-0.5,7.5,200,0.,2000.);
      hBraggR->SetOption("COLZ");
      hBraggC = new TH2F("hBraggC","hBraggC",8,-0.5,7.5,2000,0.,2000.);
      hBraggC->SetOption("COLZ");
      hPadXC = new TH2F("hPadXC","hPadXC",8,-0.5,7.5,200,-10.,10.);
      hPadXC->SetOption("COLZ");
      
      hPadXB = new TH2F("hPadXB","Beam X",48,-0.5,47.5,80,-5.,5.);
      hPadXB->SetOption("COLZ");
      hPadXB->GetXaxis()->SetTitle("Pad No.");
      hPadXB->GetXaxis()->CenterTitle(true);
      hPadXB->GetYaxis()->SetTitle("X position (cm)");
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
        if (i<48){
          sprintf(name,"hBraggB_ch%d",i);
          hBraggB_ch[i] = new TH1F(name,name,2000,0,2000);
        }
      }
      
    
    } // end if: flag_detail
  
  } // end if: flag_tpc

  if (flag_ppac && flag_tpc && flag_detail){
    
    hBraggB_30s = new TH2F("hBraggB_30s","Beam Bragg Curve Gated on ^{30}S",48,-0.5,47.5,200,0.,2000.);
    hBraggB_30s->SetOption("COLZ");
    hBraggB_30s->GetXaxis()->SetTitle("Pad No.");
    hBraggB_30s->GetXaxis()->CenterTitle(true);
    hBraggB_30s->GetYaxis()->SetTitle("#DeltaE (arb.)");
    hBraggB_30s->GetYaxis()->CenterTitle(true);
    
    hBraggB_29p = new TH2F("hBraggB_29p","hBraggB_29p",48,-0.5,47.5,200,0.,2000.);
    hBraggB_29p->SetOption("COLZ");
    hBraggB_31s = new TH2F("hBraggB_31s","hBraggB_31s",48,-0.5,47.5,200,0.,2000.);
    hBraggB_31s->SetOption("COLZ");

    hPadXB_30s = new TH2F("hPadXB_30s","Beam X Gated on ^{30}S",48,-0.5,47.5,80,-5.,5.);
    hPadXB_30s->SetOption("COLZ");
    hPadXB_30s->GetXaxis()->SetTitle("Pad No.");
    hPadXB_30s->GetXaxis()->CenterTitle(true);
    hPadXB_30s->GetYaxis()->SetTitle("X position (cm)");
    hPadXB_30s->GetYaxis()->CenterTitle(true);
    
    hPadYB_30s = new TH2F("hPadYB_30s","Beam Y Gated on ^{30}S",48,-0.5,47.5,300,0.,300.);
    hPadYB_30s->SetOption("COLZ");
    hPadYB_30s->GetXaxis()->SetTitle("Pad No.");
    hPadYB_30s->GetXaxis()->CenterTitle(true);
    hPadYB_30s->GetYaxis()->SetTitle("Y position (cm)");
    hPadYB_30s->GetYaxis()->CenterTitle(true);
    
    for (UShort_t j=0;j<15;j++){
      sprintf(name,"hBraggB_30s_jch%d",j);
      hBraggB_30s_jch[j] = new TH2F(name,name,48,-0.5,47.5,200,0.,2000.);
      hBraggB_30s_jch[j]->SetOption("COLZ");
    }
    
    hBraggB_3D = new TH3F("hBraggB_3D","hBraggB_3D",48,-0.5,47.5,4,-0.5,3.5,400,0.,2000.);
    hBraggB_3D->SetOption("GLCOL");
    hBraggB_3D_ds = new TH3F("hBraggB_3D_ds","hBraggB_3D_ds",48,-0.5,47.5,4,-0.5,3.5,400,0.,2000.);
    hBraggB_3D_ds->SetOption("GLCOL");
    hBraggB_3D_ssd = new TH3F("hBraggB_3D_ssd","hBraggB_3D_ssd",48,-0.5,47.5,4,-0.5,3.5,400,0.,2000.);
    hBraggB_3D_ssd->SetOption("GLCOL");
      
    hBraggB_30s_ds = new TH2F("hBraggB_30s_ds","Beam Bragg Curve Gated on ^{30}S DS",48,-0.5,47.5,200,0.,2000.);
    hBraggB_30s_ds->SetOption("COLZ");
    hBraggB_30s_ds->GetXaxis()->SetTitle("Pad No.");
    hBraggB_30s_ds->GetXaxis()->CenterTitle(true);
    hBraggB_30s_ds->GetYaxis()->SetTitle("#DeltaE (arb.)");
    hBraggB_30s_ds->GetYaxis()->CenterTitle(true);
    
    hBraggB_30s_ssd = new TH2F("hBraggB_30s_ssd","Beam Bragg Curve Gated on ^{30}S SSD-OR",48,-0.5,47.5,200,0.,2000.);
    hBraggB_30s_ssd->SetOption("COLZ");
    hBraggB_30s_ssd->GetXaxis()->SetTitle("Pad No.");
    hBraggB_30s_ssd->GetXaxis()->CenterTitle(true);
    hBraggB_30s_ssd->GetYaxis()->SetTitle("#DeltaE (arb.)");
    hBraggB_30s_ssd->GetYaxis()->CenterTitle(true);
    
    hBraggB_rp = new TH2F("hBraggB_rp","Beam Scattering Location for ^{30}S",48,-0.5,47.5,45,6.,50.);
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
         pad_ch[i]->Write();
      } // end for : i
    } // end if: flag_raw
    
    if (flag_detail) { // write the SSD pad detailed histograms
      hSiTcal->Write(); 
      hSiPadEcal->Write();
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
  } // end if: flag_detail & flag_ssd & flag_strip

  if (flag_ppac){ // process PPAC?
    if (flag_raw) { // write the PPAC raw histograms
       for (UShort_t i=0;i<2;i++) hRf[i]->Write();
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
       hRF0Tof->Write();
       hRF1Tof->Write();
       hPpac0XY->Write();
       hPpac1XY->Write();  
       hPpac1XYcut->Write();  
       hPpac0XRF0->Write(); 
       hPpac1XRF0->Write(); 
       hPpac0XRF1->Write(); 
       hPpac0XRF1cut->Write(); 
       hPpac1XRF1cut->Write(); 
       hPpac1XRF1->Write();
       hPpac0XRF0_29p->Write(); 
       hPpac0XRF1_29p->Write(); 
       hPpac0XRF0_30s->Write(); 
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
    } // end if: flag_detail
  } // end if: flag_ppac

  if (flag_ssd && flag_ppac){ // process PPAC & SSD?
    if (flag_detail) { // write PPAC / SSD detailed histograms
       for (UShort_t i=0;i<18;i++){
         hPpac0TpadT_30s_ch[i]->Write();
         hPpac0TpadT_29p_ch[i]->Write();
       }
    } // end if: flag_detail
  } // end if: flag_ssd & flag_ppac

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
      //hPadC->Write();
      hBraggB->Write();
      hBraggL->Write();
      hBraggR->Write();
      hBraggC->Write();
      hPadXB->Write();
      hPadYB->Write();
      hPadXC->Write();
      hPadB_3D->Write();
      for (UShort_t i=0;i<144;i++) {
         hTrough[i]->Write();
         if (i<48){
           hBraggB_ch[i]->Write();
         }
      }
    } // end if: flag_detail
  } // end if: flag_tpc
  
  if (flag_ppac && flag_tpc && flag_detail){
    hPadXB_30s->Write();
    hPadYB_30s->Write();
    hBraggB_30s->Write();
    //hBraggB_30s2->Write();
    hBraggB_29p->Write();
    //hBraggB_31s->Write();
    for (UShort_t i=0;i<15;i++) hBraggB_30s_jch[i]->Write();
    hBraggB_3D->Write();

    hBraggB_30s_ds->Write();
    hBraggB_30s_ssd->Write();
    hBraggB_rp->Write();
    hBraggB_3D_ds->Write();
    hBraggB_3D_ssd->Write();
    for (UShort_t i=0;i<48;i++){
      hBraggB_ds_ch[i]->Write();
      hBraggB_ssd_ch[i]->Write();
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
  delete hPpac0XRF0 ;
  delete hPpac0XRF0_29p ;
  delete hPpac0XRF1_29p ;
  delete hPpac0XRF0_30s ;
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
  delete hGemMstpcAll ;
  delete hGemMstpcSide ;
  delete hBraggB_30s ;
  delete hBraggB_30s_ds ;
  delete hBraggB_30s_ssd ;
  delete hPadXB ;
  delete hPadXB_30s ;
  delete hPadYB ;
  delete hPadYB_30s ;
  delete hBraggB_3D ;
  delete hBraggB_3D_ds ;
  delete hBraggB_3D_ssd ;
  delete hPadB_3D;
  delete hBraggB_29p ;
  delete hBraggB_31s ;
  delete hBraggB ;
  delete hBraggL ;
  delete hBraggR ;
  delete hBraggC ;
  delete hBraggB_rp;
  for (int i=0;i<2;i++) delete hRf[i]; 
  for (int j=0;j<15;j++) delete hBraggB_30s_jch[j];
  for (int i=0;i<97;i++){
    delete strip_ch[i];
    delete strip_cal_ch[i];
    delete strip_ch_gated30s[i];
    delete strip_ch_gated29p[i];
  }
  delete hPpacToF_30s;
  delete hPpacToF_29p;
  for (int i=0;i<18;i++){
    delete pad_ch[i];
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
  for (int i=0;i<144;i++){
    delete hfSample[i]; 
    delete hTrough[i];
    delete hTpc_dE_ch[i];
  }
  for (int i=0;i<48;i++){
    delete hBraggB_ch[i];
    delete hBraggB_ds_ch[i];
    delete hBraggB_ssd_ch[i];
  }
} // close HistClean    
