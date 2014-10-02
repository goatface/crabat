/*
 Copyright: 
 2010 daid KAHL, OTA shinsuke, HASHIMOTO takashi, YAMAGUCHI hidetoshi 
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

TF1 *finter_f1, *finter_f2;
double finter(double *x, double*par) {
   return TMath::Abs(finter_f1->EvalPar(x,par) - finter_f2->EvalPar(x,par));
}

extern "C" { // removed for now since we don't use it.  Uncomment Makefile to enable
        void enewzsub_(int *z1, float *m1, float *e, char matter1[33], int *unit_pressure, float *pressure, float *temperature, int *unit_thick, float *thick1, float *aft_ene);
	}


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
bool event=false;
int eventcounter=0;
// peak finding rebin size
const int rebin=4;

double hg_peak_left[8]={};
int hg_peak_left_counter[8]={};
double hg_peak_right[8]={};
int hg_peak_right_counter[8]={};
int hg_peak_holder[2]={};

	    double centerprojsum=0.;
	    int ytestcounter=0;
	    double ypointscounter=0;

double XpadBcalib[2];
double XpadCcalib[2];

const float padsize = 4.; // pad size mm
//const float padsize = (196./48.); // pad size mm

//scattering location declarations
bool scatter_event = false;
int nPointsB=2, nPointsP=0;
double padsHitC_double;
double zppac=0.;
double zmstpc=0.;
int padcpoints=0;
double chi2min[4], chi2nu[4];
double xmin=83.,xmed,xmax,xmedbest;
double chi2minGlobal[2]; 
int pxmedlow=3;
int pxmedhigh;
double tempx[2],tempy[2];
int temppoint=0;
int junk=0;
double dE[2];
int temppe;
double tempxe[2],tempee[2];
float pad_scatter_tmp;
int pad_scatter;
double dEpadBsum = 36.9; // 17.3 um mylar; 6.8 um kapton; 194 torr, 213 K 82.75+4 = 86.75 mm HeCO2
int gate_30s_bragg = 0;
double dEpadB_gate_30s[2][48]={{
0.,
8154.85,
7544.33,
7591.57,
7518.22,
7561.07,
7463.77,
7507.3,
7681.71,
7951.29,
7950.4,
8094.47,
8084.42,
7825.75,
7731.6	,
7643.41,
7515.75,
0.,
6848.76,
6728.66,
6117.14,
5803.94,
5509.24,
5142.19,
0.,
0.,
0.,
3426.48,
2366.12,
1697.44,
1329.89,
1117.92,
0.,
0.,
0.,
0.,
0.,
0.,
0.,
0.,
927.907,
974.823,
893.893,
787.807,
0.,
0.,
0.,
0.
},{
0.,
13216.7,
12614.3,
12901.2,
13172.6,
13674.1,
13593.8,
13767.1,
14184.2,
14541.5,
14449.2,
14702.9,
14776.1,
14665.1,
14692.5,
14711.2,
14462.7,
0.,
13436.9,
13598.2,
13208.9,
13114.2,
13186.8,
13316.5,
0.,
0.,
0.,
14397.6,
13594,
13771.6,
13264.5,
12836.2,
0.,
0.,
0.,
0.,
0.,
0.,
0.,
0.,
5477.22,
5169.21,
4925.92,
4841.45,
0.,
0.,
0.,
0.
}};

KVNucleus *beam = new KVNucleus;
KVNucleus *target = new KVNucleus;
KVNucleus *recoil = new KVNucleus;
KV2Body reaction;

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
  
  beam->Set("30S");
  target->Set("4He");
  recoil->Set("4He");
 
  KV2Body reaction;
  reaction.SetProjectile(beam);
  reaction.SetTarget(target);
  reaction.SetOutgoing(recoil);

  // initalize the variables for matter and energy loss enewzsub
 
  //SOLID TARGET PROTOTYPE - A mylar foil (so-called beam stopper)
  char matter_mylar[34] = "mylar"; // Matter the beam passes through, as defined in SNKE_MATTER.INC -- case insensitive
  int unit_pressure_mylar = 0 ; // Not used for solid target; set as 0
  float pressure_mylar = 0 ; // Not used for solid target; set as 0
  float temperature_mylar = 0 ; // Not used for solid target; set as 0
  int unit_thick_mylar = 1; // Define the units of thickness: 1=mm ; 2=mg/cm^2
  float thick_mylar = 0.0025 ; // Define the value of the thickness.

 
  //GAS TARGET PROTOTYPE - Active target fill gas
  char * matter_gas;
  matter_gas = new char[34];
  strcpy (matter_gas, "heco2");
  int unit_pressure_gas = 1 ; // Define the units of pressure: 1=Torr ; 2=mbar ; 3=atm.  If you change here, the program handles the rest.
  float pressure_gas=194; // Gas pressure (Torr with above line as I set it) -- This is input by user.
  float temperature_gas =  280 ; // Define temperature of the target in Kelvin
  int unit_thick_gas = 1; // Define the units of thickness: 1=mm ; 2=mg/cm^2
  // CHANGE FOR YOUR EXPERIMENTAL CASE THE DISTANCE FROM THE TPC ENTRANCE WINDOW TO 0 DEGREE SSD!
  //float thick_gas = 393.5; // This variable is case-by-case. It is currently set to mm units. // used this to get p 116 fits
  // this is the proper distance, but no account for the foil bowing etc 27 Aug 2012 15:39:29 
  float thick_gas = 0; // This variable is case-by-case. It is currently set to mm units. 

  //This part just fixes the name length of the matter to 32 characters for FORTRAN
  for(int j=0;j<33;j++)
  if(matter_mylar[j]=='\0')
  {
    matter_mylar[j]=' ';
    matter_mylar[j+1]='\0';
  }
  for(int i=0;i<33;i++)
    if(matter_gas[i]=='\0')
    {
      matter_gas[i]=' ';
      matter_gas[i+1]='\0';
    }
  
  int z_alpha = recoil->GetZ();
  float m_alpha = recoil->GetA();
  float e_alpha_send, e_alpha_return;
  double energy_alpha_calc;

  bool geo_verbose = false;
  bool geo_draw = false;
  // variable parameters for plotting
  // always define each time you plot
  float geo_xmin=-200.; // minimum input for plotting
  float geo_xmax=200.; // maximum input for plotting
  float Z0=0.; // Z geo_intercept
  
  TF1 *fint = new TF1("fint",finter,geo_xmin,geo_xmax,0);
  
  double p[2]; // temp to extract fit parameters
  
  // basic constants of field cage
  const float xCageSize=278.;
  const float zCageSize=295.;
  const float SsdFootDepth=30.;
  const float SsdMountWidth=252.;
  const int npCageFoot=4;
 
  const int handMarkerStyle=24;
  const float handMarkerSize=0.75;
  const int handMarkerColor=4;

  const int stripMarkerStyle=33;
  const float stripMarkerSize=1.5;
  const int stripMarkerColor=1;

  // data from Vol 4, p. 88
  
  // 45 is back of foot to front of field cage
  // thus these data are in the cage's own reference
  Float_t zCageFoot[npCageFoot]={45.,85.,145.,35.};
  Float_t zeCageFoot[npCageFoot]={}; 
  // beam axis reference
  Float_t xCageFoot[npCageFoot]={6.,9.,15.,16.5};
  // cage's own reference
  //Float_t x[npCageFoot]={-6.,-9.,-15.,-16.5};
  Float_t xeCageFoot[npCageFoot]={}; 
  

  // SET Z REFERENCE FRAME FOR ENTIRE MACRO
  // zReference:
  // 27.75 beginning of pad 0 Vol 5 p.4
  // 0 beginning of field cage is Z=0
  // -35 inner wall
  // -55 entrance window
  const float zReference=-55.;
  zCageFoot[0]-=zReference;
  TCanvas *c1; 
  if (geo_draw) c1 = new TCanvas("c1","c1",700,800);

  for (int i=0;i<npCageFoot;i++){
    if (i!=0) zCageFoot[i]+=zCageFoot[i-1]; // get the full depth, not the foot spacing
    zeCageFoot[i]=0.1; // assume the production error is 100 um
    xeCageFoot[i]=0.5; // assume ability to use a ruler is 0.5 mm precision
    if (geo_verbose) cout << "zCageFoot[" << i << "]: " << zCageFoot[i] << endl;
  }
  TGraphErrors *gCageCenter;
  gCageCenter = new TGraphErrors(npCageFoot,xCageFoot,zCageFoot,xeCageFoot,zeCageFoot);
  gCageCenter->SetMarkerColor(kBlue);
  gCageCenter->SetMarkerStyle(20);
  gCageCenter->SetMarkerSize(1);
  
  //YOUR TITLE & LABELS GO HERE!!!
  gCageCenter->SetTitle("Active Target Geometry By Laser Alignment");
  //gCageCenter->SetTitle("Field Cage Rotation By Laser Alignment");
  gCageCenter->GetYaxis()->SetTitle("Z position (mm)");
  gCageCenter->GetYaxis()->SetTitleOffset(1.35);
  gCageCenter->GetXaxis()->CenterTitle(true);
  gCageCenter->GetXaxis()->SetTitle("X position (mm)");
  gCageCenter->GetYaxis()->CenterTitle(true);
  
  gCageCenter->GetXaxis()->SetLimits(geo_xmin,geo_xmax);
  gCageCenter->GetYaxis()->SetRangeUser(-10,400);
  gCageCenter->SetMarkerStyle(handMarkerStyle);
  gCageCenter->SetMarkerSize(handMarkerSize);
  gCageCenter->SetMarkerColor(handMarkerColor);
  if (geo_draw) gCageCenter->Draw("AP9");

  TF1 *fLinear;
  fLinear = new TF1("fLinear","pol1");
  fLinear->SetLineWidth(1);
  fLinear->SetLineColor(3);
  fLinear->SetNpx(1000);
  gCageCenter->Fit("fLinear","Q");
  
  fLinear->GetParameters(&p[0]);
  const Float_t CageSlope=p[1];
  const float CageSlopeInv = -1/CageSlope;
  const Float_t CageX0=-p[0]/p[1];
  const Float_t CageZ0=p[0];
  if (geo_verbose) cout << "Field Cage Parameters\nCageZ0: " << CageZ0 << "\nCageX0: " << CageX0 << endl;
  if (geo_verbose) cout << "\nCageSlope: " << CageSlope << "CageSlopeInv: " << CageSlopeInv << endl;
  TF1 *fCageCenter = new TF1("fCageCenter","[0]+x*[1]",-1000.,1000.);
  fCageCenter->SetParameters(CageZ0,CageSlope);
  fCageCenter->SetNpx(1000);
  fCageCenter->SetLineWidth(2);
  fCageCenter->SetLineColor(2);
  //draw at the end to put overtop of bridge
  //fCageCenter->Draw("SAME9");
  TF1 *fCageCenterInv = new TF1("fCageCenterInv","(x-[0])/[1]");
  fCageCenterInv->SetParameters(CageZ0,CageSlope);

  float xCageCenter=fCageCenterInv->Eval((zCageSize/2.)-zReference); // center of FC in X
  if (geo_verbose) cout << "center of field cage (X):" <<  xCageCenter << endl;
  float zCageCenter=fCageCenter->Eval(xCageCenter);
  Z0 = zCageCenter-(CageSlopeInv*xCageCenter);
  TF1 *fCageCenterOrthogonal = new TF1("fCageCenterOrthogonal","[0]+x*[1]",-1000.,1000.);
  fCageCenterOrthogonal->SetParameters(Z0,CageSlopeInv);
  //fCageCenterOrthogonal->Draw("SAME9");
  
  // SSD LEFT POSITION BY INTERCEPTION
  float xSsdLOffsetTpc[3]={-22.,-20.,-18.}; // -13 should be -18 as shown here.  Otherwise the data point cannot fit at all
  float xCageSsdLFoot[3];
  float zCageSsdLFoot[3];
  TF1 *fCageFoot[npCageFoot];
  if (geo_verbose) cout << "xSsdL\tzSsdL" << endl;
  for (int i=0;i<npCageFoot-1;i++){
    //cout << "Iteration: " << i << endl;
    fCageFoot[i] = new TF1("fCageFoot","[0]+x*[1]",geo_xmin,geo_xmax);
    fCageFoot[i]->SetLineWidth(1);
    Z0 = zCageFoot[i]-(CageSlopeInv*xCageFoot[i]);
    // arbitrary propagation of error >_<
    //cout << "fLinear->Eval(x[i]): " << fLinear->Eval(x[i]) << endl;
    //float Z0 = fLinear->Eval(x[i])+((1/CageSlope)*x[i]);
    fCageFoot[i]->SetParameters(Z0,CageSlopeInv);
    fCageFoot[i]->SetLineColor(3);
    if (geo_draw) fCageFoot[i]->Draw("9LSAME");
    xCageSsdLFoot[i]=-(((xCageSize+30)/2)-xCageFoot[i]-xSsdLOffsetTpc[i]);
    //zSsdL[i]=
    zCageSsdLFoot[i]=fCageFoot[i]->Eval(xCageSsdLFoot[i]);
    //cout << xCageSsdLFoot[i] << "\t" << zCageSsdLFoot[i] << endl ;
  }
  
  // FRONT OF CHAMBER WALL BY INTERCEPTION
  float wall_rotation=CageSlopeInv+(2./308.); // chamber wall is slightly less rotated 
  //cout << "X CageSlope wall by beam axis: " << wall_rotation << endl;
  TF1 *fChamberWallIn = new TF1("fChamberWallIn","[0]+x*[1]",geo_xmin,geo_xmax);
  fChamberWallIn->SetParameters(-zReference-35.,wall_rotation);
  fChamberWallIn->SetLineWidth(2);
  fChamberWallIn->SetLineColor(1);
  fChamberWallIn->SetLineStyle(2);
  if (geo_draw) fChamberWallIn->Draw("9LSAME");
  TF1 *fChamberWallOut = new TF1("fChamberWallOut","[0]+x*[1]",geo_xmin,geo_xmax);
  fChamberWallOut->SetParameters(-zReference-55.,wall_rotation);
  fChamberWallOut->SetLineWidth(2);
  fChamberWallOut->SetLineColor(1);
  fChamberWallOut->SetLineStyle(2);
  if (geo_draw) fChamberWallOut->Draw("9LSAME");

  const int npSsdLFoot=6;
  float xSsdLFoot[npSsdLFoot],zSsdLFoot[npSsdLFoot],xeSsdLFoot[npSsdLFoot],zeSsdLFoot[npSsdLFoot];
  for (int i=0;i<npSsdLFoot;i++){
    xeSsdLFoot[i]=0.5; // assume ability to use a ruler is 0.5 mm precision
    zeSsdLFoot[i]=0.5; // assume the production error is 100 um
    //zSsdLErr[i]+=0.1; // assume the production error is 100 um
  }

  //cout << "Foot 0 left " << fCageFoot[0]->Eval(-148) << endl;
  //cout << "Wall left " << fChamberWall->Eval(-148) << endl;
  //cout << "Foot 0 left to wall: " << fCageFoot[0]->Eval(-148.)-fChamberWall->Eval(-148.) << endl;
  xSsdLFoot[0]=-173;
  zSsdLFoot[0]=fChamberWallIn->Eval(xSsdLFoot[0])+21.;
  //cout << "SSD L foot X\tZ" << endl;
  //cout << xSsdLFoot[0] <<  "\t" << zSsdLFoot[0]  << endl;
  
  float thetaSsdLFoot=(TMath::Pi()/2)-asin(6./SsdFootDepth);
  
  xSsdLFoot[2]=-168.5;
  zSsdLFoot[2]=(zCageSize+36.-SsdMountWidth)*sin(thetaSsdLFoot)+zSsdLFoot[0];

  xSsdLFoot[4]=-157.;
  zSsdLFoot[4]=zCageSize*sin(thetaSsdLFoot)+zSsdLFoot[0];
  //cout << xSsdLFoot[4] <<  "\t" << zSsdLFoot[4]  << endl;
  // SSD Left foot determination 
  for (int i=0;i<npCageFoot-1;i++){
    zSsdLFoot[i*2+1]=zCageSsdLFoot[i];
    xSsdLFoot[i*2+1]=xCageSsdLFoot[i];
  }
  if (geo_verbose) cout << "Data for plotting: " << endl;
  for (int i=0;i<npSsdLFoot;i++){ if (geo_verbose) cout << xSsdLFoot[i] << "\t" << zSsdLFoot[i] << endl;}
  TGraphErrors *gSsdLFoot = new TGraphErrors(npSsdLFoot,xSsdLFoot,zSsdLFoot,xeSsdLFoot,zeSsdLFoot);
  gSsdLFoot->SetMarkerColor(handMarkerColor);
  gSsdLFoot->SetMarkerStyle(handMarkerStyle);
  gSsdLFoot->SetMarkerSize(handMarkerSize);
  if (geo_draw) gSsdLFoot->Draw("9PSAME");
  
  gSsdLFoot->Fit("fLinear","Q");
  fLinear->GetParameters(&p[0]);
  
  // SSD 3a (left, d/s) strip determination
  float SsdLDetOffset=p[0]+(9.4/(sin(TMath::Pi()/2-atan(p[1])))); // 13 - 3.6 = 9.4 (spacer - frame) online II p. 160
  float SsdLDetSlope=p[1];

  //cout << "SsdLDetOffset: " << SsdLDetOffset << endl;
  //this is the downstream edge of the mount
  //xSsdL[2]=-168.5;
  float zSsdLCenter=(zCageSize+36.-(SsdMountWidth/2.))*sin(thetaSsdLFoot)+zSsdLFoot[0];
  //cout << "zSsdLCenterL " << zSsdLCenter << endl;

  // it really should be something like this, fix later for presentation
      //geo_xmin=xPadR[i]+xCageCenter+CageSlopeInv*(234./2.);
      //geo_xmax=xPadR[i]+xCageCenter-CageSlopeInv*(234./2.);
  //geo_xmin=fCageCenterInv->Eval(zSsdLCenter)+CageSlopeInv*(SsdMountWidth/2.);
  //geo_xmax=fCageCenterInv->Eval(zSsdLCenter)-CageSlopeInv*(SsdMountWidth/2.);
  geo_xmin = -1000;
  geo_xmax = 1000;
  TF1 *fSsdL = new TF1("fSsdL","[0]+x*[1]",geo_xmin,geo_xmax);
  TF1 *fSsdLInv = new TF1("fSsdLInv","(x-[0])/[1]",geo_xmin,geo_xmax);
  fSsdL->SetParameters(SsdLDetOffset,SsdLDetSlope);
  fSsdLInv->SetParameters(SsdLDetOffset,SsdLDetSlope);
  fSsdL->SetLineColor(1);
  fSsdL->SetLineWidth(1);
  fSsdL->SetNpx(1000);
  //fSsdL->Draw("L9SAME");
  
  //cout << "X of center: " << fSsdLInv->Eval(zSsdLCenter) << endl;

  const float SsdBridgeWidth=12.4; // logbook 4 p 152 has 6.2 from SSD frame 
  const float SsdStripWidth=(90.6/8.);
  const int nStrips=8;

  // depth in its own frame
  float StripSsdL[8]; // have never used notation like "Strip 0" always starts "Strip 1"
  // absolute positions
  float xStripSsdL[8];
  float zStripSsdL[8];
  //cout << endl << "LEFT SSD" << endl; 
  StripSsdL[7]=zSsdLCenter+((SsdBridgeWidth/2.)+0.5*SsdStripWidth);
  //cout << "Strip 7: "  << StripSsdL[7] << endl;
  for (int i=6;i>-1;i--){ // looping backward is strange but this is how we count the strips going out from center...
    StripSsdL[i]=StripSsdL[i+1]+SsdStripWidth;
    //cout << "Strip " << i << ": "  << StripSsdL[i] << endl;
  }
  for (int i=0;i<8;i++){
    zStripSsdL[i]=StripSsdL[i]*sin(atan(SsdLDetSlope));
    //cout << "Strip Z" << i << ": "  << zStripSsdL[i] << endl;
    xStripSsdL[i]=fSsdLInv->Eval(zStripSsdL[i]);
  }
  TGraph *gStripsL = new TGraph(nStrips,xStripSsdL,zStripSsdL);
  gStripsL->SetMarkerColor(stripMarkerColor);
  gStripsL->SetMarkerStyle(stripMarkerStyle);
  gStripsL->SetMarkerSize(stripMarkerSize);
  if (geo_draw) gStripsL->Draw("PSAME9");


  // SSD 1a (center, right) strip determination
  const int npSsdCFoot=1;
  float xSsdCFoot[npSsdCFoot],zSsdCFoot[npSsdCFoot],xeSsdCFoot[npSsdCFoot],zeSsdCFoot[npSsdCFoot];
  xSsdCFoot[0]=(146.5-131.5)/2.; // vol 4 p 152
  zSsdCFoot[0]=zCageFoot[3]+14.; // online vol 2 p 160 // 

  TGraphErrors *gSsdCFoot = new TGraphErrors(npSsdCFoot,xSsdCFoot,zSsdCFoot,xeSsdCFoot,zeSsdCFoot);
  gSsdCFoot->SetMarkerColor(handMarkerColor);
  gSsdCFoot->SetMarkerStyle(handMarkerStyle);
  gSsdCFoot->SetMarkerSize(handMarkerSize);
  if (geo_draw) gSsdCFoot->Draw("9PSAME");
  
  geo_xmin=-1*(131.5-((xCageSize-SsdMountWidth)/2));
  geo_xmax=(146.5-((xCageSize-SsdMountWidth)/2));
  TF1 *fSsdCFoot = new TF1("fSsdCFoot","[0]+x*[1]",geo_xmin,geo_xmax);
  Z0 = zSsdCFoot[0]+((1/CageSlope)*xSsdCFoot[0]);
  fSsdCFoot->SetParameters(Z0,CageSlopeInv);
  fSsdCFoot->SetLineColor(3);
  fSsdCFoot->SetLineWidth(1);
  if (geo_draw) fSsdCFoot->Draw("L9SAME");

  float zSsdC=zSsdCFoot[0]+13.-3.6; // online vol 2 p 160 nominally 368.4 from the wall
  float xSsdC=fCageCenterInv->Eval(zSsdC)-fCageCenterInv->Eval(zSsdCFoot[0])+xSsdCFoot[0]; // roughly 8, see vol 4 p 152
  //cout << endl << "CENTER SSD" << endl; 
  
  geo_xmin=-1*(131.5-((xCageSize-SsdMountWidth)/2));
  geo_xmax=(146.5-((xCageSize-SsdMountWidth)/2));
  TF1 *fSsdC = new TF1("fSsdC","[0]+x*[1]",geo_xmin,geo_xmax);
  Z0 = zSsdC+((1/CageSlope)*xSsdC);
  fSsdC->SetParameters(Z0,CageSlopeInv);
  fSsdC->SetLineColor(1);
  fSsdC->SetLineWidth(1);
  //fSsdC->Draw("L9SAME");
  
  float StripSsdC[8]; // have never used notation like "Strip 0" always starts "Strip 1"
  float xStripSsdC[8];
  float xErrStripSsdC[8];
  float zStripSsdC[8];
  float zErrStripSsdC[8];

  StripSsdC[7]=xSsdC+(SsdBridgeWidth/2.)+0.5*SsdStripWidth;
  //cout << "Strip 7: "  << StripSsdC[7] << endl;
  for (int i=6;i>-1;i--){ // looping backward is strange but this is how we count the strips going out from center...
    StripSsdC[i]=StripSsdC[i+1]+SsdStripWidth;
    //cout << "Strip " << i << ": "  << StripSsdC[i] << endl;
  }
  for (int i=0;i<8;i++){
    xStripSsdC[i]=StripSsdC[i]*cos(atan(CageSlopeInv));
    xErrStripSsdC[i]=xStripSsdC[i]-((90.6/8.)*0.5);
    //cout << "Strip X" << i << ": "  << xStripSsdC[i] << endl;
    zStripSsdC[i]=fSsdC->Eval(xStripSsdC[i]);
    zErrStripSsdC[i]=fSsdC->Eval(xErrStripSsdC[i]);
  }

  TGraph *gStripsC = new TGraph(nStrips,xStripSsdC,zStripSsdC);
  gStripsC->SetMarkerColor(stripMarkerColor);
  gStripsC->SetMarkerStyle(stripMarkerStyle);
  gStripsC->SetMarkerSize(stripMarkerSize);
  if (geo_draw) gStripsC->Draw("P9SAME9");

  // Beam GEM
  // start: 35 + 27.75
  // end : + 192
    float PadB[48],zPadB[48],xPadB[48];
    PadB[0]=-zReference+27.75+2.;
    for (int i=1;i<48;i++){
      PadB[i]=PadB[i-1]+4.;
    }
    TF1 *fPadB[48];
    for (int i=0;i<48;i++){
      zPadB[i]=PadB[i]*sin(atan(CageSlope));
      xPadB[i]=fCageCenterInv->Eval(zPadB[i]);
      geo_xmin=xPadB[i]-(110./2.);
      geo_xmax=xPadB[i]+(110./2.);
      fPadB[i] = new TF1("fPadB","[0]+x*[1]",geo_xmin,geo_xmax);
      Z0 = zPadB[i]-(CageSlopeInv*xPadB[i]);
      //Z0 = fc_center->Eval(x[i])+((1/CageSlope)*x[i]);
      fPadB[i]->SetParameters(Z0,CageSlopeInv);
      fPadB[i]->SetLineColor(1);
      fPadB[i]->SetLineWidth(1);
      fPadB[i]->SetNpx(1000);
      if (geo_draw) fPadB[i]->Draw("L9SAME");
    }
  // Central GEM
  // start : +15.5
  // end : +32
    float PadC[8],zPadC[8],xPadC[8];
    PadC[0]=-zReference+27.75+192.+15.5+2.;
    for (int i=1;i<8;i++){
      PadC[i]=PadC[i-1]+4.;
    }
    TF1 *fPadC[8];
    TF1 *fPadCBridge[8];
    for (int i=0;i<8;i++){
      zPadC[i]=PadC[i]*sin(atan(CageSlope));
      xPadC[i]=fCageCenterInv->Eval(zPadC[i]);
      geo_xmin=xPadC[i]-(110./2.);
      geo_xmax=xPadC[i]+(110./2.);
      fPadC[i] = new TF1("fPadC","[0]+x*[1]",geo_xmin,geo_xmax);
      Z0 = zPadC[i]-(CageSlopeInv*xPadC[i]);
      //Z0 = fc_center->Eval(x[i])+((1/CageSlope)*x[i]);
      fPadC[i]->SetParameters(Z0,CageSlopeInv);
      fPadC[i]->SetLineColor(1);
      fPadC[i]->SetLineWidth(1);
      fPadC[i]->SetNpx(1000);
      if (geo_draw) fPadC[i]->Draw("L9SAME");
      
      geo_xmin=xPadC[i]-(40./2.);
      geo_xmax=xPadC[i]+(40./2.);
      fPadCBridge[i] = new TF1("fPadCBridge","[0]+x*[1]",geo_xmin,geo_xmax);
      fPadCBridge[i]->SetParameters(Z0,CageSlopeInv);
      fPadCBridge[i]->SetLineColor(0);
      fPadCBridge[i]->SetLineWidth(7);
      fPadCBridge[i]->SetNpx(1000);
      if (geo_draw) fPadCBridge[i]->Draw("L9SAME");
    }
    // need to draw this here to put it over the bridge
    if (geo_draw) fCageCenter->Draw("SAMEL9");


    //Side GEM
    //
    geo_xmin=-1000.;
    geo_xmax=1000;
    TF1 *fPadSideMax = new TF1("fPadSideMax","[0]+x*[1]",geo_xmin,geo_xmax);
    Z0=-zReference+zCageSize/2.+234./2.;
    fPadSideMax->SetParameters(Z0,CageSlopeInv);
    fPadSideMax->SetNpx(1000);
    //fPadSideMax->Draw("L9SAME");
    TF1 *fPadSideMin = new TF1("fPadSideMin","[0]+x*[1]",geo_xmin,geo_xmax);
    Z0=-zReference+zCageSize/2.-234./2.;
    fPadSideMin->SetParameters(Z0,CageSlopeInv);
    fPadSideMin->SetNpx(1000);
    //fPadSideMin->Draw("L9SAME");

    // Left GEM
    // 110/2 + 21.5 + 2 from center of FC to center of inner most pad

    float PadL[8],zPadL[8],xPadL[8];

    // pad distance X from center in its own frame
    PadL[0]=-1.*((110./2.)+21.5+2.);
    // add the field cage offset -- by this method we don't need it
    //padL[0]+=xCageCenter; //
    //cout << "padL[0]: " << PadL[0] << endl;
    for (int i=1;i<8;i++){
      PadL[i]=PadL[i-1]-4.;
      //cout << "PadL[" << i << "]: " << PadL[i] << endl;
    }
    
    TF1 *fPadL[8];
    for (int i=0;i<8;i++){
    //for (int i=0;i<1;i++){      
      // pad distance X from center in rotated frame
      xPadL[i] = PadL[i];
      // don't need it?!?!
      xPadL[i] = PadL[i]*cos(atan(CageSlopeInv));
      //cout << "xPadL[" << i << "]: " << xPadL[i] << endl;
      geo_xmin=-500;
      geo_xmax=500;
      TF1 *fTemp = new TF1("fTemp","[0]+x*[1]",geo_xmin,geo_xmax);
      Z0 = CageZ0-(xPadL[i]/(sin(TMath::Pi()/2-atan(CageSlope))));
      fTemp->SetParameters(Z0,CageSlope);
      finter_f1=fTemp;
      finter_f2=fPadSideMin;
      geo_xmin=fint->GetMinimumX();
      finter_f2=fPadSideMax;
      geo_xmax=fint->GetMinimumX();
      //geo_xmin=xPadL[i]+xCageCenter+CageSlopeInv*(234./2.);
      //geo_xmax=xPadL[i]+xCageCenter-CageSlopeInv*(234./2.);
      //cout << "geo_xmin: " << geo_xmin << endl;
      //cout << "geo_xmax: " << geo_xmax << endl;
      fPadL[i] = new TF1("fPadL","[0]+x*[1]",geo_xmin,geo_xmax);
      Z0 = CageZ0-(xPadL[i]/(sin(TMath::Pi()/2-atan(CageSlope))));
      fPadL[i]->SetParameters(Z0,CageSlope);
      fPadL[i]->SetLineColor(1);
      fPadL[i]->SetLineWidth(1);
      fPadL[i]->SetNpx(1000);
      if (geo_draw) fPadL[i]->Draw("L9SAME");
      if (geo_verbose) cout << "Left pad " << i << " minimum Z: " << fPadL[i]->Eval(geo_xmin) << endl;
    }

    // Right GEM
    float PadR[8],zPadR[8],xPadR[8];

    // pad distance X from center in its own frame
    PadR[0]=((110./2.)+21.5+2.);
    if (geo_verbose) cout << "padR[0]: " << PadR[0] << endl;
    for (int i=1;i<8;i++){
      PadR[i]=PadR[i-1]+4.;
      if (geo_verbose) cout << "PadR[" << i << "]: " << PadR[i] << endl;
    }
    
    TF1 *fPadR[8];
    for (int i=0;i<8;i++){
    //for (int i=0;i<1;i++){      
      // pad distance X from center in rotated frame
      xPadR[i] = PadR[i];
      if (geo_verbose) cout << "xPadR[" << i << "]: " << xPadR[i] << endl;
      geo_xmin=-500;
      geo_xmax=500;
      TF1 *fTemp = new TF1("fTemp","[0]+x*[1]",geo_xmin,geo_xmax);
      Z0 = CageZ0-(xPadR[i]/(sin(TMath::Pi()/2-atan(CageSlope))));
      fTemp->SetParameters(Z0,CageSlope);
      finter_f1=fTemp;
      finter_f2=fPadSideMin;
      geo_xmin=fint->GetMinimumX();
      finter_f2=fPadSideMax;
      geo_xmax=fint->GetMinimumX();
      //geo_xmin=xPadR[i]+xCageCenter+CageSlopeInv*(234./2.);
      //geo_xmax=xPadR[i]+xCageCenter-CageSlopeInv*(234./2.);
      //geo_xmin=-1000.;geo_xmax=1000.;
      if (geo_verbose) cout << "geo_xmin: " << geo_xmin << endl;
      if (geo_verbose) cout << "geo_xmax: " << geo_xmax << endl;
      fPadR[i] = new TF1("fPadR","[0]+x*[1]",geo_xmin,geo_xmax);
      Z0 = CageZ0-(xPadR[i]/(sin(TMath::Pi()/2-atan(CageSlope))));
      fPadR[i]->SetParameters(Z0,CageSlope);
      fPadR[i]->SetLineColor(1);
      fPadR[i]->SetLineWidth(1);
      fPadR[i]->SetNpx(1000);
      if (geo_draw) fPadR[i]->Draw("L9SAME");
    }

  // check various geo_interceptions
  ofstream fileout_c;
  fileout_c.open("geometry_hg_c.dat", ios::out);
  //fileout_c.open("geometry_hg_c.dat", ios::out | ios::app);
  ofstream fileout_l;
  fileout_l.open("geometry_hg_l.dat", ios::out);
  //fileout_l.open("geometry_hg_l.dat", ios::out | ios::app);
  //CHANGE THESE
  // choose which GEM
  for (int PadIntercept=0;PadIntercept<2;PadIntercept++){ // 0->Central; 1->Left
    // strip number 1 to 8
    //int StripGate=7;
    for (int StripGate=1;StripGate<9;StripGate++){
      if (geo_verbose) cout << "Strip " << StripGate << endl;
      float xAlpha[2],xErrAlpha[2],zAlpha[2],zErrAlpha[2];
      xAlpha[0]=0.;
      zAlpha[0]=-zReference-(35.+14.);


      if (PadIntercept==0){
        xAlpha[1]=xStripSsdC[StripGate-1];
        xErrAlpha[1]=xErrStripSsdC[StripGate-1];
        zAlpha[1]=zStripSsdC[StripGate-1];
        zErrAlpha[1]=zErrStripSsdC[StripGate-1];
      }
      if (PadIntercept==1){
        xAlpha[1]=xStripSsdL[StripGate-1];
        zAlpha[1]=zStripSsdL[StripGate-1];
      }

      TGraph *gAlpha = new TGraph(2,xAlpha,zAlpha);
      gAlpha->SetLineColor(kRed);
      gAlpha->SetLineWidth(2);
      //gAlpha->Draw();
      
      gAlpha->Fit("fLinear","Q");
      fLinear->GetParameters(&p[0]);

      geo_xmin=-250;
      geo_xmax=150;

      TF1 *fAlpha = new TF1("fAlpha","[0]+x*[1]",geo_xmin,geo_xmax);
      fAlpha->SetParameters(p[0],p[1]);
      fAlpha->SetLineColor(kRed);
      fAlpha->SetLineWidth(2);
      fAlpha->SetNpx(1000);
      //fAlpha->Draw("SAME");
      
      /*
      if (StripGate==4  && PadIntercept==0){
        //fAlpha->Draw("SAME");
	float length;
	length = sqrt((zAlpha[1]-zAlpha[0])**2 + (xAlpha[1]-zAlpha[0])**2);
	cout << "LENGTH: " << length << endl;
      }
      */
      
      TGraph *gAlphaErr = new TGraph(2,xErrAlpha,zErrAlpha);
      gAlphaErr->SetLineColor(kRed);
      gAlphaErr->SetLineWidth(2);
      //gAlphaErr->Draw();
      
      gAlphaErr->Fit("fLinear","Q");
      fLinear->GetParameters(&p[0]);

      geo_xmin=-250;
      geo_xmax=150;

      TF1 *fAlphaErr = new TF1("fAlphaErr","[0]+x*[1]",geo_xmin,geo_xmax);
      fAlphaErr->SetParameters(p[0],p[1]);
      fAlphaErr->SetLineColor(kRed);
      fAlphaErr->SetLineWidth(2);
      fAlphaErr->SetNpx(1000);
      //fAlphaErr->Draw("SAME");
      
      float geo_intercept;
      float xErrEval[8];
      for (int i=0;i<8;i++){
        //double xint = fint->GetMinimumX();
        if (PadIntercept==0) {
          finter_f1=fAlpha;
          finter_f2=fPadC[i];
          geo_intercept = fint->GetMinimumX();
	  if (geo_verbose) cout << "Pad " << i << " " << geo_intercept << endl;
	  fileout_c << geo_intercept << "\t" ;
	  finter_f1=fAlphaErr;
	  xErrEval[i]=geo_intercept-(fint->GetMinimumX());
        }
        if (PadIntercept==1){
          finter_f1=fAlpha;
          finter_f2=fPadL[i];
          geo_intercept=fPadL[i]->Eval(fint->GetMinimumX());
	  if (geo_verbose) cout << "Pad " << i << " " << geo_intercept << endl;
          fileout_l << geo_intercept << "\t" ;
        }
      }
      if (PadIntercept==0) fileout_c << endl;
      for (int i=0;i<8;i++){
        if (PadIntercept==0) {
          fileout_c << xErrEval[i] << "\t" ;
	}
      }
      if (PadIntercept==0) fileout_c << endl;
      if (PadIntercept==1) fileout_l << endl;
    

    } //StripGate
  } // PadIntercept
  fileout_c.close();
  fileout_l.close();
  
  cout << "Analyzer::Loop processing run number " << run << endl;

  if (fChain == 0) return;
  //Fast function returns kBigNumber when it is TChain
  //Long64_t nentries = fChain->GetEntriesFast(); // I think this caused segfault -- daid
  Long64_t nentries = fChain->GetEntries();
  cout << "Entries: " << nentries << endl;

  Long64_t nbytes = 0, nb = 0;

  //Fancy color palette
  gStyle->SetPalette(1,0);

TCutG *cutg = new TCutG("mycut",6);
cutg->SetPoint(0,177.011,33.2405);
cutg->SetPoint(1,1733.76,33.1139);
cutg->SetPoint(2,4286.21,1.97468);
cutg->SetPoint(3,3709.34,2.10127);
cutg->SetPoint(4,153.305,29.443);
cutg->SetPoint(5,153.305,33.1139);
cutg->SetPoint(6,177.011,33.2405);
  //Calibration params

  Calibration *ssd_strip_calib=new Calibration("ssd_strip_calib.dat",97);//12x8
  Calibration *ssd_padE_calib=new Calibration("ssd_padE_calib.dat",19);
  Calibration *ssd_pad1a_calib=new Calibration("ssd_1a_calib.dat",9);
  Calibration *ssd_padT_calib=new Calibration("ssd_padT_calib.dat",19);
  Calibration *ppac_calib=new Calibration("ppac_calib.dat",10);
  Calibration *rf_calib=new Calibration("rf_calib.dat",2);
  Calibration *rf_ds_calib=new Calibration("rf_delay_calib.dat",2);
  Calibration *beam_x_left_calib=new Calibration("beam_x_left_calib.dat",48);
  Calibration *beam_x_right_calib=new Calibration("beam_x_right_calib.dat",48);
  Calibration *hg_xc_calib=new Calibration("hg_xc_calib.dat",48);
  Calibration *hg_x_left_calib=new Calibration("hg_x_left_calib.dat",48);
  Calibration *hg_x_right_calib=new Calibration("hg_x_right_calib.dat",48);
  Calibration *hg_center_calib=new Calibration("hg_center_calib.dat",8);
  Calibration *beam_y_delay_calib1=new Calibration("beam_y_delay_calib1.dat",48);
  Calibration *beam_y_delay_calib2=new Calibration("beam_y_delay_calib2.dat",48);
  Calibration *beam_y_calib1=new Calibration("beam_y_calib1.dat",48);
  Calibration *beam_y_calib2=new Calibration("beam_y_calib2.dat",48);
 
  TDetector MSTPC;
  vector<vector<Short_t> > tpc_ch=MSTPC.TDetector::SetTpcMap();
  vector<vector<Double_t> > padBgain=MSTPC.TDetector::SetTpcBgain();
  vector<Double_t> padXBgeo=MSTPC.TDetector::SetGeoXB();
  vector<Double_t> padXCgeo=MSTPC.TDetector::SetGeoXC();

  if (flag_detail){
    if (flag_strip){
      ssd_strip_calib->ScanFile();
      ssd_strip_calib->ShowEntries();
    }
    if (flag_ssd){
      ssd_padE_calib->ScanFile();
      ssd_padE_calib->ShowEntries();
      ssd_pad1a_calib->ScanFile();
      ssd_pad1a_calib->ShowEntries();
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

      beam_y_delay_calib1->ScanFile();
      beam_y_delay_calib1->ShowEntries();
      beam_y_delay_calib2->ScanFile();
      beam_y_delay_calib2->ShowEntries();
      beam_y_calib1->ScanFile();
      beam_y_calib1->ShowEntries();
      beam_y_calib2->ScanFile();
      beam_y_calib2->ShowEntries();
      
      hg_xc_calib->ScanFile();
      hg_xc_calib->ShowEntries();
      hg_x_right_calib->ScanFile();
      hg_x_right_calib->ShowEntries();
      hg_x_left_calib->ScanFile();
      hg_x_left_calib->ShowEntries();
      
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
         //cout << " " << jentry << "/" << nentries-1 << "\r"<< flush;
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
            //OH MY GOD CODE THIS AUTOMAGICALLY PLZ 07 May 2014 16:24:06 
	    //if ((PpacX[0]>=-14 && PpacX[0]<=0) && (PpacX[1]>=-14 && PpacX[1]<=0)) {goodPpacX_30s=true;} //30S gate // 1027
            if ((PpacX[0]>=-9 && PpacX[0]<=8) && (PpacX[1]>=-10 && PpacX[1]<=12)) {goodPpacX_30s=true;} //30S gate  // ALL OTHER PRODUCTION
            //if (goodRF_30s) {gate_30s=true;}
            if (goodRF_30s && goodPpacX_30s) {gate_30s=true;}
            // PPAC X was redone 20 Sep 2011 17:23:07 
	    // 29P redone with RFcal 16 Jan 2012 17:35:10 
	    if (  ((fRFcal[0]>=6) && (fRFcal[0]<=15)) || ((fRFcal[1]>=39)&&(fRFcal[1]<=48)) ) 
	       {goodRF_29p=true;} //29P gate
            if ((PpacX[0]>=-4 && PpacX[0]<=9) && (PpacX[1]>=-5 && PpacX[1]<=14)) {goodPpacX_29p=true;} //29P gate
            //if (goodRF_29p) {gate_29p=true;}
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
          // mystery ion...it should be some light ion at 0deg SSD in coincidence with 30S or 29P...check run 1027
	  // mystery ion RFcal 16 Jan 2012 17:42:46
	  // PPAC X redone 16 Jan 2012 17:45:43 
	  if (((fRFcal[0]>=34) && (fRFcal[0]<=43)) || ((fRFcal[1]>=6)&&(fRFcal[1]<=15))) {goodRF_17f=true;} 
          if ((PpacX[0]>=-9 && PpacX[0]<=8) && (PpacX[1]>=-11 && PpacX[1]<=2)) {goodPpacX_17f=true;} 
          if (goodRF_17f && goodPpacX_17f) {gate_17f=true;}
*/
          
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
	      //padE_ch[i]->Fill(fSiE[i][0]);
	      //if (fSiE[0][4]>200 || fSiE[0][5]>200 || fSiE[0][3]>200) padE_ch[i]->Fill(fSiE[i][0]);
	      if (fSiE[0][1]>200 || fSiE[0][2]>200 || fSiE[0][3]>200 || fSiE[0][4]) padE_ch[i]->Fill(fSiE[i][0]);
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
                //hSiPadEcal->Fill(i,fSiEcal[i][0]);	
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
              if (fSiEcal[i][j]>fSiEmax[i] && fSiE[i][j]<4000.) {//find maximum strip
              //if (fSiEcal[i][j]>fSiEmax[i] && fSiE[i][j]<4000.) {//find maximum strip
                //fSiEmax[i]=fSiEcal[i][j];
                fSiEmax[i]=fSiE[i][j];
		fSiEStripmax[i]=fSiE[i][j]; // goatface 29 Aug 2012 16:44:27 
		fSiEStripmax_ch[i]=j; // goatface 29 Aug 2012 16:44:27 
              }
              //Fill each strip data 
              //if(fSiE[i][j]>10) hSiStripE->Fill(i*8+j,fSiE[i][j]);
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
	       //strip_ch_gated29p[i*8+fSiEStripmax_ch[i]]->Fill(fSiEStripmax[i]);
	       if (PpacIsHit[0] && PpacIsHit[1] && fSiE[0][7]>1000 && gate_29p && fSiE[0][8]<1000) (hTargetXY_29p->Fill(TargetX,TargetY)); // PPACa center to window
	       //if (gate_29p && fSiE[0][7]>1000 && fSiEStripmax_ch[i]==7) (hTargetXY_29p->Fill(TargetX,TargetY)); // PPACa center to window
	      //if (WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)))  hTargetXYcut->Fill(TargetX,TargetY); // PPACa center to window
	      //if (gate_29p && fSiE[0][7]>1000 && fSiEStripmax_ch[i]==7 && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)))  hTargetXYcut->Fill(TargetX,TargetY); // PPACa center to window
	      if (gate_29p && fSiE[0][8]>1000 && fSiE[0][7]<1000)  hTargetXYcut->Fill(TargetX,TargetY); // PPACa center to window
            //}
            //if (gate_30s) {
            //if (gate_30s && fSiEStripmax[i]>2000) {
            //if ((fSiTcal[0]-tof[0])>292 && (fSiTcal[0]-tof[0])<296 && fSiEStripmax[i]>2000 && gate_30s) {
	       //strip_ch_gated30s[i*8+fSiEStripmax_ch[i]]->Fill(fSiEStripmax[i]);
	        if (fSiE[0][8]>1000)hPpac0TpadT_30s_ch[i]->Fill(tof[0]);
	        //if (fSiE[0][8]>1000)hPpac0TpadT_30s_ch[i]->Fill(fSiTcal[i]-tof[0]);
	      
	        if (fSiE[0][8]>1000 && fSiE[0][7]<1000 && PpacIsHit[0] && PpacIsHit[1] && gate_30s) hTargetXY_30s->Fill(TargetX,TargetY); // PPACa center to window
	        //if (gate_30s && fSiE[0][8]>1000 && fSiEStripmax_ch[i]==8) hTargetXY_30s->Fill(TargetX,TargetY); // PPACa center to window
	        //if (gate_30s && fSiE[0][8]>1000 && fSiEStripmax_ch[i]==8 && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.))) hTargetXYcut_30s->Fill(TargetX,TargetY); // PPACa center to window
	        if (gate_30s && fSiE[0][8]<1000 && fSiE[0][7]>1000 && fSiEStripmax_ch[i]==8) hTargetXYcut_30s->Fill(TargetX,TargetY); // PPACa center to window
	    }
          
	  for (UShort_t j=1;j<9;j++){
            if (gate_29p) {
	       strip_ch_gated29p[i*8+j]->Fill(fSiE[i][j]);
            }
            if (gate_30s) {
          //     if (j==8 && fSiE[i][j]>2000 && fSiE[i][7] > 2000) cout << "Multihit!"  << endl;
               strip_ch_gated30s[i*8+j]->Fill(fSiE[i][j]);
            }
            if (gate_17f) {
               strip_ch_gated17f[i*8+j]->Fill(fSiE[i][j]);
            }
	  } // end for: j
	  
	} // end for: i
      } // end if: flag_detail & flag_ssd & flag_strip & flag_ppac

      //GEM-MSTPC
      if (flag_tpc){  // analyze TPC data?
        if (!gate_30s) continue;
	//if (! (SiIsHit[2]) ) continue;
	if (! (SiIsHit[0] || SiIsHit[1] || SiIsHit[2] )) continue;
	//if (! (SiIsHit[0] || SiIsHit[1] )) continue;
	//if (! (SiIsHit[0] || SiIsHit[1] )) continue;
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
       
	   
	    // this is for viewing 500 sample tracks one-by-one for QA 
	    //if (id==29 && pulsetrack<500){
	    if (id==37 && pulsetrack<500){ //beam check
	    //if (id==115 && pulsetrack<500){
	      if ((fadc->fSample[0]!=0.) && (fadc->fSample[1]!=0.) && (fadc->fSample[2]!=0.) && (fadc->fSample[3]!=0.) && (fadc->fSample[4]!=0.)){
	        //if (flag_detail && flag_ppac) { // high gain left check
	        //if (flag_detail && SiIsHit[0] && fSiE[0][4]>120 ) { // high gain left check
	        //if (flag_detail && SiIsHit[1] && fSiE[0][3]>200 ) { // high gain left check // wtf...si 1 is 2a and fSiE 0 is 1a
	        //if (flag_detail && SiIsHit[2]) { // high gain left check
	        //if (flag_detail && flag_ppac && gate_29p && !SsdOR) { // beam check
	        //if (flag_detail && flag_ppac && gate_30s && !SsdOR) { // beam check
	        if (flag_detail && flag_ppac && gate_30s) { // beam check
	          
	          //for (Int_t j=0; j<fadc->fNSample; j++) {
                  //  hTpcPulse[pulsetrack]->Fill(fadc->fClock[j], fadc->fSample[j]);
	          //  specialtrack=true;
	          //}
	          for (Int_t j=0; j<fadc->fNSample; j++) {
                    hTpcPulse[pulsetrack]->Fill(fadc->fClock[j], fadc->fSample[j]-baseline[id][fPulseNo[id]]);
	            specialtrack=true;
	          }
	        pulsetrack++;
	        }
	      }
	    else continue; // what is this actually doing....25 Feb 2014 20:54:12 
	    }
            




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
            truepeak=false; // reset the true peak event
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
                        (fadc->fSample[j]   >= fadc->fSample[j+1]) &&
                        (fadc->fSample[j+1] > fadc->fSample[j+2]) &&
                        (fadc->fSample[j+2] > fadc->fSample[j+3]) &&
			((fadc->fSample[j]-baseline[id][fPulseNo[id]])>40 )
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
                      if (0) {
                      //if (specialtrack) {
                        cout << "pulsetrack: " << pulsetrack-1 << endl;
                        cout << "fSampleMax[" << id << "][" << fPulseNo[id] << "][" << MultiHit[id][fPulseNo[id]] << "]: " << fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] << endl;
                      }
                      MultiHit[id][fPulseNo[id]]++; // valid peak: increment the hit number
		      truepeak=true;
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
		
	        	
		//if (specialtrack && id==122 ) {
		//if (specialtrack && id==115 ) {
		if (0) {
		// if (specialtrack) {
		  cout << "pulsetrack: " << pulsetrack-1 << endl; 
		  cout << "fSampleMax[" << id << "][" << fPulseNo[id] << "][" << MultiHit[id][fPulseNo[id]] << "]: " << fSampleMax[id][fPulseNo[id]][MultiHit[id][fPulseNo[id]]] << endl;
		  //cout << "fSampleRebin[" << j << "]: " << fSampleRebin[j] << endl;
		}
	        

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
	      PadIsHitB[i][j]=false;
	      QratioB[i][j]=0.;
              XpadB[i][j]=-1000.;
              XpadBraw[i][j]=-1000.;
              YpadB[i][j]=-1000.;
              ZpadB[i][j]=0.;
              dEpadB[i][j]=0.;
              dEpadBcalib[i][j]=-10.;
              TpadB[i][j]=0.;
            }
          }
          for (UShort_t i=0;i<8;i++){
            TracksC[i]=0; 
            TracksL[i]=0; 
            TracksR[i]=0; 
            for (UShort_t j=0;j<20;j++) {
              PadIsHitC[i][j]=false;
	      QratioC[i][j]=0.;
	      XpadC[i][j]=-1000.;
              XpadR[i][j]=0.;
              XpadL[i][j]=0.;
              YpadC[i][j]=-1000.;
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
          double XpadC_average=0.;
          double YpadC_average=0.;
	  // BEAM
	   //double paddistance_new;
	   //const double theta_new=.01980872946019476473; //measured one -- seems correct 17 Apr 2013 17:32:51 atan((10.5*0.5)/265)
	   //Int_t rotpoint_new = 286 ; //pure geometry
	   //for(int i=0;i<48;i++) {
	     //paddistance_new=((i+0.5)*padsize)+rotpoint_new;
	   //}
          //cout << "new evt" << endl;
	  for(UShort_t i=0;i<96;i++){ //beam_itop: // all right Beam pads
            for(UShort_t j=0;j<96;j++){ //beam_jtop: // all left Beam pads
              if((tpc_ch[0][i]==-1)||(tpc_ch[1][j]==-1) || // tpc_ch will get us the pad number (0 to 47)
              (tpc_ch[0][i] != tpc_ch[1][j]) || 
              (PadIsHit[i] == false ) || (PadIsHit[j] == false )) continue; // only run for the same pads; make sure both pads are hit
              pad[i] = tpc_ch[0][i]; // place holder for the loop
              for(UShort_t k=1;k<=fPulseNo[i];k++){ //beam_ktop: // loop on left side pulse number
	        for(UShort_t l=1;l<2;l++){ // loop on left side hit number
	        //for(UShort_t l=1;l<=MultiHit[i][k];l++){ // loop on left side hit number
                  for(UShort_t m=1;m<=fPulseNo[j];m++){ //beam_mtop: // loop on right side pulse number
                    for(UShort_t n=1;n<2; n++){ // loop on right side hit number
                    //for(UShort_t n=1;n<=MultiHit[j][m]; n++){ // loop on right side hit number
                      if (TracksB[pad[i]] >= 20 ) continue; // don't screw up memory; TracksB[i] cannot access more than 19
                      if(!(fSampleMax[i][k][l] > 50) || !(fSampleMax[j][m][n] > 50)  || // both samples are nonzero?
                        //!(fSampleMax[i][k][l] < 375) || !(fSampleMax[j][m][n] < 375)  || // test against large 29P noise 26 Feb 2014 17:53:28 
                         !(fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	              if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
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
		      //
		     
		      if (run < 1006 && (pad[i]==1 || pad[i]==7 )) {
                        if (pad[i]==1){
		          XpadBcalib[0]=fSampleInt[i][k][l]*1.018;
		          XpadBcalib[1]=fSampleInt[j][m][n]*0.98202;
			}
                        if (pad[i]==7){
		          XpadBcalib[0]=fSampleInt[i][k][l]*1.0936;
		          XpadBcalib[1]=fSampleInt[j][m][n]*0.90641;
			}
		      }
                      else {
		        XpadBcalib[0]=beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l]);
		        XpadBcalib[1]=beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n]);
		      }
		      if (run < 1006){
		        dEpadBcalib[pad[i]][TracksB[pad[i]]] = padBgain[0][pad[i]]*(XpadBcalib[0]+XpadBcalib[1]); // integrated
		      }
		      else{
			dEpadBcalib[pad[i]][TracksB[pad[i]]] = padBgain[1][pad[i]]*(XpadBcalib[0]+XpadBcalib[1]); // integrated
		      }
		      // calibrated
		      dEpadB[pad[i]][TracksB[pad[i]]] = (XpadBcalib[0]+XpadBcalib[1]); // integrated
		      //raw
		      //dEpadB[pad[i]][TracksB[pad[i]]] = (fSampleInt[j][m][n]); // integrated
		      //dEpadB[pad[i]][TracksB[pad[i]]] = (fSampleInt[i][k][l]+fSampleInt[j][m][n]); // integrated
		      //dEpadB[pad[i]][TracksB[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]); // max
		      
		       //calibrated
		     
		     //   cout << "right" << endl;
		//	XpadBcalib[0]=beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l]);
		       // cout << "left" << endl;
		  //      XpadBcalib[1]=beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n]);
		     /* 
		      if (fSampleInt[i][k][l]>100 && fSampleInt[j][m][n]>100){
		        if (XpadBcalib[0]!=0){
		          XpadB[pad[i]][TracksB[pad[i]]] = 54.5*((XpadBcalib[0]-XpadBcalib[1])/
	               // (XpadBcalib[0]+XpadBcalib[1]));
	                  (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
		        }
	              }
		      */
		    
		    // some utterly hacked attempt at charge correction 10 Mar 2014 16:02:16 
		    // not better than anything else, fit quickly taken from pad 2 in alpha calib post
		    /*
		    if (fSampleInt[i][k][l]>fSampleInt[j][m][n] ) {
		      QratioB[pad[i]][TracksB[pad[i]]] = fSampleInt[i][k][l]/fSampleInt[j][m][n];
		      fSampleInt[i][k][l]+=343*QratioB[pad[i]][TracksB[pad[i]]]-312;

		    }

		    else {
		      QratioB[pad[i]][TracksB[pad[i]]] = fSampleInt[j][m][n]/fSampleInt[i][k][l];
		      fSampleInt[j][m][n]+=343*QratioB[pad[i]][TracksB[pad[i]]]-312;
		    }
		    */
			 //raw data
		      //integrated
		    XpadBraw[pad[i]][TracksB[pad[i]]] = 54.5*((XpadBcalib[0]-XpadBcalib[1])/
		    //XpadBraw[pad[i]][TracksB[pad[i]]] = 54.5*((fSampleInt[i][k][l]-fSampleInt[j][m][n])/
	                (fSampleInt[i][k][l]+fSampleInt[j][m][n])); 
		      
		      //if (fSampleInt[i][k][l]>100 && fSampleInt[j][m][n]>100){
		        if (XpadBcalib[0]!=0){
		          XpadB[pad[i]][TracksB[pad[i]]] = 54.5*((XpadBcalib[0]-XpadBcalib[1])/
	                 (fSampleInt[i][k][l]+fSampleInt[j][m][n]));
		        //outsourced rotation
			XpadB[pad[i]][TracksB[pad[i]]]+=padXBgeo[pad[i]];
		        }
	              
		      
		      //cout << XpadB[pad[i]][TracksB[pad[i]]] << endl; ;
		      // better above
		      //XpadB[pad[i]][TracksB[pad[i]]] = (55.*((beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l])-beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n]))/
	              //  (fSampleInt[i][k][l]+fSampleInt[j][m][n]))); 
		      //XpadB[pad[i]][TracksB[pad[i]]] = (55.*((beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l])-beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n]))/
	              //  (fSampleInt[i][k][l]+fSampleInt[j][m][n])))+shift_new[pad[i]]; 
		      //XpadB[pad[i]][TracksB[pad[i]]] = (55.*((beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l])-beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n]))/
	              //  (beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l])+beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n]))))+shift_new[pad[i]]; 
		      //XpadB[pad[i]][TracksB[pad[i]]] = (55.*((beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l])-beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n]))/
	              //  (beam_x_right_calib->Calib(pad[i],fSampleInt[i][k][l])+beam_x_left_calib->Calib(pad[i],fSampleInt[j][m][n])))); 
		      YpadB[pad[i]][TracksB[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	                + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	              // trigger time correction
		      if (run < 1006) {
		        if (!SsdOR) YpadB[pad[i]][TracksB[pad[i]]]=beam_y_delay_calib1->Calib(pad[i],YpadB[pad[i]][TracksB[pad[i]]]);
		        YpadB[pad[i]][TracksB[pad[i]]]=beam_y_calib1->Calib(pad[i],YpadB[pad[i]][TracksB[pad[i]]]);
		      }
		      else {
		        if (!SsdOR) YpadB[pad[i]][TracksB[pad[i]]]=beam_y_delay_calib2->Calib(pad[i],YpadB[pad[i]][TracksB[pad[i]]]);
		        YpadB[pad[i]][TracksB[pad[i]]]=beam_y_calib2->Calib(pad[i],YpadB[pad[i]][TracksB[pad[i]]]);
		      }

/*
1       0.0456732
2       0.0456235
3       0.0459355
4       0.0460176
5       0.0462915
6       0.0468342
7       0.0472294
8       0.0468313
9       0.0472931
10      0.0475924
11      0.0476542
12      0.0482263
13      0.0483719
14      0.0486252
15      0.048794
16      0.0490997
18      0.0497526
19      0.0500645
20      0.0500668
21      0.0504985
22      0.0505274
23      0.0508696
27      0.0518706
28      0.0521629
29      0.0524739
30      0.052824
31      0.0531443
40      0.05589
41      0.05614
42      0.0562514
43      0.057234
*/

		      //playing to get it near 0
		      //YpadB[pad[i]][TracksB[pad[i]]]*=(10./25.);
		      //YpadB[pad[i]][TracksB[pad[i]]]-=32.;
	              
		      ZpadB[pad[i]][TracksB[pad[i]]] = padsize*(pad[i])+2+(27.75+55.); // entrance window at Z=0
		      //ZpadB[pad[i]][TracksB[pad[i]]] = padsize*(pad[i]-23)-0.21; // no clear clue what this tries to do 08 Apr 2014 12:36:40 
                      //dEpadB[pad[i]][TracksB[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]); // peak 
                      
		      //if (fSampleInt[i][k][l]+fSampleInt[j][m][n]>10) // integrated
		      if (gate_30s) hPadQB_30s[pad[i]]->Fill(fSampleMax[i][k][l],fSampleMax[j][m][n]);
		      if (gate_29p) hPadQB_29p[pad[i]]->Fill(fSampleMax[i][k][l],fSampleMax[j][m][n]);
		      if (dEpadB[pad[i]][TracksB[pad[i]]] > 5 && XpadB[pad[i]][TracksB[pad[i]]]>-30. && XpadB[pad[i]][TracksB[pad[i]]]<30. && YpadB[pad[i]][TracksB[pad[i]]]>-30. && YpadB[pad[i]][TracksB[pad[i]]]<30. ) {
		        padsHitB[TracksB[pad[i]]]++;
			PadIsHitB[pad[i]][TracksB[pad[i]]]=true;
		        TracksB[pad[i]]++; // increase the number of tracks for that pad
                      } 
                    } // end for n: right side hit number
	          } // end m
	        } // end l
	      } // end k
            } // end for j: loop over left Beam pads 
          } // end for i: loop over left Beam pads 
	  //beam_end:

	  // Center
			hg_peak_holder[0]=0;
			hg_peak_holder[1]=0;
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
	  	        if(!(fSampleMax[i][k][l] > 50) || !(fSampleMax[j][m][n] > 50)  || // both samples are nonzero?
                           !(fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	  		// zero suppression threshold is often too low, making the clock values wrong!  14 Mar 2013 18:24:14 
			//if(TMath::Abs((fClockMax[i][k][l])-(fClockMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
			if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
			if(pad[i]==0||pad[i]==7)continue; // these are trash

			// this seems to work for the original ridf2root....but it's throwing events away?
			//if (fSampleMax[i][k][l] < 10 || fSampleMax[j][m][n]<10.) continue;
			//if (fSampleMax[i][k][l]+fSampleMax[j][m][n]<70.) continue;
	  		//if (pad[i]==0) continue; // for X checking, this pad data is way to crazy 11 Feb 2013 22:35:59 
			// X by peak
			//if (j==101) cout << "time left: " << fTimeMax[j][m][n]+fClockMax[j][m][n]  << " time right: " << fTimeMax[i][k][l]+fClockMax[i][k][l] << endl;
			//if (j==97 || j==100) cout << "time " << j << " : "  << fTimeMax[j][m][n] << endl;
			// kibun 01 Nov 2013 22:59:01  what about the integral instead of the peak here? 
			// or the value of 55 is somehow too small and could be larger?
			// hard to think about it without analying the data again
			//calib
			//XpadCcalib[0]=hg_x_right_calib->Calib(pad[i],fSampleMax[i][k][l]);
			//XpadCcalib[1]=hg_x_left_calib->Calib(pad[i],fSampleMax[j][m][n]);
		        /*
			XpadC[pad[i]][TracksC[pad[i]]] = 54.5*((XpadCcalib[0]-XpadCcalib[1])/
	                (XpadCcalib[0]+XpadCcalib[1]));
			*/
			//raw
			//peak
			XpadC[pad[i]][TracksC[pad[i]]] = 54.5*((fSampleMax[i][k][l]-fSampleMax[j][m][n])/
	  		  (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
			//int
			//XpadC[pad[i]][TracksC[pad[i]]] = 110.*((fSampleInt[i][k][l]-fSampleInt[j][m][n])/
			//XpadC[pad[i]][TracksC[pad[i]]] = 54.5*((fSampleInt[i][k][l]-fSampleInt[j][m][n])/
	  		//  (fSampleInt[i][k][l]+fSampleInt[j][m][n]));
			// first attempted scaling correction factor 26 Feb 2014 20:25:37 
			//if (XpadC[pad[i]][TracksC[pad[i]]]>0) XpadC[pad[i]][TracksC[pad[i]]]+=(0.2189*(XpadC[pad[i]][TracksC[pad[i]]])+1.94);
			//else XpadC[pad[i]][TracksC[pad[i]]]-=(-0.2189*(XpadC[pad[i]][TracksC[pad[i]]])+1.94);
			//if (XpadC[pad[i]][TracksC[pad[i]]]>0) XpadC[pad[i]][TracksC[pad[i]]]+=(0.1879*(XpadC[pad[i]][TracksC[pad[i]]])+2.168);
			//else XpadC[pad[i]][TracksC[pad[i]]]-=(-0.1879*(XpadC[pad[i]][TracksC[pad[i]]])+2.168);
			//if (XpadC[pad[i]][TracksC[pad[i]]]>-11.5) XpadC[pad[i]][TracksC[pad[i]]]+=(0.1879*(XpadC[pad[i]][TracksC[pad[i]]])+2.168);
			//else XpadC[pad[i]][TracksC[pad[i]]]-=(-0.1879*(XpadC[pad[i]][TracksC[pad[i]]])+2.168);
			//parabolic fit with left 10 Apr 2014 11:02:41 
			//XpadC[pad[i]][TracksC[pad[i]]]+=(6.05894e-03*pow(XpadC[pad[i]][TracksC[pad[i]]],2));
		        
			XpadC[pad[i]][TracksC[pad[i]]]=hg_xc_calib->Calib(pad[i],XpadC[pad[i]][TracksC[pad[i]]]);

			XpadC[pad[i]][TracksC[pad[i]]]+=padXCgeo[pad[i]];
		        
			if (pad[i]==6){
			  hg_peak_holder[0]=fSampleMax[i][k][l];
			  hg_peak_holder[1]=fSampleMax[j][m][n];
			}
			//if (j==101) cout << l << " " << n << " time left: " << fTimeMax[j][m][n]+fClockMax[j][m][n]  << " time right: " << fTimeMax[i][k][l]+fClockMax[i][k][l] << 
			//"\t" << fTimeMax[j][m][n] << " " << fClockMax[j][m][n] << " " << fTimeMax[i][k][l] << " "  << fClockMax[i][k][l] << "\t" << XpadC[pad[i]][TracksC[pad[i]]] << " \t" << fSampleMax[j][m][n] << "\t" << fSampleMax[i][k][l] << endl;
		        //if (pad[i]==5) XpadC[5][TracksC[5]]-=5.; // checking how pad 5 alters projection 19 Feb 2013 14:50:20 
		        //if (pad[i]==0) XpadC[0][TracksC[0]]+=85.; // checking how pad 0 inclusion alters projection 19 Feb 2013 14:50:20 
			// X by integration what was wrong with this way?!  kibun nice to know 01 Nov 2013 23:00:19 
			//XpadC[pad[i]][TracksC[pad[i]]] = 55.*((fSampleInt[i][k][l]-fSampleInt[j][m][n])/
	                //  (fSampleInt[i][k][l]+fSampleInt[j][m][n])); 
                        // alpha calibration
			if (fSiE[0][3]>130) {
			  hPadCPulse[0][pad[i]][0]->Fill(fSampleMax[i][k][l]);
			  hPadCPulse[1][pad[i]][0]->Fill(fSampleMax[j][m][n]);
			}
                        if (fSiE[0][4]>130) {
			  hPadCPulse[0][pad[i]][1]->Fill(fSampleMax[i][k][l]);
			  hPadCPulse[1][pad[i]][1]->Fill(fSampleMax[j][m][n]);
			}
                        if (fSiE[0][5]>130) {
			  hPadCPulse[0][pad[i]][2]->Fill(fSampleMax[i][k][l]);
			  hPadCPulse[1][pad[i]][2]->Fill(fSampleMax[j][m][n]);
			}
			// 30s production
			/*
			if(gate_29p||gate_30s){
			  if (SiIsHit[0]) {
			    hPadCPulse[0][pad[i]][0]->Fill(fSampleMax[i][k][l]);
			    hPadCPulse[1][pad[i]][0]->Fill(fSampleMax[j][m][n]);
			  }
                          if (SiIsHit[1]) {
			    hPadCPulse[0][pad[i]][1]->Fill(fSampleMax[i][k][l]);
			    hPadCPulse[1][pad[i]][1]->Fill(fSampleMax[j][m][n]);
			  }
			}
			*/
			YpadC[pad[i]][TracksC[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	  		  + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	  		YpadC[pad[i]][TracksC[pad[i]]]=(YpadC[pad[i]][TracksC[pad[i]]]-83.)*0.5;  // rough calib
			
			ZpadC[pad[i]][TracksC[pad[i]]] = ((55.+27.75+192.+15.5)+(padsize*(pad[i]))+2.); // Z=0 entrance window
			//dEpadC[pad[i]][TracksC[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
			// quick and dirty calibration for CNS interview Feb 2013
			//if (pad[i]==7 || pad[i]==8) continue;
			//double center_gain[8]={3.49,1.13,1,1,1,1.04,0,0};
			//dEpadC[pad[i]][TracksC[pad[i]]] = center_gain[pad[i]]*(fSampleMax[i][k][l]+fSampleMax[j][m][n]);
			//calib
			//if (fSampleMax[i][k][l]+fSampleMax[j][m][n]>40. ) //goatface 28 Feb 2013 20:20:31 emis
			//dEpadC[pad[i]][TracksC[pad[i]]] = hg_center_calib->Calib(pad[i],(fSampleMax[i][k][l]+fSampleMax[j][m][n]));
			//raw
			//dEpadC[pad[i]][TracksC[pad[i]]] = (fSampleInt[i][k][l]+fSampleInt[j][m][n]);
			dEpadC[pad[i]][TracksC[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
			padsHitC[TracksC[pad[i]]]++;
		        PadIsHitC[pad[i]][TracksC[pad[i]]]=true;
			XpadC_average+=XpadC[pad[i]][TracksC[pad[i]]];
		        YpadC_average+=YpadC[pad[i]][TracksC[pad[i]]];
			if (gate_30s) hPadQC[pad[i]]->Fill(fSampleMax[i][k][l],fSampleMax[j][m][n]);
                        TracksC[pad[i]]++; // increase the number of tracks for that pad
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
	  	        if(!(fSampleMax[i][k][l] > 50) || !(fSampleMax[j][m][n] > 50)  || // both samples are nonzero?
                           (!fClockMax[i][k][l] >0.1) || !(fClockMax[j][m][n] > 0.1)) continue;  // both clocks are nonzero?
	  		if(TMath::Abs((fClockMax[i][k][l]+fTimeMax[i][k][l])-(fClockMax[j][m][n]+fTimeMax[j][m][n])) >= 5 ) continue; // the clock times are similar?
			if (fSampleMax[i][k][l]+fSampleMax[j][m][n]<80.) continue;
	  		// what's the units of Z?
			ZpadL[pad[i]][TracksL[pad[i]]] = (117.*((fSampleInt[j][m][n]-fSampleInt[i][k][l])/
	  		  (fSampleInt[i][k][l]+fSampleInt[j][m][n])))+117.+89.; // 89 is from geometry3.C
			//ZpadL[pad[i]][TracksL[pad[i]]] = 117.5*((fSampleMax[j][m][n]-fSampleMax[i][k][l])/
	  		//  (fSampleMax[i][k][l]+fSampleMax[j][m][n]));
                        // alpha calibration
			for (int ii=0;ii<8;ii++){
			  if (fSiE[2][ii+1]>130 && ZpadL[pad[i]][TracksL[pad[i]]]<210.) {
			    hPadLPulseLow[0][pad[i]][ii]->Fill(fSampleMax[i][k][l]);
			    hPadLPulseLow[1][pad[i]][ii]->Fill(fSampleMax[j][m][n]);
			  }
			  if (fSiE[2][ii+1]>130 && ZpadL[pad[i]][TracksL[pad[i]]]>=210.) {
			    hPadLPulseHigh[0][pad[i]][ii]->Fill(fSampleMax[i][k][l]);
			    hPadLPulseHigh[1][pad[i]][ii]->Fill(fSampleMax[j][m][n]);
			  }
			}
                        YpadL[pad[i]][TracksL[pad[i]]] = 0.5*(fClockMax[i][k][l] 
	  		  + fTimeMax[i][k][l]+fClockMax[j][m][n] + fTimeMax[j][m][n]);
	  		XpadL[pad[i]][TracksL[pad[i]]] = -7.621 - padsize*(pad[i]); 
			XpadL[pad[i]][TracksL[pad[i]]]+=(6.05894e-03*pow(XpadL[pad[i]][TracksL[pad[i]]],2));
                        dEpadL[pad[i]][TracksL[pad[i]]] = (fSampleMax[i][k][l]+fSampleMax[j][m][n]);
			padsHitL[TracksL[pad[i]]]++;
		        hPadQL[pad[i]]->Fill(fSampleMax[i][k][l],fSampleMax[j][m][n]);
                        TracksL[pad[i]]++; // increase the number of tracks for that pad
                      } // end for n: right side hit number
            } // end for j: loop over downstream Left pads 
          
	  // right side not used in 30S
	  /*
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
			padsHitR[TracksR[pad[i]]]++;
                        TracksR[pad[i]]++; // increase the number of tracks for that pad
                    } // end for n: right side hit number
            } // end for j: loop over downstream Right pads 
            */

	  //beam
			 //PpacX[1]+=0.5;
	  
	  
	  // cleanup init
	  //int nPointsB=2, nPointsP=0;
	  //double padsHitC_double = padsHitC[0];
	  nPointsB=2; nPointsP=0;
	  padsHitC_double = padsHitC[0];
	  TargetX=TargetY=-999.;
	  XpadC_average=XpadC_average/padsHitC_double;
	  YpadC_average=YpadC_average/padsHitC_double;
	  event=0;
	  scatter_event = false;
	  //if (SsdOR && (SiIsHit[2]) && gate_30s && (padsHitB[0]>5)){
	  //if (SsdOR && (SiIsHit[2]) && gate_30s && (padsHitL[0]>2) && (padsHitB[0]>5) && (eventcounter<400)){
	  if (SsdOR && (SiIsHit[0]||SiIsHit[1]||SiIsHit[2]) && gate_30s && padsHitB[0] > 5){
	  //if (SsdOR && (SiIsHit[0]||SiIsHit[1]) && gate_30s && padsHitB[0] > 5){
	  //if (SsdOR && (SiIsHit[0]||SiIsHit[1]) && gate_30s && (padsHitC[0]>2) && (padsHitB[0]>5)){
	  //if (SsdOR && (SiIsHit[0]||SiIsHit[1]) && gate_30s && (padsHitC[0]<2) && (padsHitB[0]>5)){
	  //if (SsdOR && (SiIsHit[0]||SiIsHit[1]) && gate_30s && (padsHitB[0]>5) && (eventcounter<400)){
	    //cleanup init
	    //double zppac=0.;
	    for (int i=0;i<gBraggB_evt[eventcounter]->GetN();i++) gBraggB_evt[eventcounter]->RemovePoint(i);
	    for (int i=0;i<gPadXB_evt[eventcounter]->GetN();i++) gPadXB_evt[eventcounter]->RemovePoint(i);
	    for (int i=0;i<gPadYB_evt[eventcounter]->GetN();i++) gPadYB_evt[eventcounter]->RemovePoint(i);
	    for (int i=0;i<gPadXC_evt[eventcounter]->GetN();i++) gPadXC_evt[eventcounter]->RemovePoint(i);
	    for (int i=0;i<gPadYC_evt[eventcounter]->GetN();i++) gPadYC_evt[eventcounter]->RemovePoint(i);
	    zppac=0.;
	    for (int i=0;i<21;i++){
	          	 zppac=(i*padsize)+(padsize/2);
	          	 TargetX=fTargetX->Eval((zppac+452.));
	          	 gPpacX_evt[eventcounter]->SetPoint(nPointsP,zppac,TargetX);
			 TargetY=fTargetY->Eval((zppac+452.)); // PPACa center to window, inactive region, plus pad
	          	 gPpacY_evt[eventcounter]->SetPoint(nPointsP,zppac,TargetY);
	          	 nPointsP++;
            }
	    // cleanup init
	    //double zmstpc=0.;
	    zmstpc=0.;
	    // give dE back to pads that don't have data.
	    for(UShort_t i=0;i<48;i++){ // loop over all pads (now by pad number)
	      for(UShort_t j=0;j<1;j++){ // loop for first track only
		if (i==17 && dEpadBcalib[i-1][j]>-1 && dEpadBcalib[i+1][j]>-1) dEpadBcalib[i][j]=(dEpadBcalib[i-1][j]+dEpadBcalib[i+1][j])*0.5;
		if (((i==24)||(i==25)||(i==26)) && dEpadBcalib[23][j]>-1 && dEpadBcalib[27][j]>-1) dEpadBcalib[i][j]=(dEpadBcalib[23][j]+dEpadBcalib[27][j])*0.5;
		if ((i>0) && (i!=17) && (i!=24) && (i!=25) && (i!=26) && (i<32)){
		  if ((dEpadBcalib[i][j]<0) && (dEpadBcalib[i-1][j]>-1) && (dEpadBcalib[i+1][j]>-1)){
		    dEpadBcalib[i][j]=(dEpadBcalib[i-1][j]+dEpadBcalib[i+1][j])*0.5;
		    //dEpadBcalib[i][j]=(dEpadBcalib[i-1][j]+dEpadBcalib[i+1][j])*0.5;
		  }
		}
	      }
	    }
	    for(UShort_t i=0;i<48;i++){ // loop over all pads (now by pad number)
	      if (!PadIsHitB[i][0]) continue;
	      for(UShort_t j=0;j<1;j++){ // loop for first track only
		zmstpc=82.75+(i*padsize)+(padsize/2);
	        TargetX=fTargetX->Eval((zmstpc+452.));
	        gPpacX_evt[eventcounter]->SetPoint(nPointsP,zmstpc,TargetX);
                TargetY=fTargetY->Eval((zmstpc+452.)); // PPACa center to window, inactive region, plus pad
	        gPpacY_evt[eventcounter]->SetPoint(nPointsP,zmstpc,TargetY);
	        nPointsP++;
		if (dEpadBcalib[i][j]>-1.) gBraggB_evt[eventcounter]->SetPoint(nPointsB,zmstpc,dEpadBcalib[i][j]);
	        if (XpadB[i][j] > -999.) { 
		  gPadXB_evt[eventcounter]->SetPoint(nPointsB,zmstpc,XpadB[i][j]);
		  gPadXB_evt[eventcounter]->SetPointError(nPointsB,2.,(-0.000322*dEpadB[i][j]+6.48));
		}
	        if (YpadB[i][j] > -999.) { 
		  gPadYB_evt[eventcounter]->SetPoint(nPointsB,zmstpc,YpadB[i][j]);
		  gPadYB_evt[eventcounter]->SetPointError(nPointsB,2.,0.5);
		}
	        nPointsB++;
              }//j
	    } //i
	    gPadXB_evt[eventcounter]->SetPoint(0,-452.,PpacX[0]);
	    gPadXB_evt[eventcounter]->SetPointError(0.,1.,0.9);
	    gPadXB_evt[eventcounter]->SetPoint(1,-296.,PpacX[1]);
	    gPadXB_evt[eventcounter]->SetPointError(1,1.,0.9);
	    
	    gPadYB_evt[eventcounter]->SetPoint(0,-452.,PpacY[0]);
	    gPadYB_evt[eventcounter]->SetPointError(0.,1.,0.9);
	    gPadYB_evt[eventcounter]->SetPoint(1,-296.,PpacY[1]);
	    gPadYB_evt[eventcounter]->SetPointError(1,1.,0.9);
	    
	    // cleanup init
	    //int padcpoints=0;
	    padcpoints=0;
	    if(padsHitC[0]>2){ 
	      gPadXC_evt[eventcounter]->SetPoint(padcpoints,306.25,XpadC_average);
	      gPadYC_evt[eventcounter]->SetPoint(padcpoints,306.25,YpadC_average);
            }


	    //find the track here! 27 Jul 2014 13:13:44 
            
	    //global declarations we can move later
	    // cleanup init
	    //double chi2min[4], chi2nu[4];
            //double xmin=83.,xmed,xmax,xmedbest;
            //double chi2minGlobal[2]; 
            
            //reset chi-squared values
            for (int j=0;j<4;j++){
              chi2min[j]=100.;
              chi2nu[j]=100.;
            }
            xmin=84.75;
            xmax=260.75;
            chi2minGlobal[0]=200.;
            chi2minGlobal[1]=200.;

            //upstream scattering
            TF1 *f1b_up = new TF1("f1b_up","[0]+x*[1]",xmin,xmax);
            gPadXB_evt[eventcounter]->Fit("f1b_up","RQ");
            TF1 *f2b_up = new TF1("f2b_up","[0]+x*[1]",xmin,xmax);
            gPadYB_evt[eventcounter]->Fit("f2b_up","RQ");
            if (f1b_up->GetNDF()>0) {
              chi2nu[0]=f1b_up->GetChisquare()/f1b_up->GetNDF();
            }
            else {
              chi2nu[0]=100.;
            }
            if (f2b_up->GetNDF()>0) {
              chi2nu[1]=f2b_up->GetChisquare()/f2b_up->GetNDF();
            }
            else {
              chi2nu[1]=100.;
            }
	    if (chi2minGlobal[0]>(chi2nu[0]+chi2nu[1])){
              chi2minGlobal[0]=(chi2nu[0]+chi2nu[1]);
              for (int j=0;j<2;j++){
                chi2min[j]=chi2nu[j];
              }
            }
	    delete f1b_up;
	    delete f2b_up;
            xmedbest=-100.;
            // cleanup init
	    //int pxmedlow=3,pxmedhigh=(gPadYB_evt[eventcounter]->GetN()-1);
            //double tempx[2],tempy[2];
            pxmedlow=3;
            //pxmedlow=4;
	    pxmedhigh=(gPadYB_evt[eventcounter]->GetN()-1);
            // cleanup init
            //double tempx[2],tempy[2]; // wtf we declared this 4 lines earlier!
            pxmedlow=gBraggB_evt[eventcounter]->GetPoint(pxmedlow,tempx[0],tempy[0]);
            pxmedhigh=gBraggB_evt[eventcounter]->GetPoint(pxmedhigh,tempx[1],tempy[1]);
            //cout << "point low " << pxmedlow << "\tx: " << tempx[0]  << "\ty: " << tempy[0] << endl;
            //cout << "point high " << pxmedhigh << "\tx: " << tempx[1]  << "\ty: " << tempy[1] << endl;
            // cleanup init
            //int temppoint=0;
            temppoint=0;
            //cout << "Npoints: " << gPadYB_evt[eventcounter]->GetN() << endl;
            //xmin=80.; // this is w/o PPAC data
            xmin=-500.; // this is with PPAC data
            //scattering in active region
            for (int point=pxmedlow;point<pxmedhigh;point++){
              temppoint=gBraggB_evt[eventcounter]->GetPoint(point,tempx[0],tempy[0]);
              xmed=tempx[0];
              TF1 *f1b = new TF1("f1b","[0]+x*[1]",xmin,xmed);
              gPadXB_evt[eventcounter]->Fit("f1b","RQB");

              // cleanup init
              //int junk=0;
              junk=0;
              TF1 *f1B = new TF1("f1B","[0]+x*[1]",xmed,xmax);
              
	      if(padsHitC[0]>2){ 
	        junk = gPadXC_evt[eventcounter]->GetPoint(0,tempx[0],tempy[0]);
                if (gPadXB_evt[eventcounter]->Eval(xmed) < tempy[0]){
	          f1B->SetParLimits(1,-10.,0.);
                }else
                {
                  f1B->SetParLimits(1,0.,10.);
                }
              }
	      gPadXB_evt[eventcounter]->Fit("f1B","RQB");
              TF1 *f2b = new TF1("f2b","[0]+x*[1]",xmin,xmed);
              
	      gPadYB_evt[eventcounter]->Fit("f2b","RQB");

              TF1 *f2B = new TF1("f2B","[0]+x*[1]",xmed,xmax);
	      if(padsHitC[0]>2){ 
                junk = gPadYC_evt[eventcounter]->GetPoint(0,tempx[0],tempy[0]);
                if (gPadYB_evt[eventcounter]->Eval(xmed) < tempy[0]){
	          f2B->SetParLimits(1,-10.,0.);
                }else
                {
                  f2B->SetParLimits(1,0.,10.);
                }
	      }
              gPadYB_evt[eventcounter]->Fit("f2B","RQ");

              // cleanup init
              //double dE[2];
              //int temppe;
              //double tempxe[2],tempee[2];
              temppe=gBraggB_evt[eventcounter]->GetPoint(point-1,tempxe[0],tempee[0]);
              temppe=gBraggB_evt[eventcounter]->GetPoint(point,tempxe[1],tempee[1]);
              dE[0]=(tempee[0]+tempee[1])*0.5;
              temppe=gBraggB_evt[eventcounter]->GetPoint(point+1,tempxe[0],tempee[0]);
              temppe=gBraggB_evt[eventcounter]->GetPoint(point+2,tempxe[1],tempee[1]);
              dE[1]=(tempee[0]+tempee[1])*0.5;
	     /* 
              if (dE[0]>(dE[1])){
	        delete f1b;
	        delete f1B;
	        delete f2b;
	        delete f2B;
	        continue;
              }
	      */
              //if ((dE[0]/(dE[0]+dE[1]))>(dE[1])/(dE[0]+dE[1])) continue;
              //double dEratio=dE[0]/dE[1];
              if ((f1b->GetNDF()>0) && (f1B->GetNDF()>0)) {
                chi2nu[2]=(f1b->GetChisquare()+f1B->GetChisquare())/(f1b->GetNDF()+f1B->GetNDF()+0);
              }
              else {
                chi2nu[2]=100.;
              }

              if ((f2b->GetNDF()>0)&& (f2B->GetNDF()>0)) {
                chi2nu[3]=(f2b->GetChisquare()+f2B->GetChisquare())/(f2b->GetNDF()+f2B->GetNDF()+0);
              }
              else {
                chi2nu[3]=100.;
              }


              if (chi2minGlobal[1]>(chi2nu[2]+chi2nu[3])){
                chi2minGlobal[1]=(chi2nu[2]+chi2nu[3]);
                for (int j=2;j<4;j++){
                  chi2min[j]=chi2nu[j];
                }
                xmedbest=xmed;
              }
	      delete f1b;
	      delete f1B;
	      delete f2b;
	      delete f2B;
            }
            if ((chi2minGlobal[0] < chi2minGlobal[1]) && (chi2min[2]>1 || chi2min[3]>1) ){
              //cout << endl << "UPSTREAM SCATTERING" << endl << endl;
              xmedbest=-100.;
            }
            else if (xmedbest < 212 && xmedbest>0) { // ACTIVE TARGET SCATTERING!! For now omit pad 32 and higher
              //cout << endl << "ACTIVE TARGET SCATTERING" << endl << endl;
              //cout << xmedbest ;
	      // how we converted pad to Z so we can reverse it
	      //zmstpc=82.75+(i*padsize)+(padsize/2);
              // cleanup init
	      //float pad_scatter_tmp = (xmedbest - 82.75 - (padsize/2)) / padsize;
	      //int pad_scatter = pad_scatter_tmp;
	      //double dEpadBsum = 36.9; // 17.3 um mylar; 6.8 um kapton; 194 torr, 213 K 82.75+4 = 86.75 mm HeCO2
	      pad_scatter_tmp = (xmedbest - 82.75 - (padsize/2)) / padsize;
	      pad_scatter = pad_scatter_tmp;
	      // "randomizer" for missing pads
	      if (pad_scatter==18) pad_scatter-=(jentry%2);
	      if (pad_scatter==27) pad_scatter-=(jentry%4);
	      gate_30s_bragg = 0;
	      for (int i=1;i<=pad_scatter;i++){
	        if (i==17||i==24||i==25||i==26||i>31) continue;
		if ( (PadIsHitB[i][0]) && !(dEpadB[i][0]>dEpadB_gate_30s[0][i] && dEpadB[i][0]<dEpadB_gate_30s[1][i])) {
		  gate_30s_bragg++; 
		  //cout << pad_scatter << "\t" << i << "\t" << dEpadB[i][0] << endl;
		}
	      }
	      //cout << "EVENT: " << pad_scatter << "\t" << gate_30s_bragg << endl;
	      //if (gate_30s_bragg<=1){
	      if (gate_30s_bragg<1){
	      //if (1){
	      //if (pad_scatter<5) cout << dEpadB[1][0] << "\t" << dEpadB[2][0] << "\t" << dEpadB[3][0] << "\t" << dEpadB[4][0] << endl;
	      //cout << dEpadB[1][0] << endl;
	      //if ((pad_scatter-1)/4 >= gate_30s_bragg){
	      //if ((fSiE[0][0]< 910) || (fSiE[0][0]>980)){
	      //if ((gate_30s_bragg<1) && (fSiE[0][0]< 910) || (fSiE[0][0]>980)){
		//cout << "EVENT: " << pad_scatter << "\t" << gate_30s_bragg << endl;
		dEpadBsum = 36.9; // 17.3 um mylar; 6.8 um kapton; 194 torr, 213 K 82.75+4 = 86.75 mm HeCO2
	        for (int i=1;i<pad_scatter;i++){
	        //for (int i=1;i<=pad_scatter;i++){
	          if (dEpadBcalib[i][0]>-9) dEpadBsum-=dEpadBcalib[i][0];
                }
		dEpadBsum-=(0.5*dEpadBcalib[pad_scatter][0]);
	        //beam->SetEnergy(dEpadBsum);
		//cout << "Energy: " << dEpadBsum << endl;
		//reaction.CalculateKinematics();
		//double e1,e2;
		//double theta=10.5;
		//cout << "testing: " << reaction.GetELab(3,theta,3,e1,e2) << endl;
		//cout << "energy: " << e1 << endl;
		dEpadBsum=(dEpadBsum*(4000/34)); // convert from E_lab to E_cm and MeV to keV
	        if (cutg->IsInside(dEpadBsum,pad_scatter)){
	          hExcitation->Fill(dEpadBsum);
	          hExPad->Fill(dEpadBsum,pad_scatter);
                  if (!SiIsHit[1] && !SiIsHit[2] && SiIsHit[0]) {
		    hSiPadEcal->Fill(1.,fSiE[0][0]);	
		    scatter_event = true;
		    hEx_Si1a->Fill(dEpadBsum,fSiE[0][0]);
		  }
	        }
                //cout << "\t" << pad_scatter << endl;
	        //
              }
            }
	    
	    
	    //eventcounter++;
	  }
 

// KINEMATIC SOLUTION CALCULATOR

	  nPointsB=2; nPointsP=0;
	  padsHitC_double = padsHitC[0];
	  TargetX=TargetY=-999.;
	  //XpadC_average=XpadC_average/padsHitC_double;
	  //YpadC_average=YpadC_average/padsHitC_double;
	  event=0;
	  scatter_event = false;
	  //if (SsdOR && (SiIsHit[2]) && gate_30s && (padsHitB[0]>5)){
	  //if (SsdOR && (SiIsHit[2]) && gate_30s && (padsHitL[0]>2) && (padsHitB[0]>5) && (eventcounter<400)){
	  //if (SsdOR && (SiIsHit[0]||SiIsHit[1]||SiIsHit[2]) && gate_30s && padsHitB[0] > 5){
	  //if (SsdOR && (SiIsHit[0]||SiIsHit[1]) && gate_30s && padsHitB[0] > 5){
	  if (SsdOR && SiIsHit[0] && gate_30s && (padsHitC[0]>2) && (padsHitB[0]>5)){
	  //if (SsdOR && (SiIsHit[0]||SiIsHit[1]) && gate_30s && (padsHitC[0]<2) && (padsHitB[0]>5)){
	  //if (SsdOR && (SiIsHit[0]||SiIsHit[1]) && gate_30s && (padsHitB[0]>5) && (eventcounter<400)){
	    //cleanup init
	    //double zppac=0.;
	    for (int i=0;i<gBraggB_evt[eventcounter]->GetN();i++) gBraggB_evt[eventcounter]->RemovePoint(i);
	    for (int i=0;i<gPadXB_evt[eventcounter]->GetN();i++) gPadXB_evt[eventcounter]->RemovePoint(i);
	    for (int i=0;i<gPadYB_evt[eventcounter]->GetN();i++) gPadYB_evt[eventcounter]->RemovePoint(i);
	    for (int i=0;i<gPadXC_evt[eventcounter]->GetN();i++) gPadXC_evt[eventcounter]->RemovePoint(i);
	    for (int i=0;i<gPadYC_evt[eventcounter]->GetN();i++) gPadYC_evt[eventcounter]->RemovePoint(i);
	    zppac=0.;
	    // cleanup init
	    //double zmstpc=0.;
	    zmstpc=0.;
	    // give dE back to pads that don't have data.
	    for(UShort_t i=0;i<48;i++){ // loop over all pads (now by pad number)
	      for(UShort_t j=0;j<1;j++){ // loop for first track only
		if (i==17 && dEpadBcalib[i-1][j]>-1 && dEpadBcalib[i+1][j]>-1) dEpadBcalib[i][j]=(dEpadBcalib[i-1][j]+dEpadBcalib[i+1][j])*0.5;
		if (((i==24)||(i==25)||(i==26)) && dEpadBcalib[23][j]>-1 && dEpadBcalib[27][j]>-1) dEpadBcalib[i][j]=(dEpadBcalib[23][j]+dEpadBcalib[27][j])*0.5;
		if ((i>0) && (i!=17) && (i!=24) && (i!=25) && (i!=26) && (i<32)){
		  if ((dEpadBcalib[i][j]<0) && (dEpadBcalib[i-1][j]>-1) && (dEpadBcalib[i+1][j]>-1)){
		    dEpadBcalib[i][j]=(dEpadBcalib[i-1][j]+dEpadBcalib[i+1][j])*0.5;
		    //dEpadBcalib[i][j]=(dEpadBcalib[i-1][j]+dEpadBcalib[i+1][j])*0.5;
		  }
		}
	      }
	    }
	    for(UShort_t i=0;i<48;i++){ // loop over all pads (now by pad number)
	      if (!PadIsHitB[i][0]) continue;
	      for(UShort_t j=0;j<1;j++){ // loop for first track only
		zmstpc=82.75+(i*padsize)+(padsize/2);
		if (dEpadBcalib[i][j]>-1.) gBraggB_evt[eventcounter]->SetPoint(nPointsB,zmstpc,dEpadBcalib[i][j]);
	        if (XpadB[i][j] > -999.) { 
		  gPadXB_evt[eventcounter]->SetPoint(nPointsB,zmstpc,XpadB[i][j]);
		  gPadXB_evt[eventcounter]->SetPointError(nPointsB,2.,(-0.000322*dEpadB[i][j]+6.48));
		}
	        if (YpadB[i][j] > -999.) { 
		  gPadYB_evt[eventcounter]->SetPoint(nPointsB,zmstpc,YpadB[i][j]);
		  gPadYB_evt[eventcounter]->SetPointError(nPointsB,2.,0.5);
		}
	        nPointsB++;
              }//j
	    } //i
	    gPadXB_evt[eventcounter]->SetPoint(0,-452.,PpacX[0]);
	    gPadXB_evt[eventcounter]->SetPointError(0.,1.,0.9);
	    gPadXB_evt[eventcounter]->SetPoint(1,-296.,PpacX[1]);
	    gPadXB_evt[eventcounter]->SetPointError(1,1.,0.9);
	    
	    gPadYB_evt[eventcounter]->SetPoint(0,-452.,PpacY[0]);
	    gPadYB_evt[eventcounter]->SetPointError(0.,1.,0.9);
	    gPadYB_evt[eventcounter]->SetPoint(1,-296.,PpacY[1]);
	    gPadYB_evt[eventcounter]->SetPointError(1,1.,0.9);
	    
	    padcpoints=0;
	    if(padsHitC[0]>2){ 
	      gPadXC_evt[eventcounter]->SetPoint(padcpoints,306.25,XpadC_average);
	      gPadYC_evt[eventcounter]->SetPoint(padcpoints,306.25,YpadC_average);
            }


	    //find the track here! 27 Jul 2014 13:13:44 
            
            xmax=260.75;

            //upstream scattering
            // cleanup init
	    //int pxmedlow=3,pxmedhigh=(gPadYB_evt[eventcounter]->GetN()-1);
            //double tempx[2],tempy[2];
            pxmedlow=1;
            //pxmedlow=3;
            //pxmedlow=4;
	    pxmedhigh=(gPadYB_evt[eventcounter]->GetN()-1);
            // cleanup init
            //double tempx[2],tempy[2]; // wtf we declared this 4 lines earlier!
            pxmedlow=gBraggB_evt[eventcounter]->GetPoint(pxmedlow,tempx[0],tempy[0]);
            pxmedhigh=gBraggB_evt[eventcounter]->GetPoint(pxmedhigh,tempx[1],tempy[1]);
            //cout << "point low " << pxmedlow << "\tx: " << tempx[0]  << "\ty: " << tempy[0] << endl;
            //cout << "point high " << pxmedhigh << "\tx: " << tempx[1]  << "\ty: " << tempy[1] << endl;
            // cleanup init
            //int temppoint=0;
            temppoint=0;
            //cout << "Npoints: " << gPadYB_evt[eventcounter]->GetN() << endl;
            //xmin=80.; // this is w/o PPAC data
            xmin=-500.; // this is with PPAC data
            
	    double E_cm_kinematics = -10;
	    double delta_E_alpha_kinematics = 0.5; // 500 keV
	    
	    //scattering in active region
            
	    for (int point=pxmedlow;point<pxmedhigh;point++){
              temppoint=gBraggB_evt[eventcounter]->GetPoint(point,tempx[0],tempy[0]);
              xmed=tempx[0];
              TF1 *f1b = new TF1("f1b","[0]+x*[1]",xmin,xmed);
              gPadXB_evt[eventcounter]->Fit("f1b","RQB");
              //f1b->GetParameters(&p[0]);

              junk=0;
              
              TF1 *f2b = new TF1("f2b","[0]+x*[1]",xmin,xmed);
              
	      gPadYB_evt[eventcounter]->Fit("f2b","RQB");

	      // points 0, 1, 2 | x, y, z
	      // where: point 0 is the beam at entrance window
	      // 	point 1 is the scattering location
	      // 	point 2 is the outgoing alpha position detected
	      // 	point 3 is the SSD hitting position (extrapolated)
	      float vertex[4][3]; 
	      vertex[0][0]=f1b->Eval(0.); //point 0 x 
	      vertex[1][0]=f1b->Eval(xmed); //point 1 x "beamx"
	      
	      vertex[0][1]=f2b->Eval(0.); // point 0 y
	      vertex[1][1]=f2b->Eval(xmed); // point 1 y "beamy"
	      
	      vertex[0][2]= 0.; // point 0 z
	      // when we considered pad 0 as the point 0
	      //vertex[0][2]= fPadB[0]->Eval(vertex[0][0]); // point 0 z
	      vertex[1][2]= fPadB[point]->Eval(vertex[1][0]); // point 1 z "beamz" SHOULD BE NEARLY xmed
	      
	      vertex[2][0]=XpadC_average; // point 2 x "alphax"
	      vertex[2][1]=YpadC_average; // point 2 y "alphay"
              vertex[2][2]=(fPadC[4]->Eval(vertex[2][0])+fPadC[3]->Eval(vertex[2][0]))/2. ; //point 2 z "alphaz"
	      
              TGraph *gAlphaTrack;
              gAlphaTrack = new TGraph;
              gAlphaTrack->SetPoint(0,vertex[1][0],vertex[1][2]);
              gAlphaTrack->SetPoint(1,vertex[2][0],vertex[2][2]);
              //gAlphaTrack->Draw("L9SAME");
              gAlphaTrack->Fit("fLinear","Q");
              
              //TF1 *fAlphaTrack = new TF1("fAlphaTrack","[0]+x*[1]",-500,500);
              finter_f1=fLinear;
              finter_f2=fSsdC;
              vertex[3][0]=fint->GetMinimumX(); // point 3 x "ssdx"
              vertex[3][2]=fSsdC->Eval(vertex[3][0]); // point 3 z "ssdz"
              //cout << "X: " << vertex[3][0] << endl;
              //cout << "Z: " << vertex[3][2] << endl;

              double xAlpha1aStripMin = 100;
              int chAlpha1aStrip = -1;
              for (int Ssd1aStripHit=0; Ssd1aStripHit<8; Ssd1aStripHit++){
                if (abs(xStripSsdC[Ssd1aStripHit]-vertex[3][0])<xAlpha1aStripMin){
                  chAlpha1aStrip = Ssd1aStripHit;
                  xAlpha1aStripMin = abs(xStripSsdC[Ssd1aStripHit]-vertex[3][0]);
                  //cout << Ssd1aStripHit << "\t" << xAlpha1aStripMin << endl;
                }
              }
              chAlpha1aStrip++; // this counter started at 0, but we want to start at 1
              //cout << "Strip Hit: " << chAlpha1aStrip << endl;
              
	      TGraph *gAlphaTrackZY; 
              gAlphaTrackZY = new TGraph;
              gAlphaTrackZY->SetPoint(0,vertex[1][2],vertex[1][1]);
              gAlphaTrackZY->SetPoint(1,vertex[2][2],vertex[2][1]);
              //gAlphaTrackZY->Draw("L9SAME");
              gAlphaTrackZY->Fit("fLinear","Q");
              vertex[3][1]=fLinear->Eval(vertex[3][2]);
              //cout << "Y: " << vertex[3][1] << endl;

              double alphaTrackLength;
              alphaTrackLength = sqrt ( pow(vertex[3][0]-vertex[1][0],2.) + pow(vertex[3][1]-vertex[1][1],2.) + pow(vertex[3][2]-vertex[1][2],2.) );
              //cout << "alphaTrackLength: " << alphaTrackLength << endl;
//end here
              
	      double vertex_leg[3];
	      vertex_leg[0] = sqrt ( pow((vertex[1][0]-vertex[0][0]),2) + pow((vertex[1][1]-vertex[0][1]),2) +  pow((vertex[1][2]-vertex[0][2]),2));
	      vertex_leg[1] = sqrt ( pow((vertex[2][0]-vertex[1][0]),2) + pow((vertex[2][1]-vertex[1][1]),2) +  pow((vertex[2][2]-vertex[1][2]),2));
	      vertex_leg[2] = sqrt ( pow((vertex[2][0]-vertex[0][0]),2) + pow((vertex[2][1]-vertex[0][1]),2) +  pow((vertex[2][2]-vertex[0][2]),2));
	      double theta;

	      theta = acos ( (pow(vertex_leg[0],2) + pow(vertex_leg[1],2) - pow(vertex_leg[2],2)) / (2*vertex_leg[0]*vertex_leg[1]) );

	      dEpadBsum = 36.9; // 17.3 um mylar; 6.8 um kapton; 194 torr, 213 K 82.75+4 = 86.75 mm HeCO2
	      for (int i=1;i<point;i++){
	      //for (int i=1;i<=pad_scatter;i++){
	        if (dEpadBcalib[i][0]>-9) dEpadBsum-=dEpadBcalib[i][0];
              }
	      dEpadBsum-=(0.5*dEpadBcalib[pad_scatter][0]);
	      beam->SetEnergy(dEpadBsum);
	      //cout << "Energy: " << dEpadBsum << endl;
	      reaction.CalculateKinematics();
	      double e1,e2;
	      //double theta=10.5;
	      int nsol = reaction.GetELab(3,theta,3,e1,e2);
	      //cout << "testing: " << reaction.GetELab(3,theta,3,e1,e2) << endl;
	      //cout << "energy: " << e1 << endl;
	      
	      //enewz for alpha goes here.  
	      //initial energy is e1
	      //distance is alphaTrackLength
	      //heco2 pressure is 194 w/ pressure 280k
	      //use the 
              e_alpha_send = e1/m_alpha;
	      thick_gas = alphaTrackLength;
	      enewzsub_(&z_alpha, &m_alpha, &e_alpha_send, matter_gas, &unit_pressure_gas,
	                      &pressure_gas, &temperature_gas, &unit_thick_gas, &thick_gas, &e_alpha_return);
              e_alpha_send = (e_alpha_return);
              enewzsub_(&z_alpha, &m_alpha, &e_alpha_send, matter_mylar, &unit_pressure_mylar,
                &pressure_mylar, &temperature_mylar, &unit_thick_mylar, &thick_mylar, &e_alpha_return);
	      energy_alpha_calc = e_alpha_return*m_alpha;
	     
              double fSiE1acal=ssd_padE_calib->CalibMeV(chAlpha1aStrip,fSiE[0][0]);

	      //now we need to do something with the best solution and check for it!!
              
              if ( fSiE1acal>0.1 && energy_alpha_calc>0.1 && (abs(fSiE1acal-energy_alpha_calc)<delta_E_alpha_kinematics)){
	        E_cm_kinematics = dEpadBsum;
	        delta_E_alpha_kinematics = abs(fSiE1acal-energy_alpha_calc);
                cout << "FOUND EVENT AT E_CM: " << dEpadBsum << " with minimum: " << delta_E_alpha_kinematics<< endl;
              }

	      delete f1b;
	      delete f2b;
	      delete gAlphaTrack;
              delete gAlphaTrackZY;
            }
            // change this! 30 Sep 2014 17:38:55 
	    if (E_cm_kinematics>-1){ // reasonably small error for kinematic solution!
	    cout << endl;
	    //if (delta_E_alpha_kinematics > -1){ // reasonably small error for kinematic solution!
		dEpadBsum=(E_cm_kinematics*(4000/34)); // convert from E_lab to E_cm and MeV to keV
	        hExcitation2->Fill(dEpadBsum);
	    } 
	    
	    
	    /*
	    if (xmedbest < 212 && xmedbest>0) { // ACTIVE TARGET SCATTERING!! For now omit pad 32 and higher
              //cout << endl << "ACTIVE TARGET SCATTERING" << endl << endl;
              //cout << xmedbest ;
	      // how we converted pad to Z so we can reverse it
	      //zmstpc=82.75+(i*padsize)+(padsize/2);
              // cleanup init
	      //float pad_scatter_tmp = (xmedbest - 82.75 - (padsize/2)) / padsize;
	      //int pad_scatter = pad_scatter_tmp;
	      //double dEpadBsum = 36.9; // 17.3 um mylar; 6.8 um kapton; 194 torr, 213 K 82.75+4 = 86.75 mm HeCO2
	      pad_scatter_tmp = (xmedbest - 82.75 - (padsize/2)) / padsize;
	      pad_scatter = pad_scatter_tmp;
	      // "randomizer" for missing pads
	      if (pad_scatter==18) pad_scatter-=(jentry%2);
	      if (pad_scatter==27) pad_scatter-=(jentry%4);
	      gate_30s_bragg = 0;
	      for (int i=1;i<=pad_scatter;i++){
	        if (i==17||i==24||i==25||i==26||i>31) continue;
		if ( (PadIsHitB[i][0]) && !(dEpadB[i][0]>dEpadB_gate_30s[0][i] && dEpadB[i][0]<dEpadB_gate_30s[1][i])) {
		  gate_30s_bragg++; 
		  //cout << pad_scatter << "\t" << i << "\t" << dEpadB[i][0] << endl;
		}
	      }
	      //cout << "EVENT: " << pad_scatter << "\t" << gate_30s_bragg << endl;
	      //if (gate_30s_bragg<=1){
	      if (gate_30s_bragg<1){
	      //if (1){
	      //if (pad_scatter<5) cout << dEpadB[1][0] << "\t" << dEpadB[2][0] << "\t" << dEpadB[3][0] << "\t" << dEpadB[4][0] << endl;
	      //cout << dEpadB[1][0] << endl;
	      //if ((pad_scatter-1)/4 >= gate_30s_bragg){
	      //if ((fSiE[0][0]< 910) || (fSiE[0][0]>980)){
	      //if ((gate_30s_bragg<1) && (fSiE[0][0]< 910) || (fSiE[0][0]>980)){
		//cout << "EVENT: " << pad_scatter << "\t" << gate_30s_bragg << endl;
		
		//stole dEpadBsum from here!
		
		dEpadBsum=(dEpadBsum*(4000/34)); // convert from E_lab to E_cm and MeV to keV
	        if (cutg->IsInside(dEpadBsum,pad_scatter)){
	          hExcitation->Fill(dEpadBsum);
	          hExPad->Fill(dEpadBsum,pad_scatter);
                  if (!SiIsHit[1] && !SiIsHit[2] && SiIsHit[0]) {
		    hSiPadEcal->Fill(1.,fSiE[0][0]);	
		    scatter_event = true;
		    hEx_Si1a->Fill(dEpadBsum,fSiE[0][0]);
		  }
	        }
                //cout << "\t" << pad_scatter << endl;
	        //
              }
            }
	    */
	    
	    //eventcounter++;
	  }
  



          if (!scatter_event && gate_30s && padsHitB[0]>5)  hSiPadEcal->Fill(2.,fSiE[0][0]);	
          if (!scatter_event && gate_30s && padsHitB[0]<=5)  hSiPadEcal->Fill(3.,fSiE[0][0]);	
	  //if (SsdOR && (SiIsHit[0]||SiIsHit[1]||SiIsHit[2]) && gate_30s && padsHitB[0] > 5){
	


	  for(UShort_t i=0;i<48;i++){ // loop over all pads (now by pad number)
                         //TargetX=fTargetX->Eval((452.+82.+(i*padsize)+2.1));
			 //hTargetXZ->Fill(TargetX,i);
	    for(UShort_t j=0;j<1;j++){ // loop for first track only
	      //for(UShort_t j=0;j<TracksB[i];j++){ // loop for each track
	      //if (i==4 && TracksB[4]>1) cout << "TpadB[4][" << j << "]: " << TpadB[4][j] << endl;
	      //if (i==6) cout << "TpadB[6][" << j << "]: " << TpadB[6][j] << endl;
	      //cout << "TracksB[" << i << "]: " << TracksB[i] << endl;
	     if (!SsdOR) hBraggB->Fill(i,dEpadBcalib[i][j]);
	     TargetX=TargetY=-999.;
	     if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && gate_30s && XpadB[i][j]>-999. && dEpadB[i][j]>5 ){
	     //if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && gate_29p){
                         //checking TargetX works as expected
			 //TargetX=PpacX[0]+(((452+83+(i*padsize)+(padsize/2))/PpacSepZ)*(PpacX[1]-PpacX[0]));
			 // testing retarded stuff: 23 Apr 2013 16:28:43 
			 //TargetX=fTargetX->Eval((452.+83.+(i*padsize)+(padsize/2)));
			 // correct one I think:
			 TargetX=fTargetX->Eval((452.+82.75+(i*padsize)+(padsize/2)));
                         TargetY=fTargetY->Eval((452.+82.75+(i*padsize)+padsize/2)); // PPACa center to window, inactive region, plus pad
			 //hPadXBgeo_30s->Fill(i,TargetX); // simple target projection by PPAC for reference
			 hPadXBgeo_30s->Fill(i,XpadB[i][j]-TargetX); // proper residual w/o shifting
			 //hPadXB_30s->Fill(i,XpadBraw[i][j]-TargetX);
			 hPadXB_30s->Fill(i,XpadB[i][j]); // calibrated track
			 hBraggBvX_30s[i]->Fill(TargetX,dEpadB[i][j]);
			 //hBraggBvX[i]->Fill(XpadB[i][j],dEpadB[i][j]);
			 
             }
	     TargetX=TargetY=-999.;
	     if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && gate_29p && XpadB[i][j]>-999. && dEpadB[i][j]>5 ){
			 TargetX=fTargetX->Eval((452.+82.75+(i*padsize)+(padsize/2)));
                         TargetY=fTargetY->Eval((452.+82.75+(i*padsize)+padsize/2)); // PPACa center to window, inactive region, plus pad
			 //hPadXBgeo_29p->Fill(i,TargetX); // 
			 hPadXBgeo_29p->Fill(i,XpadB[i][j]-TargetX); // proper residual
			 //hPadXB_29p->Fill(i,XpadBraw[i][j]-TargetX); // proper residual
			 hPadXB_29p->Fill(i,XpadB[i][j]); // calibrated track
			 //hBraggBvX[i]->Fill(XpadB[i][j],dEpadB[i][j]);
			 //hBraggBvX[i]->Fill(TargetX,dEpadB[i][j]);
			 hBraggBvX_29p[i]->Fill(TargetX,dEpadB[i][j]);
	     
	     }
	     TargetX=TargetY=-999.;
	     // THIS WAS TO GET THE DRIFT TIME ZERO OFFSET, WHICH IS NEAR 80 ns
	     //if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && (gate_30s||gate_29p)  && XpadB[i][j]>-999. && dEpadB[i][j]>50 && (PpacY[0] < 0.5) && (PpacY[0] > -1.5) && (PpacY[1] < -0.5) && (PpacY[1] > -1.5) ){
	     //if (PpacIsHit[0] && PpacIsHit[1] && (PpacY[1] > 2.5) && (PpacY[1] < 3.5) && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && (gate_30s||gate_29p)  && XpadB[i][j]>-999. && dEpadB[i][j]>50){
	     if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && (gate_30s||gate_29p)  && XpadB[i][j]>-999. && dEpadB[i][j]>50){
			 TargetX=fTargetX->Eval((452.+82.75+(i*padsize)+(padsize/2)));
                         TargetY=fTargetY->Eval((452.+82.75+(i*padsize)+padsize/2)); // PPACa center to window, inactive region, plus pad
			 hPadXBgeo->Fill(i,((XpadB[i][j])-TargetX)); // proper residual
			 hPadXB->Fill(i,(XpadB[i][j]));
			 //hPadXB->Fill(i,(XpadBraw[i][j])-TargetX);
			 hPadYB1->Fill(i,YpadB[i][j]);
			 hPadYBgeo1->Fill(i,TargetY);
			 hPadYBgeo2->Fill(i,YpadB[i][j]-TargetY);
	                 
			 /*
			 if (PpacY[0]<0.5 && PpacY[0]>-0.5 && PpacY[1]<1. && PpacY[1]>0.) {
			   hPadYB1->Fill(i,YpadB[i][j]);
			   hPadYBgeo1->Fill(i,TargetY);
			 }
	                 if (PpacY[0]<3. && PpacY[0]>2. && PpacY[1]<3.5 && PpacY[1]>2.5) {
			   hPadYB2->Fill(i,YpadB[i][j]);
			   hPadYBgeo2->Fill(i,TargetY);
			 }
	                 if (PpacY[0]<5. && PpacY[0]>4. && PpacY[0]<6. && PpacY[1]>5.) {
			   hPadYB3->Fill(i,YpadB[i][j]);
			   hPadYBgeo3->Fill(i,TargetY);
			 }
			 */
			 //hPadYBgeo->Fill(i,YpadB[i][j]-TargetY);
			 //hBraggBvX[i]->Fill(XpadB[i][j],dEpadB[i][j]);
	     
	     }
	     //hPadXB->Fill(i,XpadB[i][j]);
	     if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && gate_30s && XpadB[i][j]>-999. && dEpadB[i][j]>50){
	       hPadYB_ds->Fill(i,YpadB[i][j]);
	     }
	     if (PpacIsHit[0] && PpacIsHit[1] && SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.)) && gate_30s && XpadB[i][j]>-999. && dEpadB[i][j]>50){
	       hPadYB_ssd->Fill(i,YpadB[i][j]);
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
	       if (WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.))){
	         if (fSiE[0][0]>1000.) hBraggB_30s->Fill(i,dEpadB[i][j]);
	         //hBraggB_30s->Fill(i,dEpadBcalib[i][j]);
		 }
	         //hPadXB_30s->Fill(i,XpadB[i][j]);
	         if (!SsdOR) { //downscale beam condition
	           hBraggB_30s_ds->Fill(i,dEpadBcalib[i][j]);
	         } // end if : !SsdOR
	         // PUT PHYSICS HERE
	         if ( SsdOR) { // ssd-or condition
	           hBraggB_30s_ssd->Fill(i,dEpadBcalib[i][j]);
	           hPadYB_30s->Fill(i,YpadB[i][j]);
	           // Separate the Bragg peaks by the track number. 
	           hBraggB_30s_jch[j]->Fill(i,dEpadBcalib[i][j]);
	         } // end if : SsdOR
		 if (!SsdOR) { // d/s condition
	         } // end if : !SsdOR
	       } // end if dEpadB>10
             } // end if : gate_30s
	     if (flag_ppac && gate_29p) {
	       if (dEpadB[i][j]>10){
	     //if (PpacIsHit[0] && PpacIsHit[1] && !SsdOR  && WindowCut(fTargetX->Eval(452.),fTargetY->Eval(452.))){
	     //can check window cut
	         hBraggB_29p->Fill(i,dEpadBcalib[i][j]);
	         if (!SsdOR) { //downscale beam condition
	           hBraggB_29p_ds->Fill(i,dEpadBcalib[i][j]);
                 }	
	         if ( SsdOR) { // ssd-or condition
	           hBraggB_29p_ssd->Fill(i,dEpadBcalib[i][j]);
	         }
	       } // end if: dEpadB > 10
	     } // end if: gate_29p
             if (gate_17f){
	       if (dEpadB[i][j]>10){
	         hBraggB_17f->Fill(i,dEpadBcalib[i][j]);
                 if (SsdOR) {
	         }
	         else{ // downscale

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
		  if (SiIsHit[2] && pulsetrack3a<100 && pulsetrack3a>0 && padsHitL[0]>1 && gate_30s){
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
	        if (fSiE[2][4]>200 || fSiE[2][5]>200 || fSiE[2][6]>200 || fSiE[2][7]>200 || fSiE[2][8]>200) hBraggLvZ[i]->Fill(ZpadL[i][j],dEpadL[i][j]);
	        //if (fSiE[2][1]>200 || fSiE[2][2]>200 || fSiE[2][3]>200) hBraggLvZ[i]->Fill(ZpadL[i][j],dEpadL[i][j]);
		if (SiIsHit[2] && dEpadL[i][j]>1.){
	        //if (SiIsHit[2] && dEpadL[i][j]>1. && (gate_29p || gate_30s)){
	        //if (SiIsHit[2] && dEpadL[i][j]>1.&& fSiE[2][1]>200 && fSiEStripmax_ch[2]==1){
	          hPadYL->Fill(i,YpadL[i][j]);
	          hPadZL[0]->Fill(i,ZpadL[i][j]);
	          //if (SiIsHit[3] && ssd_mult==1){
	          if (SiIsHit[2] && pulsetrack3a<100 && pulsetrack3a>0 && padsHitL[0]>1 && gate_30s){
	            hPadZL[pulsetrack3a]->Fill(i,ZpadL[i][j]);
                  }   
                }
	        for(int k=1;k<9;k++){
		  if (SiIsHit[2] && dEpadL[i][j]>1.&& fSiE[2][k]>200 && fSiEStripmax_ch[2]==k) hPadZL3a[k]->Fill(i,ZpadL[i][j]);
		}
	      }
	    }	
	  }
	  if (SiIsHit[2] && padsHitL[0]>1 && gate_30s) pulsetrack3a++;
	  if (SiIsHit[2] && gate_30s) ssd_evts[2]++;
	 //center
	  for (int i=0;i<20;i++){
	    dEpadCTotal[i]=0; // 
	  }
	  
	  for(UShort_t i=0;i<8;i++){ // loop over all pads (now by pad number)
	    //for(UShort_t j=0;j<2;j++){ // first two tracks debugging pad 0 13 Mar 2013 22:49:18 
	    for(UShort_t j=0;j<1;j++){ // first track
	    //for(UShort_t j=0;j<TracksC[i];j++){ // loop for each track
	      if ( dEpadC[i][j]>1.){// 
	      //if ( dEpadC[i][j]>1. && proton_beam ){// proton beam (I assume?)
	      //if ( dEpadC[i][j]>1. && ((fSiE[1][1]>110.) || ( fSiE[1][2]>110.) || ( fSiE[1][3]>120.))){// proton beam (I assume?)
	      //if ( dEpadC[i][j]>1. && ((fSiE[1][1]<1100. && fSiE[1][1]>900.) || (fSiE[1][2]<1100. && fSiE[1][2]>900.) || (fSiE[1][3]<1100. && fSiE[1][3]>900.))){// proton beam (I assume?)
              //1f (dEpadC[i][j]>1. && (SiIsHit[0] == true || SiIsHit[1] == true)) { // 
	        dEpadCTotal[j] = dEpadCTotal[j] + dEpadC[i][j];
	        
		hBraggCvX[i]->Fill(XpadC[i][j],dEpadC[i][j]);
			 //hBraggBvX[i]->Fill(TargetX,dEpadB[i][j]);
		if ( padsHitC[0]>1  ) hBraggC->Fill(i,dEpadC[i][j]);
	        //if (gate_17f && SiIsHit[1] && padsHitC[0]>1  ) hBraggC->Fill(i,dEpadC[i][j]);
	        if (padsHitC[0]>1){
	        //if (SiIsHit[0] && ssd_mult==1){
	          //if (1){
		  //if (SiIsHit[0]) {
		  if (fSiE[0][3]>110) {
		    if (XpadC[i][j]>-40 && XpadC[i][j]<78)
		    hPadXC1asum3->Fill(i,XpadC[i][j]);
		    hBraggC3->Fill(i,dEpadC[i][j]);
		    if (i==6){
		    hg_peak_right[3]+=hg_peak_holder[0];
		    hg_peak_left[3]+=hg_peak_holder[1];
		    hg_peak_right_counter[3]++;
		    hg_peak_left_counter[3]++;
		    //cout << dEpadC[i][j] << endl;
		    }
		  }
		  if (fSiE[0][4]>120) {
		    if (XpadC[i][j]>-40 && XpadC[i][j]<78)
		    hPadXC1asum4->Fill(i,XpadC[i][j]);
		    hBraggC4->Fill(i,dEpadC[i][j]);
		    if (i==6){
		    hg_peak_right[4]+=hg_peak_holder[0];
		    hg_peak_left[4]+=hg_peak_holder[1];
		    hg_peak_right_counter[4]++;
		    hg_peak_left_counter[4]++;
		    }
		  }
		  if (fSiE[0][5]>135) {
		    if (XpadC[i][j]>-40 && XpadC[i][j]<78)
		    hPadXC1asum5->Fill(i,XpadC[i][j]);
		    hBraggC5->Fill(i,dEpadC[i][j]);
		    if (i==6){
		    hg_peak_right[5]+=hg_peak_holder[0];
		    hg_peak_left[5]+=hg_peak_holder[1];
		    hg_peak_right_counter[5]++;
		    hg_peak_left_counter[5]++;
		    //cout << hg_peak_holder[0] << "\t" << hg_peak_holder[1] << endl;
		    }
		  }
		  
	          //if (!SiIsHit[0]) {
	          if (SiIsHit[1] && (gate_30s||gate_29p)) {
	          //if (fSiE[1][4]>200 || fSiE[1][5]>200 || fSiE[1][6]>200 || fSiE[1][7]>200 || fSiE[1][3]>200 || fSiE[1][2]>200) {
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
		//if (( SiIsHit[0] || SiIsHit[1] )&& padsHitC[0]>1){ //gate_30s &&
		if (( SiIsHit[0] || SiIsHit[1] )&& padsHitC[0]>1 && (gate_30s||gate_29p)){ //gate_30s &&
		  hPadXC->Fill(i,XpadC[i][j]);
		}
		//if (SiIsHit[0] && padsHitC[0]>1 ){ //gate_30s &&
		if (SiIsHit[0] && padsHitC[0]>1 && (gate_30s || gate_29p)){ //gate_30s &&
		//if (SiIsHit[0] && pulsetrack1a<100 && padsHitC[0]>1 && //gate_30s &&
	        // fSiE[0][8]>0 ) {
		//(fSiE[0][1]>0 || fSiE[0][2]>0 || fSiE[0][3]>0 || fSiE[0][4]>0 || 
		 //fSiE[0][5]>0 || fSiE[0][6]>0 || fSiE[0][7]>0 || fSiE[0][8]>0 )){
		//if (SiIsHit[0] && pulsetrack2<100 && pulsetrack2>0 && padsHitC[0]>1){
		  hPadXC1a[0]->Fill(i,XpadC[i][j]);
		  hPadYC1a->Fill(i,YpadC[i][j]);
		}
		if (SiIsHit[1] && pulsetrack2a<100 && padsHitC[0]>1 && gate_17f ){//&&
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
          if (padsHitC[0]>1 && fSiE[0][3]>200  && fSiEStripmax_ch[0]==3){
	    // some while or for loop which removes points? 19 Feb 2013 13:43:52 
	    int points = padsHitC[0];
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
	  }


	  if (SiIsHit[0] &&  padsHitC[0]>1 && //gate_30s && 
	  fSiE[0][8]>0 ) {
         //  (fSiE[0][1]>0 || fSiE[0][2]>0 || fSiE[0][3]>0 || fSiE[0][4]>0 || 
         //  fSiE[0][5]>0 || fSiE[0][6]>0 || fSiE[0][7]>0 || fSiE[0][8]>0 )){
	    pulsetrack1a++;
	  }
	  if (SiIsHit[1] &&  padsHitC[0]>1 && gate_17f ){//&& 
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
	   if (padsHitB[0]>0 && PpacIsHit[0] && PpacIsHit[1]) hMultB->Fill(padsHitB[0]);
	   if (padsHitL[0]>0 &&  SiIsHit[2] && dEpadLTotal[0]>50.) hMultL_3a->Fill(padsHitL[0]);
	   //if (padsHitL[0]>0 && PpacIsHit[0] && PpacIsHit[1] && SiIsHit[2] && dEpadLTotal[0]>50.) hMultL_3a->Fill(padsHitL[0]);
	   if (flag_ssd && dEpadRTotal[0]>1 && SiIsHit[0]){ 
             hBraggR_E[0]->Fill(dEpadRTotal[0],fSiE[0][0]);
	     for (int i=1;i<9;i++){
               if (padsHitR[0]==i) {
	         hBraggR_E[i]->Fill(dEpadRTotal[0],fSiE[0][0]);
	         }
	     }
	   }
	   //if (flag_ssd && dEpadCTotal[0]>1 && fSiE[0][4]>120){ 
	   //if (flag_ssd && dEpadCTotal[0]>1 && fSiEmax[0]>50){ 
	   //if (flag_ssd && dEpadCTotal[0]>0.5 && SiIsHit[0]){ 
	   if (SsdOR && SiIsHit[0] && gate_30s && (padsHitC[0]>2) && (padsHitB[0]>5)){
             //hBraggC_E[0]->Fill(dEpadCTotal[0],fSiE[0][4]);
             //hBraggC_E[0]->Fill(dEpadCTotal[0],fSiEmax[0]); // alpha-calib-post
             hBraggC_E[0]->Fill(dEpadCTotal[0],fSiE[0][0]);  // production
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
	       if (padsHitC[0]==i) {
	         // CNS Interview Feb 2013 dirty way
		 hBraggC_E[0]->Fill((dEpadCTotal[0]*(6/i)),fSiE[0][0]);
	       }
	     }
	   }
	   if (flag_ssd && dEpadLTotal[0]>1 && SiIsHit[2]){ 
	   //if (flag_ssd && dEpadCTotal[0]>1 && SiIsHit[0]){ 
             hBraggL_E[0]->Fill(dEpadLTotal[0],fSiE[2][0]);
	     for (int i=1;i<9;i++){
               if (padsHitL[0]==i) {
	         hBraggL_E[i]->Fill(dEpadLTotal[0],fSiE[2][0]);
	         }
	     }
	   }
	   if (gate_30s && SsdOR){
	     if (padsHitB[0]>0 ) hMultB_30s->Fill(padsHitB[0]);
	     
	   }
	   if (1){
	     if (dEpadCTotal[0]>10.){
	     //if (!(dEpadC[5][1]>5) && !(padsHitC[0]==1)){
	       //for (int i=0;i<8;i++){
	       
	         //if (padsHitC[0]>0 &&  SiIsHit[0] && gate_30s) hMultC_30s_1a->Fill(padsHitC[0]);
	         if (padsHitC[0]>0 && SiIsHit[0] && fSiE[0][4]>190) hMultC_30s_1a->Fill(padsHitC[0]);
	         if (padsHitC[0]>0 && SiIsHit[1] && gate_30s) hMultC_30s_2a->Fill(padsHitC[0]);
	       //}
	     }
	     if (dEpadLTotal[0]>50.){
	       if (padsHitL[0]>0 && SiIsHit[2] && gate_30s) hMultL_30s_3a->Fill(padsHitL[0]);
	     }
	   }
	 }

// end multiplicity

        } // end if: flag_detail

      } // end if: flag_tpc
          //if (!scatter_event && !gate_30s)  hSiPadEcal->Fill(4.,fSiE[0][0]);	

  // arr->Clear(); 
   fTargetX->ReleaseParameter(0);
   fTargetX->ReleaseParameter(1);
   fTargetX->ReleaseParameter(2);
   fTargetY->ReleaseParameter(0);
   fTargetY->ReleaseParameter(1);
   fTargetY->ReleaseParameter(2);
   
   } // close loop on jentry
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
