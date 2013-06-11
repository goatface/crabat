/*
 Copyright: 
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
@file Analyzer_config.cxx
@author daid
@version 1.0
 * @brief Sets configurations which vary experiment-by-experiment
 * @details If we share analysis code across experiments,\n
 * it is impractical to get new updates which over-ride experiment\n
 * specific data.  But if we put these details all only in one file\n
 * it is much easier for different experiments to use the crabat code.\n
 * Handles MSTPC pad mapping.\n
 * Handles MSTPC pad gain.\n
 * Does detailed PPAC calibration.\n
 * Defines rf frequency.
@date 05 Dec 2011 16:03:48  
*/
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

#define Analyzer_config_cxx


/**
 * @brief Sets the MSTPC pad mapping
 * @details Creates a 2D vector for left and right channel mapping.  
 * TPC mapping is non-trivial because there is not an invertible function,
 * i.e. each pad has two sides.
 * Thus we create two sets of 144 channels, half of which are not used
 * on either side.  When the vector is queried with a raw channel map,
 * it returns the pad number.
*/
vector<vector<Short_t> > TDetector::SetTpcMap()
{
  
  //basic vector usage from http://www.cplusplus.com/forum/articles/7459/
  static const Short_t HEIGHT=2;
  static const Short_t WIDTH=144;
  static const Short_t tpc_ch[HEIGHT][WIDTH]={

// low gain mapping daid 26 Oct 2011 03:36:43 22Mg, now 30S 
// rough high gain mapping daid 11 Sep 2012 00:09:38 (SET 0,15 might be inverted?)
// should be correct now I think!
{ 
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,00:OCCUPIED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,01:OCCUPIED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,02:OCCUPIED
  7,6,5,4,3,2,1,0,          // SET 0,03:Beam right upstream
  15,14,13,12,11,10,9,8,    // SET 0,04:Beam right upstream
  23,22,21,20,19,18,17,16,  // SET 0,05:Beam right upstream
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,06:OCCUPIED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,07:OCCUPIED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,08:OCCUPIED
  31,30,29,28,27,26,25,24,  // SET 0,09:Beam right downstream
  39,38,37,36,35,34,33,32,  // SET 0,10:Beam right downstream
  47,46,45,44,43,42,41,40,  // SET 0,11:Beam right downstream
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,12:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,13:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,14:DISABLED
  0,1,2,3,4,5,6,7,           // SET 0,15:high gain right center
  //7,6,5,4,3,2,1,0,          // SET 1,15:high gain right center hacking
  7,6,5,4,3,2,1,0,          // SET 1,15:high gain right downstream 
  0,1,2,3,4,5,6,7,           // SET 0,15:high gain right upstream
  //-1,-1,-1,-1,-1,-1,-1,-1,  // SET 0,17:DISABLED
   },

 { 
  16,17,18,19,20,21,22,23,  // SET 1,00:Beam left upstream
  8,9,10,11,12,13,14,15,    // SET 1,01:Beam left upstream
  0,1,2,3,4,5,6,7,           // SET 1,02:Beam left upstream
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,03:OCCUPIED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,04:OCCUPIED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,05:OCCUPIED 
  40,41,42,43,44,45,46,47,  // SET 1,06:Beam left downstream
  32,33,34,35,36,37,38,39,  // SET 1,07:Beam left downstream
  24,25,26,27,28,29,30,31,  // SET 1,08:Beam left downstream
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,09:OCCUPIED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,10:OCCUPIED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,11:OCCUPIED
  //0,1,2,3,4,5,6,7,           // SET 1,13:high gain left down
  //6,5,-1,4,3,1,1,1,          // SET 1,12:high gain left center		//mismatching
  7,6,5,4,3,2,1,0,          // SET 1,12:high gain left center
  ///0,1,2,3,7,6,5,4,          // SET 1,12:high gain left center
  7,6,5,4,3,2,1,0,          // SET 1,14:high gain left up
  0,1,2,3,4,5,6,7,           // SET 1,13:high gain left down
  ////7,6,5,4,0,1,2,3,           // SET 1,14:high gain left up
  //7,6,5,4,3,2,1,0,          // SET 1,12:high gain left center
  //0,-1,2,3,4,5,6,7,           // SET 1,14:high gain left up		//mismatching
  //-1,7,6,5,4,2,2,0,          // SET 1,12:high gain left center		//mismatching
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,15:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,16:DISABLED
  -1,-1,-1,-1,-1,-1,-1,-1,  // SET 1,17:DISABLED
  }
};

/*
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
    };*/
  
  vector<vector<Short_t> > array2D;
  array2D.resize(HEIGHT);
  for (int i = 0; i < HEIGHT; ++i)
    array2D[i].resize(WIDTH);
  for (int i = 0; i < HEIGHT; i++)
    for (int j = 0; j < WIDTH; j++)
      array2D[i][j]=tpc_ch[i][j];
  return array2D;
}

/**
 * @brief Sets the MSTPC Beam pad gain
 * @details Creates a vector for gain mapping.
 * As a vector, we can create as many different gain 
 * vectors as required, since vector sizes are dynamic.
 * Here for 30S I create a 2D vector, since the gain was
 * changed after some early production runs.
 */
vector<vector<Double_t> > TDetector::SetTpcBgain()
{
  //basic vector use suggested at http://www.cplusplus.com/forum/articles/7459/
  //internal calibration for beam pad
  static const Short_t HEIGHT=2;
  static const Short_t WIDTH=48;
  static const Double_t padBgain[HEIGHT][WIDTH]={
/*  { // runs 1001-1005
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
  }};*/
  { // runs 1001-1005
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,// 0-7
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 8-15
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 8-15
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 8-15
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 8-15
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  }, // 8-15
  
  { //other runs
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 0-7
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 8-15
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 0-7
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 8-15
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  , // 0-7
    1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.  ,1.   // 8-15
  }};
  vector<vector<Double_t> > array2D;
  array2D.resize(HEIGHT);
  for (int i = 0; i < HEIGHT; ++i)
    array2D[i].resize(WIDTH);
  for (int i = 0; i < HEIGHT; i++)
    for (int j = 0; j < WIDTH; j++)
      array2D[i][j]=padBgain[i][j];
  return array2D;
}

void inline Analyzer::SetRF()
{
   // Define RF
   // Note in class structure we cannot use static const style
   // as it is forbidden by ISO C++ for non-integers
   #define        rf 59.70; // real RF value
}

UShort_t inline Analyzer::SetPpacSepZ()
{
  return(156) ; // mm 
  //#define PpacSepZ  156; // mm 
  //static const UShort_t PpacSepZ = 156; // mm 
}

// change this to pass by address
Bool_t inline Analyzer::WindowCut(Double_t x, Double_t y)
{
// modified from http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=1145
  double radius = 10; // real value 20 mm...practically speaking 22 is the cut off for 29P bragg appearing 
  const Int_t n = 300; 
  Double_t winx[n+1],winy[n+1]; 
  Double_t rcut = 1.00*radius; 
  Double_t dphi = TMath::TwoPi()/n; 
  for (Int_t i=0;i<n;i++) { 
    winx[i] = rcut*TMath::Cos(i*dphi); 
    winy[i] = rcut*TMath::Sin(i*dphi); 
  } 
  winx[n] = winx[0]; winy[n] = winy[0]; 
  TCutG *window = new TCutG("window",n+1,winx,winy); 
  if (window->IsInside(x,y)) {
    delete window;
    return(true);
  }
  delete window;
  return (false);
}

UShort_t inline Analyzer::PpacXYCalib(Bool_t flag_detail, Bool_t flag_ppac){ 
  if (!flag_detail || !flag_ppac) return(EXIT_FAILURE);
  // PPAC calibration that is not handled in ppac_calib.dat
  Float_t dxy[2]; // 1 is X, 2 is Y; tmp variable for calibration loops
  Float_t PpacOffset[2][2];
  Float_t PpacOffsetLine[2][2];
  Float_t PpacPositionGain[2][2];
  Float_t PpacOffsetGeometry[2][2];
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
  //PpacOffsetLine[1][0]=-1.4; //05 Sep 2012 18:07:28 
  PpacOffsetLine[1][0]=-4.3; //18 Jul 2011 16:58:03 checking again 17 Apr 2013 15:57:19, and see Vol 4 p. 153...it has to be of such a value or projections wrong 
  PpacOffsetLine[1][1]=-1.0; //18 Jul 2011 16:58:04 
  PpacPositionGain[1][0]=0.6205; //18 Jul 2011 16:55:04 
  PpacPositionGain[1][1]=0.6125; //18 Jul 2011 16:55:06 
  PpacOffsetGeometry[1][0]=0.0;
  PpacOffsetGeometry[1][1]=0.0; // checked 08 Jul 2011 16:00:22 
  //PpacOffsetGeometry[1][1]=-3.5;
  // PPAC calibration for loop
  for (UShort_t i=0;i<2;i++){
    //average time
    tof[i] = 0.25*(fPPACcal[i][0]+fPPACcal[i][1]+fPPACcal[i][2]+fPPACcal[i][3]);
    // J1 results are consistent with calculations
    // But in AT data there is some offset
    // e.g. easily may have slightly different cable lengths
    //if (i==1) tof[i]-=2.;
    // we need to compare with SSDT data to see which direction
    // i.e. tof0+2 vs tof1-2
    dxy[0]=(fPPACcal[i][0]-fPPACcal[i][1]);
    //if (i==1) dxy[0]=(fPPACcal[i][1]-fPPACcal[i][0]); // inverted PPACbX 03 Sep 2012 20:23:40 goatface
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
    PpacY[i]=dxy[1]; 
    // Previously I thought the PPACbY cable was inverted...must be false from projection data! 28 Aug 2012 19:48:04 
    // PPACaY is normal
    //if (i==1) PpacY[i]=(dxy[1]*-1.0); // PPACbY cables were inverted
  } //end PPAC calibration for loop
  return(EXIT_SUCCESS);
}



