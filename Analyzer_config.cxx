/*
 Copyright: 
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
@file Analyzer_config.cxx
@author daid
@version 1.0
@brief Analyzer configurations which vary experiment-by-experiment
@details Does detailed PPAC calibration.\n
Defines rf frequency.
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
 * @brief Sets configurations which vary experiment-by-experiment
 * @details If we share analysis code across experiments,\n
 * it is impractical to get new updates which over-ride experiment\n
 * specific data.  But if we put these details all only in one file\n
 * it is much easier for different experiments to use the crabat code.
 */
void Analyzer::SetRF()
{
   // Define RF
   // Note in class structure we cannot use static const style
   // as it is forbidden by ISO C++ for non-integers
   #define        rf 59.70; // real RF value
}

UShort_t Analyzer::SetPpacSepZ()
{
  return(156) ; // mm 
  //#define PpacSepZ  156; // mm 
  //static const UShort_t PpacSepZ = 156; // mm 
}

// change this to pass by address
Bool_t Analyzer::WindowCut(Double_t x, Double_t y)
{
// modified from http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=1145
  double radius = 20; // real value 20 mm 
  const Int_t n = 30; 
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
  return (false);
}

UShort_t Analyzer::PpacXYCalib(Bool_t flag_detail, Bool_t flag_ppac){ 
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
  PpacOffsetLine[1][0]=-4.3; //18 Jul 2011 16:58:03 
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
    if (i==0) PpacY[i]=dxy[1]; // PPACaY is normal
    if (i==1) PpacY[i]=(dxy[1]*-1.0); // PPACbY cables were inverted
  } //end PPAC calibration for loop
  return(EXIT_SUCCESS);
}



