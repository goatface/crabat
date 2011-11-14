/*
 Copyright: 
 2010 daid KAHL, YAMAGUCHI hidetoshi 
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
//#include <iostream>
//#include <stdio.h>
#include <Rtypes.h>
#include "Calibration.h"

Calibration::Calibration() {
  Init("calibration.dat",200);
}

void Calibration::Init(const char* calibfile, int nmax_in) {//nmax...number of channels
  nch=0;
  nmax=nmax_in; 
  offset=new Double_t [nmax]; 
  gain=new Double_t [nmax]; 
  SetFileName(calibfile);

}

Calibration::Calibration(const char* mapfile,int nmax_in) {
  Init(mapfile,nmax_in);
}

Calibration::~Calibration() {

  delete[] offset;
  delete[] gain;
}

void Calibration::AddEntry(Double_t offset_in, Double_t gain_in) {
  if (nch<=nmax) {
    offset[nch]=offset_in;
    gain[nch]=gain_in;
    nch++;
  } else {
    cerr << "Error: Too many calibration parameters! " << endl;
    cerr << "Offset=" << offset_in;
    cerr << " Gain=" << gain_in << endl; 
  }
  
}

void Calibration::ScanFile() {

  //read calibration params from a file.
  //Format:
  //# Lines starting form "#" are comments.
  //# One channel/one line 
  //# offset gain  
  //20.20   0.0152
  //30.12   0.0145
  //...

  float g,o;
  //  Double_t g,o;

  //static allocation[200] -> better to modify
  char line[200];

  fstr.open(filename);
  
  if (!fstr.is_open()) {
    cerr << "Could not open file: " << filename << endl;
    return;
  }

  while (!fstr.eof()) {
    
    fstr.getline(line,200);

    if (line[0]=='#')  continue;

    //Line-> Name/Module/Ch/Id
    int nmatch=sscanf(line,"%g %g",&o,&g);//want to avoid using sscanf, if possible, 
    
    if (nmatch<2) {
      continue;
    }
    //ch
    AddEntry(o,g);
  }
  if (fstr.is_open()) {fstr.close();}
}

void Calibration::ShowEntries() {
  cout << "Calibration file: " << filename << "\n";
  cout << "[Ch: Offset, Gain] " << "\n";
  cout << "-------------------------------" << "\n";

  for (int i=0;i<nch;i++) {
    cout << i << ": " << offset[i] <<", "<< gain[i] << "\n";
  }
  cout << "-------------------------------" << "\n";
  cout << flush;
}

