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
  
 In order to publish results obtained from crabat or any
 future derivative software, the authors must make reference
 to the paper describing it (in preparation): D. Kahl et al.
 ``Algorithms for Active Target Data Analysis'' 2013.  This is an 
 additional term to the license as stipulated in section 7b.
 Publication includes, but is not limited to: books, peer-
 reviewed journals, conference proceedings, annual reports, etc.
 
*/
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#define MAX_FILENAME 200   //filename ..up to 200 letters

using namespace std;


class Calibration {
private:

  char filename[MAX_FILENAME];
  int nch;
  int nmax;

  Double_t *offset;
  Double_t *gain;

  ifstream fstr;


public:
  Calibration();
  Calibration(const char* calibfile,int nmax_in);

  virtual ~Calibration();

  void Init(const char* calibfile,int nmax_in); 
  void AddEntry(Double_t offset_in, Double_t gain_in);

  void SetFileName(const char *str) {strncpy(filename,str,MAX_FILENAME);}
  void ScanFile();
  void ShowEntries();

  Double_t Offset(int n) {if (n<nch) return offset[n]; else return 0;}
  Double_t Gain(int n) {if (n<nch) return gain[n]; else return 1;}
  //Get calibration parameter, (x-offset)*gain 
  Double_t Calib(int n,Double_t x) {return ((x-Offset(n))*Gain(n));}
  Double_t CalibMeV(int n,Double_t x) {return ((x-Offset(n))*Gain(n));}
  //Double_t CalibMeV(int n,Double_t x) {return ((x-Offset(n))*Gain(n)/1000);}

};

