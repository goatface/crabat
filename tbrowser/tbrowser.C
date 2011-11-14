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

// Quick way to open TBrowser from the command line, and nice colors

#include "TROOT.h"
#include "TRint.h"
#include "TBrowser.h"
#include "TStyle.h"
#include "TColor.h"
#include <iostream>
using namespace std;

int main(int argc, char **argv)
{
    // Create interactive interface
    TRint *theApp = new TRint("ROOT CInt", &argc, argv, NULL, 0);
    new TBrowser();
    
    // set nice plot colors
    // from http://ultrahigh.org/2007/08/20/making-pretty-root-color-palettes
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont); 
    
    cout << "" << endl;
    // Run interactive interface
    theApp->Run();
    return(0);
}
