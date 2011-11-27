/**
@example example_enewzsub.cxx 
*/
{
  // Testing KaliVeda 
  KVNucleus *beam = new KVNucleus;
  beam->Set("30S");
  beam->SetEnergy(64);
  
  int z1 = beam->GetZ() ; // Beam nuclear charge Z in integer units of the elementary charge e
  float m1 = beam->GetA() ; // Beam mass M in AMU
  float e = beam->GetEnergyPerNucleon(); // Beam energy E in MeV/u
  
  // testing enewzsub
  //
  //GAS TARGET PROTOTYPE - Uses cryogenic He gas like the CRIB production target
  char matter1[34] = "heco2"; // Matter the beam passes through, as defined in SNKE_MATTER.INC -- Case-insensitive.
  int unit_pressure = 1 ; // Define the units of pressure: 1=Torr ; 2=mbar ; 3=atm
  float pressure = 194 ; // Define the value of the pressure.
  float temperature = 300 ; // Define temperature of the target in Kelvin
  int unit_thick = 1; // Define the units of thickness: 1=mm ; 2=mg/cm^2
  float thick1 = 272 ; // Define the value of the thickness.
  
  //SOLID TARGET PROTOTYPE - A single havar foil
  char matter1[34] = "havar"; // Matter the beam passes through, as defined in SNKE_MATTER.INC -- case insensitive
  int unit_pressure = 0 ; // Not used for solid target; set as 0
  float pressure = 0 ; // Not used for solid target; set as 0
  float temperature = 0 ; // Not used for solid target; set as 0
  int unit_thick = 1; // Define the units of thickness: 1=mm ; 2=mg/cm^2
  float thick1 = 0.0025 ; // Define the value of the thickness.  With this and the above line, it is 2.5 um
  
  float aft_ene; // Beam energy E in MeV/u after energy loss in matter
  
  //The matter name must be 32 characters plus the termination charcater \0 or Fortran will not read it correctly
  //The for loop takes the input matter name and appends the approrpiate number of spaces.
  for(int i=0;i<33;i++)
    if(matter1[i]=='\0')
  {
  	matter1[i]=' ';
  	matter1[i+1]='\0';
  }
  
  enewzsub_(&z1, &m1, &e, matter1, &unit_pressure, &pressure, &temperature, &unit_thick, &thick1, &aft_ene);
  cout << "The energy of the beam after the active target region is: " << aft_ene << endl;
  delete beam;
}
