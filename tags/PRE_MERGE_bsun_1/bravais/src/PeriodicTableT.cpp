// DEVELOPMENT
/* $Id: PeriodicTableT.cpp,v 1.7 2002-11-14 01:47:33 saubry Exp $ */

#include "PeriodicTableT.h"

#include "dArrayT.h"
#include "dArray2DT.h"

#include "ExceptionCodes.h"

#include "StringT.h"
#include "ArrayT.h"
#include "PerTabEntryT.h"

using namespace Tahoe;

PeriodicTableT::PeriodicTableT():PT(110) {
	Initialize();
}

void PeriodicTableT::Initialize() 
{
  PT[1].SetName("hydrogen");
  PT[1].SetSymbol("H");
  PT[1].SetMass(1.008);
  PT[1].SetLattParams(3.77,6.16);
  PT[1].SetLattType("hcp");
  
  PT[2].SetName("helium");
  PT[2].SetSymbol("He");
  PT[2].SetMass(4.003);
  PT[2].SetLattParams(3.53,4.24);
  PT[2].SetLattType("hcp");
  
  PT[5].SetName("boron");
  PT[5].SetSymbol("B");
  PT[5].SetMass(10.81);
  PT[5].SetLattParams(8.74,5.06);
  PT[5].SetLattType("tetragonal");
  
  PT[6].SetName("carbon");
  PT[6].SetSymbol("C");
  PT[6].SetMass(12.01);
  PT[6].SetLattParams(3.57);
  PT[6].SetLattType("diamond");
  
  PT[7].SetName("nitrogen");
  PT[7].SetSymbol("N");
  PT[7].SetMass(14.01);
  PT[7].SetLattParams(5.64);
  PT[7].SetLattType("sc");
  
  PT[8].SetName("oxygen");
  PT[8].SetSymbol("O");
  PT[8].SetMass(16.00);
  PT[8].SetLattParams(6.83);
  PT[8].SetLattType("sc");
  
  PT[9].SetName("fluorine");
  PT[9].SetSymbol("Fl");
  PT[9].SetMass(19.00);
  PT[9].SetLattParams(6.67);
  PT[9].SetLattType("sc");
  
  PT[10].SetName("neon");
  PT[10].SetSymbol("Ne");
  PT[10].SetMass(20.18);
  PT[10].SetLattParams(4.45);
  PT[10].SetLattType("fcc");
  
  PT[11].SetName("sodium");
  PT[11].SetSymbol("Na");
  PT[11].SetMass(22.99);
  PT[11].SetLattParams(3.77,6.15);
  PT[11].SetLattType("hcp");
  
  PT[12].SetName("magnesium");
  PT[12].SetSymbol("Mg");
  PT[12].SetMass(24.31);
  PT[12].SetLattParams(3.21,5.21);
  PT[12].SetLattType("hcp");
  
  PT[13].SetName("aluminum");
  PT[13].SetSymbol("Al");
  PT[13].SetMass(26.98);
  PT[13].SetLattParams(4.05);
  PT[13].SetLattType("fcc");
  
  PT[14].SetName("silicon");
  PT[14].SetSymbol("Si");
  PT[14].SetMass(28.09);
  PT[14].SetLattParams(5.43);
  PT[14].SetLattType("diamond");
  
  PT[15].SetName("phosphorus");
  PT[15].SetSymbol("P");
  PT[15].SetMass(30.97);
  PT[15].SetLattParams(18.51);
  PT[15].SetLattType("sc");
  
  PT[16].SetName("sulfur");
  PT[16].SetSymbol("S");
  PT[16].SetMass(32.07);
  PT[16].SetLattParams(10.46,12.87,24.49);
  PT[16].SetLattType("orthorhombic");
  
  PT[17].SetName("chlorine");
  PT[17].SetSymbol("Cl");
  PT[17].SetMass(35.45);
  PT[17].SetLattParams(6.24,4.48,8.26);
  PT[17].SetLattType("orthorhombic");
  
  PT[18].SetName("argon");
  PT[18].SetSymbol("Ar");
  PT[18].SetMass(39.95);
  PT[18].SetLattParams(5.31);
  PT[18].SetLattType("fcc");
  
  PT[19].SetName("potassium");
  PT[19].SetSymbol("K");
  PT[19].SetMass(39.10);
  PT[19].SetLattParams(5.33);
  PT[19].SetLattType("bcc");
  
  PT[20].SetName("calcium");
  PT[20].SetSymbol("Ca");
  PT[20].SetMass(40.08);
  PT[20].SetLattParams(5.59);
  PT[20].SetLattType("fcc");
  
  PT[22].SetName("titanium");
  PT[22].SetSymbol("Ti");
  PT[22].SetMass(47.88);
  PT[22].SetLattParams(2.95,4.68);
  PT[22].SetLattType("hcp");
  
  PT[23].SetName("vandium");
  PT[23].SetSymbol("V");
  PT[23].SetMass(50.94);
  PT[23].SetLattParams(3.02);
  PT[23].SetLattType("bcc");
  
  PT[24].SetName("chromium");
  PT[24].SetSymbol("Cr");
  PT[24].SetMass(52.00);
  PT[24].SetLattParams(2.88);
  PT[24].SetLattType("bcc");
  
  PT[25].SetName("manganese");
  PT[25].SetSymbol("Mn");
  PT[25].SetMass(54.94);
  PT[25].SetLattParams(8.91);
  PT[25].SetLattType("bcc");
  
  PT[26].SetName("iron");
  PT[26].SetSymbol("Fe");
  PT[26].SetMass(55.85);
  PT[26].SetLattParams(2.87);
  PT[26].SetLattType("bcc");
  
  PT[27].SetName("cobalt");
  PT[27].SetSymbol("Co");
  PT[27].SetMass(58.93);
  PT[27].SetLattParams(3.54);
  PT[27].SetLattType("fcc");
  
  PT[28].SetName("nickel");
  PT[28].SetSymbol("Ni");
  PT[28].SetMass(58.69);
  PT[28].SetLattParams(3.52);
  PT[28].SetLattType("fcc");

  PT[29].SetName("copper");
  PT[29].SetSymbol("Cu");
  PT[29].SetMass(63.55);
  PT[29].SetLattParams(3.615);
  PT[29].SetLattType("fcc");

  PT[30].SetName("zinc");
  PT[30].SetSymbol("Zn");
  PT[30].SetMass(65.39);
  PT[30].SetLattParams(2.66,4.95);
  PT[30].SetLattType("hcp");

  PT[31].SetName("gallium");
  PT[31].SetSymbol("Ga");
  PT[31].SetMass(69.72);
  PT[31].SetLattParams(4.52,7.66,4.53);
  PT[31].SetLattType("orthorhombic");
  
  PT[32].SetName("germanium");
  PT[32].SetSymbol("Ge");
  PT[32].SetMass(72.61);
  PT[32].SetLattParams(5.66);
  PT[32].SetLattType("diamond");
  
  PT[33].SetName("arsenic");
  PT[33].SetSymbol("As");
  PT[33].SetMass(74.92);
  PT[33].SetLattParams(4.13,54.1667);
  PT[33].SetLattType("rhombohedral");
  
  PT[36].SetName("krypton");
  PT[36].SetSymbol("Kr");
  PT[36].SetMass(83.80);
  PT[36].SetLattParams(5.72);
  PT[36].SetLattType("fcc");

  PT[40].SetName("zirconium");
  PT[40].SetSymbol("Zr");
  PT[40].SetMass(91.22);
  PT[40].SetLattParams(3.23,5.15);
  PT[40].SetLattType("hcp");
  
  PT[41].SetName("niobium");
  PT[41].SetSymbol("Nb");
  PT[41].SetMass(92.91);
  PT[41].SetLattParams(3.30);
  PT[41].SetLattType("bcc");
  
  PT[42].SetName("molybdenum");
  PT[42].SetSymbol("Mo");
  PT[42].SetMass(95.94);
  PT[42].SetLattParams(3.15);
  PT[42].SetLattType("bcc");
  
  PT[44].SetName("ruthenium");
  PT[44].SetSymbol("Ru");
  PT[44].SetMass(101.07);
  PT[44].SetLattParams(2.71,4.28);
  PT[44].SetLattType("hcp");
  
  PT[45].SetName("rhodium");
  PT[45].SetSymbol("Rh");
  PT[45].SetMass(102.91);
  PT[45].SetLattParams(3.80);
  PT[45].SetLattType("fcc");
  
  PT[46].SetName("palladium");
  PT[46].SetSymbol("Pd");
  PT[46].SetMass(106.42);
  PT[46].SetLattParams(3.89);
  PT[46].SetLattType("fcc");
  
  PT[47].SetName("silver");
  PT[47].SetSymbol("Ag");
  PT[47].SetMass(107.87);
  PT[47].SetLattParams(4.09);
  PT[47].SetLattType("fcc");

  PT[50].SetName("tin");
  PT[50].SetSymbol("Sn");
  PT[50].SetMass(118.71);
  PT[50].SetLattParams(5.83,3.18);
  PT[50].SetLattType("tetragonal");

  PT[51].SetName("antimony");
  PT[51].SetSymbol("Sb");
  PT[51].SetMass(121.75);
  PT[51].SetLattParams(4.51,57.1167);
  PT[51].SetLattType("rhombohedral");

  PT[54].SetName("xenon");
  PT[54].SetSymbol("Xe");
  PT[54].SetMass(131.29);
  PT[54].SetLattParams(6.19);
  PT[54].SetLattType("fcc");

  PT[55].SetName("cesium");
  PT[55].SetSymbol("Cs");
  PT[55].SetMass(132.91);
  PT[55].SetLattParams(6.14);
  PT[55].SetLattType("bcc");

  PT[73].SetName("tantulum");
  PT[73].SetSymbol("Ta");
  PT[73].SetMass(180.95);
  PT[73].SetLattParams(3.30);
  PT[73].SetLattType("bcc");
  
  PT[74].SetName("tungsten");
  PT[74].SetSymbol("W");
  PT[74].SetMass(183.85);
  PT[74].SetLattParams(3.17);
  PT[74].SetLattType("bcc");
  
  PT[78].SetName("platinum");
  PT[78].SetSymbol("Pt");
  PT[78].SetMass(195.08);
  PT[78].SetLattParams(3.92);
  PT[78].SetLattType("fcc");
  
  PT[79].SetName("gold");
  PT[79].SetSymbol("Au");
  PT[79].SetMass(196.97);
  PT[79].SetLattParams(4.08);
  PT[79].SetLattType("fcc");
  
  PT[82].SetName("lead");
  PT[82].SetSymbol("Pb");
  PT[82].SetMass(207.2);
  PT[82].SetLattParams(4.95);
  PT[82].SetLattType("fcc");
  
  PT[83].SetName("bismuth");
  PT[83].SetSymbol("Bi");
  PT[83].SetMass(208.98);
  PT[83].SetLattParams(4.75,57.2333);
  PT[83].SetLattType("rhombohedral");
}

PerTabEntryT PeriodicTableT::operator[] (const char * s) 
{
  int i;
  int loc = 0;
  /* StringT = dummy; */
  
  for(i=1;i<110;i++) 
    {
      if (s == PT[i].GetName()) loc = i;
      if (s == PT[i].GetSymbol()) loc = i;
    }

  if (loc>0) 
    {
      return PT[loc];
    }
  else 
    {
      throw eBadInputValue;
    }
}

PerTabEntryT PeriodicTableT::operator[] (const int m) 
{
  if (m>0 && m<110) 
    {
      return PT[m];
    }
  else 
    {
      throw eBadInputValue;
    }
}

