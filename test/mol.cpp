/**********************************************************************
mol.cpp - Unit tests for Open Babel OBMol class

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <stdio.h>
#include <iostream>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string d2file = testdatadir + "test2d.xyz";
  string d3file = testdatadir + "test3d.xyz";
#else
  string d2file = "files/test2d.xyz";
  string d3file = "files/test3d.xyz";
#endif

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      cout << "Usage: mol" << endl;
      cout << " Unit tests for OBMol " << endl;
      return(-1);
    }

  cout << "# Unit tests for OBMol \n";

  cout << "ok 1\n"; // for loading tests

  OBMol emptyMol, testMol1;
  cout << "ok 2\n"; // ctor works

  testMol1.ReserveAtoms(-1);
  testMol1.ReserveAtoms(0);
  testMol1.ReserveAtoms(2);
  cout << "ok 3\n";

  // atom component tests

  if (testMol1.NumAtoms() == 0) {
    cout << "ok 4\n";
  } else {
    cout << "not ok 4\n";
  }

  testMol1.NewAtom();
  if (testMol1.NumAtoms() == 1) {
    cout << "ok 5\n";
  } else {
    cout << "not ok 5\n";
  }

  testMol1.NewAtom();
  testMol1.AddBond(1, 2, 1);
  if (testMol1.NumBonds() == 1) {
    cout << "ok 6\n";
  } else {
    cout << "not ok 6\n";
  }

  testMol1.Clear();
  if (testMol1.NumAtoms() == 0) {
    cout << "ok 7\n";
  } else {
    cout << "not ok 7\n";
  }

  ifstream ifs1(d3file.c_str());
  if (!ifs1)
    {
      cout << "Bail out! Cannot read input file!" << endl;
      return(-1);
    }
  OBConversion conv(&ifs1, &cout);
  OBFormat* pFormat;
  
  pFormat = conv.FindFormat("XYZ");
  if ( pFormat == NULL )
    {
      cout << "Bail out! Cannot read file format!" << endl;
      return(-1);
    }
  if (! conv.SetInAndOutFormats(pFormat, pFormat))
    {
      cout << "Bail out! File format isn't loaded" << endl;
      return (-1);
    }

  OBMol testMol2D, testMol3D;
  if (conv.Read(&testMol3D))
    cout << "ok 8\n";
  else
    cout << "not ok 8\n";
  testMol3D.Center();
  
  // test bond insertion (PR#1665649)
  OBMol doubleBondMol;
  OBAtom *a1, *a2, *a3;
  OBBond *b;
  doubleBondMol.BeginModify();
  a1 = doubleBondMol.NewAtom();
  a1->SetVector(0.0, 0.0, 0.0);
  a1->SetAtomicNum(6);
  a2 = doubleBondMol.NewAtom();
  a2->SetVector(1.6, 0.0, 0.0);
  a2->SetAtomicNum(6);
  b = doubleBondMol.NewBond();
  b->SetBegin(a1);
  b->SetEnd(a2);
  a1->AddBond(b);
  a2->AddBond(b);
  doubleBondMol.EndModify();
  cout << "ok 9" << endl;

  // test AddHydrogens
  OBMol testMolH;
  testMolH.BeginModify();
  OBAtom *testAtom = testMolH.NewAtom();
  testAtom->SetVector(0.5f, 0.5f, 0.5f);
  testAtom->SetAtomicNum(6);
  testMolH.EndModify();
  testMolH.AddHydrogens();
  if (testMolH.NumAtoms() == 5) {
    cout << "ok 10" << endl;
  } else {
    cout << "not ok 10" << endl;
  }

  // test AddHydrogens (pr #1665519)
  OBMol testMolH2;
  OBAtom *testAtom2 = testMolH2.NewAtom();
  testAtom2->SetVector(0.5f, 0.5f, 0.5f);
  testAtom2->SetAtomicNum(6);
  testMolH2.AddHydrogens();
  if (testMolH2.NumAtoms() == 5) {
    cout << "ok 11" << endl;
  } else {
    cout << "not ok 11" << endl;
  }
  
  // test "h3" remove hydrogens case -- PR#1704551
  OBMol testMolH3;
  OBAtom *testAtom3 = testMolH3.NewAtom();
  testAtom3->SetAtomicNum(1);
  testAtom3 = testMolH3.NewAtom();
  testAtom3->SetAtomicNum(1);
  testAtom3 = testMolH3.NewAtom();
  testAtom3->SetAtomicNum(1);
  testMolH3.AddBond(1,2,1);
  testMolH3.AddBond(2,3,1);
  // OK, we now have 3 hydrogen atoms, 2 bonds
  testMolH3.DeleteHydrogens();
  // should have 0 atoms, 0 bonds
  if (testMolH3.NumAtoms() == 0) {
    cout << "ok 12" << endl;
  } else {
    cout << "not ok 12" << endl;
  }
  if (testMolH3.NumBonds() == 0) {
    cout << "ok 13" << endl;
  } else {
    cout << "not ok 13 # should delete all H-H bonds" << endl;
  }  

  // the total number of tests for "prove"
  // update when you add more tests!
  cout << "1..13\n";

  return(0);
}
