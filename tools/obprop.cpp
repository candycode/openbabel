/**********************************************************************
obprop - Open Babel properties calculation

Copyright (C) 2003 Fabien Fontaine
Some portions Copyright (C) 2004-2007 Geoffrey R. Hutchison
 
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
#include <openbabel/groupcontrib.h>
#include <openbabel/obconversion.h>
#ifndef _MSC_VER
	#include <unistd.h>
#endif

using namespace std;
using namespace OpenBabel;
void some_tests(OBMol &mol);
// PROTOTYPES /////////////////////////////////////////////////////////////////
int nrings(OBMol &mol);
string sequence(OBMol &mol);

///////////////////////////////////////////////////////////////////////////////
//! \brief Compute some properties easy to access from open babel
//


int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  char *FileIn = NULL;

  if (argc != 2)
    {
      string err = "Usage: ";
      err += program_name;
      err += " <filename>\n"
      "Output format:\n"
        "name NAME\n"
        "formula  FORMULA\n"
        "mol_weight MOLECULAR_WEIGHT\n"
        "exact_mass ISOTOPIC MASS\n"
        "canonical_SMILES STRING\n"
        "InChI  STRING\n"
        "num_atoms  NUM\n"
        "num_bonds  NUM\n"
        "num_residues  NUM\n"
        "sequence RESIDUE_SEQUENCE\n"
        "num_rings NUMBER_OF_RING_(SSSR)\n"
        "logP   NUM\n"
        "PSA    POLAR_SURFACE_AREA\n";
      err += "$$$$";
      ThrowError(err);
      exit(-1);
    }
  else
    {
      FileIn  = argv[1];
    }

  // Find Input filetype
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(FileIn);
    
  if (!format || !conv.SetInFormat(format))
    {
      cerr << program_name << ": cannot read input format!" << endl;
      exit (-1);
    }

  ifstream ifs;

  // Read the file
  ifs.open(FileIn);
  if (!ifs)
    {
      cerr << program_name << ": cannot read input file!" << endl;
      exit (-1);
    }
  
  OBMol mol;
  OBLogP logP;
  OBPSA psa;
  OBMR mr;
  OBFormat *canSMIFormat = conv.FindFormat("can");
  OBFormat *inchiFormat = conv.FindFormat("inchi");


  ////////////////////////////////////////////////////////////////////////////
  // List of properties
  // Name
  // Molecular weight (Standard molar mass given by IUPAC atomic masses)
  // Number of rings : the size of the smallest set of smallest rings (SSSR)
  
  //.....ADD YOURS HERE.....
  
  for (c = 1;; ++c)
    {
      mol.Clear();
      conv.Read(&mol, &ifs);
      if (mol.Empty())
        break;
      
      if (!mol.HasHydrogensAdded())
        mol.AddHydrogens();
      // Print the properties
      if (strlen(mol.GetTitle()) != 0)
        cout << "name             " << mol.GetTitle() << endl;
      else 
        cout << "name             " << FileIn << " " << c << endl;

      cout << "formula          " << mol.GetFormula() << endl;
      cout << "mol_weight       " << mol.GetMolWt() << endl;
      cout << "exact_mass       " << mol.GetExactMass() << endl;

      if (canSMIFormat) {
        conv.SetOutFormat(canSMIFormat);
        cout << "canonical_SMILES " << conv.WriteString(&mol);
      }
      else
        cout << "canonical_SMILES " << "-" << endl;

      if (inchiFormat) {
        conv.SetOutFormat(inchiFormat);
        cout << "InChI            " << conv.WriteString(&mol);
      }
      else
        cout << "InChI            " << "-" << endl;

      cout << "num_atoms        " << mol.NumAtoms() << endl;
      cout << "num_bonds        " << mol.NumBonds() << endl;
      cout << "num_residues     " << mol.NumResidues() << endl;
      if (mol.NumResidues() > 0)
        cout << "sequence         " << sequence(mol) << endl;
      else
        cout << "sequence         " << "-" << endl;

      cout << "num_rings        " << nrings(mol) << endl;
      cout << "logP             " << logP.Predict(mol) << endl;
      cout << "PSA              " << psa.Predict(mol) << endl;
      cout << "MR               " << mr.Predict(mol) << endl;
      cout << "$$$$" << endl; // SDF like end of compound descriptor list
      
    } // end for loop
  
  return(1);
}



///////////////////////////////////////////////////////////////////////////////
//! \return the number of size of the set of smallest rings (SSSR)
int nrings(OBMol &mol)
{
  int nr;
  vector<OBRing*> vr;
  
  vr = mol.GetSSSR();
  nr = vr.size();
  return (nr);
}

//! \return the sequence of residues ordered by chain
string sequence(OBMol &mol)
{
  unsigned int currentChain = 0;
  string residueSequence;
  FOR_RESIDUES_OF_MOL(r, mol)
    {
      if (r->GetName().find("HOH") != string::npos)
        continue;
      
      if (r->GetChainNum() != currentChain) {
        if (residueSequence.size() != 0) { // remove the trailing "-"
          residueSequence.erase(residueSequence.size() - 1);
          residueSequence += ", "; // separate different chains
        }
        
        currentChain = r->GetChainNum();
      }
      residueSequence += r->GetName();
      residueSequence += "-";
    }
  if (residueSequence.size() != 0) // remove the trailing "-"
    residueSequence.erase(residueSequence.size() - 1);

  return residueSequence;
}

/* obprop man page*/
/** \page obprop print standard molecular properties
*
* \n
* \par SYNOPSIS
*
* \b obprop \<filename\>
*
* \par DESCRIPTION
*
* The obprop program is a tool to print a set of standard molecular properties
* for all molecules in a file. It also serves as example code for using the
* Open Babel library (libopenbabel).
* 
* \par EXAMPLES
*
*   obprop pyridines.sdf
*
* \par AUTHORS
*
* The obprop program was contributed by \b Fabien \b Fontaine.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.sourceforge.net/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
*  Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison \n \n
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation version 2 of the License.\n \n
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
* \par SEE ALSO
*   The web pages for Open Babel can be found at: http://openbabel.sourceforge.net/ \n
**/
