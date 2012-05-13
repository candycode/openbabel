/**********************************************************************
Copyright (C) 1998-2003 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/math/matrix3x3.h>
#include <openbabel/obmolecformat.h>

using namespace std;

namespace OpenBabel
{
    	OBElementTable etab;
}

namespace OpenBabel
{

  class ShelXFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ShelXFormat()
    {
      OBConversion::RegisterFormat("res",this, "chemical/x-shelx");
      OBConversion::RegisterFormat("ins",this, "chemical/x-shelx");
    }

    virtual const char* Description() //required
    {
      return
        "ShelX format\n \
       Read Options e.g. -as\n\
        s  Output single bonds only\n\
        b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://shelx.uni-ac.gwdg.de/SHELX/";}; //optional

    virtual const char* GetMIMEType() 
    { return "chemical/x-shelx"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTWRITABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  ShelXFormat theShelXFormat;

  vector3 ParseVector( const std::string s, std::string::value_type v, double& trans )
  {
	  std::vector< std::string > tokens;
	  istringstream is( s );
	  int vpos = std::string::npos;
	  int i = 0;
	  while( is )
	  {
		  std::string buf;
		  getline( is, buf, ' ' );
		  tokens.push_back( buf );
		  ++i;
	  }
	  trans = 0.;
	  return vector3( 0, 0, 0 );

  }

  /////////////////////////////////////////////////////////////////
  bool ShelXFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream& iss = *pConv->GetInStream();
	ostringstream oss;
	std::string buf;
	while( iss )
	{
		std::getline( iss, buf );
		oss << buf << '\n';
	}
	istringstream ifs( oss.str() );// = *pConv->GetInStream();
     
	OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename

    char buffer[BUFF_SIZE];
    //  int natoms; CM
    double A,B,C,Alpha,Beta,Gamma;
    matrix3x3 m;

    ifs.getline(buffer,BUFF_SIZE);
    mol.SetTitle(buffer);
    while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"CELL",4));

    if (!EQn(buffer,"CELL",4))
      return(false);
    vector<string> vs;
    tokenize(vs,buffer," \n\t,");
    if (vs.size() != 8)
      return(false);

    //parse cell values
    A = atof((char*)vs[2].c_str());
    B = atof((char*)vs[3].c_str());
    C = atof((char*)vs[4].c_str());
    Alpha = atof((char*)vs[5].c_str());
    Beta  = atof((char*)vs[6].c_str());
    Gamma = atof((char*)vs[7].c_str());
    OBUnitCell *uc = new OBUnitCell;
    uc->SetOrigin(fileformatInput);
    uc->SetData(A, B, C, Alpha, Beta, Gamma);
    mol.SetData(uc);
    m = uc->GetOrthoMatrix();

    //  int i; CM
    double x,y,z;
    char type[16], *j;
    OBAtom *atom;
    vector3 v;

	//read symmetry
	
	std::vector< OpenBabel::matrix3x3 > matrices;
	std::vector< OpenBabel::vector3   > translation;
	while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"SYMM",4))
	{
		vector3 row1( 1, 0, 0 );
		vector3 row2( 0, 1, 0 );
		vector3 row3( 0, 0, 1 );
		
		if( ifs )
		{
			istringstream ss( buffer );
			std::string t;
			ss >> t;
			if( std::find( buffer, buffer + 128, 'X' ) != buffer + 128 &&
				std::find( buffer, buffer + 128, 'Y' ) != buffer + 128 &&
				std::find( buffer, buffer + 128, 'Z' ) != buffer + 128 )
			{
				double t1 = double();
				double t2 = double();
				double t3 = double();
				std::getline( ss, t, ',' );
				row1 = ParseVector( t, 'X', t1 );
				t = "";
				getline( ss, t, ',' );
				row2 = ParseVector( t, 'Y', t2 );
				t = "";
				getline( ss, t );
				row3 = ParseVector( t, 'Z', t3 );
				matrix3x3 m( row1, row2, row3 );
				vector3 xt( t1, t2, t3 );
				matrices.push_back( m );
				translation.push_back( xt );
			}
		}
	}
	
    //skip until FVAR or SFAC found
	
	std::vector< double > freeVars;
	const istream::pos_type pos = ifs.tellg();
    while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"FVAR",4));
	if( !ifs )
	{
		ifs.clear();
		ifs.str( oss.str() );
		ifs.seekg( 0, ios::beg );
		while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"SFAC",4));
	}
	else 
	{
		double v = 0.;
		std::istringstream is( buffer );
		std::string fvar;
		is >> fvar;
		while( is )
		{
			is >> v;
			freeVars.push_back( v );
		}		
	}

	mol.BeginModify();
	
	buffer[0] = '\0';
    //read atom records
    while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"HKLF",4))
      {
        tokenize(vs,buffer," \n\t,");

        //skip AFIX and PART instructions
        if (vs.size() < 6 + freeVars.size() ) continue;
        atom = mol.NewAtom();

        x = atof((char*)vs[2].c_str());
        y = atof((char*)vs[3].c_str());
        z = atof((char*)vs[4].c_str());
        v.Set(x,y,z);
        v *= m;

        strncpy(type,vs[0].c_str(), sizeof(type));
        type[sizeof(type) - 1] = '\0';
        j = strpbrk(type, "0123456789");
        j[0] = '\0';
        atom->SetAtomicNum( etab.GetAtomicNum(type) );
		std::cout << type << ' ' << etab.GetName( etab.GetAtomicNum(type) ) << ' ' << etab.GetName( 300 ) << std::endl;
        atom->SetVector(v);
		
        //skip next line if anisotropic atoms.
        if (vs.size() == 9)
          ifs.getline(buffer,BUFF_SIZE);
      } //while

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    return(true);
  }

} //namespace OpenBabel
