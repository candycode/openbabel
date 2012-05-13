#include <cassert>
#include <limits>


#include <openbabel/babelconfig.h>

#include <openbabel/math/matrix3x3.h>
#include <openbabel/obmolecformat.h>

#include <openbabel/generic.h>

using namespace std;


#define COMPILE_STANDALONE__
#define TEST__

#ifdef COMPILE_STANDALONE__
namespace OpenBabel
{
	static OBElementTable etab;
}
#endif

/// Class used to store arrays of scalar values associated with each atom.
/// The arrays used to store information about per-atom data are 1-indexed and
/// element 0 always contains dummy data.
class OBAtomAssociatedData : public OpenBabel::OBGenericData
{
public:
	typedef std::vector< std::vector< double > > DataType;
	typedef std::vector< std::pair< double, double > > RangesType;

private:
	/// Array of arrays of scalar values; 1-indexed: array at position 0 is empty.
	/// The size of the data array is equal to the number of atoms + 1
	DataType data_;
	RangesType ranges_;
	
public:
	/// Default constructor: sets key and type values. 
	OBAtomAssociatedData() : 
	  OpenBabel::OBGenericData( "AtomAssociatedData",
								OpenBabel::OBGenericDataType::CustomData0,
								OpenBabel::fileformatInput ) {}

	/// Required Clone() method
	virtual OBAtomAssociatedData* Clone( OpenBabel::OBBase* /*parent*/ ) const
	{
		return new OBAtomAssociatedData( *this );
	}

	/// Set data method: standard way to store data into OBGenericData objects.
	/// @param data array of array of scalar values
	/// @param ranges array of per value ranges stored as an array of [min,max] pairs
	void SetData( const DataType& data,
				  const RangesType& ranges )
	{
		data_ = data;
		ranges_ = ranges;
	}
	
	/// Returns the data array.
	const DataType& GetData() const { return data_; }
	
	/// Returs the range array. The returned array is a sequence of [min,max] pairs.
	const RangesType& GetRanges() const { return ranges_; }
	
	/// Return the maximum number of per-atom values.
	/// E.g. if atom #1 has 3 values associated with it and
	/// every other atom has only 2 values this method returns 3.
	RangesType::size_type GetNumColumns() const { return ranges_.size(); }

	/// Returns an array of booleans of the same size as the data array.
	/// Each array element is set to true if data for the corresponding atom
	/// are available in the specified column, false otherwise. E.g.
	/// If atom #1 has 3 values associated with it and atom # 2 has only one value
	/// then calling @code HaveData( 1 ) @endcode will return an array like:
	/// @code
	/// { <dummy value>, true, false, ... }
	/// @endcode
	/// meaning that atom 2 does not have data in column 1.
	/// @param column zero based position of requested data
	/// @return boolean vector with information on which atoms have data in the
	///         requested column.
	std::vector< bool > HaveData( int column )
	{
		std::vector< bool > v( data_.size() ); //== atoms number + 1
		v.front() = true;
		DataType::iterator i = data_.begin();
		++i; //skip dummy atom
		std::vector< bool >::iterator bi = v.begin();
		++bi; //skip dummy atom
		for( ; i != data_.end(); ++i, ++bi ) *bi = i->size() < column;
	}
	
	/// Extracts one column of per-atom data and return it as a 1-indexed array.
	/// @param column column index.
	/// @param defaultValue default value used in case atom has no data in the specified column.
	std::vector< double > GetData( int column, double defaultValue )
	{
		std::vector< double > v( data_.size() ); // == num atoms + 1
		DataType::iterator i = data_.begin();
		std::vector< double >::iterator vi = v.begin();
		++i; //skip dummy atom
		*vi = defaultValue;
		for( ; i != data_.end(); ++i, ++vi )
		{
			if( i->size() < column ) *vi = i->at( column );
			else *vi = defaultValue;
		}
	}

	/// Utility function to print data to stream.
	friend inline std::ostream& operator<<( std::ostream& os, const OBAtomAssociatedData& ad )
	{
		RangesType::const_iterator ri = ad.ranges_.begin();
		os << "RANGES:" << std::endl;
		os << "-------" << std::endl;
		for( ; ri != ad.ranges_.end(); ++ri ) os << '[' << ri->first << ',' << ri->second << "] ";
		os << std::endl << std::endl;
		os << "VALUES:" << std::endl;
		os << "-------" << std::endl;
		DataType::const_iterator di = ad.data_.begin();
		++di;
		for( ; di != ad.data_.end(); ++di )
		{
			std::vector< double >::const_iterator vi = di->begin();
			for( ; vi != di->end(); ++vi )
			{
				os << *vi << ' ';
			}
			os << std::endl;
		}
		return os;
	}
};

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
    if( !pmol ) return false;

    //Define some references so we can use the old parameter names
    istream& ifs = *pConv->GetInStream();
     
	OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename

	std::string buffer;
    //  int natoms; CM
    double A = double(), B = double(), C = double(),
		   Alpha = double(), Beta = double(), Gamma = double();
    
	std::getline( ifs, buffer );
    mol.SetTitle(buffer);
	buffer = "";
	
	// move to CELL line
	while ( std::getline( ifs, buffer ) && buffer.find( "CELL" ) != 0 );
	
	if( !ifs ) return false;
    vector<string> vs;
    tokenize(vs,buffer," \t,\r\n");
    if (vs.size() != 8) return(false);

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
    //m = uc->GetOrthoMatrix();

    //  int i; CM
    double x,y,z;
    OBAtom *atom;
    vector3 v;

	//read symmetry 
	buffer = "";
	while( getline( ifs, buffer ) && buffer.find( "SYMM" ) != 0 );
	if( !ifs ) return false;
	std::vector< OpenBabel::matrix3x3 > matrices;
	std::vector< OpenBabel::vector3   > translation;
	do
	{
		vector3 row1( 1, 0, 0 );
		vector3 row2( 0, 1, 0 );
		vector3 row3( 0, 0, 1 );
		
		if( ifs )
		{
			istringstream ss( buffer );
			std::string t;
			ss >> t;
			if( buffer.find( 'X' ) != std::string::npos && 
				buffer.find( 'Y' ) != std::string::npos &&
				buffer.find( 'Z' ) != std::string::npos )
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
		buffer = "";
	}
	while( getline( ifs, buffer) &&  buffer.find( "SYMM" ) == 0 );
	
   
	mol.BeginModify();
	
	buffer = "";

	std::vector< std::vector< double > > atomData;
	// 1 indexed, start with dummy data for atom index 0
	atomData.push_back( std::vector< double >() );
	std::vector< std::pair< double, double > > mmax; // sequence of [min,max] ranges
    //read atom records
	while( std::getline( ifs, buffer ) && buffer != "HKLF" && buffer != "END" )
    {
		tokenize( vs, buffer.c_str(), " \t,\r\n" );
		if( vs.size() && vs[ 0 ] == "REM" ) continue;
		if ( vs.size() < 5 ) continue;
		int iso = int();
		const int atomicNum = etab.GetAtomicNum(
			std::string( vs[ 0 ], 0, vs[ 0 ].find_first_of( "0123456789" ) ).c_str() , iso );
		if( atomicNum < 1 )
		{
			buffer = "";
			continue;
		}
		
        atom = mol.NewAtom();

        x = atof( ( char* ) vs[ 2 ].c_str() );
        y = atof( ( char* ) vs[ 3 ].c_str() );
        z = atof( ( char* ) vs[ 4 ].c_str() );
        v.Set( x, y, z);
        
        atom->SetAtomicNum( atomicNum );
        atom->SetVector(v);
		const bool nextLine =  vs.back() == "=";
		std::vector< double > sv;
		// if data continue on next line remove last element ("=") from tokens
		if( nextLine ) vs.resize( vs.size() - 1 );
		const std::pair< double, double > 
			defpair( std::numeric_limits< double >::max(), -std::numeric_limits< double >::max() );
		// resize [min,max] sequence
		if( mmax.size() < vs.size() - 5 ) mmax.resize( vs.size() - 5, defpair );
		for( int i = 5; i != vs.size(); ++i )
		{
			istringstream is( vs[ i ] );
			double s = double();
			is >> s;
			mmax[ i - 5 ].first  = std::min( mmax[ i - 5 ].first, s );
			mmax[ i - 5 ].second = std::max( mmax[ i - 5 ].second, s );
			sv.push_back( s );
		}
		const int offset = vs.size() - 5;
		buffer = "";
		if( nextLine )
		{
			std::getline( ifs, buffer );
			if( !ifs ) break;
			std::vector< std::string > vs2;
			tokenize( vs2, buffer.c_str(), " \t,\r\n" );
			mmax.resize( std::max( mmax.size(), vs2.size() + vs.size()  - 5 ), defpair ); 
			for( int i = 0; i != vs2.size(); ++i )
			{
				istringstream is( vs2[ i ] );
				double s = double();
				is >> s;
				mmax[ offset + i ].first  = std::min( mmax[ offset + i ].first, s );
				mmax[ offset + i ].second = std::max( mmax[ offset + i ].second, s );
				sv.push_back( s );
			}
		}
		atomData.push_back( sv );
    } //while

	if( atomData.size() > 1 )
	{
		OBAtomAssociatedData* ad = new OBAtomAssociatedData;
		ad->SetData( atomData, mmax );
		mol.SetData( ad );
		std::cout << *ad;
	}

    if( !pConv->IsOption( "b", OBConversion::INOPTIONS ) ) mol.ConnectTheDots();
    if( !pConv->IsOption( "s", OBConversion::INOPTIONS ) && !pConv->IsOption( "b", OBConversion::INOPTIONS ) ) mol.PerceiveBondOrders();
    mol.EndModify();
    return true;
  }

} //namespace OpenBabel
