#INCHI FORMAT
SET ( INCHI_SOURCES
        ./src/formats/getinchi.cpp
        ./src/formats/inchiformat.cpp )
    
#XML FORMATS
SET ( XML_SOURCES
        ./src/formats/xml/cmlreactlformat.cpp
        ./src/formats/xml/xml.cpp
        ./include/openbabel/xml.h
        ./src/formats/xml/cdxmlformat.cpp
        ./src/formats/xml/cmlformat.cpp
        ./src/formats/xml/pubchem.cpp
        ./src/formats/xml/xmlformat.cpp )

#DEFAULT: common classes + all formats except xml and inchi

SET( DL_HANDLER ./src/dlhandler_unix.cpp )

IF(WIN32)
  SET( DL_HANDLER ./src/dlhandler_win32.cpp )
ENDIF( WIN32 )

SET ( SOURCES
        ./src/atom.cpp
	./src/base.cpp
	./src/bitvec.cpp
	./src/bond.cpp
	./src/bondtyper.cpp
	./src/canon.cpp
	./src/chains.cpp
	./src/chiral.cpp
	./src/formats/cifformat.cpp	
	./src/data.cpp
	#./src/dlhandler_win32.cpp
        ${DL_HANDLER}
	./src/fingerprints/finger2.cpp
	./src/fingerprints/finger3.cpp
	./src/fingerprint.cpp
	./src/forcefield.cpp
	./src/generic.cpp	
	./src/grid.cpp
	./src/groupcontrib.cpp
	./src/kekulize.cpp
	./src/matrix.cpp
	./src/math/matrix3x3.cpp
	./src/mol.cpp
	./src/molchrg.cpp
	./src/obconversion.cpp
	./src/oberror.cpp
	./src/obiter.cpp
	./src/obutil.cpp
	./src/parsmart.cpp
	./src/patty.cpp
	./src/phmodel.cpp
	./src/rand.cpp
	./src/residue.cpp
	./src/ring.cpp
	./src/rotamer.cpp
	./src/rotor.cpp
	./src/tokenst.cpp
	./src/transform.cpp
	./src/typer.cpp
	./src/math/vector3.cpp	
	./include/openbabel/atom.h
	./include/openbabel/babelconfig.h
	./include/openbabel/base.h
	./include/openbabel/bitvec.h
	./include/openbabel/bond.h
	./include/openbabel/bondtyper.h
	./include/openbabel/canon.h
	./include/openbabel/chains.h
	./include/openbabel/chiral.h
	./include/openbabel/data.h
	./include/openbabel/dlhandler.h
	./include/openbabel/fingerprint.h
	./include/openbabel/forcefield.h
	./src/forcefields/forcefieldghemical.h
	./include/openbabel/generic.h
	./include/openbabel/grid.h
	./include/openbabel/groupcontrib.h
	./include/openbabel/internalcoord.h
	./include/openbabel/kinetics.h
	./include/openbabel/matrix.h
	./include/openbabel/math/matrix3x3.h
	./include/openbabel/mol.h
	./include/openbabel/molchrg.h
	./include/openbabel/obconversion.h
	./include/openbabel/oberror.h
	./include/openbabel/obiter.h
	./include/openbabel/obmolecformat.h
	./include/openbabel/obutil.h
	./include/openbabel/parsmart.h
	./include/openbabel/patty.h
	./include/openbabel/phmodel.h
	./include/openbabel/pluginiter.h
	./include/openbabel/rand.h
	./include/openbabel/reaction.h
	./include/openbabel/residue.h
	./include/openbabel/ring.h
	./include/openbabel/rotamer.h
	./include/openbabel/rotor.h
	./include/openbabel/typer.h
	./include/openbabel/math/vector3.h	
	./src/formats/acrformat.cpp
	./src/formats/alchemyformat.cpp
	./src/formats/amberformat.cpp
	./src/formats/APIInterface.cpp
	./src/formats/balstformat.cpp
	./src/formats/bgfformat.cpp
	./src/formats/boxformat.cpp
	./src/formats/cacaoformat.cpp
	./src/formats/cacheformat.cpp
	./src/formats/cansmilesformat.cpp
	./src/formats/carformat.cpp
	./src/formats/cccformat.cpp	
	./src/formats/chem3dformat.cpp
	./src/formats/chemdrawcdx.cpp
	./src/formats/chemdrawct.cpp
	./src/formats/chemtoolformat.cpp	
	./src/formats/copyformat.cpp
	./src/formats/crkformat.cpp
	./src/formats/CSRformat.cpp
	./src/formats/cssrformat.cpp
	./src/formats/dmolformat.cpp
	./src/formats/fastaformat.cpp
	./src/formats/fastsearchformat.cpp
	./src/formats/fchkformat.cpp
	./src/formats/featformat.cpp
	./src/formats/fhformat.cpp
	./src/formats/fingerprintformat.cpp
	./src/forcefields/forcefieldghemical.cpp
	./src/formats/freefracformat.cpp
	./src/formats/gamessformat.cpp
	./src/formats/gaussformat.cpp
	./src/formats/ghemicalformat.cpp
	./src/formats/gromos96format.cpp
	./src/formats/hinformat.cpp	
	./src/formats/jaguarformat.cpp
	./src/formats/mdlformat.cpp
	./src/formats/mmodformat.cpp
	./src/formats/mol2format.cpp
	./src/formats/molreport.cpp
	./src/formats/mopacformat.cpp
	./src/formats/mpdformat.cpp
	./src/formats/mpqcformat.cpp
	./src/formats/nwchemformat.cpp
	./src/formats/obmolecformat.cpp
	./src/formats/pcmodelformat.cpp
	./src/formats/pdbformat.cpp
	./src/formats/povrayformat.cpp
	./src/formats/PQSformat.cpp	
	./src/formats/qchemformat.cpp
	./src/formats/reportformat.cpp
	./src/formats/rxnformat.cpp
	./src/formats/shelxformat.cpp
	./src/formats/smilesformat.cpp
	./src/formats/thermoformat.cpp
	./src/formats/tinkerformat.cpp
	./src/formats/titleformat.cpp
	./src/formats/turbomoleformat.cpp
	./src/formats/unichemformat.cpp
	./src/formats/viewmolformat.cpp
	./src/formats/xedformat.cpp	
	./src/formats/xyzformat.cpp
	./src/formats/yasaraformat.cpp
	./src/formats/zindoformat.cpp )
