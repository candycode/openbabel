PROJECT( OpenBabel )

INCLUDE_DIRECTORIES( ./include ./data )

INCLUDE( ./sourcefiles.cmake )

IF( MSVC )
  ADD_DEFINITIONS( -DUSING_DYNAMIC_LIBS
                   -DOBDLL_EXPORTS
                   -DOBCONV_EXPORTS
                   -DOBFPRT_EXPORTS
                   -DOBCOMMON_EXPORTS
                   -DOBERROR_EXPORTS )

  IF( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
    SET( LIBNAME openbabeld ) 
  ELSE( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
    SET( LIBNAME openbabel )
  ENDIF( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )   
ELSE( MSVC )
  SET( LIBNAME openbabel )
ENDIF( MSVC )  

IF( MSVC )
  SET( INCHILIB libinchi )
  SET( XMLLIB libxml2 )
ELSE( MSVC )
  SET( INCHILIB inchi )
  SET( XMLLIB xml2 )
ENDIF( MSVC )  

#Always build inchi
SET( INCHI ON CACHE BOOL "Build inchi support" ) 

#Inchi
IF( INCHI )
  INCLUDE( ./libinchi_sourcefiles.cmake )
  INCLUDE_DIRECTORIES( ./include/inchi )
  ADD_LIBRARY( ${INCHILIB} ${LIBINCHI_SOURCES} )
ENDIF( INCHI )

#Always build XML
SET( XML ON CACHE BOOL "Build XML support" )

#XML - Needs libxml2
IF( XML )
  IF( WIN32 )
    IF( MSVC )
      SET( LIBXML_INCLUDE_DIR c:/vclibs/libxml2-2.6.30.win32/include CACHE PATH "libxml2 include dir (where libxml/ is)" )
      SET( LIBXML_LINK_DIR c:/vclibs/libxml2-2.6.30.win32/lib CACHE PATH "libxml2 lib dir" )
    ELSE( MSVC )
      SET( LIBXML_INCLUDE_DIR C:/mingw/devlib/release/include/libxml2 CACHE PATH "libxml2 include dir (where libxml/ is)" )
      SET( LIBXML_LINK_DIR C:/mingw/devlib/release/lib CACHE PATH "libxml2 lib dir" )
    ENDIF( MSVC )
  ELSE( WIN32 )
      SET( LIBXML_INCLUDE_DIR /usr/local/include CACHE PATH "libxml2 include dir (where libxml/ is)" )
      SET( LIBXML_LINK_DIR /usr/local/lib CACHE PATH "libxml2 lib dir" )
  ENDIF( WIN32 )
  INCLUDE_DIRECTORIES( ${LIBXML_INCLUDE_DIR} )
  LINK_DIRECTORIES( ${LIBXML_LINK_DIR} )
  IF( MSVC )
      SET( LIBICONV_INCLUDE_DIR C:/vclibs/iconv-1.9.2.win32/include CACHE PATH "iconv include dir" )
      SET( LIBICONV_LINK_DIR C:/vclibs/iconv-1.9.2.win32/lib CACHE PATH "iconv lib dir" )
      INCLUDE_DIRECTORIES( ${LIBICONV_INCLUDE_DIR} )
      LINK_DIRECTORIES( ${LIBICONV_LINK_DIR} )
  ENDIF( MSVC ) 	
ENDIF( XML )

##### LIBRARY
ADD_LIBRARY( ${LIBNAME} SHARED ${SOURCES} ${INCHI_SOURCES} ${XML_SOURCES} )

IF( INCHI )
  TARGET_LINK_LIBRARIES( ${LIBNAME} ${INCHILIB} )
ENDIF( INCHI )

IF( XML )
  TARGET_LINK_LIBRARIES( ${LIBNAME} ${XMLLIB} )
ENDIF( XML )

IF( MSVC AND XML )
  TARGET_LINK_LIBRARIES( ${LIBNAME} iconv )
ENDIF( MSVC AND XML )  
