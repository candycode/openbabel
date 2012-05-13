PROJECT( OpenBabelTools )


SET( OPENBABEL_INCLUDE_DIR /usr/local/include CACHE PATH "OpenBabel include dir" )
SET( OPENBABEL_LINK_DIR /usr/local/lib CACHE PATH "OpenBabel link dir" )

INCLUDE_DIRECTORIES( ${OPENBABEL_INCLUDE_DIR} )
LINK_DIRECTORIES( ${OPENBABEL_LINK_DIR} )


IF( MSVC )
  IF( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
    SET( LIBNAME openbabeld ) 
  ELSE( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
    SET( LIBNAME openbabel )
  ENDIF( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )   
ELSE( MSVC )
  SET( LIBNAME openbabel )
ENDIF( MSVC )  


##### EXECUTABLES
ADD_EXECUTABLE( babel babel.cpp )
TARGET_LINK_LIBRARIES( babel ${LIBNAME} )

ADD_EXECUTABLE( obchiral obchiral.cpp )
TARGET_LINK_LIBRARIES( obchiral ${LIBNAME} )

ADD_EXECUTABLE( obenergy obenergy.cpp )
TARGET_LINK_LIBRARIES( obenergy ${LIBNAME} )

ADD_EXECUTABLE( obfit obfit.cpp )
TARGET_LINK_LIBRARIES( obfit ${LIBNAME} )

#WARNING: requires XGetOpt on MSVC/Windows environments
#ADD_EXECUTABLE( obgrep obgrep.cpp )
#TARGET_LINK_LIBRARIES( obgrep ${LIBNAME} )

ADD_EXECUTABLE( obminimize obminimize.cpp )
TARGET_LINK_LIBRARIES( obminimize ${LIBNAME} )

ADD_EXECUTABLE( obprop obprop.cpp )
TARGET_LINK_LIBRARIES( obprop ${LIBNAME} )

ADD_EXECUTABLE( obrotamer obrotamer.cpp )
TARGET_LINK_LIBRARIES( obrotamer ${LIBNAME} )

ADD_EXECUTABLE( obrotate obrotate.cpp )
TARGET_LINK_LIBRARIES( obrotate ${LIBNAME} )