# The YADE has the following parameters to configure:
#  INSTALL_PREFIX: path, where Yade will be installed (/usr/local by default)
#  LIBRARY_OUTPUT_PATH: path to install libraries (lib by default)
#  DEBUG: compile in debug-mode (OFF by default)
#  CMAKE_VERBOSE_MAKEFILE: output additional information during compiling (OFF by default)
#  SUFFIX: suffix, added after binary-names (version number by default)
#  NOSUFFIX: do not add a suffix after binary-name (OFF by default)
#  YADE_VERSION: explicitely set version number (is defined from git-directory by default)
#  ENABLE_GUI: enable GUI option (ON by default)
#  ENABLE_CGAL: enable CGAL option (ON by default)
#  ENABLE_VTK: enable VTK-export option (ON by default)
#  ENABLE_OPENMP: enable OpenMP-parallelizing option (ON by default)
#  ENABLE_GTS: enable GTS-option (ON by default)
#  ENABLE_GL2PS: enable GL2PS-option (ON by default)
#  ENABLE_LINSOLV: enable LINSOLV-option (ON by default)
#  ENABLE_PFVFLOW: enable PFVFLOW-option, FlowEngine (ON by default)
#  ENABLE_SPH: enable SPH-option, Smoothed Particle Hydrodynamics (OFF by default, experimental)
#  ENABLE_LBMFLOW: enable LBMFLOW-option, LBM_ENGINE (OFF by default)
#  ENABLE_LIQCONTROL: enable LIQCONTROL-option, (OFF by default)
#  runtimePREFIX: used for packaging, when install directory is not the same is runtime directory (/usr/local by default)
#  CHUNKSIZE: set >1, if you want several sources to be compiled at once. Increases compilation speed and RAM-consumption during it (1 by default).


project(Yade C CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cMake")

INCLUDE(FindPythonLibs)
INCLUDE(FindPythonInterp)
INCLUDE(FindOpenMP)
INCLUDE(FindQt4)
INCLUDE(FindVTK)
INCLUDE(FindPkgConfig)
INCLUDE(GetVersion)
INCLUDE(FindOpenGL)
INCLUDE(FindGTS)
INCLUDE(FindGL2PS)
INCLUDE(FindCGAL)
INCLUDE(FindNumPy)
INCLUDE(FindLoki)
INCLUDE(FindPythonModule)
INCLUDE(GNUInstallDirs)
#===========================================================
# HACK!!! If the version of gcc is 4.8 or greater, we add -ftrack-macro-expansion=0 
# and -save-temps into compiler to reduce the memory consumption during compilation.
# See http://bugs.debian.org/726009 for more information
# Can be removed later, if gcc fixes its regression
# Taken from http://stackoverflow.com/questions/4058565/check-gcc-minor-in-cmake
 
EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
IF (GCC_VERSION VERSION_GREATER 4.8 OR GCC_VERSION VERSION_EQUAL 4.8)
  MESSAGE(STATUS "GCC Version >= 4.8. Adding -ftrack-macro-expansion=0 and -save-temps")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -ftrack-macro-expansion=0 -save-temps")
ENDIF()

#===========================================================

IF ("${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" MATCHES ".*clang")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -ftemplate-depth-512")
ENDIF()

#===========================================================
ADD_DEFINITIONS(" -DYADE_PTR_CAST=static_pointer_cast -DYADE_CAST=static_cast ")
IF (CMAKE_CXX_FLAGS)
  #If flags are set, add only neccessary flags
  IF (DEBUG)
    SET(CMAKE_VERBOSE_MAKEFILE 1)
    SET(CMAKE_BUILD_TYPE Debug)
    ADD_DEFINITIONS("-DYADE_DEBUG")
  ELSE (DEBUG)
    SET(CMAKE_BUILD_TYPE Release)
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fPIC")
  ENDIF (DEBUG)
ELSE (CMAKE_CXX_FLAGS)
  #If flags are not set, add all useful flags
  IF (DEBUG)
    SET(CMAKE_VERBOSE_MAKEFILE 1)
    SET(CMAKE_BUILD_TYPE Debug)
    ADD_DEFINITIONS("-DYADE_DEBUG")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fPIC -g -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security")
  ELSE (DEBUG)
    SET(CMAKE_BUILD_TYPE Release)
    SET(CMAKE_CXX_FLAGS  "-Wall -fPIC -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -s")
  ENDIF (DEBUG)
ENDIF (CMAKE_CXX_FLAGS)
#===========================================================
# Add possibility to use local boost installation (e.g. -DLocalBoost=1.46.1)

FIND_PACKAGE(Boost ${LocalBoost} COMPONENTS python thread date_time filesystem iostreams regex serialization system REQUIRED)
INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})
# for checking purpose
MESSAGE("--   Boost_VERSION: " ${Boost_VERSION})
MESSAGE("--   Boost_LIB_VERSION: " ${Boost_LIB_VERSION})
MESSAGE("--   Boost_INCLUDE_DIRS: " ${Boost_INCLUDE_DIRS})
MESSAGE("--   Boost_LIBRARIES: " ${Boost_LIBRARIES})

#===========================================================
FIND_PACKAGE(NumPy REQUIRED)
  INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_DIRS})

FIND_PACKAGE(Loki REQUIRED)
  INCLUDE_DIRECTORIES(${LOKI_INCLUDE_DIR})
FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(BZip2 REQUIRED)
FIND_PACKAGE(ZLIB REQUIRED)
#===========================================================

SET(DEFAULT ON CACHE INTERNAL "Default value for enabled by default options")
SET(LINKLIBS "")
SET(CONFIGURED_FEATS "")
SET(DISABLED_FEATS "")

OPTION(ENABLE_VTK "Enable VTK" ${DEFAULT})
OPTION(ENABLE_OPENMP "Enable OpenMP" ${DEFAULT})
OPTION(ENABLE_GTS "Enable GTS" ${DEFAULT})
OPTION(ENABLE_GUI "Enable GUI" ${DEFAULT})
OPTION(ENABLE_CGAL "Enable CGAL" ${DEFAULT})
OPTION(ENABLE_GL2PS "Enable GL2PS" ${DEFAULT})
OPTION(ENABLE_LINSOLV "Enable direct solver for the flow engines (experimental)" ${DEFAULT})
OPTION(ENABLE_PFVFLOW "Enable flow engine (experimental)" ${DEFAULT})
OPTION(ENABLE_SPH "Enable SPH" OFF)
OPTION(ENABLE_LBMFLOW "Enable LBM engine (very experimental)" OFF)
OPTION(ENABLE_LIQCONTROL "Enable liquid control (very experimental)" OFF)

#===========================================================
# Use Eigen3 by default
IF (EIGEN3_FOUND)
  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
  MESSAGE(STATUS "Found Eigen3")
  SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} Eigen3")
ENDIF(EIGEN3_FOUND)
#===========================================================
INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
SET(LINKLIBS  "${LINKLIBS};${BZIP2_LIBRARIES};${ZLIB_LIBRARIES};")
#===========================================================
IF((Boost_MAJOR_VERSION EQUAL 1) OR (Boost_MAJOR_VERSION GREATER 1) AND
  ((Boost_MINOR_VERSION EQUAL 53) OR (Boost_MINOR_VERSION GREATER 53)))
  SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} Odeint")
  ADD_DEFINITIONS("-DYADE_ODEINT")
ELSE((Boost_MAJOR_VERSION EQUAL 1) OR (Boost_MAJOR_VERSION GREATER 1) AND
    ((Boost_MINOR_VERSION EQUAL 53) OR (Boost_MINOR_VERSION GREATER 53)))
  SET(DISABLED_FEATS "${DISABLED_FEATS} Odeint")
  MESSAGE(STATUS "Boost Odeint can be enabled, only if Boost>=1.53 is used")
  MESSAGE(STATUS ${Boost_MAJOR_VERSION})
  MESSAGE(STATUS ${Boost_MINOR_VERSION})
ENDIF((Boost_MAJOR_VERSION EQUAL 1) OR (Boost_MAJOR_VERSION GREATER 1) AND
     ((Boost_MINOR_VERSION EQUAL 53) OR (Boost_MINOR_VERSION GREATER 53)))
#===========================================================
IF(ENABLE_VTK)
  FIND_PACKAGE(VTK) 
  IF(VTK_FOUND)
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
    LINK_DIRECTORIES( ${VTK_LIBRARY_DIRS} )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DYADE_VTK")
    SET(features "${features} vtk")
    MESSAGE(STATUS "Found VTK")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} VTK")
  ELSE(VTK_FOUND)
    MESSAGE(STATUS "VTK NOT found")
    SET(ENABLE_VTK OFF)
    SET(DISABLED_FEATS "${DISABLED_FEATS} VTK")
  ENDIF(VTK_FOUND)
ELSE(ENABLE_VTK)
  SET(DISABLED_FEATS "${DISABLED_FEATS} VTK")
ENDIF(ENABLE_VTK)
#===========================================================
IF(ENABLE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DYADE_OPENMP ${OpenMP_CXX_FLAGS}")
    MESSAGE(STATUS "Found OpenMP")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} OpenMP")
    SET(features "${features} openmp")
  ELSE(OPENMP_FOUND)
    MESSAGE(STATUS "OpenMP NOT found")
    SET(ENABLE_OPENMP OFF)
    SET(DISABLED_FEATS "${DISABLED_FEATS} OPENMP")
  ENDIF(OPENMP_FOUND)
ELSE(ENABLE_OPENMP)
  SET(DISABLED_FEATS "${DISABLED_FEATS} OPENMP")
ENDIF(ENABLE_OPENMP)
#===========================================================
IF(ENABLE_GTS)
  FIND_PACKAGE(GTS)
  FIND_PACKAGE(glib2)
  IF(GTS_FOUND AND GLIB2_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DYADE_GTS ${CMAKE_GTS_CXX_FLAGS}")
    SET(CMAKE_LD_FLAGS  "${CMAKE_LD_FLAGS} ${GTS_LIBRARIES}")
    INCLUDE_DIRECTORIES(${GTS_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${GLIB2_INCLUDE_DIRS})
    MESSAGE(STATUS "Found GTS")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} GTS")
    SET(features "${features} gts")
  ELSE(GTS_FOUND AND GLIB2_FOUND)
    MESSAGE(STATUS "GTS NOT found")
    SET(DISABLED_FEATS "${DISABLED_FEATS} GTS")
    SET(ENABLE_GTS OFF)
  ENDIF(GTS_FOUND AND GLIB2_FOUND)
ELSE(ENABLE_GTS)
  SET(DISABLED_FEATS "${DISABLED_FEATS} GTS")
ENDIF(ENABLE_GTS)
#===========================================================
IF(ENABLE_GUI)
  FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui QtOpenGL)
  FIND_PACKAGE(OpenGL)
  FIND_PACKAGE(GLUT)
  FIND_PACKAGE(glib2)
  FIND_PACKAGE(QGLVIEWER)
  IF(QT4_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND GLIB2_FOUND AND QGLVIEWER_FOUND)
    SET(GUI_LIBS ${GLUT_LIBRARY} ${OPENGL_LIBRARY} ${QGLVIEWER_LIBRARIES})
    SET(GUI_SRC_LIB "lib/opengl/GLUtils.cpp")
    SET(features "${features} qt4 opengl")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DYADE_OPENGL")
    INCLUDE_DIRECTORIES(${GLIB2_INCLUDE_DIRS})
    INCLUDE_DIRECTORIES(${QT_INCLUDES})
    
    MESSAGE(STATUS "Found GUI-LIBS")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} GUI")
  ELSE(QT4_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND GLIB2_FOUND AND QGLVIEWER_FOUND)
    MESSAGE(STATUS "GUI-LIBS NOT found")
    SET(DISABLED_FEATS "${DISABLED_FEATS} GUI")
    SET(ENABLE_GUI OFF)
  ENDIF(QT4_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND GLIB2_FOUND AND QGLVIEWER_FOUND)
ELSE(ENABLE_GUI)
  SET(DISABLED_FEATS "${DISABLED_FEATS} GUI")
ENDIF(ENABLE_GUI)
#===========================================================
# This one will automatically enable CGAL
IF(ENABLE_PFVFLOW)
  IF (NOT ENABLE_CGAL)
    MESSAGE(STATUS "PFVFLOW depends on CGAL, attempting to turn ENABLE_CGAL ON")
    SET(ENABLE_CGAL ON)
  ENDIF (NOT ENABLE_CGAL)
ENDIF(ENABLE_PFVFLOW)
#===========================================================
IF(ENABLE_CGAL)  
  INCLUDE(FindGMP)
  FIND_PACKAGE(CGAL)
  FIND_PACKAGE(GMP)
  IF(CGAL_FOUND AND GMP_FOUND AND (NOT("${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" MATCHES ".*clang"))) #Check for clang should be removed, when CGAL will be compilable by clang
    #INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR})
    #INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
    #INCLUDE_DIRECTORIES(${GMPXX_INCLUDE_DIR})
    
    FILE(GLOB CGAL_SRC_LIB "lib/triangulation/*.cpp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_DEFINITIONS} -DYADE_CGAL")
    SET(LINKLIBS  "${LINKLIBS};${CGAL_LIBRARIES};${GMP_LIBRARIES};${GMPXX_LIBRARIES};")
    MESSAGE(STATUS "Found CGAL")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} CGAL")

    IF(ENABLE_PFVFLOW)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFLOW_ENGINE")
      SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} PFVflow")
    ELSE(ENABLE_PFVFLOW)
      SET(DISABLED_FEATS "${DISABLED_FEATS} PFVflow")
    ENDIF(ENABLE_PFVFLOW)

  ELSE(CGAL_FOUND AND GMP_FOUND AND (NOT("${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" MATCHES ".*clang")))
    MESSAGE(STATUS "CGAL NOT found")
    SET(DISABLED_FEATS "${DISABLED_FEATS} CGAL")
    SET(ENABLE_CGAL OFF)
    IF(ENABLE_PFVFLOW)
      SET(DISABLED_FEATS "${DISABLED_FEATS} PFVflow")
      MESSAGE(STATUS "CGAL NOT found: PFVFLOW disabled")
    ENDIF(ENABLE_PFVFLOW)

  ENDIF(CGAL_FOUND AND GMP_FOUND AND (NOT("${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" MATCHES ".*clang")))
ELSE(ENABLE_CGAL)
  SET(DISABLED_FEATS "${DISABLED_FEATS} CGAL")
ENDIF(ENABLE_CGAL)
#===========================================================
IF(ENABLE_LINSOLV)
  FIND_PACKAGE(Cholmod)
  FIND_PACKAGE(OpenBlas)
  FIND_PACKAGE(Metis)
  IF(CHOLMOD_FOUND AND OPENBLAS_FOUND AND METIS_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_DEFINITIONS} -DLINSOLV -DFLOW_ENGINE")
    SET(LINKLIBS  "${LINKLIBS};${CHOLMOD_LIBRARIES};${AMD_LIBRARY};${CAMD_LIBRARY};${COLAMD_LIBRARY};${CCOLAMD_LIBRARY};${OPENBLAS_LIBRARY};${METIS_LIBRARY}")
    INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR} ${CHOLMOD_INCLUDE_DIR})
    MESSAGE(STATUS "Found Cholmod")
    MESSAGE(STATUS "Found OpenBlas")
    MESSAGE(STATUS "Found Metis")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} LinSolv")
  ELSE(CHOLMOD_FOUND AND OPENBLAS_FOUND AND METIS_FOUND)
    MESSAGE(STATUS "Missing dependency for LINSOLV, disabled")
    SET(DISABLED_FEATS "${DISABLED_FEATS} LinSolv")
    SET(ENABLE_LINSOLV OFF)
  ENDIF(CHOLMOD_FOUND AND OPENBLAS_FOUND AND METIS_FOUND)
ELSE(ENABLE_LINSOLV)
  SET(DISABLED_FEATS "${DISABLED_FEATS} LinSolv")
ENDIF(ENABLE_LINSOLV)
#===============================================
IF(ENABLE_SPH)
  ADD_DEFINITIONS("-DYADE_SPH")
  SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} SPH")
ELSE(ENABLE_SPH)
  SET(DISABLED_FEATS "${DISABLED_FEATS} SPH")
ENDIF(ENABLE_SPH)
#===============================================
IF(ENABLE_LIQCONTROL)
  ADD_DEFINITIONS("-DYADE_LIQCONTROL")
  SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} LIQCONTROL")
ELSE(ENABLE_LIQCONTROL)
  SET(DISABLED_FEATS "${DISABLED_FEATS} LIQCONTROL")
ENDIF(ENABLE_LIQCONTROL)
#===============================================
IF(ENABLE_GL2PS)
  FIND_PACKAGE(GL2PS)
  IF(GL2PS_FOUND)
    INCLUDE_DIRECTORIES(${GL2PS_INCLUDE_DIR})
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DYADE_GL2PS")
    SET(LINKLIBS  "${LINKLIBS};${GL2PS_LIBRARIES};")
    MESSAGE(STATUS "Found GL2PS")
    SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} GL2PS")
  ELSE(GL2PS_FOUND)
    MESSAGE(STATUS "GL2PS NOT found")
    SET(DISABLED_FEATS "${DISABLED_FEATS} GL2PS")
    SET(ENABLE_GL2PS OFF)
  ENDIF(GL2PS_FOUND)
ELSE(ENABLE_GL2PS)
  SET(DISABLED_FEATS "${DISABLED_FEATS} GL2PS")
ENDIF(ENABLE_GL2PS)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/extra/floating_point_utilities_v3)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

#===========================================================
IF(ENABLE_LBMFLOW)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLBM_ENGINE")
  SET(CONFIGURED_FEATS "${CONFIGURED_FEATS} LBMFLOW")
  MESSAGE("LBMFLOW is still experimental, building and running LBM engine are at your own risk!")
ELSE(ENABLE_LBMFLOW)
  SET(DISABLED_FEATS "${DISABLED_FEATS} LBMFLOW")
ENDIF(ENABLE_LBMFLOW)

#===========================================================

IF (NOT INSTALL_PREFIX)
 SET(CMAKE_INSTALL_PREFIX "/usr/local")
 MESSAGE("Yade will be installed to default path ${CMAKE_INSTALL_PREFIX}, if you want to override it use -DINSTALL_PREFIX option.")
ELSE (NOT INSTALL_PREFIX)
	GET_FILENAME_COMPONENT(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} ABSOLUTE)
	MESSAGE("Yade will be installed to ${CMAKE_INSTALL_PREFIX}")
ENDIF (NOT INSTALL_PREFIX)

IF (NOT SUFFIX)
  SET (SUFFIX "-${YADE_VERSION}")
ENDIF (NOT SUFFIX)

IF(NOSUFFIX)   #For packaging
  SET (SUFFIX "")
ENDIF(NOSUFFIX)   #For packaging

IF(NOT LIBRARY_OUTPUT_PATH)   #For packaging
  SET (LIBRARY_OUTPUT_PATH ${CMAKE_INSTALL_LIBDIR})
ENDIF(NOT LIBRARY_OUTPUT_PATH)   #For packaging

IF (NOT runtimePREFIX)
  SET (runtimePREFIX ${CMAKE_INSTALL_PREFIX})
ENDIF (NOT runtimePREFIX)

MESSAGE (STATUS "Suffix is set to " ${SUFFIX})
MESSAGE (STATUS "LIBRARY_OUTPUT_PATH is set to " ${LIBRARY_OUTPUT_PATH})
MESSAGE (STATUS "runtimePREFIX is set to " ${runtimePREFIX})
#===========================================================

SET(YADE_LIB_PATH ${CMAKE_INSTALL_PREFIX}/${LIBRARY_OUTPUT_PATH}/yade${SUFFIX})
SET(YADE_EXEC_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
SET(YADE_PY_PATH ${YADE_LIB_PATH}/py)
SET(YADE_DOC_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/doc/yade${SUFFIX})
SET(YADE_MAN_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_MANDIR})


SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${YADE_LIB_PATH};${YADE_PY_PATH};${YADE_PY_PATH}/yade/;${YADE_PY_PATH}/yade/qt;${YADE_PY_PATH}/gts")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#===========================================================
IF(ENABLE_GUI)
  ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/gui")
ENDIF(ENABLE_GUI)
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/py")
#===========================================================

# this is only a temporary hack, headers should be installed in the build-dir directly
# perhaps there is a proper way, have to ask
IF(NOT EXISTS "${CMAKE_BINARY_DIR}/yade")
  EXECUTE_PROCESS(COMMAND ln -s ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}/yade)
ENDIF ()
#===========================================================

FILE(GLOB SRC_CORE "core/*.cpp")
FILE(GLOB_RECURSE SRC_PKG  "pkg/*.cpp")
FILE(GLOB SRC_LIB  "lib/*.cpp")

SET(SRC_LIB "${SRC_LIB};lib/base/Math.cpp;lib/factory/ClassFactory.cpp;lib/factory/DynLibManager.cpp;lib/multimethods/Indexable.cpp;lib/serialization/Serializable.cpp;lib/pyutil/gil.cpp;core/main/pyboot.cpp;${GUI_SRC_LIB};${CGAL_SRC_LIB}")

#===========================================================

IF (CHUNKSIZE)
  INCLUDE(CombineSources)
  COMBINE_SOURCES(${CMAKE_BINARY_DIR}/core "${SRC_CORE}" ${CHUNKSIZE})
  FILE(GLOB SRC_CORE_COMBINED "${CMAKE_BINARY_DIR}/core.*.cpp")
  COMBINE_SOURCES(${CMAKE_BINARY_DIR}/pkg "${SRC_PKG}" ${CHUNKSIZE})
  FILE(GLOB SRC_PKG_COMBINED "${CMAKE_BINARY_DIR}/pkg.*.cpp")
  COMBINE_SOURCES(${CMAKE_BINARY_DIR}/lib "${SRC_LIB}" ${CHUNKSIZE})
  FILE(GLOB SRC_LIB_COMBINED "${CMAKE_BINARY_DIR}/lib.*.cpp")
  ADD_LIBRARY(yade SHARED ${SRC_LIB_COMBINED} ${SRC_CORE_COMBINED} ${SRC_PKG_COMBINED})
ELSE (CHUNKSIZE)
  ADD_LIBRARY(yade SHARED ${SRC_CORE} ${SRC_PKG} ${SRC_LIB})
ENDIF (CHUNKSIZE)
#===========================================================
find_python_module(minieigen)
IF (PY_minieigen)
  MESSAGE(STATUS "Use system minieigen version")
ELSE (PY_minieigen)
  MESSAGE(STATUS "Use embedded version of minieigen. Please, consider installing the corresponding package")
ENDIF (PY_minieigen)
find_python_module(Tkinter REQUIRED)
#===========================================================

ADD_LIBRARY(boot SHARED ${CMAKE_CURRENT_SOURCE_DIR}/core/main/pyboot.cpp)
SET_TARGET_PROPERTIES(boot PROPERTIES PREFIX "")
TARGET_LINK_LIBRARIES(yade ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${LINKLIBS} ${LOKI_LIBRARY} -lrt)
TARGET_LINK_LIBRARIES(boot yade)

IF(ENABLE_VTK)
  IF(${VTK_MAJOR_VERSION} EQUAL 6)
    TARGET_LINK_LIBRARIES(yade ${VTK_LIBRARIES})
    ADD_DEFINITIONS("-DYADE_VTK6")
    MESSAGE(STATUS "VTK version 6 is found")
  ELSE(${VTK_MAJOR_VERSION} EQUAL 6)
    TARGET_LINK_LIBRARIES(yade vtkHybrid)
  ENDIF(${VTK_MAJOR_VERSION} EQUAL 6)
ENDIF(ENABLE_VTK)
IF(ENABLE_GUI)
  TARGET_LINK_LIBRARIES(yade _GLViewer ${GUI_LIBS})
ENDIF(ENABLE_GUI)

#====================================
IF (NOT (PY_minieigen))
  ADD_LIBRARY(miniEigen SHARED py/mathWrap/miniEigen.cpp)
  SET_TARGET_PROPERTIES(miniEigen PROPERTIES PREFIX "")
  TARGET_LINK_LIBRARIES(miniEigen ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
  INSTALL(TARGETS miniEigen DESTINATION ${YADE_PY_PATH})
ENDIF (NOT (PY_minieigen))
#====================================
#Back compatibility with scons
SET (realVersion ${YADE_VERSION})
SET (version ${YADE_VERSION})
SET (pyExecutable ${PYTHON_EXECUTABLE})
SET (profile "default")
SET (sourceRoot "${CMAKE_CURRENT_SOURCE_DIR}")
#====================================

CONFIGURE_FILE(core/main/yade-batch.in "${CMAKE_BINARY_DIR}/bins/yade${SUFFIX}-batch")
CONFIGURE_FILE(core/main/main.py.in "${CMAKE_BINARY_DIR}/bins/yade${SUFFIX}")
CONFIGURE_FILE(py/config.py.in "${CMAKE_BINARY_DIR}/config.py")
CONFIGURE_FILE(py/__init__.py.in "${CMAKE_BINARY_DIR}/__init__.py")

#===========================================================

INSTALL(PROGRAMS "${CMAKE_BINARY_DIR}/bins/yade${SUFFIX}-batch" DESTINATION ${YADE_EXEC_PATH}/)
INSTALL(PROGRAMS "${CMAKE_BINARY_DIR}/bins/yade${SUFFIX}" DESTINATION ${YADE_EXEC_PATH}/)
INSTALL(FILES "${CMAKE_BINARY_DIR}/config.py" DESTINATION ${YADE_PY_PATH}/yade/)
INSTALL(FILES "${CMAKE_BINARY_DIR}/__init__.py" DESTINATION ${YADE_PY_PATH}/yade/)
FILE(GLOB filesPYChecks "${CMAKE_CURRENT_SOURCE_DIR}/scripts/checks-and-tests/checks/*.py")
INSTALL(FILES ${filesPYChecks} DESTINATION ${YADE_PY_PATH}/yade/tests/checks)
FILE(GLOB filesPYChecksData "${CMAKE_CURRENT_SOURCE_DIR}/scripts/checks-and-tests/checks/data/*")
INSTALL(FILES ${filesPYChecksData} DESTINATION ${YADE_PY_PATH}/yade/tests/checks/data)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/doc/yade-logo-note.png" DESTINATION "${YADE_DOC_PATH}/img")

INSTALL(TARGETS boot DESTINATION "${YADE_PY_PATH}/yade/")
INSTALL(TARGETS yade DESTINATION ${YADE_LIB_PATH})

#===========================================================
MESSAGE(STATUS "===========================================================")
MESSAGE(STATUS "Yade configured with following features:${CONFIGURED_FEATS}")
MESSAGE(STATUS "Disabled features:${DISABLED_FEATS}")
IF (DEBUG)
  MESSAGE(STATUS "Debug build")
  SET (debugbuild " (debug build)")
ELSE (DEBUG)
  MESSAGE(STATUS "Optimized build")
ENDIF (DEBUG)
IF (CHUNKSIZE)
  MESSAGE(STATUS "CHUNKSIZE is set to " ${CHUNKSIZE})
ENDIF (CHUNKSIZE)
MESSAGE(STATUS "===========================================================")
#===========================================================
#Building doc
ADD_CUSTOM_TARGET(doc)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/sphinx)
ADD_CUSTOM_COMMAND( 
                  TARGET doc PRE_BUILD
                  COMMAND rm -rf ${CMAKE_BINARY_DIR}/doc/sphinx/_build
                  COMMAND cp -r ${CMAKE_CURRENT_SOURCE_DIR}/doc/* ${CMAKE_BINARY_DIR}/doc
                  COMMAND PYTHONPATH=${CMAKE_BINARY_DIR}/doc/sphinx ${YADE_EXEC_PATH}/yade${SUFFIX} yadeSphinx.py
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc/sphinx
                  DEPENDS ${YADE_EXEC_PATH}/yade${SUFFIX}
                  )
ADD_CUSTOM_COMMAND( 
                  TARGET doc POST_BUILD
                  COMMAND xelatex Yade.tex
                  COMMAND xelatex Yade.tex
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc/sphinx/_build/latex
                  DEPENDS ${YADE_EXEC_PATH}/yade${SUFFIX} ${CMAKE_BINARY_DIR}/doc/sphinx/_build/latex/Yade.tex
                  )
ADD_CUSTOM_COMMAND( 
                  TARGET doc POST_BUILD
                  COMMAND rm -rf ${YADE_DOC_PATH}/html
                  COMMAND mv ${CMAKE_BINARY_DIR}/doc/sphinx/_build/html ${YADE_DOC_PATH}/html
                  COMMAND mv ${CMAKE_BINARY_DIR}/doc/sphinx/_build/latex/Yade.pdf ${YADE_DOC_PATH}/
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc/sphinx/_build/latex
                  DEPENDS ${YADE_EXEC_PATH}/yade${SUFFIX} ${CMAKE_BINARY_DIR}/doc/sphinx/_build/latex/Yade.tex
                  )
ADD_CUSTOM_COMMAND( 
                  TARGET doc POST_BUILD
                  COMMAND mv ${CMAKE_BINARY_DIR}/doc/sphinx/_build/epub/Yade.epub ${YADE_DOC_PATH}/ || true
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc/sphinx/_build/latex
                  DEPENDS ${YADE_EXEC_PATH}/yade${SUFFIX} ${CMAKE_BINARY_DIR}/doc/sphinx/_build/latex/Yade.tex
                  )
#===========================================================
#Building manpage
ADD_CUSTOM_TARGET(manpage)
ADD_CUSTOM_COMMAND( 
                  TARGET manpage POST_BUILD
                  COMMAND help2man ${YADE_EXEC_PATH}/yade${SUFFIX} > yade${SUFFIX}.1
                  COMMAND help2man ${YADE_EXEC_PATH}/yade${SUFFIX}-batch > yade${SUFFIX}-batch.1
                  COMMAND mkdir -p ${YADE_MAN_PATH}
                  COMMAND mv *.1 ${YADE_MAN_PATH}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                  DEPENDS ${YADE_EXEC_PATH}/yade${SUFFIX} /usr/bin/help2man
                  )
#===========================================================
#Execute standard checks
ADD_CUSTOM_TARGET(check)
ADD_CUSTOM_COMMAND( 
                  TARGET check POST_BUILD
                  COMMAND ${YADE_EXEC_PATH}/yade${SUFFIX} --test
                  COMMAND ${YADE_EXEC_PATH}/yade${SUFFIX} --check
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                  DEPENDS ${YADE_EXEC_PATH}/yade${SUFFIX}
                  )
#===========================================================
