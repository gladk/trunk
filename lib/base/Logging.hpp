// 2006-2008 © Václav Šmilauer <eudoxos@arcig.cz> 
#pragma once
/*
 * This file defines various useful logging-related macros - userspace stuff is
 * - LOG_* for actual logging,
 * - DECLARE_LOGGER; that should be used in class header to create separate logger for that class,
 * - CREATE_LOGGER(className); that must be used in class implementation file to create the static variable.
 *
 * Note that the latter 2 may change their name to something like LOG_DECLARE and LOG_CREATE, to be consistent.
 * Some other macros will be very likely added, to allow for easy variable tracing etc. Suggestions welcome.
 *
 * All of user macros should come in 2 flavors, depending on whether we use log4cxx or not (backward compatibility).
 * The default is not to use it, unless the preprocessor macro YADE_LOG4CXX is defined. In that case, you want to #include
 * yade-core/logging.h and link with log4cxx.
 *
 * TODO:
 * 1. [done] for optimized builds, at least debugging macros should become no-ops
 * 2. add the TRACE level; this has to be asked about on log4cxx mailing list perhaps. Other levels may follow easily. [will be done with log4cxx 0.10, once debian ships packages and we can safely migrate]
 *
 * For more information, see http://logging.apache.org/log4cxx/, especially the part on configuration files, that allow
 * very flexibe, runtime and fine-grained output redirections, filtering etc.
 *
 * Yade has the logging config file by default in ~/.yade-$VERSION/logging.conf.
 *
 */

#ifdef YADE_LOG4CXX

#	include<log4cxx/logger.h>
#	include<log4cxx/basicconfigurator.h>
#	include<log4cxx/propertyconfigurator.h>
#	include<log4cxx/helpers/exception.h>

#	define _LOG_HEAD __FILE__":"<<__LINE__<<" "<<__FUNCTION__<<": "
	// logger is local for every class, but if it is missing, we will use the parent's class logger automagically.
// TRACE doesn't really exist in log4cxx 0.9.7 (does in 0.10), otput through DEBUG then
#ifdef NDEBUG
#	define LOG_TRACE(msg){}
#	define LOG_DEBUG(msg){}
#else
#	ifdef LOG4CXX_TRACE
#		define LOG_TRACE(msg) {LOG4CXX_TRACE(logger, _LOG_HEAD<<msg);}
#	else
#		define LOG_TRACE(msg) {LOG4CXX_DEBUG(logger, _LOG_HEAD<<msg);}
#	endif
#	define LOG_DEBUG(msg) {LOG4CXX_DEBUG(logger, _LOG_HEAD<<msg);}
#endif
#	define LOG_INFO(msg)  {LOG4CXX_INFO(logger,  _LOG_HEAD<<msg);}
#	define LOG_WARN(msg)  {LOG4CXX_WARN(logger,  _LOG_HEAD<<msg);}
#	define LOG_ERROR(msg) {LOG4CXX_ERROR(logger, _LOG_HEAD<<msg);}
#	define LOG_FATAL(msg) {LOG4CXX_FATAL(logger, _LOG_HEAD<<msg);}

#	define DECLARE_LOGGER public: static log4cxx::LoggerPtr logger
#	define CREATE_LOGGER(classname) log4cxx::LoggerPtr classname::logger = log4cxx::Logger::getLogger("yade." #classname)

#else

#include<iostream>

#	define _POOR_MANS_LOG(level,msg) {std::cerr<<level " "<<_LOG_HEAD<<msg<<std::endl;}
#	define _LOG_HEAD __FILE__ ":"<<__LINE__<<" "<<__FUNCTION__<<": "
#	define LOG_TRACE(msg) // _POOR_MANS_LOG("TRACE",msg)
#	define LOG_DEBUG(msg) // _POOR_MANS_LOG("DEBUG",msg)
#	define LOG_INFO(msg) _POOR_MANS_LOG("INFO ",msg)
#	define LOG_WARN(msg) _POOR_MANS_LOG("WARN ",msg)
#	define LOG_ERROR(msg) _POOR_MANS_LOG("ERROR",msg)
#	define LOG_FATAL(msg) _POOR_MANS_LOG("FATAL",msg)

#	define DECLARE_LOGGER
#	define CREATE_LOGGER(classname)

#endif

// these macros are temporary
#define TRACE LOG_TRACE("Been here")
#define _TRVHEAD cerr<<__FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<": "
#define _TRV(x) #x"="<<x<<"; "
#define _TRVTAIL "\n"
#define TRVAR1(a) LOG_TRACE( _TRV(a) );
#define TRVAR2(a,b) LOG_TRACE( _TRV(a) << _TRV(b) );
#define TRVAR3(a,b,c) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) );
#define TRVAR4(a,b,c,d) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) << _TRV(d) );
#define TRVAR5(a,b,c,d,e) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) << _TRV(d) << _TRV(e) );
#define TRVAR6(a,b,c,d,e,f) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) << _TRV(d) << _TRV(e) << _TRV(f));
// prints boost matrix
#define TRMAT(MAT) _TRVHEAD<< #MAT "=(";for(unsigned i=0; i<MAT.size1(); i++){cerr<<"(";for(unsigned j=0; j<MAT.size2(); j++){ cerr<<MAT(i,j)<<" "; } cerr<<")"; } cerr<<")"; cerr<<_TRVTAIL
// dtto, but for matrix of vectors; maybe the previos macro could handle that also.
#define TRMATVEC(MAT) _TRVHEAD<< #MAT "=(";for(unsigned i=0; i<MAT.size1(); i++){cerr<<"(";for(unsigned j=0; j<MAT.size2(); j++){ cerr<<"["<<static_cast<Vector3r>(MAT(i,j))<<"]"; } cerr<<")"; } cerr<<")"; cerr<<_TRVTAIL
// show Matrix3 from the wm3 library 
#define TRWM3MAT(_M)	LOG_TRACE(#_M "=(("<<_M(0,0)<<" "<<_M(0,1)<<" "<<_M(0,2)<<")("<<_M(1,0)<<" "<<_M(1,1)<<" "<<_M(1,2)<<")("<<_M(2,0)<<" "<<_M(2,1)<<" "<<_M(2,2)<<"))");
#define TRWM3VEC(_V) LOG_TRACE(#_V "=("<<_V[0]<<" "<<_V[1]<<" "<<_V[2]<<")")
#define TRWM3QUAT(_Q) LOG_TRACE(#_Q "=("<<_Q[0]<<" "<<_Q[1]<<" "<<_Q[2]<<" "<<_Q[3]<<")")

