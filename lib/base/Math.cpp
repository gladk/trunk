#include<yade/lib-base/Math.hpp>
#ifndef YADE_WM3
	template<> const Real Math<Real>::EPSILON = DBL_EPSILON;
	template<> const Real Math<Real>::ZERO_TOLERANCE = 1e-20;
	template<> const Real Math<Real>::MAX_REAL = DBL_MAX;
	template<> const Real Math<Real>::PI = 4.0*atan(1.0);
	template<> const Real Math<Real>::TWO_PI = 2.0*Math<Real>::PI;
	template<> const Real Math<Real>::HALF_PI = 0.5*Math<Real>::PI;
	template<> const Real Math<Real>::DEG_TO_RAD = Math<Real>::PI/180.0;
	template<> const Real Math<Real>::RAD_TO_DEG = 180.0/Math<Real>::PI;
#endif
