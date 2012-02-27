// © 2010 Václav Šmilauer <eudoxos@arcig.cz>
#pragma once

#if 0 // broken, do not use
// optimize as much as possible even in the debug mode (effective?)
#if defined(__GNUG__) && __GNUC__ >= 4 && __GNUC_MINOR__ >=4
	#pragma GCC push_options
	#pragma GCC optimize "2"
#endif
#endif

#ifdef QUAD_PRECISION
	typedef long double quad;
	typedef quad Real;
#else
	typedef double Real;
#endif

#include<limits>
#include<cstdlib>

/*
 * use Eigen http://eigen.tuxfamily.org
 */
#ifndef YADE_WM3
	// different macros for different versions of eigen:
	//  http://bitbucket.org/eigen/eigen/issue/96/eigen_dont_align-doesnt-exist-in-205-but-appears-in-web
	#define EIGEN_DONT_VECTORIZE
	#define EIGEN_DONT_ALIGN
	#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
	#define EIGEN_NO_DEBUG
		#include<Eigen/Core>
		#include<Eigen/Geometry>
		#include<Eigen/Array>
		#include<Eigen/QR>
		#include<Eigen/LU>
		#include<float.h>
	// USING_PART_OF_NAMESPACE_EIGEN
	//using namespace eigen; // for eigen3
	// 
	// templates of those types with single parameter are not possible (for compat with Wm3), use macros for now
	#define VECTOR2_TEMPLATE(Scalar) Eigen::Matrix<Scalar,2,1>
	#define VECTOR3_TEMPLATE(Scalar) Eigen::Matrix<Scalar,3,1>
	#define MATRIX3_TEMPLATE(Scalar) Eigen::Matrix<Scalar,3,3>
	// this would be the proper way, but only works in c++-0x (not yet supported by gcc (4.5))
	#if 0
		template<typename Scalar> using Vector2=Eigen::Matrix<Scalar,2,1>;
		template<typename Scalar> using Vector3=Eigen::Matrix<Scalar,3,1>;
		template<typename Scalar> using Matrix3=Eigen::Matrix<Scalar,3,3>;
		typedef Vector2<int> Vector2i;
		typedef Vector2<Real> Vector2r;
		// etc
	#endif

	typedef VECTOR2_TEMPLATE(int) Vector2i;
	typedef VECTOR2_TEMPLATE(Real) Vector2r;
	typedef VECTOR3_TEMPLATE(int) Vector3i;
	typedef VECTOR3_TEMPLATE(Real) Vector3r;
	typedef MATRIX3_TEMPLATE(Real) Matrix3r;

	typedef Eigen::Quaternion<Real> Quaternionr;
	typedef Eigen::AngleAxis<Real> AngleAxisr;
	using Eigen::AngleAxis; using Eigen::Quaternion;

	// io
	template<class Scalar> std::ostream & operator<<(std::ostream &os, const VECTOR2_TEMPLATE(Scalar)& v){ os << v.x()<<" "<<v.y(); return os; };
	template<class Scalar> std::ostream & operator<<(std::ostream &os, const VECTOR3_TEMPLATE(Scalar)& v){ os << v.x()<<" "<<v.y()<<" "<<v.z(); return os; };
	template<class Scalar> std::ostream & operator<<(std::ostream &os, const Eigen::Quaternion<Scalar>& q){ os<<q.w()<<" "<<q.x()<<" "<<q.y()<<" "<<q.z(); return os; };
	// operators
	//template<class Scalar> VECTOR3_TEMPLATE(Scalar) operator*(Scalar s, const VECTOR3_TEMPLATE(Scalar)& v) {return v*s;}
	//template<class Scalar> MATRIX3_TEMPLATE(Scalar) operator*(Scalar s, const MATRIX3_TEMPLATE(Scalar)& m) { return m*s; }
	//template<class Scalar> Quaternion<Scalar> operator*(Scalar s, const Quaternion<Scalar>& q) { return q*s; }
	template<typename Scalar> void matrixEigenDecomposition(const MATRIX3_TEMPLATE(Scalar) m, MATRIX3_TEMPLATE(Scalar)& mRot, MATRIX3_TEMPLATE(Scalar)& mDiag){ Eigen::SelfAdjointEigenSolver<MATRIX3_TEMPLATE(Scalar)> a(m); mRot=a.eigenvectors(); mDiag=a.eigenvalues().asDiagonal(); }
	//__attribute__((warning("Replace this function with direct AngleAxis constructor from Quaternion")))
	template<typename Scalar> AngleAxis<Scalar> angleAxisFromQuat(const Quaternion<Scalar>& q){ return AngleAxis<Scalar>(q); }
	// http://eigen.tuxfamily.org/dox/TutorialGeometry.html
	template<typename Scalar> MATRIX3_TEMPLATE(Scalar) matrixFromEulerAnglesXYZ(Scalar x, Scalar y, Scalar z){ MATRIX3_TEMPLATE(Scalar) m; m=AngleAxis<Scalar>(x,VECTOR3_TEMPLATE(Scalar)::UnitX())*AngleAxis<Scalar>(y,VECTOR3_TEMPLATE(Scalar)::UnitY())*AngleAxis<Scalar>(z,VECTOR3_TEMPLATE(Scalar)::UnitZ()); return m;}
	template<typename Scalar> bool operator==(const Quaternion<Scalar>& u, const Quaternion<Scalar>& v){ return u.x()==v.x() && u.y()==v.y() && u.z()==v.z() && u.w()==v.w(); }
	template<typename Scalar> bool operator!=(const Quaternion<Scalar>& u, const Quaternion<Scalar>& v){ return !(u==v); }
	template<typename Scalar> bool operator==(const MATRIX3_TEMPLATE(Scalar)& m, const MATRIX3_TEMPLATE(Scalar)& n){ for(int i=0;i<3;i++)for(int j=0;j<3;j++)if(m(i,j)!=n(i,j)) return false; return true; }
	template<typename Scalar> bool operator!=(const MATRIX3_TEMPLATE(Scalar)& m, const MATRIX3_TEMPLATE(Scalar)& n){ return !(m==n); }
	template<typename Scalar> bool operator==(const VECTOR3_TEMPLATE(Scalar)& u, const VECTOR3_TEMPLATE(Scalar)& v){ return u.x()==v.x() && u.y()==v.y() && u.z()==v.z(); }
	template<typename Scalar> bool operator!=(const VECTOR3_TEMPLATE(Scalar)& u, const VECTOR3_TEMPLATE(Scalar)& v){ return !(u==v); }
	template<typename Scalar> bool operator==(const VECTOR2_TEMPLATE(Scalar)& u, const VECTOR2_TEMPLATE(Scalar)& v){ return u.x()==v.x() && u.y()==v.y(); }
	template<typename Scalar> bool operator!=(const VECTOR2_TEMPLATE(Scalar)& u, const VECTOR2_TEMPLATE(Scalar)& v){ return !(u==v); }
	template<typename Scalar> Quaternion<Scalar> operator*(Scalar f, const Quaternion<Scalar>& q){ return Quaternion<Scalar>(q.coeffs()*f); }
	template<typename Scalar> Quaternion<Scalar> operator+(Quaternion<Scalar> q1, const Quaternion<Scalar>& q2){ return Quaternion<Scalar>(q1.coeffs()+q2.coeffs()); }	/* replace all those by standard math functions
		this is a non-templated version, to avoid compilation because of static constants;
	*/
	template<typename Scalar>
	struct Math{
		static const Scalar PI;
		static const Scalar HALF_PI;
		static const Scalar TWO_PI;
		static const Scalar MAX_REAL;
		static const Scalar DEG_TO_RAD;
		static const Scalar RAD_TO_DEG;
		static const Scalar EPSILON;
		static const Scalar ZERO_TOLERANCE;
		static Scalar Sign(Scalar f){ if(f<0) return -1; if(f>0) return 1; return 0; }
		static Scalar FAbs(Scalar f){ return abs(f); }
		static Scalar Sqrt(Scalar f){ return sqrt(f); }
		static Scalar Log(Scalar f){ return log(f); }
		static Scalar Exp(Scalar f){ return exp(f); }
		static Scalar ATan(Scalar f){ return atan(f); }
		static Scalar Tan(Scalar f){ return tan(f); }
		static Scalar Pow(Scalar base,Scalar exponent){ return pow(base,exponent); }

		static Scalar UnitRandom(){ return ((double)rand()/((double)(RAND_MAX))); }
		static Scalar SymmetricRandom(){ return 2.*(((double)rand())/((double)(RAND_MAX)))-1.; }
		static Scalar FastInvCos0(Scalar fValue){ Scalar fRoot = sqrt(((Scalar)1.0)-fValue); Scalar fResult = -(Scalar)0.0187293; fResult *= fValue; fResult += (Scalar)0.0742610; fResult *= fValue; fResult -= (Scalar)0.2121144; fResult *= fValue; fResult += (Scalar)1.5707288; fResult *= fRoot; return fResult; }
	};
	typedef Math<Real> Mathr;
#else
	// wm3 fallback

	#include<Wm3Vector2.h>
	#include<Wm3Vector3.h>
	#include<Wm3Matrix3.h>
	#include<Wm3Quaternion.h>
	#include<Wm3Math.h>
#if 0
	namespace Wm3 {
		template<class T> class Math; 		template<class T> class Quaternion; 
		template<class T> class Matrix3; //typedef Matrix3<Real> Matrix3r;
		template<class T> class Vector2; //typedef Vector2<Real> Vector2r; typedef Vector2<int> Vector2i;
		template<class T> class Vector3; //typedef Vector3<Real> Vector3r; typedef Vector3<int> Vector3i;
		template<class T> class AngleAxis;  // added to the Wm3Quaternion header
	}
#endif
	using Wm3::Quaternion;
	using Wm3::AngleAxis;
	typedef Wm3::Math<Real> Mathr;
	typedef Quaternion<Real> Quaternionr;
	typedef AngleAxis<Real> AngleAxisr;
	// workaround for template aliases, to be compatible with Eigen
	#define VECTOR2_TEMPLATE(Scalar) Wm3::Vector2<Scalar>
	#define VECTOR3_TEMPLATE(Scalar) Wm3::Vector3<Scalar>
	#define MATRIX3_TEMPLATE(Scalar) Wm3::Matrix3<Scalar>
	// typedefs
	typedef VECTOR2_TEMPLATE(int) Vector2i;
	typedef VECTOR2_TEMPLATE(Real) Vector2r;
	typedef VECTOR3_TEMPLATE(int) Vector3i;
	typedef VECTOR3_TEMPLATE(Real) Vector3r;
	typedef MATRIX3_TEMPLATE(Real) Matrix3r;
	//
	template<typename Scalar> AngleAxis<Scalar> angleAxisFromQuat(const Quaternion<Scalar>& q){ AngleAxis<Scalar> aa; q.ToAxisAngle_(aa.axis(),aa.angle()); return aa; }
	template<typename Scalar> void matrixEigenDecomposition(const MATRIX3_TEMPLATE(Scalar)& m, MATRIX3_TEMPLATE(Scalar)& mRot, MATRIX3_TEMPLATE(Scalar)& mDiag){ m.EigenDecomposition_(mRot,mDiag); }
	template<typename Scalar> MATRIX3_TEMPLATE(Scalar) matrixFromEulerAnglesXYZ(Scalar x, Scalar y, Scalar z){ return MATRIX3_TEMPLATE(Scalar)().FromEulerAnglesXYZ_(x,y,z); }
	// only needed for Wm3, Eigen defines already
	template<class Scalar> VECTOR3_TEMPLATE(Scalar) operator*(const Quaternion<Scalar>& q, const VECTOR3_TEMPLATE(Scalar)& v){ return q.rotate(v); /* check whether it converts to matrix internally, which is supposedly faster*/ }
#endif

/*
 * Compatibility bridge for Wm3 and eigen (will be removed once Wm3 is dropped and replaced by respective eigen constructs;
 * see https://www.yade-dem.org/wiki/Wm3→Eigen
 */

// replace by outer product of 2 vectors: v1*v2.transpose();
template<typename Scalar>
MATRIX3_TEMPLATE(Scalar) makeTensorProduct (const VECTOR3_TEMPLATE(Scalar)& rkU,  const VECTOR3_TEMPLATE(Scalar)& rkV)
{
	MATRIX3_TEMPLATE(Scalar) ret;
   ret(0,0) = rkU[0]*rkV[0];
   ret(0,1) = rkU[0]*rkV[1];
   ret(0,2) = rkU[0]*rkV[2];
   ret(1,0) = rkU[1]*rkV[0];
   ret(1,1) = rkU[1]*rkV[1];
   ret(1,2) = rkU[1]*rkV[2];
   ret(2,0) = rkU[2]*rkV[0];
   ret(2,1) = rkU[2]*rkV[1];
   ret(2,2) = rkU[2]*rkV[2];
   return ret;
}
// eigen2: a.cwise().min(b) / max(b)
// eigen3: a.array().min(b) / max(b)
template<typename Scalar> VECTOR2_TEMPLATE(Scalar) componentMaxVector(const VECTOR2_TEMPLATE(Scalar)& a, const VECTOR2_TEMPLATE(Scalar)& b){ return VECTOR2_TEMPLATE(Scalar)(std::max(a.x(),b.x()),std::max(a.y(),b.y()));}
template<typename Scalar> VECTOR2_TEMPLATE(Scalar) componentMinVector(const VECTOR2_TEMPLATE(Scalar)& a, const VECTOR2_TEMPLATE(Scalar)& b){ return VECTOR2_TEMPLATE(Scalar)(std::min(a.x(),b.x()),std::min(a.y(),b.y()));}
template<typename Scalar> VECTOR3_TEMPLATE(Scalar) componentMaxVector(const VECTOR3_TEMPLATE(Scalar)& a, const VECTOR3_TEMPLATE(Scalar)& b){ return VECTOR3_TEMPLATE(Scalar)(std::max(a.x(),b.x()),std::max(a.y(),b.y()),std::max(a.z(),b.z()));}
template<typename Scalar> VECTOR3_TEMPLATE(Scalar) componentMinVector(const VECTOR3_TEMPLATE(Scalar)& a, const VECTOR3_TEMPLATE(Scalar)& b){ return VECTOR3_TEMPLATE(Scalar)(std::min(a.x(),b.x()),std::min(a.y(),b.y()),std::min(a.z(),b.z()));}
// eigen2: v1.cwise()*v2;
// eigen3: v1.array()*v2.array()
template<typename Scalar> VECTOR3_TEMPLATE(Scalar) diagMult(const VECTOR3_TEMPLATE(Scalar)& a, const VECTOR3_TEMPLATE(Scalar)& b){return VECTOR3_TEMPLATE(Scalar)(a.x()*b.x(),a.y()*b.y(),a.z()*b.z());}
template<typename Scalar> VECTOR3_TEMPLATE(Scalar) diagDiv(const VECTOR3_TEMPLATE(Scalar)& a, const VECTOR3_TEMPLATE(Scalar)& b){return VECTOR3_TEMPLATE(Scalar)(a.x()/b.x(),a.y()/b.y(),a.z()/b.z());}
// eigen: m << m00,m01,m02,m10,m11,m12,m20,m21,m22;
template<typename Scalar> MATRIX3_TEMPLATE(Scalar) matrixFromElements(Scalar m00, Scalar m01, Scalar m02, Scalar m10, Scalar m11, Scalar m12, Scalar m20, Scalar m21, Scalar m22){ MATRIX3_TEMPLATE(Scalar) m; m(0,0)=m00; m(0,1)=m01; m(0,2)=m02; m(1,0)=m10; m(1,1)=m11; m(1,2)=m12; m(2,0)=m20; m(2,1)=m21; m(2,2)=m22; return m; }


/*
 * Extra yade math functions and classes
 */

__attribute__((unused))
const Real NaN(std::numeric_limits<Real>::signaling_NaN());

// void quaternionToEulerAngles (const Quaternionr& q, Vector3r& eulerAngles,Real threshold=1e-6f);
template<typename Scalar> void quaterniontoGLMatrix(const Quaternion<Scalar>& q, Scalar m[16]){
	Scalar w2=2.*q.w(), x2=2.*q.x(), y2=2.*q.y(), z2=2.*q.z();
	Scalar x2w=w2*q.w(), y2w=y2*q.w(), z2w=z2*q.w();
	Scalar x2x=x2*q.x(), y2x=y2*q.x(), z2x=z2*q.x();
	Scalar x2y=y2*q.y(), y2y=y2*q.y(), z2y=z2*q.y();
	Scalar x2z=z2*q.z(), y2z=y2*q.z(), z2z=z2*q.z();
	m[0]=1.-(y2y+z2z); m[4]=y2x-z2w;      m[8]=z2x+y2w;       m[12]=0;
	m[1]=y2x+z2w;      m[5]=1.-(x2x+z2z); m[9]=z2y-x2w;       m[13]=0;
	m[2]=z2x-y2w;      m[6]=z2y+x2w;      m[10]=1.-(x2x+y2y); m[14]=0;
	m[3]=0.;           m[7]=0.;           m[11]=0.;           m[15]=1.;
}



// se3
template <class Scalar>
class Se3
{
	public :
		VECTOR3_TEMPLATE(Scalar)	position;
		Quaternion<Scalar>	orientation;
		Se3(){};
		Se3(VECTOR3_TEMPLATE(Scalar) rkP, Quaternion<Scalar> qR){ position = rkP; orientation = qR; }
		Se3(const Se3<Scalar>& s){position = s.position;orientation = s.orientation;}
		Se3(Se3<Scalar>& a,Se3<Scalar>& b){
			position  = b.orientation.inverse()*(a.position - b.position);  
			orientation = b.orientation.inverse()*a.orientation;
		}
		Se3<Scalar> inverse(){ return Se3(-(orientation.inverse()*position), orientation.inverse());}
		void toGLMatrix(float m[16]){ orientation.toGLMatrix(m); m[12] = position[0]; m[13] = position[1]; m[14] = position[2];}
		VECTOR3_TEMPLATE(Scalar) operator * (const VECTOR3_TEMPLATE(Scalar)& b ){return orientation*b+position;}
		Se3<Scalar> operator * (const Quaternion<Scalar>& b ){return Se3<Scalar>(position , orientation*b);}
		Se3<Scalar> operator * (const Se3<Scalar>& b ){return Se3<Scalar>(orientation*b.position+position,orientation*b.orientation);}
};

// functions
template<typename Scalar> Scalar unitVectorsAngle(const VECTOR3_TEMPLATE(Scalar)& a, const VECTOR3_TEMPLATE(Scalar)& b){ return acos(a.dot(b)); }
// operators

/*
 * typedefs
 */
typedef Se3<Real> Se3r;


/*
 * Serialization of math classes
 */


// gccxml chokes on the boost::serialization code; this part is not needed if wrapping miniWm3 anyway
#if !defined(__GCCXML__) and defined(YADE_BOOST_SERIALIZATION)

#include<boost/serialization/nvp.hpp>
#include<boost/serialization/is_bitwise_serializable.hpp>

// fast serialization (no version infor and no tracking) for basic math types
// http://www.boost.org/doc/libs/1_42_0/libs/serialization/doc/traits.html#bitwise
BOOST_IS_BITWISE_SERIALIZABLE(Vector2r);
BOOST_IS_BITWISE_SERIALIZABLE(Vector2i);
BOOST_IS_BITWISE_SERIALIZABLE(Vector3r);
BOOST_IS_BITWISE_SERIALIZABLE(Vector3i);
BOOST_IS_BITWISE_SERIALIZABLE(Quaternionr);
BOOST_IS_BITWISE_SERIALIZABLE(Se3r);
BOOST_IS_BITWISE_SERIALIZABLE(Matrix3r);

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, Vector2r & g, const unsigned int version){
	Real &x=g[0], &y=g[1];
	ar & BOOST_SERIALIZATION_NVP(x) & BOOST_SERIALIZATION_NVP(y);
}

template<class Archive>
void serialize(Archive & ar, Vector2i & g, const unsigned int version){
	int &x=g[0], &y=g[1];
	ar & BOOST_SERIALIZATION_NVP(x) & BOOST_SERIALIZATION_NVP(y);
}

template<class Archive>
void serialize(Archive & ar, Vector3r & g, const unsigned int version)
{
	Real &x=g[0], &y=g[1], &z=g[2];
	ar & BOOST_SERIALIZATION_NVP(x) & BOOST_SERIALIZATION_NVP(y) & BOOST_SERIALIZATION_NVP(z);
}

template<class Archive>
void serialize(Archive & ar, Vector3i & g, const unsigned int version){
	int &x=g[0], &y=g[1], &z=g[2];
	ar & BOOST_SERIALIZATION_NVP(x) & BOOST_SERIALIZATION_NVP(y) & BOOST_SERIALIZATION_NVP(z);
}

template<class Archive>
void serialize(Archive & ar, Quaternionr & g, const unsigned int version)
{
	Real &w=g.w(), &x=g.x(), &y=g.y(), &z=g.z();
	ar & BOOST_SERIALIZATION_NVP(w) & BOOST_SERIALIZATION_NVP(x) & BOOST_SERIALIZATION_NVP(y) & BOOST_SERIALIZATION_NVP(z);
}

template<class Archive>
void serialize(Archive & ar, Se3r & g, const unsigned int version){
	Vector3r& position=g.position; Quaternionr& orientation=g.orientation;
	ar & BOOST_SERIALIZATION_NVP(position) & BOOST_SERIALIZATION_NVP(orientation);
}

template<class Archive>
void serialize(Archive & ar, Matrix3r & m, const unsigned int version){
	Real &m00=m(0,0), &m01=m(0,1), &m02=m(0,2), &m10=m(1,0), &m11=m(1,1), &m12=m(1,2), &m20=m(2,0), &m21=m(2,1), &m22=m(2,2);
	ar & BOOST_SERIALIZATION_NVP(m00) & BOOST_SERIALIZATION_NVP(m01) & BOOST_SERIALIZATION_NVP(m02) &
		BOOST_SERIALIZATION_NVP(m10) & BOOST_SERIALIZATION_NVP(m11) & BOOST_SERIALIZATION_NVP(m12) &
		BOOST_SERIALIZATION_NVP(m20) & BOOST_SERIALIZATION_NVP(m21) & BOOST_SERIALIZATION_NVP(m22);
}

} // namespace serialization
} // namespace boost

#endif /* __GCCXML */

#if 0
// revert optimization options back
#if defined(__GNUG__) && __GNUC__ >= 4 && __GNUC_MINOR__ >=4
	#pragma GCC pop_options
#endif
#endif
