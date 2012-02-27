	// Geometric Tools, Inc.
	// http://www.geometrictools.com
	// Copyright (c) 1998-2006.  All Rights Reserved
	//
	// The Wild Magic Library (WM3) source code is supplied under the terms of
	// the license agreement
	//     http://www.geometrictools.com/License/WildMagic3License.pdf
	// and may not be copied or disclosed except in accordance with the terms
	// of that agreement.

#ifndef WM3QUATERNION_H
#define WM3QUATERNION_H

#include "EigenCompat.h"
#include "Wm3FoundationLIB.h"
#include "Wm3Matrix3.h"
#include "Wm3Vector3.h"


	namespace Wm3
	{

	template<class Real>
	class AngleAxis{
		public:
			AngleAxis(){}
			AngleAxis(const Real& angle, const Vector3<Real> axis): _angle(angle),_axis(axis){}
			Real& angle(){ return _angle;}
			Real angle() const { return _angle;}
			Vector3<Real>& axis(){ return _axis;}
			Vector3<Real> axis() const {return _axis;}
		private:
			Real _angle;
			Vector3<Real> _axis;
	};

	template <class Real>
	class Quaternion
	{
	public:
		 // A quaternion is q = w + x*i + y*j + z*k where (w,x,y,z) is not
		 // necessarily a unit length vector in 4D.

		 // construction
		 Quaternion ();  // uninitialized
		 // this ordering is used by eigen as well
		 Quaternion (Real fW, Real fX, Real fY, Real fZ);
		 Quaternion (const Quaternion& rkQ);

		 // quaternion for the input rotation matrix
		 Quaternion (const Matrix3<Real>& rkRot);
		
		 EIG_FUN Quaternion(const AngleAxis<Real>& aa){ FromAxisAngle_(aa.axis(),aa.angle()); }

		 // quaternion for the rotation of the axis-angle pair
		 WM3_OLD Quaternion (const Vector3<Real>& rkAxis, Real fAngle);

		 // quaternion for the rotation matrix with specified columns
		 WM3_FUN Quaternion (const Vector3<Real> akRotColumn[3]);

		 // member access:  0 = w, 1 = x, 2 = y, 3 = z
		 operator const Real* () const;
		 operator Real* ();
		 // different ordering between eigen and wm3!!
		 __attribute__((warning("Quaternion_index_access: Wm3 and Eigen uses different storage order for quaternions, better use x(), y(), z(), w()!")))
		 	WM3_FUN Real operator[] (int i) const;
		 __attribute__((warning("Quaternion_index_access: Wm3 and Eigen uses different storage order for quaternions, better use x(), y(), z(), w()!")))
		 WM3_FUN Real& operator[] (int i);

		 EIG_FUN Real w () const;
		 EIG_FUN Real& w ();
		 EIG_FUN Real x () const;
		 EIG_FUN Real& x ();
		 EIG_FUN Real y () const;
		 EIG_FUN Real& y ();
		 EIG_FUN Real z () const;
		 EIG_FUN Real& z ();

		 WM3_OLD Real X() const{return x();}
		 WM3_OLD Real& X(){return x();}
		 WM3_OLD Real Y() const{return y();}
		 WM3_OLD Real& Y(){return y();}
		 WM3_OLD Real Z() const{return z();}
		 WM3_OLD Real& Z(){return z();}
		 WM3_OLD Real W() const{return w();}
		 WM3_OLD Real& W(){return w();}


		 // assignment
		 Quaternion& operator= (const Quaternion& rkQ);

    // comparison
    bool operator== (const Quaternion& rkQ) const;
    bool operator!= (const Quaternion& rkQ) const;
    WM3_OLD bool operator<  (const Quaternion& rkQ) const;
    WM3_OLD bool operator<= (const Quaternion& rkQ) const;
    WM3_OLD bool operator>  (const Quaternion& rkQ) const;
    WM3_OLD bool operator>= (const Quaternion& rkQ) const;

    // arithmetic operations
    Quaternion operator+ (const Quaternion& rkQ) const;
    Quaternion operator- (const Quaternion& rkQ) const;
    Quaternion operator* (const Quaternion& rkQ) const;
    Quaternion operator* (Real fScalar) const;
    Quaternion operator/ (Real fScalar) const;
    Quaternion operator- () const;

    // arithmetic updates
    Quaternion& operator+= (const Quaternion& rkQ);
    Quaternion& operator-= (const Quaternion& rkQ);
    Quaternion& operator*= (Real fScalar);
    Quaternion& operator/= (Real fScalar);

    // conversion between quaternions, matrices, and axis-angle
    WM3_FUN Quaternion& FromRotationMatrix (const Matrix3<Real>& rkRot) { return FromRotationMatrix_(rkRot); }
            Quaternion& FromRotationMatrix_ (const Matrix3<Real>& rkRot); // hidden function, do not call from outside
    EIG_FUN Matrix3<Real> toRotationMatrix () const;
    WM3_FUN Quaternion& FromRotationMatrix (const Vector3<Real> akRotColumn[3]);
    WM3_OLD void ToRotationMatrix (Vector3<Real> akRotColumn[3]) const;
            Quaternion& FromAxisAngle_ (const Vector3<Real>& rkAxis, Real fAngle); // hidden function, do not call from outside
    WM3_OLD Quaternion& FromAxisAngle (const Vector3<Real>& rkAxis, Real fAngle){return FromAxisAngle_(rkAxis,fAngle);}
            void ToAxisAngle_ (Vector3<Real>& rkAxis, Real& rfAngle) const; // hidden function, do not call from outside
    WM3_OLD void ToAxisAngle (Vector3<Real>& rkAxis, Real& rfAngle) const {return ToAxisAngle_(rkAxis,rfAngle);}

	 WM3_OLD Matrix3<Real> ToRotationMatrix(Matrix3<Real>& rkRot) const { Matrix3<Real> ret=toRotationMatrix(); return ret;}

    // functions of a quaternion
    WM3_OLD Real Dot (const Quaternion& rkQ) const;  // dot product of 4-tuples
    EIG_FUN Real normalize ();  // make the 4-tuple unit length
    EIG_FUN Real norm () const;
    EIG_FUN Real squaredNorm () const;
    EIG_FUN Quaternion conjugate () const;
    WM3_OLD Quaternion Inverse () const;  // apply to non-zero quaternion
    WM3_OLD Quaternion Exp () const;  // apply to quaternion with w = 0
    WM3_OLD Quaternion Log () const;  // apply to unit-length quaternion
	 WM3_OLD void Normalize() { (void)normalize(); } 
	 WM3_OLD Quaternion Conjugate() const { return conjugate(); } 


    // rotation of a vector by a quaternion
    EIG_FUN Vector3<Real> rotate (const Vector3<Real>& rkVector) const;
    WM3_OLD Vector3<Real> Rotate (const Vector3<Real>& v) const {return rotate(v);}

    // spherical linear interpolation
    WM3_OLD Quaternion& Slerp (Real fT, const Quaternion& rkP, const Quaternion& rkQ);

    WM3_OLD Quaternion& SlerpExtraSpins (Real fT, const Quaternion& rkP,
        const Quaternion& rkQ, int iExtraSpins);

    // intermediate terms for spherical quadratic interpolation
    WM3_OLD Quaternion& Intermediate (const Quaternion& rkQ0,
        const Quaternion& rkQ1, const Quaternion& rkQ2);

    // spherical quadratic interpolation
    WM3_OLD Quaternion& Squad (Real fT, const Quaternion& rkQ0,
        const Quaternion& rkA0, const Quaternion& rkA1,
        const Quaternion& rkQ1);

    // Compute a quaternion that rotates unit-length vector V1 to unit-length
    // vector V2.  The rotation is about the axis perpendicular to both V1 and
    // V2, with angle of that between V1 and V2.  If V1 and V2 are parallel,
    // any axis of rotation will do, such as the permutation (z2,x2,y2), where
    // V2 = (x2,y2,z2).
    EIG_FUN Quaternion& setFromTwoVectors (const Vector3<Real>& rkV1, const Vector3<Real>& rkV2);
	 WM3_OLD Quaternion& Align(const Vector3<Real>& a, const Vector3<Real>& b){ return setFromTwoVectors(a,b); }

    // Decompose a quaternion into q = q_twist * q_swing, where q is 'this'
    // quaternion.  If V1 is the input axis and V2 is the rotation of V1 by
    // q, q_swing represents the rotation about the axis perpendicular to
    // V1 and V2 (see Quaternion::Align), and q_twist is a rotation about V1.
    WM3_FUN void DecomposeTwistTimesSwing (const Vector3<Real>& rkV1,
        Quaternion& rkTwist, Quaternion& rkSwing);

    // Decompose a quaternion into q = q_swing * q_twist, where q is 'this'
    // quaternion.  If V1 is the input axis and V2 is the rotation of V1 by
    // q, q_swing represents the rotation about the axis perpendicular to
    // V1 and V2 (see Quaternion::Align), and q_twist is a rotation about V1.
    WM3_FUN void DecomposeSwingTimesTwist (const Vector3<Real>& rkV1,
        Quaternion& rkSwing, Quaternion& rkTwist);

    // special values
    WM3_FUN WM3_ITEM static const Quaternion IDENTITY;  // the identity rotation
    WM3_FUN WM3_ITEM static const Quaternion ZERO;
	 EIG_FUN static const Quaternion Identity(){ return IDENTITY; }

private:
    // support for comparisons
    int CompareArrays (const Quaternion& rkQ) const;

    // support for FromRotationMatrix
    static int ms_iNext[3];

    Real m_afTuple[4];
};

template <class Real>
Quaternion<Real> operator* (Real fScalar, const Quaternion<Real>& rkQ);

#include "Wm3Quaternion.inl"

typedef Quaternion<float> Quaternionf;
typedef Quaternion<double> Quaterniond;

}

#endif

