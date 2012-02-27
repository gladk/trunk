// Geometric Tools, Inc.
// http://www.geometrictools.com
// Copyright (c) 1998-2006.  All Rights Reserved
//
// The Wild Magic Library (WM3) source code is supplied under the terms of
// the license agreement
//     http://www.geometrictools.com/License/WildMagic3License.pdf
// and may not be copied or disclosed except in accordance with the terms
// of that agreement.

#ifndef WM3MATRIX3_H
#define WM3MATRIX3_H

// Matrix operations are applied on the left.  For example, given a matrix M
// and a vector V, matrix-times-vector is M*V.  That is, V is treated as a
// column vector.  Some graphics APIs use V*M where V is treated as a row
// vector.  In this context the "M" matrix is really a transpose of the M as
// represented in Wild Magic.  Similarly, to apply two matrix operations M0
// and M1, in that order, you compute M1*M0 so that the transform of a vector
// is (M1*M0)*V = M1*(M0*V).  Some graphics APIs use M0*M1, but again these
// matrices are the transpose of those as represented in Wild Magic.  You
// must therefore be careful about how you interface the transformation code
// with graphics APIS.
//
// For memory organization it might seem natural to chose Real[N][N] for the
// matrix storage, but this can be a problem on a platform/console that
// chooses to store the data in column-major rather than row-major format.
// To avoid potential portability problems, the matrix is stored as Real[N*N]
// and organized in row-major order.  That is, the entry of the matrix in row
// r (0 <= r < N) and column c (0 <= c < N) is stored at index i = c+N*r
// (0 <= i < N*N).

// The (x,y,z) coordinate system is assumed to be right-handed.  Coordinate
// axis rotation matrices are of the form
//   RX =    1       0       0
//           0     cos(t) -sin(t)
//           0     sin(t)  cos(t)
// where t > 0 indicates a counterclockwise rotation in the yz-plane
//   RY =  cos(t)    0     sin(t)
//           0       1       0
//        -sin(t)    0     cos(t)
// where t > 0 indicates a counterclockwise rotation in the zx-plane
//   RZ =  cos(t) -sin(t)    0
//         sin(t)  cos(t)    0
//           0       0       1
// where t > 0 indicates a counterclockwise rotation in the xy-plane.

#include "EigenCompat.h"
#include "Wm3FoundationLIB.h"
#include "Wm3Vector3.h"

namespace Wm3
{

template <class Real>
class Matrix3
{
public:
    // If bZero is true, create the zero matrix.  Otherwise, create the
    // identity matrix.
    Matrix3 (){ setIdentity(); }
	 WM3_OLD Matrix3 (bool setZero);

    // copy constructor
    Matrix3 (const Matrix3& rkM);

    // input Mrc is in row r, column c.
    Matrix3 (Real fM00, Real fM01, Real fM02,
             Real fM10, Real fM11, Real fM12,
             Real fM20, Real fM21, Real fM22);

    // Create a matrix from an array of numbers.  The input array is
    // interpreted based on the Boolean input as
    //   true:  entry[0..8]={m00,m01,m02,m10,m11,m12,m20,m21,m22} [row major]
    //   false: entry[0..8]={m00,m10,m20,m01,m11,m21,m02,m12,m22} [col major]
    WM3_FUN Matrix3 (const Real afEntry[9], bool bRowMajor);

    // Create matrices based on vector input.  The Boolean is interpreted as
    //   true: vectors are columns of the matrix
    //   false: vectors are rows of the matrix
    WM3_FUN Matrix3 (const Vector3<Real>& rkU, const Vector3<Real>& rkV,
        const Vector3<Real>& rkW, bool bColumns);
    WM3_FUN Matrix3 (const Vector3<Real>* akV, bool bColumns);

    // create a diagonal matrix
    WM3_FUN Matrix3 (Real fM00, Real fM11, Real fM22);

    // Create rotation matrices (positive angle - counterclockwise).  The
    // angle must be in radians, not degrees.
    WM3_FUN Matrix3 (const Vector3<Real>& rkAxis, Real fAngle);

    // create a tensor product U*V^T
    WM3_FUN Matrix3 (const Vector3<Real>& rkU, const Vector3<Real>& rkV);

    // create various matrices
    EIG_FUN Matrix3& setZero ();
    EIG_FUN Matrix3& setIdentity ();
    WM3_OLD Matrix3& MakeDiagonal (Real fM00, Real fM11, Real fM22) { return MakeDiagonal_(fM00,fM11,fM22); }
            Matrix3& MakeDiagonal_ (Real fM00, Real fM11, Real fM22); // hidden function, do not use from outside
    WM3_FUN Matrix3& FromAxisAngle (const Vector3<Real>& rkAxis, Real fAngle);
    WM3_FUN Matrix3& MakeTensorProduct (const Vector3<Real>& rkU,
        const Vector3<Real>& rkV);

    // member access
    WM3_FUN operator const Real* () const;
    WM3_FUN operator Real* ();
    WM3_FUN const Real* operator[] (int iRow) const;
    WM3_FUN Real* operator[] (int iRow);
    Real operator() (int iRow, int iCol) const;
    Real& operator() (int iRow, int iCol);
    WM3_FUN void SetRow (int iRow, const Vector3<Real>& rkV);
    WM3_FUN Vector3<Real> GetRow (int iRow) const;
    WM3_FUN void SetColumn (int iCol, const Vector3<Real>& rkV);
    EIG_FUN Vector3<Real> col (int iCol) const;
	 WM3_OLD Vector3<Real> GetColumn(int i) const { return col(i); }
    WM3_FUN void GetColumnMajor (Real* afCMajor) const;

    // assignment
    Matrix3& operator= (const Matrix3& rkM);

    // comparison
    bool operator== (const Matrix3& rkM) const;
    bool operator!= (const Matrix3& rkM) const;
    WM3_FUN bool operator<  (const Matrix3& rkM) const;
    WM3_FUN bool operator<= (const Matrix3& rkM) const;
    WM3_FUN bool operator>  (const Matrix3& rkM) const;
    WM3_FUN bool operator>= (const Matrix3& rkM) const;

    // arithmetic operations
    Matrix3 operator+ (const Matrix3& rkM) const;
    Matrix3 operator- (const Matrix3& rkM) const;
    Matrix3 operator* (const Matrix3& rkM) const;
    Matrix3 operator* (Real fScalar) const;
    Matrix3 operator/ (Real fScalar) const;
    Matrix3 operator- () const;

    // arithmetic updates
    Matrix3& operator+= (const Matrix3& rkM);
    Matrix3& operator-= (const Matrix3& rkM);
    Matrix3& operator*= (Real fScalar);
    Matrix3& operator/= (Real fScalar);

    // matrix times vector
    Vector3<Real> operator* (const Vector3<Real>& rkV) const;  // M * v

    // other operations
    EIG_FUN Matrix3 transpose () const;  // M^T
    WM3_FUN Matrix3 TransposeTimes (const Matrix3& rkM) const;  // this^T * M
    WM3_FUN Matrix3 TimesTranspose (const Matrix3& rkM) const;  // this * M^T
    EIG_FUN Matrix3 inverse () const;
    EIG_FUN Matrix3 adjoint () const;
    EIG_FUN Real determinant () const;
    WM3_FUN Real QForm (const Vector3<Real>& rkU,
        const Vector3<Real>& rkV) const;  // u^T*M*v
    WM3_FUN Matrix3 TimesDiagonal (const Vector3<Real>& rkDiag) const;  // M*D
    WM3_FUN Matrix3 DiagonalTimes (const Vector3<Real>& rkDiag) const;  // D*M

	 WM3_OLD Matrix3 Transpose() const {return transpose();}
	 WM3_OLD Matrix3 Adjoint() const {return adjoint();}
	 WM3_OLD Matrix3 Determinant() const {return determinant();}
	 WM3_OLD Matrix3 Inverse() const {return inverse();}

    // The matrix must be a rotation for these functions to be valid.  The
    // last function uses Gram-Schmidt orthonormalization applied to the
    // columns of the rotation matrix.  The angle must be in radians, not
    // degrees.
    WM3_FUN void ToAxisAngle (Vector3<Real>& rkAxis, Real& rfAngle) const;
    WM3_FUN void Orthonormalize ();

    // The matrix must be symmetric.  Factor M = R * D * R^T where
    // R = [u0|u1|u2] is a rotation matrix with columns u0, u1, and u2 and
    // D = diag(d0,d1,d2) is a diagonal matrix whose diagonal entries are d0,
    // d1, and d2.  The eigenvector u[i] corresponds to eigenvector d[i].
    // The eigenvalues are ordered as d0 <= d1 <= d2.
    WM3_OLD void EigenDecomposition(Matrix3& rkRot, Matrix3& rkDiag) const { EigenDecomposition_(rkRot,rkDiag);}
    void EigenDecomposition_ (Matrix3& rkRot, Matrix3& rkDiag) const; // hidden function, do not call from outside

    // The matrix must be orthonormal.  The decomposition is yaw*pitch*roll
    // where yaw is rotation about the Up vector, pitch is rotation about the
    // Right axis, and roll is rotation about the Direction axis.
            Matrix3& FromEulerAnglesXYZ_ (Real fYAngle, Real fPAngle, Real fRAngle); // hidden function, do not call from outside
    WM3_OLD Matrix3& FromEulerAnglesXZY (Real fYAngle, Real fPAngle, Real fRAngle);
    WM3_OLD Matrix3& FromEulerAnglesYXZ (Real fYAngle, Real fPAngle, Real fRAngle);
    WM3_OLD Matrix3& FromEulerAnglesYZX (Real fYAngle, Real fPAngle, Real fRAngle);
    WM3_OLD Matrix3& FromEulerAnglesZXY (Real fYAngle, Real fPAngle, Real fRAngle);
    WM3_OLD Matrix3& FromEulerAnglesZYX (Real fYAngle, Real fPAngle, Real fRAngle);
    WM3_OLD bool ToEulerAnglesXYZ (Real& rfYAngle, Real& rfPAngle,
        Real& rfRAngle) const;
    WM3_OLD bool ToEulerAnglesXZY (Real& rfYAngle, Real& rfPAngle,
        Real& rfRAngle) const;
    WM3_OLD bool ToEulerAnglesYXZ (Real& rfYAngle, Real& rfPAngle,
        Real& rfRAngle) const;
    WM3_OLD bool ToEulerAnglesYZX (Real& rfYAngle, Real& rfPAngle,
        Real& rfRAngle) const;
    WM3_OLD bool ToEulerAnglesZXY (Real& rfYAngle, Real& rfPAngle,
        Real& rfRAngle) const;
    WM3_OLD bool ToEulerAnglesZYX (Real& rfYAngle, Real& rfPAngle,
        Real& rfRAngle) const;

    // SLERP (spherical linear interpolation) without quaternions.  Computes
    // R(t) = R0*(Transpose(R0)*R1)^t.  If Q is a rotation matrix with
    // unit-length axis U and angle A, then Q^t is a rotation matrix with
    // unit-length axis U and rotation angle t*A.
    WM3_FUN static Matrix3 Slerp (Real fT, const Matrix3& rkR0,
        const Matrix3& rkR1);

    // Singular value decomposition, M = L*S*R, where L and R are orthogonal
    // and S is a diagonal matrix whose diagonal entries are nonnegative.
    WM3_FUN void SingularValueDecomposition (Matrix3& rkL, Matrix3& rkS,
        Matrix3& rkR) const;
    WM3_FUN void SingularValueComposition (const Matrix3& rkL, const Matrix3& rkS,
        const Matrix3& rkR);

    // factor M = Q*D*U with orthogonal Q, diagonal D, upper triangular U
    WM3_FUN void QDUDecomposition (Matrix3& rkQ, Matrix3& rkD, Matrix3& rkU) const;

    // special matrices
    WM3_FUN WM3_ITEM static const Matrix3 ZERO;
    WM3_FUN WM3_ITEM static const Matrix3 IDENTITY;
	 EIG_FUN static const Matrix3 Identity(){ return IDENTITY; }
	 EIG_FUN static const Matrix3 Zero(){ return ZERO; }


private:
    // Support for eigendecomposition.  The Tridiagonalize function applies
    // a Householder transformation to the matrix.  If that transformation
    // is the identity (the matrix is already tridiagonal), then the return
    // value is 'false'.  Otherwise, the transformation is a reflection and
    // the return value is 'true'.  The QLAlgorithm returns 'true' iff the
    // QL iteration scheme converged.
    bool Tridiagonalize (Real afDiag[3], Real afSubd[2]);
    bool QLAlgorithm (Real afDiag[3], Real afSubd[2]);

    // support for singular value decomposition
    static void Bidiagonalize (Matrix3& rkA, Matrix3& rkL, Matrix3& rkR);
    static void GolubKahanStep (Matrix3& rkA, Matrix3& rkL, Matrix3& rkR);

    // for indexing into the 1D array of the matrix, iCol+N*iRow
    static int I (int iRow, int iCol);

    // support for comparisons
    int CompareArrays (const Matrix3& rkM) const;

    Real m_afEntry[9];
};

// c * M
template <class Real>
Matrix3<Real> operator* (Real fScalar, const Matrix3<Real>& rkM);

// v^T * M
template <class Real>
Vector3<Real> operator* (const Vector3<Real>& rkV, const Matrix3<Real>& rkM);

#include "Wm3Matrix3.inl"

typedef Matrix3<float> Matrix3f;
typedef Matrix3<double> Matrix3d;

}

#endif

