// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#include<boost/python.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/trim.hpp>
#include<string>
#include<stdexcept>
#include<sstream>
#include<iostream>

#include<yade/lib-base/Math.hpp>
#include<yade/lib-pyutil/doc_opts.hpp>


namespace bp=boost::python;

#define IDX_CHECK(i,MAX){ if(i<0 || i>=MAX) { PyErr_SetString(PyExc_IndexError, ("Index out of range 0.." + boost::lexical_cast<std::string>(MAX-1)).c_str()); bp::throw_error_already_set(); } }
#define IDX2_CHECKED_TUPLE_INTS(tuple,max2,arr2) {int l=bp::len(tuple); if(l!=2) { PyErr_SetString(PyExc_IndexError,"Index must be integer or a 2-tuple"); bp::throw_error_already_set(); } for(int _i=0; _i<2; _i++) { bp::extract<int> val(tuple[_i]); if(!val.check()) throw std::runtime_error("Unable to convert "+boost::lexical_cast<std::string>(_i)+"-th index to int."); int v=val(); IDX_CHECK(v,max2[_i]); arr2[_i]=v; }  }

void Vector3r_set_item(Vector3r & self, int idx, Real value){ IDX_CHECK(idx,3); self[idx]=value; }
void Vector3i_set_item(Vector3i & self, int idx, int  value){ IDX_CHECK(idx,3); self[idx]=value; }
void Vector2r_set_item(Vector2r & self, int idx, Real value){ IDX_CHECK(idx,2); self[idx]=value; }
void Vector2i_set_item(Vector2i & self, int idx, int  value){ IDX_CHECK(idx,2); self[idx]=value; }

void Quaternionr_set_item(Quaternionr & self, int idx, Real value){ IDX_CHECK(idx,4);  if(idx==0) self.x()=value; else if(idx==1) self.y()=value; else if(idx==2) self.z()=value; else if(idx==3) self.w()=value; }
void Matrix3r_set_item(Matrix3r & self, bp::tuple _idx, Real value){ int idx[2]; int mx[2]={3,3}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); self(idx[0],idx[1])=value; }
void Matrix3r_set_item_linear(Matrix3r & self, int idx, Real value){ IDX_CHECK(idx,9); self(idx/3,idx%3)=value; }

Real Vector3r_get_item(const Vector3r & self, int idx){ IDX_CHECK(idx,3); return self[idx]; }
int  Vector3i_get_item(const Vector3r & self, int idx){ IDX_CHECK(idx,3); return self[idx]; }
Real Vector2r_get_item(const Vector2r & self, int idx){ IDX_CHECK(idx,2); return self[idx]; }
int  Vector2i_get_item(const Vector2r & self, int idx){ IDX_CHECK(idx,2); return self[idx]; }

Real Quaternionr_get_item(const Quaternionr & self, int idx){ IDX_CHECK(idx,4); if(idx==0) return self.x(); if(idx==1) return self.y(); if(idx==2) return self.z(); return self.w(); }
Real Matrix3r_get_item(Matrix3r & self, bp::tuple _idx){ int idx[2]; int mx[2]={3,3}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); return self(idx[0],idx[1]); }
Real Matrix3r_get_item_linear(Matrix3r & self, int idx){ IDX_CHECK(idx,9); return self(idx/3,idx%3); }

std::string Vector3r_str(const Vector3r & self){ return std::string("Vector3(")+boost::lexical_cast<std::string>(self[0])+","+boost::lexical_cast<std::string>(self[1])+","+boost::lexical_cast<std::string>(self[2])+")";}
std::string Vector3i_str(const Vector3i & self){ return std::string("Vector3i(")+boost::lexical_cast<std::string>(self[0])+","+boost::lexical_cast<std::string>(self[1])+","+boost::lexical_cast<std::string>(self[2])+")";}
std::string Vector2r_str(const Vector2r & self){ return std::string("Vector2(")+boost::lexical_cast<std::string>(self[0])+","+boost::lexical_cast<std::string>(self[1])+")";}
std::string Vector2i_str(const Vector2i & self){ return std::string("Vector2i(")+boost::lexical_cast<std::string>(self[0])+","+boost::lexical_cast<std::string>(self[1])+")";}
std::string Quaternionr_str(const Quaternionr & self){ AngleAxisr aa(angleAxisFromQuat(self)); return std::string("Quaternion((")+boost::lexical_cast<std::string>(aa.axis()[0])+","+boost::lexical_cast<std::string>(aa.axis()[1])+","+boost::lexical_cast<std::string>(aa.axis()[2])+"),"+boost::lexical_cast<std::string>(aa.angle())+")";}
std::string Matrix3r_str(const Matrix3r & self){ std::ostringstream oss; oss<<"Matrix3("; for(int i=0; i<3; i++) for(int j=0; j<3; j++) oss<<self(i,j)<<((i==2 && j==2)?")":",")<<((i<2 && j==2)?" ":""); return oss.str(); }

int Vector3r_len(){return 3;}
int Vector3i_len(){return 3;}
int Vector2r_len(){return 2;}
int Vector2i_len(){return 2;}
int Quaternionr_len(){return 4;}
int Matrix3r_len(){return 9;}

// pickling support
struct Matrix3r_pickle: bp::pickle_suite{static bp::tuple getinitargs(const Matrix3r& x){ return bp::make_tuple(x(0,0),x(0,1),x(0,2),x(1,0),x(1,1),x(1,2),x(2,0),x(2,1),x(2,2));} };
struct Quaternionr_pickle: bp::pickle_suite{static bp::tuple getinitargs(const Quaternionr& x){ return bp::make_tuple(x.w(),x.x(),x.y(),x.z());} };
struct Vector3r_pickle: bp::pickle_suite{static bp::tuple getinitargs(const Vector3r& x){ return bp::make_tuple(x[0],x[1],x[2]);} };
struct Vector3i_pickle: bp::pickle_suite{static bp::tuple getinitargs(const Vector3i& x){ return bp::make_tuple(x[0],x[1],x[2]);} };
struct Vector2r_pickle: bp::pickle_suite{static bp::tuple getinitargs(const Vector2r& x){ return bp::make_tuple(x[0],x[1]);} };
struct Vector2i_pickle: bp::pickle_suite{static bp::tuple getinitargs(const Vector2i& x){ return bp::make_tuple(x[0],x[1]);} };

#undef IDX_CHECK

// automagic converters from sequence (list, tuple, …) to Vector{2,3}{r,i}
struct custom_Vector3r_from_sequence{
	custom_Vector3r_from_sequence(){	bp::converter::registry::push_back(&convertible,&construct,bp::type_id<Vector3r>()); }
	static void* convertible(PyObject* obj_ptr){ if(!PySequence_Check(obj_ptr) || PySequence_Size(obj_ptr)!=3) return 0;	return obj_ptr; }
	static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
		void* storage=((bp::converter::rvalue_from_python_storage<Vector3r>*)(data))->storage.bytes;
		new (storage) Vector3r(bp::extract<Real>(PySequence_GetItem(obj_ptr,0)),bp::extract<Real>(PySequence_GetItem(obj_ptr,1)),bp::extract<Real>(PySequence_GetItem(obj_ptr,2)));
		data->convertible=storage;
	}
};
struct custom_Vector3i_from_sequence{
	custom_Vector3i_from_sequence(){	bp::converter::registry::push_back(&convertible,&construct,bp::type_id<Vector3i>()); }
	static void* convertible(PyObject* obj_ptr){ if(!PySequence_Check(obj_ptr) || PySequence_Size(obj_ptr)!=3) return 0;	return obj_ptr; }
	static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
		void* storage=((bp::converter::rvalue_from_python_storage<Vector3i>*)(data))->storage.bytes;
		new (storage) Vector3i(bp::extract<int>(PySequence_GetItem(obj_ptr,0)),bp::extract<int>(PySequence_GetItem(obj_ptr,1)),bp::extract<int>(PySequence_GetItem(obj_ptr,2)));
		data->convertible=storage;
	}
};
struct custom_Vector2r_from_sequence{
	custom_Vector2r_from_sequence(){	bp::converter::registry::push_back(&convertible,&construct,bp::type_id<Vector2r>()); }
	static void* convertible(PyObject* obj_ptr){ if(!PySequence_Check(obj_ptr) || PySequence_Size(obj_ptr)!=2) return 0;	return obj_ptr; }
	static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
		void* storage=((bp::converter::rvalue_from_python_storage<Vector2r>*)(data))->storage.bytes;
		new (storage) Vector2r(bp::extract<Real>(PySequence_GetItem(obj_ptr,0)),bp::extract<Real>(PySequence_GetItem(obj_ptr,1)));
		data->convertible=storage;
	}
};
struct custom_Vector2i_from_sequence{
	custom_Vector2i_from_sequence(){	bp::converter::registry::push_back(&convertible,&construct,bp::type_id<Vector2i>()); }
	static void* convertible(PyObject* obj_ptr){ if(!PySequence_Check(obj_ptr) || PySequence_Size(obj_ptr)!=2) return 0;	return obj_ptr; }
	static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
		void* storage=((bp::converter::rvalue_from_python_storage<Vector2i>*)(data))->storage.bytes;
		new (storage) Vector2i(bp::extract<int>(PySequence_GetItem(obj_ptr,0)),bp::extract<int>(PySequence_GetItem(obj_ptr,1)));
		data->convertible=storage;
	}
};

static Matrix3r* Matrix3r_fromElements(Real m00, Real m01, Real m02, Real m10, Real m11, Real m12, Real m20, Real m21, Real m22){ Matrix3r* m(new Matrix3r); (*m)<<m00,m01,m02,m10,m11,m12,m20,m21,m22; return m; }
static Quaternionr Quaternionr_setFromTwoVectors(Quaternionr& q, const Vector3r& u, const Vector3r& v){ return q.setFromTwoVectors(u,v); }
static Vector3r Quaternionr_Rotate(Quaternionr& q, const Vector3r& u){ return q*u; }
static Quaternionr* Quaternionr_fromAxisAngle(const Vector3r& axis, const Real angle){ return new Quaternionr(AngleAxisr(angle,axis)); }
static Quaternionr* Quaternionr_fromAngleAxis(const Real angle, const Vector3r& axis){ return new Quaternionr(AngleAxisr(angle,axis)); }
static bp::tuple Quaternionr_toAxisAngle(const Quaternionr& self){ AngleAxisr aa(self); return bp::make_tuple(aa.axis(),aa.angle());}
static bp::tuple Quaternionr_toAngleAxis(const Quaternionr& self){ AngleAxisr aa(self); return bp::make_tuple(aa.angle(),aa.axis());}

static Real Vector3r_dot(const Vector3r& self, const Vector3r& v){ return self.dot(v); }
static Real Vector3i_dot(const Vector3i& self, const Vector3i& v){ return self.dot(v); }
static Real Vector2r_dot(const Vector2r& self, const Vector2r& v){ return self.dot(v); }
static Real Vector2i_dot(const Vector2i& self, const Vector2i& v){ return self.dot(v); }
static Vector3r Vector3r_cross(const Vector3r& self, const Vector3r& v){ return self.cross(v); }
static Vector3i Vector3i_cross(const Vector3i& self, const Vector3i& v){ return self.cross(v); }
static bool Quaternionr__eq__(const Quaternionr& q1, const Quaternionr& q2){ return q1==q2; }
static bool Quaternionr__neq__(const Quaternionr& q1, const Quaternionr& q2){ return q1!=q2; }
#ifndef YADE_WM3
	#include<Eigen/SVD>
	static bp::tuple Matrix3r_polarDecomposition(const Matrix3r& self){ Matrix3r unitary,positive; Eigen::SVD<Matrix3r>(self).computeUnitaryPositive(&unitary,&positive); return bp::make_tuple(unitary,positive); }
#endif

#define WM3_COMPAT

#ifdef WM3_COMPAT
	#define _PYCLASS(klass) boost::algorithm::trim_right_copy_if(std::string(klass),boost::algorithm::is_any_of("r"))
	#define WM3_OLD_METH0(klass,old,neww) static klass klass##_##old(){ std::cerr<<"WARN: "<<_PYCLASS(#klass)<<"."<<#old<<" is deprecated, use "<<_PYCLASS(#klass)<<"."<<#neww<<" instead"<<std::endl; return klass().neww(); }
	#define WM3_OLD_METH1(klass,old,neww,ret)  static ret klass##_##old(klass& self){ std::cerr<<"WARN: "<<_PYCLASS(#klass)<<"."<<#old<<" is deprecated, use "<<_PYCLASS(#klass)<<"."<<#neww<<" instead"<<std::endl; if(typeid(ret)!=typeid(void)) return self.neww(); return ret(); }
	#define WM3_OLD_METH2(klass,klass2,old,neww,ret)  static ret klass##_##old(klass& self,const klass2& arg){ std::cerr<<"WARN: "<<_PYCLASS(#klass)<<"."<<#old<<" is deprecated, use "<<_PYCLASS(#klass)<<"."<<#neww<<" instead"<<std::endl; if(typeid(ret)!=typeid(void)) return self.neww(arg); return ret(); }
	#define WM3_OLD_METH3(klass,klass2,klass3,old,neww,ret)  static ret klass##_##old(klass& self,const klass2& arg1, const klass3& arg2){ std::cerr<<"WARN: "<<_PYCLASS(#klass)<<"."<<#old<<" is deprecated, use "<<_PYCLASS(#klass)<<"."<<#neww<<" instead"<<std::endl; if(typeid(ret)!=typeid(void)) return self.neww(arg1,arg2); return ret(); }
	WM3_OLD_METH0(Matrix3r,IDENTITY,Identity)
	WM3_OLD_METH0(Matrix3r,ZERO,Zero)
	WM3_OLD_METH1(Matrix3r,Determinant,determinant,Real)
	WM3_OLD_METH1(Matrix3r,Inverse,inverse,Matrix3r)
	WM3_OLD_METH1(Matrix3r,Transpose,transpose,Matrix3r)

	WM3_OLD_METH0(Quaternionr,IDENTITY,Identity)
	bp::tuple Quaternionr_ToAxisAngle(const Quaternionr& self) { std::cerr<<"WARN: Quaternion.ToAxisAngle is deprecated, use Quaternion.toAxisAngle instead"<<std::endl; return Quaternionr_toAxisAngle(self); }
	WM3_OLD_METH1(Quaternionr,Conjugate,conjugate,Quaternionr)
	WM3_OLD_METH1(Quaternionr,Inverse,inverse,Quaternionr)
	WM3_OLD_METH1(Quaternionr,Normalize,normalize,void)
	WM3_OLD_METH1(Quaternionr,Length,norm,Real)
	WM3_OLD_METH3(Quaternionr,Vector3r,Vector3r,Align,setFromTwoVectors,Quaternionr)


	WM3_OLD_METH0(Vector3r,ZERO,Zero)
	WM3_OLD_METH0(Vector3r,ONE,Ones)
	WM3_OLD_METH0(Vector3r,UNIT_X,UnitX)
	WM3_OLD_METH0(Vector3r,UNIT_Y,UnitY)
	WM3_OLD_METH0(Vector3r,UNIT_Z,UnitZ)
	WM3_OLD_METH1(Vector3r,Length,norm,Real)
	WM3_OLD_METH1(Vector3r,Normalize,normalize,void)
	WM3_OLD_METH1(Vector3r,SquaredLength,squaredNorm,Real)
	WM3_OLD_METH2(Vector3r,Vector3r,Dot,dot,Real)
	WM3_OLD_METH2(Vector3r,Vector3r,Cross,cross,Vector3r)

	WM3_OLD_METH0(Vector2r,ZERO,Zero)
	WM3_OLD_METH0(Vector2r,ONE,Ones)
	WM3_OLD_METH0(Vector2r,UNIT_X,UnitX)
	WM3_OLD_METH0(Vector2r,UNIT_Y,UnitY)
	WM3_OLD_METH1(Vector2r,Length,norm,Real)
	WM3_OLD_METH1(Vector2r,Normalize,normalize,void)
	WM3_OLD_METH1(Vector2r,SquaredLength,squaredNorm,Real)
	WM3_OLD_METH2(Vector2r,Vector2r,Dot,dot,Real)

	static void Matrix3r_fromAxisAngle(Matrix3r& self, const Vector3r& axis, const Real angle){ std::cerr<<"Matrix3.fromAxisAngle is deprecated, use Quaternion.toRotationMatrix instead"<<std::endl; self=AngleAxisr(angle,axis).toRotationMatrix(); }
#endif

#define EIG_WRAP_METH1(klass,meth) static klass klass##_##meth(const klass& self){ return self.meth(); }
//#define EIG_WRAP_METH0(klass,meth) static klass klass##_##meth(){ return klass().meth(); }
#define EIG_WRAP_METH0(klass,meth) static const klass klass##_##meth=klass().meth();
EIG_WRAP_METH1(Matrix3r,transpose);
EIG_WRAP_METH1(Matrix3r,inverse);

EIG_WRAP_METH0(Matrix3r,Zero);
EIG_WRAP_METH0(Matrix3r,Identity);
EIG_WRAP_METH0(Vector3r,Zero); EIG_WRAP_METH0(Vector3r,UnitX); EIG_WRAP_METH0(Vector3r,UnitY); EIG_WRAP_METH0(Vector3r,UnitZ); EIG_WRAP_METH0(Vector3r,Ones);
EIG_WRAP_METH0(Vector3i,Zero); EIG_WRAP_METH0(Vector3i,UnitX); EIG_WRAP_METH0(Vector3i,UnitY); EIG_WRAP_METH0(Vector3i,UnitZ); EIG_WRAP_METH0(Vector3i,Ones);
EIG_WRAP_METH0(Vector2r,Zero); EIG_WRAP_METH0(Vector2r,UnitX); EIG_WRAP_METH0(Vector2r,UnitY); EIG_WRAP_METH0(Vector2r,Ones);
EIG_WRAP_METH0(Vector2i,Zero); EIG_WRAP_METH0(Vector2i,UnitX); EIG_WRAP_METH0(Vector2i,UnitY); EIG_WRAP_METH0(Vector2i,Ones);
EIG_WRAP_METH0(Quaternionr,Identity);

#define EIG_OP1(klass,op,sym) typeof((sym klass()).eval()) klass##op(const klass& self){ return (sym self).eval();}
#define EIG_OP2(klass,op,sym,klass2) typeof((klass() sym klass2()).eval()) klass##op##klass2(const klass& self, const klass2& other){ return (self sym other).eval(); }
#define EIG_OP2_INPLACE(klass,op,sym,klass2) klass klass##op##klass2(klass& self, const klass2& other){ self sym other; return self; }


EIG_OP1(Matrix3r,__neg__,-)
EIG_OP2(Matrix3r,__add__,+,Matrix3r) EIG_OP2_INPLACE(Matrix3r,__iadd__,+=,Matrix3r)
EIG_OP2(Matrix3r,__sub__,-,Matrix3r) EIG_OP2_INPLACE(Matrix3r,__isub__,-=,Matrix3r)
EIG_OP2(Matrix3r,__mul__,*,Real) EIG_OP2(Matrix3r,__rmul__,*,Real) EIG_OP2_INPLACE(Matrix3r,__imul__,*=,Real)
EIG_OP2(Matrix3r,__mul__,*,int) EIG_OP2(Matrix3r,__rmul__,*,int) EIG_OP2_INPLACE(Matrix3r,__imul__,*=,int)
EIG_OP2(Matrix3r,__mul__,*,Vector3r) EIG_OP2(Matrix3r,__rmul__,*,Vector3r)
EIG_OP2(Matrix3r,__mul__,*,Matrix3r) EIG_OP2_INPLACE(Matrix3r,__imul__,*=,Matrix3r)
EIG_OP2(Matrix3r,__div__,/,Real) EIG_OP2_INPLACE(Matrix3r,__idiv__,/=,Real)
EIG_OP2(Matrix3r,__div__,/,int) EIG_OP2_INPLACE(Matrix3r,__idiv__,/=,int)

EIG_OP1(Vector3r,__neg__,-);
EIG_OP2(Vector3r,__add__,+,Vector3r); EIG_OP2_INPLACE(Vector3r,__iadd__,+=,Vector3r)
EIG_OP2(Vector3r,__sub__,-,Vector3r); EIG_OP2_INPLACE(Vector3r,__isub__,-=,Vector3r)
EIG_OP2(Vector3r,__mul__,*,Real) EIG_OP2(Vector3r,__rmul__,*,Real) EIG_OP2_INPLACE(Vector3r,__imul__,*=,Real) EIG_OP2(Vector3r,__div__,/,Real) EIG_OP2_INPLACE(Vector3r,__idiv__,/=,Real)
EIG_OP2(
Vector3r,__mul__,*,int) EIG_OP2(Vector3r,__rmul__,*,int) EIG_OP2_INPLACE(Vector3r,__imul__,*=,int) EIG_OP2(Vector3r,__div__,/,int) EIG_OP2_INPLACE(Vector3r,__idiv__,/=,int)

EIG_OP1(Vector3i,__neg__,-);
EIG_OP2(Vector3i,__add__,+,Vector3i); EIG_OP2_INPLACE(Vector3i,__iadd__,+=,Vector3i)
EIG_OP2(Vector3i,__sub__,-,Vector3i); EIG_OP2_INPLACE(Vector3i,__isub__,-=,Vector3i)
EIG_OP2(Vector3i,__mul__,*,int) EIG_OP2(Vector3i,__rmul__,*,int)  EIG_OP2_INPLACE(Vector3i,__imul__,*=,int)

EIG_OP1(Vector2r,__neg__,-);
EIG_OP2(Vector2r,__add__,+,Vector2r); EIG_OP2_INPLACE(Vector2r,__iadd__,+=,Vector2r)
EIG_OP2(Vector2r,__sub__,-,Vector2r); EIG_OP2_INPLACE(Vector2r,__isub__,-=,Vector2r)
EIG_OP2(Vector2r,__mul__,*,Real) EIG_OP2(Vector2r,__rmul__,*,Real) EIG_OP2_INPLACE(Vector2r,__imul__,*=,Real) EIG_OP2(Vector2r,__div__,/,Real) EIG_OP2_INPLACE(Vector2r,__idiv__,/=,Real)
EIG_OP2(Vector2r,__mul__,*,int) EIG_OP2(Vector2r,__rmul__,*,int) EIG_OP2_INPLACE(Vector2r,__imul__,*=,int) EIG_OP2(Vector2r,__div__,/,int) EIG_OP2_INPLACE(Vector2r,__idiv__,/=,int)

EIG_OP1(Vector2i,__neg__,-);
EIG_OP2(Vector2i,__add__,+,Vector2i); EIG_OP2_INPLACE(Vector2i,__iadd__,+=,Vector2i)
EIG_OP2(Vector2i,__sub__,-,Vector2i); EIG_OP2_INPLACE(Vector2i,__isub__,-=,Vector2i)
EIG_OP2(Vector2i,__mul__,*,int)  EIG_OP2_INPLACE(Vector2i,__imul__,*=,int) EIG_OP2(Vector2i,__rmul__,*,int)


BOOST_PYTHON_MODULE(miniEigen){
	bp::scope().attr("__doc__")="Basic math functions for Yade: small matrix, vector and quaternion classes. This module internally wraps small parts of the `Eigen <http://eigen.tuxfamily.org>`_ library. Refer to its documentation for details. All classes in this module support pickling.";

	YADE_SET_DOCSTRING_OPTS;

	custom_Vector3r_from_sequence();
	custom_Vector3i_from_sequence();
	custom_Vector2r_from_sequence();
	custom_Vector2i_from_sequence();

	bp::class_<Matrix3r>("Matrix3","3x3 float matrix.\n\nSupported operations (``m`` is a Matrix3, ``f`` if a float/int, ``v`` is a Vector3): ``-m``, ``m+m``, ``m+=m``, ``m-m``, ``m-=m``, ``m*f``, ``f*m``, ``m*=f``, ``m/f``, ``m/=f``, ``m*m``, ``m*=m``, ``m*v``, ``v*m``, ``m==m``, ``m!=m``.",bp::init<>())
		.def(bp::init<Matrix3r const &>((bp::arg("m"))))
		.def("__init__",bp::make_constructor(&Matrix3r_fromElements,bp::default_call_policies(),(bp::arg("m00"),bp::arg("m01"),bp::arg("m02"),bp::arg("m10"),bp::arg("m11"),bp::arg("m12"),bp::arg("m20"),bp::arg("m21"),bp::arg("m22"))))
		.def_pickle(Matrix3r_pickle())
		//
		.def("determinant",&Matrix3r::determinant)
		.def("inverse",&Matrix3r_inverse)
		.def("transpose",&Matrix3r_transpose)
		#ifndef YADE_WM3
			.def("polarDecomposition",&Matrix3r_polarDecomposition)
		#endif

		//
		.def("__neg__",&Matrix3r__neg__)
		.def("__add__",&Matrix3r__add__Matrix3r).def("__iadd__",&Matrix3r__iadd__Matrix3r)
		.def("__sub__",&Matrix3r__sub__Matrix3r).def("__isub__",&Matrix3r__isub__Matrix3r)
		.def("__mul__",&Matrix3r__mul__Real).def("__rmul__",&Matrix3r__rmul__Real).def("__imul__",&Matrix3r__imul__Real)
		.def("__mul__",&Matrix3r__mul__int).def("__rmul__",&Matrix3r__rmul__int).def("__imul__",&Matrix3r__imul__int)
		.def("__mul__",&Matrix3r__mul__Vector3r).def("__rmul__",&Matrix3r__rmul__Vector3r)
		.def("__mul__",&Matrix3r__mul__Matrix3r).def("__imul__",&Matrix3r__imul__Matrix3r)
		.def("__div__",&Matrix3r__div__Real).def("__idiv__",&Matrix3r__idiv__Real)
		.def("__div__",&Matrix3r__div__int).def("__idiv__",&Matrix3r__idiv__int)
		.def(bp::self == bp::self)
		.def(bp::self != bp::self)
		//
 		.def("__len__",&::Matrix3r_len).staticmethod("__len__").def("__setitem__",&::Matrix3r_set_item).def("__getitem__",&::Matrix3r_get_item).def("__str__",&::Matrix3r_str).def("__repr__",&::Matrix3r_str)
		/* extras for matrices */
		.def("__setitem__",&::Matrix3r_set_item_linear).def("__getitem__",&::Matrix3r_get_item_linear)
		.def_readonly("Identity",&Matrix3r_Identity).def_readonly("Zero",&Matrix3r_Zero)
		#ifdef WM3_COMPAT
			// wm3 compat
			.add_static_property("IDENTITY",&Matrix3r_IDENTITY,"|ydeprecated|").add_static_property("ZERO",&Matrix3r_ZERO,"|ydeprecated|")
			.def("Determinant",&Matrix3r_Determinant,"|ydeprecated|").def("Inverse",&Matrix3r_Inverse,"|ydeprecated|").def("Transpose",&Matrix3r_Transpose,"|ydeprecated|")
			.def("FromAxisAngle",&Matrix3r_fromAxisAngle,"|ydeprecated|")
		#endif
	;
	bp::class_<Quaternionr>("Quaternion","Quaternion representing rotation.\n\nSupported operations (``q`` is a Quaternion, ``v`` is a Vector3): ``q*q`` (rotation composition), ``q*=q``, ``q*v`` (rotating ``v`` by ``q``), ``q==q``, ``q!=q``.",bp::init<>())
		.def("__init__",bp::make_constructor(&Quaternionr_fromAxisAngle,bp::default_call_policies(),(bp::arg("axis"),bp::arg("angle"))))
		.def("__init__",bp::make_constructor(&Quaternionr_fromAngleAxis,bp::default_call_policies(),(bp::arg("angle"),bp::arg("axis"))))
		.def(bp::init<Real,Real,Real,Real>((bp::arg("w"),bp::arg("x"),bp::arg("y"),bp::arg("z")),"Initialize from coefficients.\n\n.. note:: The order of coefficients is *w*, *x*, *y*, *z*. The [] operator numbers them differently, 0…4 for *x* *y* *z* *w*!"))
		.def(bp::init<Quaternionr>((bp::arg("other"))))
		.def_pickle(Quaternionr_pickle())
		// properties
		.def_readonly("Identity",&Quaternionr_Identity)
		// methods
		.def("setFromTwoVectors",&Quaternionr_setFromTwoVectors,((bp::arg("u"),bp::arg("v"))))
		.def("conjugate",&Quaternionr::conjugate)
		.def("toAxisAngle",&Quaternionr_toAxisAngle).def("toAngleAxis",&Quaternionr_toAngleAxis)
		.def("toRotationMatrix",&Quaternionr::toRotationMatrix)
		.def("Rotate",&Quaternionr_Rotate,((bp::arg("v"))))
		.def("inverse",&Quaternionr::inverse)
		.def("norm",&Quaternionr::norm)
		.def("normalize",&Quaternionr::normalize)
		// operators
		.def(bp::self * bp::self)
		.def(bp::self *= bp::self)
		.def(bp::self * bp::other<Vector3r>())
		//.def(bp::self != bp::self).def(bp::self == bp::self) // these don't work... (?)
		.def("__eq__",&Quaternionr__eq__).def("__neq__",&Quaternionr__neq__)
		// specials
		.def("__len__",&Quaternionr_len).staticmethod("__len__")
		.def("__setitem__",&Quaternionr_set_item).def("__getitem__",&Quaternionr_get_item)
		.def("__str__",&Quaternionr_str).def("__repr__",&Quaternionr_str)
		#ifdef WM3_COMPAT
			.def("Align",&Quaternionr_Align,"|ydeprecated|").def("Conjugate",&Quaternionr_Conjugate,"|ydeprecated|").def("Inverse",&Quaternionr_Inverse,"|ydeprecated|").def("Length",&Quaternionr_Length,"|ydeprecated|").def("Normalize",&Quaternionr_Normalize,"|ydeprecated|").add_static_property("IDENTITY",&Quaternionr_IDENTITY)
			.def("ToAxisAngle",&Quaternionr_ToAxisAngle,"|ydeprecated|")
		#endif
	;
	bp::class_<Vector3r>("Vector3","3-dimensional float vector.\n\nSupported operatrions (``f`` if a float/int, ``v`` is a Vector3): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*f``, ``f*v``, ``v*=f``, ``v/f``, ``v/=f``, ``v==v``, ``v!=v``, plus operations with ``Matrix3`` and ``Quaternion``.\n\nImplicit conversion from sequence (list,tuple, …) of 3 floats.",bp::init<>())
		.def(bp::init<Vector3r>((bp::arg("other"))))
		.def(bp::init<Real,Real,Real>((bp::arg("x"),bp::arg("y"),bp::arg("z"))))
		.def_pickle(Vector3r_pickle())
		// properties
		.def_readonly("Ones",&Vector3r_Ones).def_readonly("Zero",&Vector3r_Zero)
		.def_readonly("UnitX",&Vector3r_UnitX).def_readonly("UnitY",&Vector3r_UnitY).def_readonly("UnitZ",&Vector3r_UnitZ)
		// methods
		.def("dot",&Vector3r_dot).def("cross",&Vector3r_cross)
		.def("norm",&Vector3r::norm).def("squaredNorm",&Vector3r::squaredNorm).def("normalize",&Vector3r::normalize)
		// operators
		.def("__neg__",&Vector3r__neg__) // -v
		.def("__add__",&Vector3r__add__Vector3r).def("__iadd__",&Vector3r__iadd__Vector3r) // +, +=
		.def("__sub__",&Vector3r__sub__Vector3r).def("__isub__",&Vector3r__isub__Vector3r) // -, -=
		.def("__mul__",&Vector3r__mul__Real).def("__rmul__",&Vector3r__rmul__Real) // f*v, v*f
		.def("__div__",&Vector3r__div__Real).def("__idiv__",&Vector3r__idiv__Real) // v/f, v/=f
		.def("__mul__",&Vector3r__mul__int).def("__rmul__",&Vector3r__rmul__int) // f*v, v*f
		.def("__div__",&Vector3r__div__int).def("__idiv__",&Vector3r__idiv__int) // v/f, v/=f
		.def(bp::self != bp::self).def(bp::self == bp::self)
		// specials
		.def("__len__",&::Vector3r_len).staticmethod("__len__")
		.def("__setitem__",&::Vector3r_set_item).def("__getitem__",&::Vector3r_get_item)
		.def("__str__",&::Vector3r_str).def("__repr__",&::Vector3r_str)
		#ifdef WM3_COMPAT
			// wm3 compat
			.def("Dot",&Vector3r_Dot,"|ydeprecated|").def("Cross",&Vector3r_Cross,"|ydeprecated|").def("Length",&Vector3r_Length,"|ydeprecated|").def("SquaredLength",&Vector3r_SquaredLength,"|ydeprecated|").def("Normalize",&Vector3r_Normalize,"|ydeprecated|")
			.add_static_property("ONE",&Vector3r_ONE,"|ydeprecated|").add_static_property("ZERO",&Vector3r_ZERO,"|ydeprecated|").add_static_property("UNIT_X",&Vector3r_UNIT_X,"|ydeprecated|").add_static_property("UNIT_Y",&Vector3r_UNIT_Y,"|ydeprecated|").add_static_property("UNIT_Z",&Vector3r_UNIT_Z,"|ydeprecated|")
		#endif
	;	
	bp::class_<Vector3i>("Vector3i","3-dimensional integer vector.\n\nSupported operatrions (``i`` if an int, ``v`` is a Vector3i): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*i``, ``i*v``, ``v*=i``, ``v==v``, ``v!=v``.\n\nImplicit conversion from sequence  (list,tuple, …) of 3 integers.",bp::init<>())
		.def(bp::init<Vector3i>((bp::arg("other"))))
		.def(bp::init<int,int,int>((bp::arg("x"),bp::arg("y"),bp::arg("z"))))
		.def_pickle(Vector3i_pickle())
		// properties
		.def_readonly("Ones",&Vector3i_Ones).def_readonly("Zero",&Vector3i_Zero)
		.def_readonly("UnitX",&Vector3i_UnitX).def_readonly("UnitY",&Vector3i_UnitY).def_readonly("UnitZ",&Vector3i_UnitZ)
		// methods
		.def("dot",&Vector3i_dot).def("cross",&Vector3i_cross)
		.def("norm",&Vector3i::norm).def("squaredNorm",&Vector3i::squaredNorm)
		// operators
		.def("__neg__",&Vector3i__neg__) // -v
		.def("__add__",&Vector3i__add__Vector3i).def("__iadd__",&Vector3i__iadd__Vector3i) // +, +=
		.def("__sub__",&Vector3i__sub__Vector3i).def("__isub__",&Vector3i__isub__Vector3i) // -, -=
		.def("__mul__",&Vector3i__mul__int).def("__rmul__",&Vector3i__rmul__int) // f*v, v*f
		.def(bp::self != bp::self).def(bp::self == bp::self)
		// specials
		.def("__len__",&::Vector3i_len).staticmethod("__len__")
		.def("__setitem__",&::Vector3i_set_item).def("__getitem__",&::Vector3i_get_item)
		.def("__str__",&::Vector3i_str).def("__repr__",&::Vector3i_str)
	;	
	bp::class_<Vector2r>("Vector2","3-dimensional float vector.\n\nSupported operatrions (``f`` if a float/int, ``v`` is a Vector3): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*f``, ``f*v``, ``v*=f``, ``v/f``, ``v/=f``, ``v==v``, ``v!=v``.\n\nImplicit conversion from sequence (list,tuple, …) of 2 floats.",bp::init<>())
		.def(bp::init<Vector2r>((bp::arg("other"))))
		.def(bp::init<Real,Real>((bp::arg("x"),bp::arg("y"))))
		.def_pickle(Vector2r_pickle())
		// properties
		.def_readonly("Ones",&Vector2r_Ones).def_readonly("Zero",&Vector2r_Zero)
		.def_readonly("UnitX",&Vector2r_UnitX).def_readonly("UnitY",&Vector2r_UnitY)
		// methods
		.def("dot",&Vector2r_dot)
		.def("norm",&Vector2r::norm).def("squaredNorm",&Vector2r::squaredNorm).def("normalize",&Vector2r::normalize)
		// operators
		.def("__neg__",&Vector2r__neg__) // -v
		.def("__add__",&Vector2r__add__Vector2r).def("__iadd__",&Vector2r__iadd__Vector2r) // +, +=
		.def("__sub__",&Vector2r__sub__Vector2r).def("__isub__",&Vector2r__isub__Vector2r) // -, -=
		.def("__mul__",&Vector2r__mul__Real).def("__rmul__",&Vector2r__rmul__Real) // f*v, v*f
		.def("__div__",&Vector2r__div__Real).def("__idiv__",&Vector2r__idiv__Real) // v/f, v/=f
		.def("__mul__",&Vector2r__mul__int).def("__rmul__",&Vector2r__rmul__int) // f*v, v*f
		.def("__div__",&Vector2r__div__int).def("__idiv__",&Vector2r__idiv__int) // v/f, v/=f
		.def(bp::self != bp::self).def(bp::self == bp::self)
		// specials
		.def("__len__",&::Vector2r_len).staticmethod("__len__")
		.def("__setitem__",&::Vector2r_set_item).def("__getitem__",&::Vector2r_get_item)
		.def("__str__",&::Vector2r_str).def("__repr__",&::Vector2r_str)
		#ifdef WM3_COMPAT
			// wm3 compat
			.def("Dot",&Vector2r_Dot,"|ydeprecated|").def("Length",&Vector2r_Length,"|ydeprecated|").def("SquaredLength",&Vector2r_SquaredLength,"|ydeprecated|").def("Normalize",&Vector2r_Normalize,"|ydeprecated|")
			.add_static_property("ONE",&Vector2r_ONE,"|ydeprecated|").add_static_property("ZERO",&Vector2r_ZERO,"|ydeprecated|").add_static_property("UNIT_X",&Vector2r_UNIT_X,"|ydeprecated|").add_static_property("UNIT_Y",&Vector2r_UNIT_Y,"|ydeprecated|")
		#endif
	;	
	bp::class_<Vector2i>("Vector2i","2-dimensional integer vector.\n\nSupported operatrions (``i`` if an int, ``v`` is a Vector2i): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*i``, ``i*v``, ``v*=i``, ``v==v``, ``v!=v``.\n\nImplicit conversion from sequence (list,tuple, …) of 2 integers.",bp::init<>())
		.def(bp::init<Vector2i>((bp::arg("other"))))
		.def(bp::init<int,int>((bp::arg("x"),bp::arg("y"))))
		.def_pickle(Vector2i_pickle())
		// properties
		.def_readonly("Ones",&Vector2i_Ones).def_readonly("Zero",&Vector2i_Zero)
		.def_readonly("UnitX",&Vector2i_UnitX).def_readonly("UnitY",&Vector2i_UnitY)
		// methods
		.def("dot",&Vector2i_dot)
		.def("norm",&Vector2i::norm).def("squaredNorm",&Vector2i::squaredNorm).def("normalize",&Vector2i::normalize)
		// operators
		.def("__neg__",&Vector2i__neg__) // -v
		.def("__add__",&Vector2i__add__Vector2i).def("__iadd__",&Vector2i__iadd__Vector2i) // +, +=
		.def("__sub__",&Vector2i__sub__Vector2i).def("__isub__",&Vector2i__isub__Vector2i) // -, -=
		.def("__mul__",&Vector2i__mul__int).def("__rmul__",&Vector2i__rmul__int) // f*v, v*f
		.def(bp::self != bp::self).def(bp::self == bp::self)
		// specials
		.def("__len__",&::Vector2i_len).staticmethod("__len__")
		.def("__setitem__",&::Vector2i_set_item).def("__getitem__",&::Vector2i_get_item)
		.def("__str__",&::Vector2i_str).def("__repr__",&::Vector2i_str)
	;	
	
};








