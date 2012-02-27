// 2008,2009 © Václav Šmilauer <eudoxos@arcig.cz>
#include<boost/lexical_cast.hpp>
#include<string>
#include<stdexcept>
#include<sstream>

#define IDX_CHECK(i,MAX){ if(i<0 || i>=MAX) { PyErr_SetString(PyExc_IndexError, ("Index out of range 0.." + boost::lexical_cast<std::string>(MAX-1)).c_str()); boost::python::throw_error_already_set(); } }
#define IDX2_CHECKED_TUPLE_INTS(tuple,max2,arr2) {int l=boost::python::len(tuple); if(l!=2) { PyErr_SetString(PyExc_IndexError,"Index must be integer or a 2-tuple"); boost::python::throw_error_already_set(); } for(int _i=0; _i<2; _i++) { boost::python::extract<int> val(tuple[_i]); if(!val.check()) throw std::runtime_error("Unable to convert "+boost::lexical_cast<std::string>(_i)+"-th index to int."); int v=val(); IDX_CHECK(v,max2[_i]); arr2[_i]=v; }  }

void Vector2_set_item(Vector2r & self, int idx, Real value){ IDX_CHECK(idx,2); self[idx]=value; }
void Vector3_set_item(Vector3r & self, int idx, Real value){ IDX_CHECK(idx,3); self[idx]=value; }
void Quaternion_set_item(Quaternionr & self, int idx, Real value){ IDX_CHECK(idx,4);  if(idx==0) self.w()=value; else if(idx==1) self.x()=value; else if(idx==2) self.y()=value; else if(idx==3) self.z()=value; }
void Matrix3_set_item(Matrix3r & self, boost::python::tuple _idx, Real value){ int idx[2]; int mx[2]={3,3}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); self(idx[0],idx[1])=value; }
void Matrix3_set_item_linear(Matrix3r & self, int idx, Real value){ IDX_CHECK(idx,9); self(idx/3,idx%3)=value; }

Real Vector2_get_item(const Vector2r & self, int idx){ IDX_CHECK(idx,2); return self[idx]; }
Real Vector3_get_item(const Vector3r & self, int idx){ IDX_CHECK(idx,3); return self[idx]; }
Real Quaternion_get_item(const Quaternionr & self, int idx){ IDX_CHECK(idx,4); if(idx==0) return self.w(); if(idx==1) return self.x(); if(idx==2) return self.y(); return self.z(); }
Real Matrix3_get_item(Matrix3r & self, boost::python::tuple _idx){ int idx[2]; int mx[2]={3,3}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); return self(idx[0],idx[1]); }
Real Matrix3_get_item_linear(Matrix3r & self, int idx){ IDX_CHECK(idx,9); return self(idx/3,idx%3); }

std::string Vector2_str(const Vector2r & self){ return std::string("Vector2(")+boost::lexical_cast<std::string>(self[0])+","+boost::lexical_cast<std::string>(self[1])+")";}
std::string Vector3_str(const Vector3r & self){ return std::string("Vector3(")+boost::lexical_cast<std::string>(self[0])+","+boost::lexical_cast<std::string>(self[1])+","+boost::lexical_cast<std::string>(self[2])+")";}
std::string Quaternion_str(const Quaternionr & self){ AngleAxisr aa(angleAxisFromQuat(self)); return std::string("Quaternion((")+boost::lexical_cast<std::string>(aa.axis()[0])+","+boost::lexical_cast<std::string>(aa.axis()[1])+","+boost::lexical_cast<std::string>(aa.axis()[2])+"),"+boost::lexical_cast<std::string>(aa.angle())+")";}
std::string Matrix3_str(const Matrix3r & self){ std::ostringstream oss; oss<<"Matrix3("; for(int i=0; i<3; i++) for(int j=0; j<3; j++) oss<<self(i,j)<<((i==2 && j==2)?")":",")<<((i<2 && j==2)?" ":""); return oss.str(); }

int Vector2_len(){return 2;}
int Vector3_len(){return 3;}
int Quaternion_len(){return 4;}
int Matrix3_len(){return 9;}
#undef IDX_CHECK

#if 1
	// workarounds for static methods; disabled for now, see comment in miniWm3Wrap-generate.py
	#define _WORKAROUND(clss,member,member2) clss clss##_##member(const clss& self){return clss::member2();}
	_WORKAROUND(Vector2r,ONE,Ones); _WORKAROUND(Vector2r,UNIT_X,UnitX); _WORKAROUND(Vector2r,UNIT_Y,UnitY); _WORKAROUND(Vector2r,ZERO,Zero);
	_WORKAROUND(Vector3r,ONE,Ones); _WORKAROUND(Vector3r,UNIT_X,UnitX); _WORKAROUND(Vector3r,UNIT_Y,UnitY); _WORKAROUND(Vector3r,UNIT_Z,UnitZ); _WORKAROUND(Vector3r,ZERO,Zero);
	_WORKAROUND(Quaternionr,IDENTITY,Identity);
	_WORKAROUND(Matrix3r,IDENTITY,Identity); _WORKAROUND(Matrix3r,ZERO,Zero);
	#undef _WORKAROUND
#endif

// automagic converter from tuple to Vector3r
struct custom_Vector3r_from_tuple{
	custom_Vector3r_from_tuple(){
		boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<Vector3r>());
	}
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr) || PySequence_Size(obj_ptr)!=3) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data){
		void* storage=((boost::python::converter::rvalue_from_python_storage<Vector3r>*)(data))->storage.bytes;
		new (storage) Vector3r(boost::python::extract<double>(PySequence_GetItem(obj_ptr,0)),boost::python::extract<double>(PySequence_GetItem(obj_ptr,1)),boost::python::extract<double>(PySequence_GetItem(obj_ptr,2)));
		data->convertible=storage;
	}
};

using namespace Wm3;
