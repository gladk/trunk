// 2009 © Václav Šmilauer <eudoxos@arcig.cz>


// this is not currently used, but can be enabled if needed
// probably breaks compilation for older (like <=1.35 or so)
// boost::python
#if 0
	#include<indexing_suite/container_suite.hpp>
	#include<indexing_suite/vector.hpp>
#endif

#include<yade/lib-pyutil/numpy.hpp>

#include<boost/python.hpp>
#include<boost/python/class.hpp>
#include<boost/python/module.hpp>
#include<boost/foreach.hpp>
#ifndef FOREACH
	#define FOREACH BOOST_FOREACH
#endif

#include<vector>
#include<string>
#include<stdexcept>
#include<iostream>
#include<map>

#include<yade/lib-base/Math.hpp>

#include<yade/core/Engine.hpp>

#include<yade/pkg-common/BoundFunctor.hpp>
#include<yade/pkg-common/InteractionGeometryFunctor.hpp>
#include<yade/pkg-common/InteractionPhysicsFunctor.hpp>
#include<yade/pkg-common/LawFunctor.hpp>
#include<yade/pkg-common/Callbacks.hpp>




using namespace boost::python;

// will be removed later
#ifdef YADE_WM3
struct custom_Vector3r_from_seq{
	custom_Vector3r_from_seq(){
		 converter::registry::push_back(&convertible,&construct,type_id<Vector3r>());
	}
	static void* convertible(PyObject* obj_ptr){
		 if(!PySequence_Check(obj_ptr) || PySequence_Size(obj_ptr)!=3) return 0;
		 return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data){
		 void* storage=((converter::rvalue_from_python_storage<Vector3r>*)(data))->storage.bytes;
		 new (storage) Vector3r(extract<Real>(PySequence_GetItem(obj_ptr,0)),extract<Real>(PySequence_GetItem(obj_ptr,1)),extract<Real>(PySequence_GetItem(obj_ptr,2)));
		 data->convertible=storage;
	}
};

struct custom_vector3i_to_seq{
	static PyObject* convert(const Vector3i v3i){
		python::tuple ret(python::make_tuple(v3i[0],v3i[1],v3i[2]));
		return incref(ret.ptr());
	}
};
#endif

// move this to the miniEigen wrapper later

/* two-way se3 handling */
struct custom_se3_to_tuple{
	static PyObject* convert(const Se3r& se3){
		python::tuple ret=python::make_tuple(se3.position,se3.orientation);
		return incref(ret.ptr());
	}
};
struct custom_Se3r_from_seq{
	custom_Se3r_from_seq(){
		 converter::registry::push_back(&convertible,&construct,type_id<Se3r>());
	}
	static void* convertible(PyObject* obj_ptr){
		 if(!PySequence_Check(obj_ptr)) return 0;
		 if(PySequence_Size(obj_ptr)!=2 && PySequence_Size(obj_ptr)!=7) return 0;
		 return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data){
		void* storage=((converter::rvalue_from_python_storage<Se3r>*)(data))->storage.bytes;
		new (storage) Se3r; Se3r* se3=(Se3r*)storage;
		if(PySequence_Size(obj_ptr)==2){ // from vector and quaternion
			se3->position=extract<Vector3r>(PySequence_GetItem(obj_ptr,0));
			se3->orientation=extract<Quaternionr>(PySequence_GetItem(obj_ptr,1));
		} else if(PySequence_Size(obj_ptr)==7){ // 3 vector components, 3 axis components, angle
			se3->position=Vector3r(extract<Real>(PySequence_GetItem(obj_ptr,0)),extract<Real>(PySequence_GetItem(obj_ptr,1)),extract<Real>(PySequence_GetItem(obj_ptr,2)));
			Vector3r axis=Vector3r(extract<Real>(PySequence_GetItem(obj_ptr,3)),extract<Real>(PySequence_GetItem(obj_ptr,4)),extract<Real>(PySequence_GetItem(obj_ptr,5)));
			Real angle=extract<Real>(PySequence_GetItem(obj_ptr,6));
			se3->orientation=Quaternionr(AngleAxisr(angle,axis));
		} else throw std::logic_error(__FILE__ ": First, the sequence size for Se3r object was 2 or 7, but now is not? (programming error, please report!");
		data->convertible=storage;
	}
};

/*** c++-vector to python-list and python-(list,sequence) to c++-vector converters ***/
template<typename containedType>
struct custom_vector_to_list{
	static PyObject* convert(const std::vector<containedType>& v){
		python::list ret; FOREACH(const containedType& e, v) ret.append(e);
		return incref(ret.ptr());
	}
};
template<typename containedType>
struct custom_vector_from_seq{
	custom_vector_from_seq(){ converter::registry::push_back(&convertible,&construct,type_id<std::vector<containedType> >()); }
	static void* convertible(PyObject* obj_ptr){
		// the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data){
		 void* storage=((converter::rvalue_from_python_storage<std::vector<containedType> >*)(data))->storage.bytes;
		 new (storage) std::vector<containedType>();
		 std::vector<containedType>* v=(std::vector<containedType>*)(storage);
		 int l=PySequence_Size(obj_ptr); if(l<0) abort(); /*std::cerr<<"l="<<l<<"; "<<typeid(containedType).name()<<std::endl;*/ v->reserve(l); for(int i=0; i<l; i++) { v->push_back(extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		 data->convertible=storage;
	}
};



template<typename numT, int dim>
struct custom_numpyBoost_to_py{
	static PyObject* convert(numpy_boost<numT, dim> nb){
		return nb.py_ptr(); // handles incref internally
	}
};




using namespace boost::python;

BOOST_PYTHON_MODULE(_customConverters){
	// class_<std::vector<int> >("vecInt").def(indexing::container_suite<std::vector<int> >());
	custom_Se3r_from_seq(); to_python_converter<Se3r,custom_se3_to_tuple>();
	// Vector3i to python (not implemented the other way around yet)
#ifdef YADE_WM3
	custom_Vector3r_from_seq(); // Vector3r is wrapped, it is returned as a Vector3 instance; no to-python converter needed
	custom_vector3i_to_seq(); to_python_converter<Vector3i,custom_vector3i_to_seq>();
#endif
	// StrArrayMap (typedef for std::map<std::string,numpy_boost>) → python dictionary
	//custom_StrArrayMap_to_dict();
	// register from-python converter and to-python converter

	// register 2-way conversion between c++ vector and python homogeneous sequence (list/tuple) of corresponding type
	#define VECTOR_SEQ_CONV(Type) custom_vector_from_seq<Type>(); to_python_converter<std::vector<Type>, custom_vector_to_list<Type> >();
		VECTOR_SEQ_CONV(int);
		VECTOR_SEQ_CONV(Real);
		VECTOR_SEQ_CONV(Se3r);
		VECTOR_SEQ_CONV(Vector3r);
		VECTOR_SEQ_CONV(std::string);
		VECTOR_SEQ_CONV(shared_ptr<Body>);
		VECTOR_SEQ_CONV(shared_ptr<Engine>);
		VECTOR_SEQ_CONV(shared_ptr<Material>);
		VECTOR_SEQ_CONV(shared_ptr<Serializable>);
		VECTOR_SEQ_CONV(shared_ptr<BoundFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<InteractionGeometryFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<InteractionPhysicsFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<LawFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<IntrCallback>);
		VECTOR_SEQ_CONV(shared_ptr<BodyCallback>);
	#undef VECTOR_SEQ_CONV

	import_array();
	to_python_converter<numpy_boost<Real,1>, custom_numpyBoost_to_py<Real,1> >();
	to_python_converter<numpy_boost<Real,2>, custom_numpyBoost_to_py<Real,2> >();
	to_python_converter<numpy_boost<int,1>, custom_numpyBoost_to_py<int,1> >();
	to_python_converter<numpy_boost<int,2>, custom_numpyBoost_to_py<int,2> >();

}





