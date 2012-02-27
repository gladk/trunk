# encoding: utf-8
# 2009 © Václav Šmilauer <eudoxos@arcig.cz>
import os, os.path
from pyplusplus import module_builder
from pyplusplus.decl_wrappers import *
from pyplusplus.module_builder import call_policies
from pyplusplus import function_transformers as FT




#Creating an instance of class that will help you to expose your declarations
mb = module_builder.module_builder_t( [os.path.abspath("miniWm3Wrap-toExpose.hpp")]
                                      , working_directory=r"."
                                      , include_paths=['../../lib','../../lib/miniWm3','../../../build-trunk-dbg/include/yade-trunk','/usr/include']
                                      , define_symbols=['USE_MINIWM3','YADE_EIGEN'] )
# exclude exerything first
mb.decls().exclude()
# include what is in Wm3
mb.decls(lambda d: ('Wm3' in d.location.file_name) if d.location else False).include()


## things we want to exclude:
# we don't need system things from wm3
mb.decls(lambda decl: 'System' in decl.name).exclude() 
# do not expose instances we created just to instantiate templates
mb.decls(lambda decl: decl.name.startswith('noexpose')).exclude() 
# exclude casting operators
mb.casting_operators().exclude()
# exclude functions&operators returning double& (not representable)
mb.member_functions(return_type='double &').exclude()
mb.global_ns.operators(return_type='double &').exclude() # global_ns is a bug workaround: http://www.mail-archive.com/cplusplus-sig@python.org/msg00730.html
## exclude everything (member functions, operators, ...) taking or returning types we do not wrap
## mb.decls(lambda d: 'Matrix2' in str(d) or 'Vector4' in str(d)).exclude()
# exclude operator[] since we implement __getitem__/__setitem__ by ourselves
mb.member_operators(lambda o: o.symbol=='[]').exclude()
# exclude From/To Euler angles, since we don't need it
mb.member_functions(lambda f: 'Euler' in str(f)).exclude()
## exclude free functions
mb.decls(lambda d: d.name in ('componentMaxVector','componentMinVector','componentSum','diagDiv','diagMult','quaternionFromAxes','quaternionToAxes','quaternionToEulerAngles','quaterniontoGLMatrix','unitVectorsAngle')).exclude()


## exclude what is not in the wm3-eigen glue
if 1:
	undef="""
Adjoint DiagonalTimes GetColumnMajor Orthonormalize QDUDecomposition QForm SingularValueComposition SingularValueDecomposition Slerp TimesDiagonal TimesTranspose TransposeTimes MakeTensorProduct
DecomposeSwingTimesTwist DecomposeTwistTimesSwing Exp FromRotationMatrix Intermediate Log SlerpExtraSpins Squad
ComputeExtremes DotPerp GenerateOrthonormalBasis GetBarycentrics Orthonormalize Perp UnitPerp
ComputeExtremes GenerateOrthonormalBasis GetBarycentrics Orthonormalize
MakeDiagonal MakeZero MakeIdentity""".split()
	mb.member_functions(lambda d: d.name in undef).exclude()


# register manual wraps
v2,v3,q4,m3=mb.class_(lambda d: d.name=='Vector2<double>'),mb.class_(lambda d: d.name=='Vector3<double>'),mb.class_(lambda d: d.name=='Quaternion<double>'),mb.class_(lambda d: d.name=='Matrix3<double>')
v2.add_registration_code('def("__len__",&::Vector2_len)   .staticmethod("__len__").def("__setitem__",&::Vector2_set_item)   .def("__getitem__",&::Vector2_get_item)   .def("__str__",&::Vector2_str)   .def("__repr__",&::Vector2_str)')
v3.add_registration_code('def("__len__",&::Vector3_len)   .staticmethod("__len__").def("__setitem__",&::Vector3_set_item)   .def("__getitem__",&::Vector3_get_item)   .def("__str__",&::Vector3_str)   .def("__repr__",&::Vector3_str)')
q4.add_registration_code('def("__len__",&::Quaternion_len).staticmethod("__len__").def("__setitem__",&::Quaternion_set_item).def("__getitem__",&::Quaternion_get_item).def("__str__",&::Quaternion_str).def("__repr__",&::Quaternion_str)')
# expose both [0…2,0…2] access for matrices as well as row-major [0…8]
m3.add_registration_code('def("__len__",&::Matrix3_len).staticmethod("__len__")   .def("__setitem__",&::Matrix3_set_item)   .def("__getitem__",&::Matrix3_get_item)   .def("__str__",&::Matrix3_str)   .def("__repr__",&::Matrix3_str)  /* extras for matrices */ .def("__setitem__",&::Matrix3_set_item_linear).def("__getitem__",&::Matrix3_get_item_linear)   ')

## workarounds for def_readonly on static members
# disabled, doesn't work on types (e.g. Vector3.ZERO), must pass instance (e.g. Vector3().ZERO)
if 0:
	mb.decls(lambda d: 'UNIT_' in str(d) or 'ZERO' in str(d) or 'IDENTITY' in str(d) or 'ONE' in str(d)).exclude()
	v2.add_registration_code('.'.join('add_property("%s",::Vector2r_%s)'%(prop,prop) for prop in ('ZERO','UNIT_X','UNIT_Y','ONE')))
	v3.add_registration_code('.'.join('add_property("%s",::Vector3r_%s)'%(prop,prop) for prop in ('ZERO','UNIT_X','UNIT_Y','UNIT_Z','ONE')))
	q4.add_registration_code('.'.join('add_property("%s",::Quaternionr_%s)'%(prop,prop) for prop in ('ZERO','IDENTITY')))
	m3.add_registration_code('.'.join('add_property("%s",::Matrix3r_%s)'%(prop,prop) for prop in ('ZERO','IDENTITY')))

# arg0 and arg1 will be returned as tuple, instead of the references ToAxisAngle takes
mb.member_functions('ToAxisAngle').add_transformation(FT.output(0),FT.output(1))
# quaternion operations returning itself
mb.member_functions(lambda mf: 'Quaternion<double> &' in str(mf.return_type)).call_policies=call_policies.return_self()

mb.add_registration_code("custom_Vector3r_from_tuple();")

mb.member_functions(lambda mf: 'Matrix3<double> &' in str(mf.return_type)).call_policies=call_policies.return_self()
mb.decls(lambda decl: decl.name.startswith('ToEulerAngles')).add_transformation(FT.output(0),FT.output(1),FT.output(2))
if 0:
	#mb.member_functions(return_type='Wm3::Matrix3<double> &').call_policies=call_policies.return_self
	mb.member_functions('ToAngle').add_transformation(FT.output(0)) # Matrix2

mb.build_code_creator(module_name='miniWm3Wrap')
mb.code_creator.add_include("miniWm3Wrap-funcs.ipp")
mb.write_module('miniWm3Wrap.cpp')

# cleaup unneeded generated files
toClean='exposed_decl.pypp.txt','miniWm3Wrap.cpp~','named_tuple.py'
for f in toClean:
	if os.path.exists(f):
		try:
			os.remove(f)
		except IOError: pass
## remove absolute path from the generated file (ugly, oh well)
os.system(r"perl -pi -e 's@^#include\s*\"/.*/(.*)\"\s*$@#include \"\1\"\n@' miniWm3Wrap.cpp")
os.system(r"perl -pi -e 's@(::)?Wm3::@@g' miniWm3Wrap.cpp")
os.system(r"perl -pi -e 's@^.*double const \*.*$@@' miniWm3Wrap.cpp")
os.system(r"perl -pi -e 's@^.*implicitly_convertible.*$@@' miniWm3Wrap.cpp")

