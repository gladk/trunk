// (C) 2004 by Olivier Galizzi <olivier.galizzi@imag.fr>

#include"PhysicalParameters.hpp"
#include<boost/foreach.hpp>
#include<stdexcept>
#ifndef FOREACH
	#define FOREACH BOOST_FOREACH
#endif

std::vector<std::string> PhysicalParameters::blockedDOFs_vec_get() const {
	std::vector<std::string> ret;
	#define _SET_DOF(DOF_ANY,str) if((blockedDOFs & PhysicalParameters::DOF_ANY)!=0) ret.push_back(str);
	_SET_DOF(DOF_X,"x"); _SET_DOF(DOF_Y,"y"); _SET_DOF(DOF_Z,"z"); _SET_DOF(DOF_RX,"rx"); _SET_DOF(DOF_RY,"ry"); _SET_DOF(DOF_RZ,"rz");
	#undef _SET_DOF
	return ret;
}
void PhysicalParameters::blockedDOFs_vec_set(const std::vector<std::string>& dofs){
	FOREACH(const std::string s, dofs){
		#define _GET_DOF(DOF_ANY,str) if(s==str) { blockedDOFs|=PhysicalParameters::DOF_ANY; continue; }
		_GET_DOF(DOF_X,"x"); _GET_DOF(DOF_Y,"y"); _GET_DOF(DOF_Z,"z"); _GET_DOF(DOF_RX,"rx"); _GET_DOF(DOF_RY,"ry"); _GET_DOF(DOF_RZ,"rz");
		#undef _GET_DOF
		throw std::invalid_argument("Invalid  DOF specification `"+s+"', must be ∈{x,y,z,rx,ry,rz}.");
	}
}

