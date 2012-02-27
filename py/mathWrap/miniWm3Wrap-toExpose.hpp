#include<base/Math.hpp>
//#include<boost/lexical_cast.hpp>
#include<string>

// define python names for these classes
// see http://language-binding.net/pyplusplus/documentation/how_to/hints.html#pyplusplus-aliases-namespace
namespace pyplusplus{ namespace aliases {
	typedef Vector3r Vector3;
	typedef Vector2r Vector2;
	typedef Quaternionr Quaternion;
	typedef Matrix3r Matrix3;
}}
// instantiate those types so that they are exposed
// but name them noexpose_*, as we ask for exclusion of such things in the wrapper script
Vector3r noexpose_v3r;
Vector2r noexpose_v2r;
Quaternionr noexpose_qr;
Matrix3r noexpose_m3r;
//Matrix2r noexpose_m2r;

