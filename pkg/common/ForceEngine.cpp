// 2004 © Janek Kozicki <cosurgi@berlios.de> 
// 2009 © Václav Šmilauer <eudoxos@arcig.cz> 


#include"ForceEngine.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/LinearInterpolate.hpp>
#include<yade/pkg-dem/Shop.hpp>

#include<yade/core/IGeom.hpp>
#include<yade/core/IPhys.hpp>

YADE_PLUGIN((ForceEngine)(InterpolatingDirectedForceEngine));

void ForceEngine::action(){
	FOREACH(Body::id_t id, ids){
		if (!(scene->bodies->exists(id))) continue;
		scene->forces.addForce(id,force);
	}
}

void InterpolatingDirectedForceEngine::action(){
	Real virtTime=wrap ? Shop::periodicWrap(scene->time,*times.begin(),*times.rbegin()) : scene->time;
	direction.normalize(); 
	force=linearInterpolate<Real,Real>(virtTime,times,magnitudes,_pos)*direction;
	ForceEngine::action();
}


