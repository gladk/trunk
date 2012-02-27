/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact.hpp"
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-snow/IstSnowLayersContact.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-snow/BssSnowGrain.hpp>

#include<yade/lib-base/yadeWm3Extra.hpp>

#ifndef MINIWM3
	#include <Wm3ApprPlaneFit3.h>
	#include <Wm3Plane3.h>
#endif
YADE_REQUIRE_FEATURE(geometricalmodel);
CREATE_LOGGER(Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact);

bool is_point_inside_box(Box* b, Vector3r P)
{
	return		std::abs(P[0]) < std::abs(b->extents[0])
		 &&	std::abs(P[1]) < std::abs(b->extents[1])
		 &&	std::abs(P[2]) < std::abs(b->extents[2]);
};


Vector3r find_cross_point(Box* m,Vector3r in,Vector3r out, Quaternionr q2, Quaternionr q1, Vector3r pos1, Vector3r pos2,int depth = 4)
{
	Vector3r mid((in+out)*0.5);
	if(depth == 0)
		return mid;
	if(is_point_inside_box(m, q2.Conjugate()*(q1 * mid + pos1-pos2)))
	{
		return find_cross_point(m,mid,out,q2,q1,pos1,pos2,depth-1);
	}
	else
	{
		return find_cross_point(m,in, mid,q2,q1,pos1,pos2,depth-1);
	}
}


bool Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact::go(
		const shared_ptr<Shape>& cm1,
		const shared_ptr<Shape>& cm2,
		const Se3r& se31,
		const Se3r& se32,
		const Vector3r& shift2,
		const bool& force,
		const shared_ptr<Interaction>& c)
{

	//Box* m1=static_cast<Box*>(cm1.get()), *m2=static_cast<BssSnowGrain*>(cm2.get());
	//std::cerr << "------------------- " << __FILE__ << "\n";

	if(cm1->getClassName() != std::string("Box") || cm2->getClassName() != std::string("BssSnowGrain"))
	{
		std::cerr << cm1->getClassName() << " " << cm2->getClassName() << "\n";
		std::cerr << "whooooooooops =66=\n";
		return false;
	}
				Box* m1=static_cast<Box*>(cm1.get());
				BssSnowGrain *m2=static_cast<BssSnowGrain*>(cm2.get());
				if(m1==0 || m2==0)
				{
					std::cerr << cm1->getClassName() << " " << cm2->getClassName() << "\n";
					std::cerr << "whooooooooops =6= " << __FILE__ << "\n"; 
					return false;
				}

	shared_ptr<IstSnowLayersContact> scm;
	if(c->interactionGeometry) scm=dynamic_pointer_cast<IstSnowLayersContact>(c->interactionGeometry);
	else { scm=shared_ptr<IstSnowLayersContact>(new IstSnowLayersContact()); c->interactionGeometry=scm; }
	
//	Vector3r normal_from_to=se32.position-se31.position;

	std::vector<Vector3r> cross_section;
	Vector3r    pos1(se31.position);
	Vector3r    pos2(se32.position);
	Quaternionr q1(se31.orientation);
	Quaternionr q2(se32.orientation);
	
	typedef std::pair<Vector3r, std::set<Vector3r> > t_edge;
//	const std::map<Vector3r, std::set<Vector3r> >& edges1 = m1->m_copy.get_edges_const_ref();
	const std::map<Vector3r, std::set<Vector3r> >& edges2 = m2->m_copy.get_edges_const_ref();
//	BOOST_FOREACH(const t_edge& e,edges1)
//	{
//			Vector3r v(e.first);
//			if(m2->m_copy.is_point_inside_polyhedron( q2.Conjugate()*(q1 * v + pos1-pos2)))
//			{
//				BOOST_FOREACH(Vector3r v2,e.second)
//				{
//					if( ! (m2->m_copy.is_point_inside_polyhedron( q2.Conjugate()*(q1 * v2 + pos1-pos2))))
//					{
//						Vector3r v_c(q1 * find_cross_point(m2,v,v2, q2,q1,pos1,pos2) + pos1);
//						cross_section.push_back(v_c);
//					}
//				}
//			}
//	}
	BOOST_FOREACH(const t_edge& e,edges2)
	{
			Vector3r v(e.first);
			if(is_point_inside_box(m1, q1.Conjugate()*(q2 * v + pos2-pos1)))
			{
				BOOST_FOREACH(Vector3r v2,e.second)
				{
					if( ! (is_point_inside_box(m1, q1.Conjugate()*(q2 * v2 + pos2-pos1))))
					{
						Vector3r v_c(q2 * find_cross_point(m1,v,v2, q1,q2,pos2,pos1) + pos2);
						cross_section.push_back(v_c);
					}
				}
			}
	}


	int id1 = c->getId1();
	//int id2 = c->getId2();

	if(!cross_section.empty())
	{
		// find the contact plane with least squre fitting from wm3 library
#ifndef MINIWM3
		Plane3<Real> plane(OrthogonalPlaneFit3 (cross_section.size(), &cross_section[0]));
		Vector3r N = plane.Normal;
		//const Vector3r N((normal_from_to.Dot(N) > 0) ? plane.Normal : 1.0*plane.Normal);
		Vector3r C = Vector3r(0,0,0);
		int i=0;
		if( std::abs(N[1]) > std::abs(N[0]) && std::abs(N[1]) > std::abs(N[2]) ) i=1;
		if( std::abs(N[2]) > std::abs(N[1]) && std::abs(N[2]) > std::abs(N[0]) ) i=2;
		C[i] = plane.Constant/N[i];

		// project all points onto C,N
		BOOST_FOREACH(Vector3r& v,cross_section) v -= N.Dot(v - C) * N;
		std::vector<Vector3r> cross_section_circumference(g.find_boundary(cross_section));
		if(cross_section_circumference.size() < 2)
		{
			std::cerr << "\n -------- WTF, wrong cross section " << __FILE__ << " --------- \n\n";
			return false;
		}
		Vector3r center_point(0,0,0);
		BOOST_FOREACH(Vector3r& v,cross_section_circumference) center_point +=v;
		center_point /= (Real)(cross_section_circumference.size());


		//Real min_dist(0),max_dist(0),count(0); // that's for penetration depth
		// find all points from layers projected on this plane
//		for(size_t i=0;i < m1->m_copy.slices.size();++i)
//		{
//			for(size_t j=0 ; j < m1->m_copy.slices[i].size() ; ++j)
//			{
//				Vector3r v(m1->m_copy.slices[i][j]);
//				v = q1 * v + pos1;
//				Real point_plane_distance = N.Dot(v - C);
//				
//				if(m1->depths.find(id2) == m1->depths.end())
//					m1->depths[id2] = std::set< depth_one >();
//
//				std::set<depth_one>::iterator it(m1->depths[id2].find(depth_one(0,i,j)));
//				if(it != m1->depths[id2].end()) it->set_current_depth(point_plane_distance);
//
//				if(std::abs(point_plane_distance) < m1->radius*0.2 
//					&& is_point_inside_cross_section(v,point_plane_distance,cross_section_circumference,center_point,N))
//				{
//					m1->depths[id2].insert(depth_one(point_plane_distance , i ,j));
//		
//					min_dist = std::min(min_dist,point_plane_distance);
//					max_dist = std::max(max_dist,point_plane_distance);
//				}
//			}
//		}
		for(size_t i=0;i < m2->m_copy.slices.size();++i)
		{
			for(size_t j=0 ; j < m2->m_copy.slices[i].size() ; ++j)
			{
				Vector3r v(m2->m_copy.slices[i][j]);
				v = q2 * v + pos2;
				Real point_plane_distance = N.Dot(v - C);
				
				if(m2->depths.find(id1) == m2->depths.end())
					m2->depths[id1] = std::set< depth_one >();

				std::set<depth_one>::iterator it(m2->depths[id1].find(depth_one(0,i,j)));
				if(it != m2->depths[id1].end()) it->set_current_depth(point_plane_distance);
				
				if(std::abs(point_plane_distance) < m2->radius*0.2
					&& g.is_point_inside_cross_section(v,point_plane_distance,cross_section_circumference,center_point,N))
				{
					m2->depths[id1].insert(depth_one(point_plane_distance , i ,j));
		
					//std::set<depth_one>::iterator it(m2->depths[id1].find(depth_one(0,i,j)));
					//min_dist = std::min(min_dist,point_plane_distance - it->original_depth);
					//max_dist = std::max(max_dist,point_plane_distance - it->original_depth);
					//count++;
				}
			}
		}

//		Real penetrationDepth = max_dist - min_dist;
//		scm->contactPoint = center_point;
//		scm->normal = N;
//		//scm->normal = (normal_from_to.Dot(N) > 0) ? N : -1.0*N; //N;
//		scm->penetrationDepth = penetrationDepth;
//
////		if(m1->sphere_depth.find(id2) == m1->sphere_depth.end())
////			m1->sphere_depth[id2] = penetrationDepth;
//
//		if(m2->sphere_depth.find(id1) == m2->sphere_depth.end())
//			m2->sphere_depth[id1] = penetrationDepth;
//
////		Real d1 = m1->sphere_depth[id2];
//		Real d2 = m2->sphere_depth[id1];
////		if(d1 != d2)
////			std::cerr << "bad initial penetration?\n";
//
//		penetrationDepth -= d2;
//		scm->penetrationDepth=penetrationDepth;
//
//		scm->radius1=m2->radius*2.0;
//		scm->radius2=m2->radius;

		//std::cerr << "id1: " << id1 << " " << m2->depths[id1].size() << " ---- " << "id2: " << id2 << " " << m1->depths[id2].size() << "\n";

//	std::cerr << __FILE__ << " " << scm->getClassName() << "\n";
		
		// FIXME: ScGeom (components that are not from IstSnowLayersContact itself) are calculated by "parent" class
		//        the penetration depth, contact point and normal. I couldn't make stimulation to be stable without this.
//FIXME//		bool old_n = c->isNew;
//FIXME//		c->isNew=false;
		//ggg.assist=true;
///////////////??????		c->init();
		bool res = ggg.go(cm1,cm2,se31,se32,c);
//FIXME//		c->isNew=old_n;
		return res;

		//return true;
#else
		LOG_FATAL("Using miniWm3; recompile with full Wm3 support to make snow fully functional.");
		throw runtime_error("full wm3 required (message above).");
#endif
	}
//	if(! m1->depths[id2].empty()) m1->depths[id2].clear();
	if(! m2->depths[id1].empty()) m2->depths[id1].clear();
	return false;
}


bool Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact::goReverse(	const shared_ptr<Shape>& cm1,
						const shared_ptr<Shape>& cm2,
						const Se3r& se31,
						const Se3r& se32,
						const Vector3r& shift2,
						const bool& force,
						const shared_ptr<Interaction>& c)
{
//	bool result = go(cm2,cm1,se32,se31,c);
//	if(result)
//	{
//		shared_ptr<ScGeom> scm;
//		if(c->interactionGeometry) scm=YADE_PTR_CAST<ScGeom>(c->interactionGeometry);
//		else { std::cerr << "whooooooooops_2!" << __FILE__ << "\n"; return false; }
//		scm->normal *= -1.0;
//		std::swap(scm->radius1,scm->radius2);
//	}
//	return result;


	//Box* m1=static_cast<Box*>(cm1.get()), *m2=static_cast<BssSnowGrain*>(cm2.get());
////	std::cerr << "----- reverse ----- " << __FILE__ << "\n";

	if(cm1->getClassName() != std::string("BssSnowGrain") || cm2->getClassName() != std::string("Box"))
	{
		std::cerr << cm1->getClassName() << " " << cm2->getClassName() << "\n";
		std::cerr << "whooooooooops =11=\n";
		return false;
	}
				Box* m2=static_cast<Box*>(cm2.get());
				BssSnowGrain *m1=static_cast<BssSnowGrain*>(cm1.get());
				if(m1==0 || m2==0)
				{
					std::cerr << cm1->getClassName() << " " << cm2->getClassName() << "\n";
					std::cerr << "\n whooooooooops =1= " << __FILE__ << "\n"; 
					return false;
				}

	shared_ptr<IstSnowLayersContact> scm;
	if(c->interactionGeometry) scm=dynamic_pointer_cast<IstSnowLayersContact>(c->interactionGeometry);
	else { scm=shared_ptr<IstSnowLayersContact>(new IstSnowLayersContact()); c->interactionGeometry=scm; }

	//Vector3r normal_from_to=se32.position-se31.position;
	
	std::vector<Vector3r> cross_section;
	Vector3r    pos1(se31.position);
	Vector3r    pos2(se32.position);
	Quaternionr q1(se31.orientation);
	Quaternionr q2(se32.orientation);
	
	typedef std::pair<Vector3r, std::set<Vector3r> > t_edge;
	const std::map<Vector3r, std::set<Vector3r> >& edges1 = m1->m_copy.get_edges_const_ref();
//	const std::map<Vector3r, std::set<Vector3r> >& edges2 = m2->m_copy.get_edges_const_ref();
	BOOST_FOREACH(const t_edge& e,edges1)
	{
			Vector3r v(e.first);
			if(is_point_inside_box(m2, q2.Conjugate()*(q1 * v + pos1-pos2)))
			{
				BOOST_FOREACH(Vector3r v2,e.second)
				{
					if( ! (is_point_inside_box(m2, q2.Conjugate()*(q1 * v2 + pos1-pos2))))
					{
						Vector3r v_c(q1 * find_cross_point(m2,v,v2, q2,q1,pos1,pos2) + pos1);
						cross_section.push_back(v_c);
					}
				}
			}
	}
//	BOOST_FOREACH(const t_edge& e,edges2)
//	{
//			Vector3r v(e.first);
//			if(is_point_inside_box(m1, q1.Conjugate()*(q2 * v + pos2-pos1)))
//			{
//				BOOST_FOREACH(Vector3r v2,e.second)
//				{
//					if( ! (is_point_inside_box(m1, q1.Conjugate()*(q2 * v2 + pos2-pos1))))
//					{
//						Vector3r v_c(q2 * find_cross_point(m1,v,v2, q1,q2,pos2,pos1) + pos2);
//						cross_section.push_back(v_c);
//					}
//				}
//			}
//	}


	//int id1 = c->getId1();
	int id2 = c->getId2();

	if(!cross_section.empty())
	{
#ifndef MINIWM3
		// find the contact plane with least squre fitting from wm3 library
		Plane3<Real> plane(OrthogonalPlaneFit3 (cross_section.size(), &cross_section[0]));
		Vector3r N = plane.Normal;
		//const Vector3r N((normal_from_to.Dot(N) > 0) ? plane.Normal : 1.0*plane.Normal);
		Vector3r C = Vector3r(0,0,0);
		int i=0;
		if( std::abs(N[1]) > std::abs(N[0]) && std::abs(N[1]) > std::abs(N[2]) ) i=1;
		if( std::abs(N[2]) > std::abs(N[1]) && std::abs(N[2]) > std::abs(N[0]) ) i=2;
		C[i] = plane.Constant/N[i];

		// project all points onto C,N
		BOOST_FOREACH(Vector3r& v,cross_section) v -= N.Dot(v - C) * N;
		std::vector<Vector3r> cross_section_circumference(g.find_boundary(cross_section));
		if(cross_section_circumference.size() < 2)
		{
			std::cerr << "\n -------- WTF, wrong cross section " << __FILE__ << " --------- \n\n";
			return false;
		}
		Vector3r center_point(0,0,0);
		BOOST_FOREACH(Vector3r& v,cross_section_circumference) center_point +=v;
		center_point /= (Real)(cross_section_circumference.size());


		//Real min_dist(0),max_dist(0),count(0); // that's for penetration depth
		// find all points from layers projected on this plane
		for(size_t i=0;i < m1->m_copy.slices.size();++i)
		{
			for(size_t j=0 ; j < m1->m_copy.slices[i].size() ; ++j)
			{
				Vector3r v(m1->m_copy.slices[i][j]);
				v = q1 * v + pos1;
				Real point_plane_distance = N.Dot(v - C);
				
				if(m1->depths.find(id2) == m1->depths.end())
					m1->depths[id2] = std::set< depth_one >();

				std::set<depth_one>::iterator it(m1->depths[id2].find(depth_one(0,i,j)));
				if(it != m1->depths[id2].end()) it->set_current_depth(point_plane_distance);

				if(std::abs(point_plane_distance) < m1->radius*0.2 
					&& g.is_point_inside_cross_section(v,point_plane_distance,cross_section_circumference,center_point,N))
				{
					m1->depths[id2].insert(depth_one(point_plane_distance , i ,j));
		
					//std::set<depth_one>::iterator it(m1->depths[id2].find(depth_one(0,i,j)));
					//min_dist = std::min(min_dist,point_plane_distance - it->original_depth);
					//max_dist = std::max(max_dist,point_plane_distance - it->original_depth);
					//count++;
				}
			}
		}
//		for(size_t i=0;i < m2->m_copy.slices.size();++i)
//		{
//			for(size_t j=0 ; j < m2->m_copy.slices[i].size() ; ++j)
//			{
//				Vector3r v(m2->m_copy.slices[i][j]);
//				v = q2 * v + pos2;
//				Real point_plane_distance = N.Dot(v - C);
//				
//				if(m2->depths.find(id1) == m2->depths.end())
//					m2->depths[id1] = std::set< depth_one >();
//
//				std::set<depth_one>::iterator it(m2->depths[id1].find(depth_one(0,i,j)));
//				if(it != m2->depths[id1].end()) it->set_current_depth(point_plane_distance);
//				
//				if(std::abs(point_plane_distance) < m2->radius*0.2
//					&& g.is_point_inside_cross_section(v,point_plane_distance,cross_section_circumference,center_point,N))
//				{
//					m2->depths[id1].insert(depth_one(point_plane_distance , i ,j));
//		
//					min_dist = std::min(min_dist,point_plane_distance);
//					max_dist = std::max(max_dist,point_plane_distance);
//				}
//			}
//		}

//		Real penetrationDepth = max_dist - min_dist;
//		scm->contactPoint = center_point;
//		scm->normal = N;
//		//scm->normal = (normal_from_to.Dot(N) > 0) ? N : -1.0*N; //N;
//		scm->penetrationDepth = penetrationDepth;
//
//		if(m1->sphere_depth.find(id2) == m1->sphere_depth.end())
//			m1->sphere_depth[id2] = penetrationDepth;
//
////		if(m2->sphere_depth.find(id1) == m2->sphere_depth.end())
////			m2->sphere_depth[id1] = penetrationDepth;
//
//		Real d1 = m1->sphere_depth[id2];
////		Real d2 = m2->sphere_depth[id1];
////		if(d1 != d2)
////			std::cerr << "bad initial penetration?\n";
//
//		penetrationDepth -= d1;
//		scm->penetrationDepth=penetrationDepth;
//
//		scm->radius1=m1->radius;
//		scm->radius2=m1->radius*2.0;
//
//		//std::cerr << "id1: " << id1 << " " << m2->depths[id1].size() << " ---- " << "id2: " << id2 << " " << m1->depths[id2].size() << "\n";
		
//	std::cerr << __FILE__ << " " << scm->getClassName() << "\n";
		
		// FIXME: ScGeom (components that are not from IstSnowLayersContact itself) are calculated by "parent" class
		//        the penetration depth, contact point and normal. I couldn't make stimulation to be stable without this.
//FIXME//		bool old_n = c->isNew;
//FIXME//		c->isNew=false;
		//ggg.assist=true;
		bool res = ggg.goReverse(cm1,cm2,se31,se32,c);
//FIXME//		c->isNew=old_n;
		return res;

		//return true;
#else
		LOG_FATAL("Using miniWm3; recompile with full Wm3 support to make snow folly functional.");
		throw runtime_error("full wm3 required (message above).");
#endif
	}
	if(! m1->depths[id2].empty()) m1->depths[id2].clear();
//	if(! m2->depths[id1].empty()) m2->depths[id1].clear();
	return false;

}

YADE_PLUGIN((Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact));

YADE_REQUIRE_FEATURE(PHYSPAR);

