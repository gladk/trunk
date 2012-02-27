// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2004 Janek Kozicki <cosurgi@berlios.de>
// © 2007 Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>

#include"Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact.hpp"
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-snow/IstSnowLayersContact.hpp>
#include<yade/pkg-snow/BssSnowGrain.hpp>
#include<yade/pkg-common/Sphere.hpp>

#include<yade/lib-base/yadeWm3Extra.hpp>
#include<yade/core/Omega.hpp>

#ifndef MINIWM3
	#include <Wm3ApprPlaneFit3.h>
	#include <Wm3Plane3.h>
#endif

YADE_REQUIRE_FEATURE(geometricalmodel);

bool is_point_orthogonally_projected_on_triangle_from_both_sides(const Vector3r& a,const Vector3r& b,const Vector3r c,Vector3r& N,Vector3r& P,Real point_plane_distance)
{
	Vector3r d(P - point_plane_distance*N);
	// now check if the point (when projected on a plane) is within triangle a,b,c
	// it could be faster with methods from http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm
	// but I don't understand them, so I prefer to use the method which I derived myself
	Vector3r c1((a - b).Cross(d - a)); // a,b,c can be clockwise or counterclockwise, I don't know here
	Vector3r c2((c - a).Cross(d - c)); // 'd' will be inside if it preserves clockwiseness
	Vector3r c3((b - c).Cross(d - b));
	if(c1.Dot(N) >= 0 && c2.Dot(N) >= 0 && c3.Dot(N) >= 0)
		return true;
	if(c1.Dot(N) <= 0 && c2.Dot(N) <= 0 && c3.Dot(N) <= 0) // so dot producs must all have the same sign
		return true;
	return false;
};


Vector3r find_cross_point(BssSnowGrain* m,Vector3r in,Vector3r out, Quaternionr q2, Quaternionr q1, Vector3r pos1, Vector3r pos2,int depth = 4)
{
	Vector3r mid((in+out)*0.5);
	if(depth == 0)
		return mid;
	if(m->m_copy.is_point_inside_polyhedron( q2.Conjugate()*(q1 * mid + pos1-pos2)))
	{
		return find_cross_point(m,mid,out,q2,q1,pos1,pos2,depth-1);
	}
	else
	{
		return find_cross_point(m,in, mid,q2,q1,pos1,pos2,depth-1);
	}
}

bool Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact::is_point_inside_cross_section(Vector3r v,Real point_plane_distance,const std::vector<Vector3r>& cross_section_circumference,Vector3r center,Vector3r N)
{
	for(size_t i = 0 ; i<cross_section_circumference.size() - 1 ; i++)
	{
		const Vector3r& a(cross_section_circumference[i]);
		const Vector3r& b(cross_section_circumference[i+1]);
		if(is_point_orthogonally_projected_on_triangle_from_both_sides(a,b,center,N,v,point_plane_distance))
			return true;
	}
	return false;
}

/// all points must share the same plane!
// points that are "inside" circumference of set of points are removed from the set of all points
std::vector<Vector3r> Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact::find_boundary(std::vector<Vector3r> points)
{
	std::vector<Vector3r> res;
	if(points.size()>3)
	{
		Vector3r V(points[0].Cross(points[1]));
		V.Normalize();
		for(size_t I=0; I<points.size()-1 ; I++)
		for(size_t J=I+1; J<points.size() ; J++)
		//BOOST_FOREACH(Vector3r& a,points)
		//BOOST_FOREACH(Vector3r& b,points) // line a-b
		{
			Vector3r& a(points[I]);
			Vector3r& b(points[J]);
			//if(a != b)
			{
				Vector3r N(V.Cross(b-a));
				N.Normalize();
				int side(0);
				bool ok(true);
				BOOST_FOREACH(Vector3r& P,points) // checking side of c
				{
					if(P !=a && P != b)
					{
						Real point_plane_distance = N.Dot(P - a);
						int sign = ((point_plane_distance > 0) ? (1) : (-1));
						if(side == 0)
							side = sign;
						else
						{
							if(side != sign)
								ok = false;
						}
						if(!ok)
							break;
					}
				}
				if(ok)
				{
					res.push_back(a);
					res.push_back(b);
				}
			}
		}
	}
	else
		res=points;
	return res;
}

Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact::Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact()
{
}


bool Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact::go(	const shared_ptr<Shape>& cm1,
							const shared_ptr<Shape>& cm2,
							const Se3r& se31,
							const Se3r& se32,
							const Vector3r& shift2,
							const bool& force,
							const shared_ptr<Interaction>& c)
{
//	bool result = g.go(cm1,cm2,se31,se32,c);
//	std::cerr << "------------------- " << __FILE__ << "\n";
//	return result;
	//if(box)
	//{
	//	std::cerr << "----------------- box -- Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact\n";
	//	return false;
	//}

	BssSnowGrain* m1=static_cast<BssSnowGrain*>(cm1.get()), *m2=static_cast<BssSnowGrain*>(cm2.get());

	shared_ptr<IstSnowLayersContact> scm;
	if(c->interactionGeometry) scm=YADE_PTR_CAST<IstSnowLayersContact>(c->interactionGeometry);
	else { scm=shared_ptr<IstSnowLayersContact>(new IstSnowLayersContact()); c->interactionGeometry=scm; }

	//std::list<Vector3r> nodes;
	//std::list<Vector3r> inside_nodes;

	//Vector3r normal_from_to=se32.position-se31.position;
	
	std::vector<Vector3r> cross_section;
	Vector3r    pos1(se31.position);
	Vector3r    pos2(se32.position);
	Quaternionr q1(se31.orientation);
	Quaternionr q2(se32.orientation);
	
	typedef std::pair<Vector3r, std::set<Vector3r> > t_edge;
	const std::map<Vector3r, std::set<Vector3r> >& edges1 = m1->m_copy.get_edges_const_ref();
	const std::map<Vector3r, std::set<Vector3r> >& edges2 = m2->m_copy.get_edges_const_ref();
	BOOST_FOREACH(const t_edge& e,edges1)
	{
			Vector3r v(e.first);
			if(m2->m_copy.is_point_inside_polyhedron( q2.Conjugate()*(q1 * v + pos1-pos2)))
			{
				BOOST_FOREACH(Vector3r v2,e.second)
				{
					if( ! (m2->m_copy.is_point_inside_polyhedron( q2.Conjugate()*(q1 * v2 + pos1-pos2))))
					{
						Vector3r v_c(q1 * find_cross_point(m2,v,v2, q2,q1,pos1,pos2) + pos1);
						cross_section.push_back(v_c);
					}
				}
//				Vector3r pos(q1 * v + pos1);
//				inside_nodes.push_back(pos);
			}
	}
	BOOST_FOREACH(const t_edge& e,edges2)
	{
			Vector3r v(e.first);
			if(m1->m_copy.is_point_inside_polyhedron( q1.Conjugate()*(q2 * v + pos2-pos1)))
			{
				BOOST_FOREACH(Vector3r v2,e.second)
				{
					if( ! (m1->m_copy.is_point_inside_polyhedron( q1.Conjugate()*(q2 * v2 + pos2-pos1))))
					{
						Vector3r v_c(q2 * find_cross_point(m1,v,v2, q1,q2,pos2,pos1) + pos2);
						cross_section.push_back(v_c);
					}
				}
//				Vector3r pos(q2 * v + pos2);
//				inside_nodes.push_back(pos);
			}
	}


	int id1 = c->getId1();
	int id2 = c->getId2();

	#ifdef MINIWM3
		throw runtime_error(__FILE__ ": Your build uses miniWm3, which doesn't have OrthogonalPlaneFit3 function. Aborting.");
	#else
	if(!cross_section.empty())
	{
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
		std::vector<Vector3r> cross_section_circumference(find_boundary(cross_section));
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
					&& is_point_inside_cross_section(v,point_plane_distance,cross_section_circumference,center_point,N))
				{
					m1->depths[id2].insert(depth_one(point_plane_distance , i ,j));

				}
				//it = m1->depths[id2].find(depth_one(0,i,j));
				//min_dist += point_plane_distance - it->original_depth;
				////max_dist = std::max(max_dist,point_plane_distance - it->original_depth);
				//count++;
			}
		}
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
					&& is_point_inside_cross_section(v,point_plane_distance,cross_section_circumference,center_point,N))
				{
					m2->depths[id1].insert(depth_one(point_plane_distance , i ,j));
		
				}
				//it=m2->depths[id1].find(depth_one(0,i,j));
				//min_dist += point_plane_distance - it->original_depth;
				////max_dist = std::max(max_dist,point_plane_distance - it->original_depth);
				//count++;
			}
		}

		//if(count == 0)	std::cerr << ".............\n";

	//	Real penetrationDepth= (m1->radius+m2->radius) - (se32.position-se31.position).Length();
	//	//Real penetrationDepth = /*max_dist -*/ -1.0* min_dist/count;
	//	scm->contactPoint = center_point;
	//	scm->normal = N;
	//	//scm->normal = (normal_from_to.Dot(N) > 0) ? N : -1.0*N; //N;
	//	scm->penetrationDepth = penetrationDepth;
        //
	//	if(m1->sphere_depth.find(id2) == m1->sphere_depth.end())
	//		m1->sphere_depth[id2] = penetrationDepth;
	//	if(m2->sphere_depth.find(id1) == m2->sphere_depth.end())
	//		m2->sphere_depth[id1] = penetrationDepth;
        //
	//	Real d1 = m1->sphere_depth[id2];
	//	Real d2 = m2->sphere_depth[id1];
	//	if(d1 != d2)
	//		std::cerr << "bad initial penetration?\n";
        //
	//	penetrationDepth -= d1;
	//	scm->penetrationDepth=penetrationDepth;
        //
	//	scm->radius1=m1->radius;
	//	scm->radius2=m2->radius;

		//std::cerr << "id1: " << id1 << " " << m2->depths[id1].size() << " ---- " << "id2: " << id2 << " " << m1->depths[id2].size() << "\n";



	//{
	//	Vector3r n2(se32.position-se31.position);
	//	Real penetrationDepth=m1->radius+m2->radius - n2.Normalize();
	//	scm->contactPoint=se31.position+(m1->radius-0.5*penetrationDepth)*n2;//0.5*(pt1+pt2);
	//	scm->normal=n2;
	//	//scm->contactPoint=center_point;
	//	//scm->normal=N;
	//	scm->penetrationDepth=penetrationDepth;

	//	if(m1->sphere_depth.find(id2) == m1->sphere_depth.end())
	//		m1->sphere_depth[id2] = penetrationDepth;
	//	if(m2->sphere_depth.find(id1) == m2->sphere_depth.end())
	//		m2->sphere_depth[id1] = penetrationDepth;

	//	Real d1 = m1->sphere_depth[id2];
	//	Real d2 = m2->sphere_depth[id1];
	//	if(d1 != d2)
	//		std::cerr << "bad initial penetration?\n";

	//	penetrationDepth -= d1;
	//	scm->penetrationDepth=penetrationDepth;

	//	scm->radius1=m1->radius;
	//	scm->radius2=m2->radius;
	//}


		// FIXME: ScGeom (components that are not from IstSnowLayersContact itself) are calculated by "parent" class
		//        the penetration depth, contact point and normal. I couldn't make stimulation to be stable without this.
		//g.assist=true;
		return g.go(cm1,cm2,se31,se32,c);
		
		////bool old_n = c->isNew;
		////c->isNew=false;
		////g.assist=true;
		/////*bool res = */g.go(cm1,cm2,se31,se32,c);
		////c->isNew=old_n;

		//return true;
	}
	#endif
	if(! m1->depths[id2].empty()) m1->depths[id2].clear();
	if(! m2->depths[id1].empty()) m2->depths[id1].clear();
	return false;






/*	
	for(size_t i=0;i < m1->m_copy.slices.size();++i)
	{
		for(size_t j=0 ; j < m1->m_copy.slices[i].size() ; ++j)
		{
			Vector3r v(m1->m_copy.slices[i][j]);
			if(m2->m_copy.is_point_inside_polyhedron( q2.Conjugate()*(q1 * v + pos1-pos2)))
			{
				nodes.push_back(q1 * v + pos1);
			}
		}
	}
	for(size_t i=0;i < m2->m_copy.slices.size();++i)
	{
		for(size_t j=0 ; j < m2->m_copy.slices[i].size() ; ++j)
		{
			Vector3r v(m2->m_copy.slices[i][j]);
			if(m1->m_copy.is_point_inside_polyhedron( q1.Conjugate()*(q2 * v + pos2-pos1)))
			{
				nodes.push_back(q2 * v + pos2);
			}
		}
	}
*/
/*
	if(!nodes.empty())
	{
		// plane of contact will be C (average pos of all points),N (normal of the contact - between two grains)
		Vector3r C(0,0,0);
		Real count(0);
		BOOST_FOREACH(Vector3r& v,nodes) C+=v,++count;
		C /= count;
		Vector3r N = pos2 - pos1;
		N.Normalize();
		// project all points onto C,N
		Real min_dist(0),max_dist(0); // that's for penetration depth
		BOOST_FOREACH(Vector3r& P,nodes)
		{
			Real point_plane_distance = N.Dot(P - C);
			min_dist = std::min(min_dist,point_plane_distance);
			max_dist = std::max(max_dist,point_plane_distance);
			P -= point_plane_distance * N;
		}
		Real penetrationDepth = max_dist - min_dist;
		scm->contactPoint = C;
		scm->normal = N;
		scm->penetrationDepth = penetrationDepth;

		int id1 = c->getId1();
		int id2 = c->getId2();

		if(m1->depth.find(id2) == m1->depth.end())
			m1->depth[id2] = penetrationDepth;
		if(m2->depth.find(id1) == m2->depth.end())
			m2->depth[id1] = penetrationDepth;

		Real d1 = m1->depth[id2];
		Real d2 = m2->depth[id1];
		if(d1 != d2)
			std::cerr << "bad initial penetration?\n";

		penetrationDepth -= d1;
		scm->penetrationDepth=penetrationDepth;

		scm->radius1=m1->radius;
		scm->radius2=m2->radius;

		return true;
	}
	return false;
*/
}


bool Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact::goReverse(	const shared_ptr<Shape>& cm1,
								const shared_ptr<Shape>& cm2,
								const Se3r& se31,
								const Se3r& se32,
								const Vector3r& shift2,
								const bool& force,
								const shared_ptr<Interaction>& c)
{
	std::cerr << "---- goReverse ---- " << __FILE__ << "\n";
	return go(cm1,cm2,se31,se32,c);
//	
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
}

YADE_PLUGIN((Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact));

YADE_REQUIRE_FEATURE(PHYSPAR);

