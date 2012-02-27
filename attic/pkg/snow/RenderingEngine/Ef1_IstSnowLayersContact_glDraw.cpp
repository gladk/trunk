// © 2009 Janek Kozicki <cosurgi@mail.berlios.de>

#include"Ef1_IstSnowLayersContact_glDraw.hpp"
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-common/NormShearPhys.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>
#include<yade/pkg-snow/IstSnowLayersContact.hpp>
#include<yade/pkg-snow/BshSnowGrain.hpp>
#include<yade/pkg-snow/BssSnowGrain.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>

#include<yade/lib-opengl/OpenGLWrapper.hpp>
#include<yade/lib-opengl/GLUtils.hpp>

#ifndef MINIWM3
	#include <Wm3ApprPlaneFit3.h>
	#include <Wm3Plane3.h>
#endif
YADE_REQUIRE_FEATURE(geometricalmodel)
YADE_PLUGIN((Ef1_IstSnowLayersContact_glDraw));
/*
/// all points must share the same plane!
std::vector<Vector3r> find_boundary(std::vector<Vector3r> points)
{
	std::vector<Vector3r> res;
	if(points.size()>3)
	{
		Vector3r V(points[0].Cross(points[1]));
		V.Normalize();
		for(int I=0; I<points.size()-1 ; I++)
		for(int J=I+1; J<points.size() ; J++)
//		BOOST_FOREACH(Vector3r& a,points)
//		BOOST_FOREACH(Vector3r& b,points) // line a-b
		{
			Vector3r& a(points[I]);
			Vector3r& b(points[J]);
			if(a != b)
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
		// now sort it..

	}
	else
		res=points;
	return res;
}



Vector3r find_cross_point(BshSnowGrain* m,Vector3r in,Vector3r out, Quaternionr q2, Quaternionr q1, Vector3r pos1, Vector3r pos2,int depth = 4)
{
	Vector3r mid((in+out)*0.5);
	if(depth == 0)
		return mid;
	if(m->is_point_inside_polyhedron( q2.Conjugate()*(q1 * mid + pos1-pos2)))
	{
		return find_cross_point(m,mid,out,q2,q1,pos1,pos2,depth-1);
	}
	else
	{
		return find_cross_point(m,in, mid,q2,q1,pos1,pos2,depth-1);
	}
}
*/

void Ef1_IstSnowLayersContact_glDraw::go(
		const shared_ptr<InteractionGeometry>& ig,
		const shared_ptr<Interaction>& ip,
		const shared_ptr<Body>& b1,
		const shared_ptr<Body>& b2,
		bool wireFrame)
{
	if(!ip->isReal())
		return;

	IstSnowLayersContact* sc = static_cast<IstSnowLayersContact*>(ig.get());

	//const Se3r& se31=b1->physicalParameters->dispSe3,se32=b2->physicalParameters->dispSe3;
	//const Vector3r& pos1=se31.position,pos2=se32.position;

	if(wireFrame)
	{
		glPushMatrix();
			glTranslatev(sc->contactPoint);
			glBegin(GL_LINES);
				glColor3(0.0,1.0,1.0);
				glVertex3(0.0,0.0,0.0);
				glVertex3v(-1.0*(sc->normal*sc->radius1));
				glVertex3(0.0,0.0,0.0);
				glVertex3v( 1.0*(sc->normal*sc->radius2));
			glEnd();
		glPopMatrix();
	}
	else
	{
	}

	int id1 = ip->getId1();
	int id2 = ip->getId2();
	Real LEN(sc->radius1+sc->radius2);
	BssSnowGrain* m1 = dynamic_cast<BssSnowGrain*>(b1->shape.get());
	BssSnowGrain* m2 = dynamic_cast<BssSnowGrain*>(b2->shape.get());
		//std::map<int,std::list<boost::tuple<Real,int,int> > >	depths; // body_id , < depth, i, j (indexes of point on a slice) >
		//std::map<int,Real>	sphere_depth;
		//BshSnowGrain	m_copy;
	Vector3r    pos1(b1->physicalParameters->se3.position);
	Vector3r    pos2(b2->physicalParameters->se3.position);
	Quaternionr q1(b1->physicalParameters->se3.orientation);
	Quaternionr q2(b2->physicalParameters->se3.orientation);

//	Real dir_f=1.0; ////////////////// sometimes this must be -1.0. The normal should decide about it, becuase it always points from id1 to od2, but still there is a problem here and the lines are drawn in wrong direction
////	//Real r1=sc->radius1;
////	//Real r2=sc->radius2;
////	//if( (pos2-pos1).Dot(sc->normal) < 0 ) dir_f=-1.0;
////	//if(!m1 &&  m2 && (pos2-pos1).Dot(sc->normal) < 0 ) dir_f=-1.0;
////	//if( m1 && !m2 && (pos2-pos1).Dot(sc->normal) < 0 ) dir_f=-1.0;
////	if(m1 && m1->depths[id2].size() > 3)
////	{
////		int count=0;
////		std::vector<Vector3r> posX;posX.resize(m1->depths[id2].size());
////		BOOST_FOREACH(const depth_one& t,m1->depths[id2])
////		{
////			Real d = t.current_depth;
////			int i  = t.i;
////			int j  = t.j;
////			Vector3r pos(q1 * m1->m_copy.slices[i][j] + pos1);
////			posX[count++]=(pos - (dir_f) * d * sc->normal );
////			//if(count == 2)
////			//	break;
////		}
////		if(((posX[0]-posX[posX.size()-1]).Cross(posX[1]-posX[posX.size()-1])).Dot(- (sc->normal)) < 0.995)
////			dir_f*=-1.0;
////	}
////	if(m2 && m2->depths[id1].size() > 3)
////	{
////		int count=0;
////		std::vector<Vector3r> posX;posX.resize(m2->depths[id1].size());
////		BOOST_FOREACH(const depth_one& t,m2->depths[id1])
////		{
////			Real d = t.current_depth;
////			int i  = t.i;
////			int j  = t.j;
////			Vector3r pos(q2 * m2->m_copy.slices[i][j] + pos2);
////			posX[count++]=(pos - (dir_f) * d * sc->normal );
////			//if(count == 2)
////			//	break;
////		}
////		if(((posX[0]-posX[posX.size()-1]).Cross(posX[1]-posX[posX.size()-1])).Dot(- (sc->normal)) < 0.995)
////			dir_f*=-1.0;
////	}
	
	glColor3(1.0,0.0,0.0);
	if(m1)
	{
		//typedef boost::tuple<Real,int,int> tt;
		BOOST_FOREACH(const depth_one& t,m1->depths[id2])
		{
			Real d = t.current_depth;
			int i  = t.i;
			int j  = t.j;
			Vector3r pos(q1 * m1->m_copy.slices[i][j] + pos1);
				glTranslatev(pos);
				glutSolidCube(LEN*0.01);
				glTranslatev(-pos);
			glDisable(GL_LIGHTING);
			glLineWidth(3.0);
			glBegin(GL_LINES);
				glVertex3v(pos);
				glVertex3v(pos - std::abs(d) * sc->normal );
			glEnd();
			glLineWidth(1.0);
			glEnable(GL_LIGHTING);
		}
	}

	glColor3(0.0,1.0,0.0);
	if(m2)
	{
		//typedef boost::tuple<Real,int,int> tt;
		BOOST_FOREACH(const depth_one& t,m2->depths[id1])
		{
			Real d = t.current_depth;
			int i  = t.i;
			int j  = t.j;
			Vector3r pos(q2 * m2->m_copy.slices[i][j] + pos2);
				glTranslatev(pos);
				glutSolidCube(LEN*0.01);
				glTranslatev(-pos);
			glDisable(GL_LIGHTING);
			glLineWidth(3.0);
			glBegin(GL_LINES);
				glVertex3v(pos);
				glVertex3v(pos + std::abs(d) * sc->normal );
			glEnd();
			glLineWidth(1.0);
			glEnable(GL_LIGHTING);
		}
	}



/*
	if(!(id1 == (int)(Omega::instance().getScene()->selectedBody) || id2 == (int)(Omega::instance().getScene()->selectedBody)))
		return;

	assert(Omega::instance().getScene()->bodies->exists(id1));
	assert(Omega::instance().getScene()->bodies->exists(id2));
	BshSnowGrain* m1 = dynamic_cast<BshSnowGrain*>((*(Omega::instance().getScene()->bodies))[id1]->geometricalModel.get());
	BshSnowGrain* m2 = dynamic_cast<BshSnowGrain*>((*(Omega::instance().getScene()->bodies))[id2]->geometricalModel.get());
	if(!m1 || !m2)
	{
		std::cerr << "not BshSnowGrain\n";
		return;
	}
	std::vector<Vector3r> inside;
	std::vector<Vector3r> cross_section;
	Vector3r    pos1((*(Omega::instance().getScene()->bodies))[id1]->physicalParameters->se3.position);
	Vector3r    pos2((*(Omega::instance().getScene()->bodies))[id2]->physicalParameters->se3.position);
	Quaternionr q1((*(Omega::instance().getScene()->bodies))[id1]->physicalParameters->se3.orientation);
	Quaternionr q2((*(Omega::instance().getScene()->bodies))[id2]->physicalParameters->se3.orientation);

	glColor3(1.0,0.0,0.0);
	for(size_t i=0;i < m1->slices.size();++i)
	{
		for(size_t j=0 ; j < m1->slices[i].size() ; ++j)
		{
			Vector3r v(m1->slices[i][j]);
			if(m2->is_point_inside_polyhedron( q2.Conjugate()*(q1 * v + pos1-pos2)))
			{
				Vector3r pos(q1 * v + pos1);
				inside.push_back(pos);
				glTranslatev(pos);
				glutSolidCube(LEN*0.01);
				glTranslatev(-pos);
			}
		}
	}
	glColor3(0.0,1.0,0.0);
	for(size_t i=0;i < m2->slices.size();++i)
	{
		for(size_t j=0 ; j < m2->slices[i].size() ; ++j)
		{
			Vector3r v(m2->slices[i][j]);
			if(m1->is_point_inside_polyhedron( q1.Conjugate()*(q2 * v + pos2-pos1)))
			{
				Vector3r pos(q2 * v + pos2);
				inside.push_back(pos);
				glTranslatev(pos);
				glutSolidCube(LEN*0.01);
				glTranslatev(-pos);
			}
		}
	}
*/
////////////////////////////////////////
/*
	typedef std::pair<Vector3r, std::set<Vector3r> > t_edge;
	const std::map<Vector3r, std::set<Vector3r> >& edges1 = m1->get_edges_const_ref();
	const std::map<Vector3r, std::set<Vector3r> >& edges2 = m2->get_edges_const_ref();
	glColor3(1.0,0.5,0.0);
	BOOST_FOREACH(const t_edge& e,edges1)
	{
			Vector3r v(e.first);
			if(m2->is_point_inside_polyhedron( q2.Conjugate()*(q1 * v + pos1-pos2)))
			{
				BOOST_FOREACH(Vector3r v2,e.second)
				{
					if( ! (m2->is_point_inside_polyhedron( q2.Conjugate()*(q1 * v2 + pos1-pos2))))
					{
						Vector3r v_c(q1 * find_cross_point(m2,v,v2, q2,q1,pos1,pos2) + pos1);
						cross_section.push_back(v_c);
						glColor3(0.3,0.3,1.0);
						glTranslatev(v_c);
						glutSolidCube(LEN*0.008);
						glTranslatev(-v_c);
						glColor3(1.0,0.5,0.0);
					}
				}
				Vector3r pos(q1 * v + pos1);
			//	inside.push_back(pos);
				glTranslatev(pos);
				glutSolidSphere(LEN*0.008,6,6);
				glTranslatev(-pos);
			}
	}
	glColor3(0.5,1.0,0.0);
	BOOST_FOREACH(const t_edge& e,edges2)
	{
			Vector3r v(e.first);
			if(m1->is_point_inside_polyhedron( q1.Conjugate()*(q2 * v + pos2-pos1)))
			{
				BOOST_FOREACH(Vector3r v2,e.second)
				{
					if( ! (m1->is_point_inside_polyhedron( q1.Conjugate()*(q2 * v2 + pos2-pos1))))
					{
						Vector3r v_c(q2 * find_cross_point(m1,v,v2, q1,q2,pos2,pos1) + pos2);
						cross_section.push_back(v_c);
						glColor3(0.3,0.3,1.0);
						glTranslatev(v_c);
						glutSolidCube(LEN*0.008);
						glTranslatev(-v_c);
						glColor3(0.5,1.0,0.0);
					}
				}
				Vector3r pos(q2 * v + pos2);
			//	inside.push_back(pos);
				glTranslatev(pos);
				glutSolidSphere(LEN*0.008,6,6);
				glTranslatev(-pos);
			}
	}
*/


/*
	if(!wireFrame)
	{
		// got a list of points.
		Vector3r C(0,0,0);
		Real count(0);
		BOOST_FOREACH(Vector3r& v,inside) C+=v,++count;
		C /= count;
		Vector3r N(sc->normal);
		N.Normalize();
		std::vector<Vector3r> ins_copy(inside);
		// plane of contact is C,N
		// project all points onto C,N
		BOOST_FOREACH(Vector3r& P,inside) P -= N.Dot(P - C) * N;
		
	//	glColor3(1.0,1.0,1.0);
	//	glBegin(GL_LINES);
	//	BOOST_FOREACH(Vector3r& a,inside)
	//	BOOST_FOREACH(Vector3r& b,inside)
	//	{
	//		glVertex3v(a);
	//		glVertex3v(b);
	//	}
	//	glEnd();






		Plane3<Real> plane(OrthogonalPlaneFit3 (ins_copy.size(), &ins_copy[0]));
		N = plane.Normal;
		C = Vector3r(0,0,0);
		int i=0;
		if( std::abs(N[1]) > std::abs(N[0]) && std::abs(N[1]) > std::abs(N[2]) ) i=1;
		if( std::abs(N[2]) > std::abs(N[1]) && std::abs(N[2]) > std::abs(N[0]) ) i=2;
		C[i] = plane.Constant/N[i];

		// project all points onto C,N
		BOOST_FOREACH(Vector3r& P,ins_copy) P -= N.Dot(P - C) * N;

	//	glColor3(1.0,1.0,0.0);
	//	glBegin(GL_LINES);
	//	BOOST_FOREACH(Vector3r& a,ins_copy)
	//	BOOST_FOREACH(Vector3r& b,ins_copy)
	//	{
	//		glVertex3v(a);
	//		glVertex3v(b);
	//	}
	//	glEnd();






		Plane3<Real> plane2(OrthogonalPlaneFit3 (cross_section.size(), &cross_section[0]));
		N = plane2.Normal;
		C = Vector3r(0,0,0);
		i=0;
		if( std::abs(N[1]) > std::abs(N[0]) && std::abs(N[1]) > std::abs(N[2]) ) i=1;
		if( std::abs(N[2]) > std::abs(N[1]) && std::abs(N[2]) > std::abs(N[0]) ) i=2;
		C[i] = plane2.Constant/N[i];

		// project all points onto C,N
		BOOST_FOREACH(Vector3r& P,cross_section) P -= N.Dot(P - C) * N;

	//	glColor3(1.0,1.0,0.0);
	//	glBegin(GL_LINES);
	//	BOOST_FOREACH(Vector3r& a,cross_section)
	//	BOOST_FOREACH(Vector3r& b,cross_section)
	//	{
	//		glVertex3v(a);
	//		glVertex3v(b);
	//	}
	//	glEnd();






	//	std::vector<Vector3r> poly1(find_boundary(inside));
	//	std::vector<Vector3r> poly2(find_boundary(ins_copy));
		std::vector<Vector3r> poly3(find_boundary(cross_section));
//
//		glColor3(1.0,1.0,1.0);
//		glBegin(GL_POLYGON);
//		BOOST_FOREACH(Vector3r& a,poly1)
//			glVertex3v(a);
//		glEnd();
//
//		glColor3(1.0,1.0,0.0);
//		glBegin(GL_POLYGON);
//		BOOST_FOREACH(Vector3r& a,poly2)
//			glVertex3v(a);
//		glEnd();
//
		glColor3(0.3,0.3,1.0);
		glBegin(GL_POLYGON);
		BOOST_FOREACH(Vector3r& a,poly3)
			glVertex3v(a);
		glEnd();
		glColor3(0.8,0.8,1.0);
		BOOST_FOREACH(Vector3r& a,poly3)
		{
			glTranslatev(a);
			glutSolidCube(LEN*0.014);
			glTranslatev(-a);
		}
	}
*/
}


YADE_REQUIRE_FEATURE(PHYSPAR);

