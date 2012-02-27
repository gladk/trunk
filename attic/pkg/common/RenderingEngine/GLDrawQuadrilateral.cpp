/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "GLDrawQuadrilateral.hpp"
#include<yade/pkg-common/Quadrilateral.hpp>
#include<yade/core/Scene.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>
YADE_REQUIRE_FEATURE(geometricalmodel);

Vector3r makeColor(double val,double min,double max)
{
	Real sc01 = ((val-min)/(max-min))*2.0-1.0;
	sc01*=-1.0;
	if(sc01 < 0) return Vector3r(0.9,0.9,1.0) + sc01 * Vector3r(0.9,0.9,0.0);
	if(sc01 > 0) return Vector3r(1.0,0.9,0.9) - sc01 * Vector3r(0.0,0.9,0.9);
	return Vector3r(0.9,0.9,0.9);
}

void GLDrawQuadrilateral::go(const shared_ptr<GeometricalModel>& gm, const shared_ptr<PhysicalParameters>&,bool wire)
{	

// FIXME - that's a hack also. We must give more power to the GUI
//
	static bool initialized=false;
	if(!initialized)
	{
		Omega::instance().isoSec=0;
		initialized=true;
	}
	int WTT = Omega::instance().isoSec%3; 

	Quadrilateral* q = static_cast<Quadrilateral*>(gm.get());
	Scene* mb = Omega::instance().getScene().get();
	
	if(   !(
		   mb->bodies->exists(q->n1)
		&& mb->bodies->exists(q->n2)
		&& mb->bodies->exists(q->n3)
		&& mb->bodies->exists(q->n4)))
		return;

	double ox1,oy1,ox2,oy2,ox3,oy3,ox4,oy4;
	double nx1,ny1,nx2,ny2,nx3,ny3,nx4,ny4;
	double e11,e22,e12,e21,x,y;
	
	#define Q(x) ( (*(mb->bodies))[q->x]->physicalParameters->se3.position )
	ox1 = q->i1[0];
	ox2 = q->i2[0];
	ox3 = q->i3[0];
	ox4 = q->i4[0];

	oy1 = q->i1[1];
	oy2 = q->i2[1];
	oy3 = q->i3[1];
	oy4 = q->i4[1];

	nx1 = Q(n1)[0];
	nx2 = Q(n2)[0];
	nx3 = Q(n3)[0];
	nx4 = Q(n4)[0];

	ny1 = Q(n1)[1];
	ny2 = Q(n2)[1];
	ny3 = Q(n3)[1];
	ny4 = Q(n4)[1];

	calculateStrainQuadrilateral(ox1,oy1,ox2,oy2,ox3,oy3,ox4,oy4,nx1,ny1,nx2,ny2,nx3,ny3,nx4,ny4,0,0,e11,e22,e12,e21,x,y);

//	static Real lastx = 0;
//	//std::cout << " x: " << x << " y: " << y << "e11: " << e11 << " e22: " << e22 << " e12: " << e12 << " e21: " << e21  << "\n";
//	if( x>=0.20 && x <=0.40 && y>0.037 && y <0.10 )
//	//if( x>=0.01 && x <=0.31 && y>0.023 && y <0.097 )
//		std::cout << x << " " << y << " " << e11 << " " << e22 << " " << e12 << " " << e21 << "\n";
//	if(x<lastx)
//		std::cout << "\n";
//	lastx=x;

	//double m=-0.00005,M = 0.0001;
	double m=-0.005,M = 0.005; // dodatni to rozciaganie, niebieski
	//double m=-1,M = 1;//GLDrawLatticeSetGeometry
	//M = std::max(e22,M);
	//m = std::min(e22,m);
	//static int zzz=0;
	//if((++zzz%1000)==0)
	//	std::cerr << m << " " << M << "\n";
	static int ZZZ=0;
	switch(WTT)
	{
		case 0: glColor3v(makeColor(e11,m,M));if(ZZZ++%10000==0)std::cerr << "e11 ";break;
		case 1: glColor3v(makeColor(e22,m,M));if(ZZZ++%10000==0)std::cerr << "e22 ";break;
		case 2: glColor3v(makeColor(e12,m,M));if(ZZZ++%10000==0)std::cerr << "e12 ";break;
		default: std::cerr << "WTF?\n";
	}
	//if(std::abs(e22) > 0.0003 /*|| std::abs(e11) > r*/) return;

	if (gm->wire || wire)
	{
		glBegin(GL_LINE_LOOP);
		glDisable(GL_LIGHTING);

		glVertex3v(Q(n1));
		glVertex3v(Q(n2));
		glVertex3v(Q(n3));
		glVertex3v(Q(n4));

		glEnd();
	}
	else
	{
		glDisable(GL_CULL_FACE);
		glDisable(GL_LIGHTING);
		glBegin(GL_QUADS);

		glVertex3v(Q(n1));
		glVertex3v(Q(n2));
		glVertex3v(Q(n3));
		glVertex3v(Q(n4));

		glEnd();
	}
	#undef Q
}

double H1(double r,double s)		{	return 0.25*(1-r)*(1-s);	}
double H2(double r,double s)		{	return 0.25*(1+r)*(1-s);	}
double H3(double r,double s)		{	return 0.25*(1+r)*(1+s);	}
double H4(double r,double s)		{	return 0.25*(1-r)*(1+s);	}

double dH1_dr(double r,double s)	{	return -0.25*(1-s);		}
double dH1_ds(double r,double s)	{	return -0.25*(1-r);		}
double dH2_dr(double r,double s)	{	return  0.25*(1-s);		}
double dH2_ds(double r,double s)	{	return -0.25*(1+r);		}
double dH3_dr(double r,double s)	{	return  0.25*(1+s);		}
double dH3_ds(double r,double s)	{	return  0.25*(1+r);		}
double dH4_dr(double r,double s)	{	return -0.25*(1+s);		}
double dH4_ds(double r,double s)	{	return  0.25*(1-r);		}


void GLDrawQuadrilateral::calculateStrainQuadrilateral(
	 double ox1 ,double oy1 ,double ox2 ,double oy2 ,double ox3 ,double oy3 ,double ox4 ,double oy4 
	,double nx1 ,double ny1 ,double nx2 ,double ny2 ,double nx3 ,double ny3 ,double nx4 ,double ny4 
	,double r ,double s
	,double& e11 ,double& e22 ,double& e12 ,double& e21,double& x,double& y)
{
	// coordinates of integration point
	Vector3d X =	  H1(r,s) * Vector3d(ox1,oy1,0)
			+ H2(r,s) * Vector3d(ox2,oy2,0)
			+ H3(r,s) * Vector3d(ox3,oy3,0)
			+ H4(r,s) * Vector3d(ox4,oy4,0);
	
	// displacements, delta u, delta v
	double Du1   =    nx1 - ox1; //C(Du1);
	double Du2   =    nx2 - ox2; //C(Du2);
	double Du3   =    nx3 - ox3; //C(Du3);
	double Du4   =    nx4 - ox4; //C(Du4);
	
	double Dv1   =    ny1 - oy1; //C(Dv1);
	double Dv2   =    ny2 - oy2; //C(Dv2);
	double Dv3   =    ny3 - oy3; //C(Dv3);
	double Dv4   =    ny4 - oy4; //C(Dv4);

	// Jacobian matrix components
	double dX_dr =    dH1_dr(r,s) * nx1
			+ dH2_dr(r,s) * nx2
			+ dH3_dr(r,s) * nx3
			+ dH4_dr(r,s) * nx4;

	double dX_ds =    dH1_ds(r,s) * nx1
			+ dH2_ds(r,s) * nx2
			+ dH3_ds(r,s) * nx3
			+ dH4_ds(r,s) * nx4;

	double dY_dr =    dH1_dr(r,s) * ny1
			+ dH2_dr(r,s) * ny2
			+ dH3_dr(r,s) * ny3
			+ dH4_dr(r,s) * ny4;

	double dY_ds =    dH1_ds(r,s) * ny1
			+ dH2_ds(r,s) * ny2
			+ dH3_ds(r,s) * ny3
			+ dH4_ds(r,s) * ny4;

	// its determinant
	double detJ  = dX_dr*dY_ds-dY_dr*dX_ds;

	// componentes of inversed matrix
	double dr_dX =  dY_ds/detJ;
	double ds_dX = -dY_dr/detJ;
	double dr_dY = -dX_ds/detJ;
	double ds_dY =  dX_dr/detJ;

	// derivatives of shape function
	double dH1_dX = dH1_dr(r,s)*dr_dX + dH1_ds(r,s)*ds_dX;
	double dH2_dX = dH2_dr(r,s)*dr_dX + dH2_ds(r,s)*ds_dX;
	double dH3_dX = dH3_dr(r,s)*dr_dX + dH3_ds(r,s)*ds_dX;
	double dH4_dX = dH4_dr(r,s)*dr_dX + dH4_ds(r,s)*ds_dX;

	double dH1_dY = dH1_dr(r,s)*dr_dY + dH1_ds(r,s)*ds_dY;
	double dH2_dY = dH2_dr(r,s)*dr_dY + dH2_ds(r,s)*ds_dY;
	double dH3_dY = dH3_dr(r,s)*dr_dY + dH3_ds(r,s)*ds_dY;
	double dH4_dY = dH4_dr(r,s)*dr_dY + dH4_ds(r,s)*ds_dY;

	// calculated strian
	e11 = dH1_dX*Du1 + dH2_dX*Du2 + dH3_dX*Du3 + dH4_dX*Du4;
	e22 = dH1_dY*Dv1 + dH2_dY*Dv2 + dH3_dY*Dv3 + dH4_dY*Dv4;
	e12 = dH1_dY*Du1 + dH2_dY*Du2 + dH3_dY*Du3 + dH4_dY*Du4;
	e21 = dH1_dX*Dv1 + dH2_dX*Dv2 + dH3_dX*Dv3 + dH4_dX*Dv4;
	
	x = X[0];
	y = X[1];
};

YADE_PLUGIN((GLDrawQuadrilateral));
YADE_REQUIRE_FEATURE(OPENGL)

YADE_REQUIRE_FEATURE(PHYSPAR);

