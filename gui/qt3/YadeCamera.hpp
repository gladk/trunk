/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#ifndef YADE_LOCAL_QGLVIEWER
	#include<QGLViewer/qglviewer.h>
#else
	#include<yade/lib-QGLViewer/qglviewer.h>
#endif
#include<iostream>


class YadeCamera : public qglviewer::Camera
{	
	Q_OBJECT 
	private:
		float cuttingDistance;
        public :
		YadeCamera():cuttingDistance(0){};
		virtual float zNear() const;
		virtual void setCuttingDistance(float s){cuttingDistance=s;};
};



