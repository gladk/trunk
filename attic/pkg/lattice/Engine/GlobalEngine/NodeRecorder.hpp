/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/DataRecorder.hpp>

#include <string>
#include <fstream>
#include <list>
#include <vector>


class NodeRecorder : public DataRecorder
{ // records average position of nodes in a given region
	private :
		std::ofstream ofile; 
		
		std::vector<std::list<unsigned int> >	subscribedBodies; // [region]<list of bodies' ids>
		bool			first;
		int			nodeGroupMask,beamGroupMask;

	public :
		std::string	 outputFile;
		unsigned int	 interval;
		std::list<std::pair<Vector3r,Vector3r> > regions; // a list of min/max pairs describing each region.
		std::vector<Vector3r> regions_min;
		std::vector<Vector3r> regions_max;
		Real	only_this_stiffness;

		NodeRecorder ();

		virtual void action();
		virtual bool isActivated(Scene*);
	REGISTER_ATTRIBUTES(DataRecorder,
		(outputFile)
		(interval)
		//(regions)
		(regions_min)
		(regions_max)
		(only_this_stiffness)
		//	(first)
		//	(subscribedBodies)
	);
	protected :
		virtual void preProcessAttributes(bool deserializing);
		virtual void postProcessAttributes(bool deserializing);
	REGISTER_CLASS_NAME(NodeRecorder);
	REGISTER_BASE_CLASS_NAME(DataRecorder);
};

REGISTER_SERIALIZABLE(NodeRecorder);


