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


class MovingSupport : public DataRecorder
{ // given a section plane it deletes all the beams that cross it, and remembers the nodes that they did connect.
  // Then it ensures that the nodes have constant distance in the 'direction' coordinate. But nodes can move in other directions
  // thus it creates a place that transfers only displacement in 'direction' and ignores movement in other directions.
  // It like a moving support.
	private :
		std::ofstream ofile; 
		std::vector<std::pair<std::list<std::pair<std::pair<unsigned int,unsigned int>,Real> >, int> >	subscribedBodies; // [region]<(list of node pairs,distance), direction>
		bool			first;
		std::list<unsigned int> futureDeletes;
		int			nodeGroupMask,beamGroupMask;

	public :
		std::list<std::pair<Vector3r,std::pair< int , Real > > > sections; // a list of section planes: <a midpoint, direction , half length >
		// int direction refers to: which coordinate to take from the Beam's.direction: x, y or z (0, 1 or 2)
		std::vector<Vector3r> sections_midpoints;
		std::vector<int> sections_directions;
		std::vector<Real> sections_halflength;

		MovingSupport ();

		virtual void action();
		virtual bool isActivated(Scene*);
	REGISTER_ATTRIBUTES(DataRecorder,/*(sections)*/(sections_midpoints)(sections_halflength)(sections_directions)/*(first)(subscribedBodies)*/);
	protected :
		virtual void preProcessAttributes(bool deserializing);
		virtual void postProcessAttributes(bool deserializing);
	REGISTER_CLASS_NAME(MovingSupport);
	REGISTER_BASE_CLASS_NAME(DataRecorder);
};

REGISTER_SERIALIZABLE(MovingSupport);


