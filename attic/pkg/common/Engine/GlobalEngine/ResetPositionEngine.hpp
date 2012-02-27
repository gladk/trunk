/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/PeriodicEngines.hpp>
#include <Wm3Vector3.h>
#include<yade/lib-base/yadeWm3.hpp>

/*! \brief Reset the spatial position of all subscribed bodies to the desired value. 
 *
 *
 * */
class ResetPositionEngine : public PeriodicEngine {
	public:
		void action();
		bool isActivated(Scene* rootBody) {if (first) return true; else return PeriodicEngine::isActivated(rootBody);}
		ResetPositionEngine();

		Real Y_min;
		/// \brief import/export desired positions from file.
		/// If fileName is set and initial_positions are not set, then positions are read from file.
		/// If fileName is not set and initial_positions are set, then positions are read bodies se3
		std::string fileName; 
		std::vector< int >   subscribedBodies;
		std::vector<Vector3r> initial_positions; // for serialization
		
		DECLARE_LOGGER;
	protected:
		virtual void postProcessAttributes(bool);
	REGISTER_ATTRIBUTES(PeriodicEngine,(Y_min)(subscribedBodies)(initial_positions)(fileName));
	REGISTER_CLASS_NAME(ResetPositionEngine);
	REGISTER_BASE_CLASS_NAME(PeriodicEngine);
	private:
		std::vector<Vector3r> ini_pos;
		bool first;
		void initialize();
};

REGISTER_SERIALIZABLE(ResetPositionEngine);


