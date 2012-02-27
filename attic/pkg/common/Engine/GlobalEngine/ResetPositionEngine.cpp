/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"ResetPositionEngine.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/ParticleParameters.hpp>
#include<boost/foreach.hpp>

CREATE_LOGGER(ResetPositionEngine);

ResetPositionEngine::ResetPositionEngine(){
		Y_min=-1;
		ini_pos.clear();
		initial_positions.clear();
		subscribedBodies.clear();
		fileName="";
		first=true;
}

void ResetPositionEngine::postProcessAttributes(bool deserializing){}


void ResetPositionEngine::action()
{
	if (first) { initialize(ncb); return; }

	shared_ptr<BodyContainer>& bodies = ncb->bodies;
	for(int i=0,e=subscribedBodies.size(); i<e; ++i)
	{
			ParticleParameters* pp = YADE_CAST<ParticleParameters*>((*bodies)[subscribedBodies[i]]->physicalParameters.get());
			if (pp->se3.position[1]<Y_min)
			{
					pp->se3.position = ini_pos[i];
					pp->velocity = Vector3r(0,0,0);
			}
	}
}

void ResetPositionEngine::initialize()
{
	first=false;
	if (fileName=="")
	{ 
		if (initial_positions.size()==0) 	
		{ // initialize positions from bodies se3
			LOG_INFO("Initialize positions from bodies se3");
			initial_positions.resize(subscribedBodies.size());
			shared_ptr<BodyContainer>& bodies = ncb->bodies;
			for(int i=0,e=subscribedBodies.size(); i<e; ++i)
				initial_positions[i]=(*bodies)[subscribedBodies[i]]->physicalParameters->se3.position;
		}
		ini_pos.assign(initial_positions.begin(),initial_positions.end());
		return;
	}
	
	std::ifstream is(fileName.c_str());
	if (is) 
	{// reading positions from file 
		LOG_INFO("Reading positions from file: " << fileName);
		ini_pos.resize(subscribedBodies.size());
		for(int i=0,e=subscribedBodies.size(); i<e && !is.eof(); ++i)
			is >> ini_pos[i][0] >> ini_pos[i][1] >> ini_pos[i][2];
		return;
	}
	
	// initialize positions form bodies se3 if need and export to file
	if (initial_positions.size()==0)
	{
		LOG_INFO("Initialize positions from bodies se3");
		ini_pos.resize(subscribedBodies.size());
		shared_ptr<BodyContainer>& bodies = ncb->bodies;
		for(int i=0,e=subscribedBodies.size(); i<e; ++i)
			ini_pos[i]=(*bodies)[subscribedBodies[i]]->physicalParameters->se3.position;
	}
	else ini_pos.swap(initial_positions);
	std::ofstream os(fileName.c_str());
	if (!os) {
		LOG_ERROR("Can't open file to export positions: "<<fileName<<"!");
		return;
	}
	LOG_INFO("Export positions to file: "<<fileName);
	for(int i=0,e=subscribedBodies.size(); i<e; ++i)
		os << ini_pos[i][0]<< '\t' << ini_pos[i][1]<< '\t' << ini_pos[i][2]<< std::endl;
}
YADE_PLUGIN((ResetPositionEngine));

YADE_REQUIRE_FEATURE(PHYSPAR);

