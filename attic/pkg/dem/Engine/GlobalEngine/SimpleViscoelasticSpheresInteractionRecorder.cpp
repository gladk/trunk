/*************************************************************************
*  Copyright (C) 2006 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include "SimpleViscoelasticSpheresInteractionRecorder.hpp"
#include <yade/core/Omega.hpp>
#include <yade/core/Scene.hpp>
#include <boost/lexical_cast.hpp>

SimpleViscoelasticSpheresInteractionRecorder::SimpleViscoelasticSpheresInteractionRecorder (): DataRecorder(), interactionSphere(new ScGeom), viscoelasticInteraction(new ViscoelasticInteraction)
{
	outputBase = "contacts";
	interval=50;
}


SimpleViscoelasticSpheresInteractionRecorder::~SimpleViscoelasticSpheresInteractionRecorder ()
{

}


void SimpleViscoelasticSpheresInteractionRecorder::postProcessAttributes(bool deserializing)
{
	if(deserializing)
	{
        //
	}
}



bool SimpleViscoelasticSpheresInteractionRecorder::isActivated(Scene*)
{
   return ((Omega::instance().getCurrentIteration() % interval == 0));
}

void SimpleViscoelasticSpheresInteractionRecorder::action()
{
    shared_ptr<InteractionContainer>& interactions = ncb->interactions;
  
    InteractionContainer::iterator ii    = interactions->begin();
    InteractionContainer::iterator iiEnd = interactions->end();
    for(  ; ii!=iiEnd ; ++ii )
    {
	const shared_ptr<Interaction>& i = *ii;
	if ( !i->isReal() ) continue;
	if ( i->interactionGeometry->getClassIndex() != interactionSphere->getClassIndex() ) continue;
	if ( i->interactionPhysics->getClassIndex() != viscoelasticInteraction->getClassIndex() ) continue;

	const ScGeom* s = static_cast<ScGeom*>(i->interactionGeometry.get());
	const ViscoelasticInteraction* p = static_cast<ViscoelasticInteraction*>(i->interactionPhysics.get());

	std::string outputFile = outputBase+lexical_cast<string>(i->getId1())+"-"+lexical_cast<string>(i->getId2())+".dat";
	
	std::ofstream ofs;

	if (interaction_ids.insert(pair<int,int>(i->getId1(),i->getId2())).second)
	{
	    ofs.open(outputFile.c_str());
	    if (!ofs)
	    {
		std::cerr << "ERROR: Can't open file: " << outputFile << std::endl;
		return;
	    }
	    ofs << "#(t - time; d - penetration depth; F - normal force; S - shear force)\n"
		<< "#t(1)\td(2)\tFx(3)\tFy(4)\tFz(5)\tSx(6)\tSy(7)\tSz(8)"
		<< std::endl;
	}
	else
	{
	    ofs.open(outputFile.c_str(),std::ios::app);
	    if (!ofs)
	    {
		std::cerr << "ERROR: Can't open file: " << outputFile << std::endl;
		return;
	    }
	}

	ofs << lexical_cast<string>( Omega::instance().getCurrentIteration() ) << '\t'
	    << s->penetrationDepth << '\t'
	    << p->normalForce[0] << '\t'
	    << p->normalForce[1] << '\t'
	    << p->normalForce[2] << '\t'
	    << p->shearForce[0] << '\t'
	    << p->shearForce[1] << '\t'
	    << p->shearForce[2] << '\t'
	    << std::endl;
	ofs.close();
    }
}

YADE_PLUGIN((SimpleViscoelasticSpheresInteractionRecorder));

YADE_REQUIRE_FEATURE(PHYSPAR);

