/*************************************************************************
*  Copyright (C) 2009 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "MGPRecorder.hpp"

#include <yade/core/Omega.hpp>
#include <yade/core/Scene.hpp>
#include <boost/lexical_cast.hpp>

#include <yade/pkg-common/Sphere.hpp>

MGPRecorder::MGPRecorder(): DataRecorder()
{
  outputBase = "mgp.out.";
  interval   = 50;
  stateId    = 0;
}


MGPRecorder::~MGPRecorder ()
{

}


void MGPRecorder::postProcessAttributes(bool deserializing)
{
  if(deserializing)
  {
  //
  }
}


bool MGPRecorder::isActivated(Scene*)
{
   return ((Omega::instance().getCurrentIteration() % interval == 0));
}

void MGPRecorder::action()
{
  shared_ptr<BodyContainer>& bodies = ncb->bodies;

  std::string outputFile = outputBase+lexical_cast<string>(stateId);
  std::ofstream ofs(outputFile.c_str());

  ofs << "<?xml version=\"1.0\"?>" << endl
      << " <mgpost mode=\"3D\">" << endl
      << "  <state id=\"" << stateId++
      << "\" time=\"" << lexical_cast<string>(Omega::instance().getSimulationTime()) << "\">" << endl;

   BodyContainer::iterator bi    = bodies->begin();
   BodyContainer::iterator biEnd = bodies->end();
   for(; bi != biEnd ; ++bi)
   {
     const shared_ptr<Body>& b = *bi;
	 if (b->isClump()) continue;
	 
     const RigidBodyParameters* p = YADE_CAST<RigidBodyParameters*>(b->physicalParameters.get());
     const Shape* gm   = YADE_CAST<Shape*>(b->shape.get());

     if ( typeid(*gm) == typeid(Sphere) )
     {
      ofs << "   <body>" << endl;
      ofs << "    <SPHER id=\"" << b->getId() << "\" r=\"" << YADE_CAST<Sphere*>(b->shape.get())->radius << "\">" << endl;
      ofs << "     <position x=\"" << p->se3.position[0] << "\" y=\"" 
          << p->se3.position[1] << "\" z=\"" << p->se3.position[2] << "\"/>" << endl;   
      ofs << "     <velocity x=\"" << p->velocity[0] << "\" y=\"" 
          << p->velocity[1] << "\" z=\"" << p->velocity[2] << "\"/>" << endl; 
      ofs << "    </SPHER>" << endl;
      ofs << "   </body>" << endl << flush;
     }



   }


  ofs << "  </state>" << endl
      << " </mgpost>" << endl << flush;


}

YADE_PLUGIN((MGPRecorder));

YADE_REQUIRE_FEATURE(PHYSPAR);

