/*************************************************************************
*  Copyright (C) 2009 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"ElawSnowLayersDeformation.hpp"
#include<yade/pkg-dem/CohFrictMat.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/CohFrictPhys.hpp>
#include<yade/pkg-dem/SDECLinkPhysics.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg-snow/BssSnowGrain.hpp>
#include<yade/pkg-snow/BshSnowGrain.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

ElawSnowLayersDeformation::ElawSnowLayersDeformation() : InteractionSolver()
{
	sdecGroupMask=1;
	creep_viscosity = 1000.0;
}


void ElawSnowLayersDeformation::action()
{
	//return;
	shared_ptr<BodyContainer>& bodies = ncb->bodies;

//	Real dt = Omega::instance().getTimeStep();
	InteractionContainer::iterator ii    = ncb->interactions->begin();
	InteractionContainer::iterator iiEnd = ncb->interactions->end();
	for (  ; ii!=iiEnd ; ++ii )
	{
		if ((*ii)->isReal())
		{
			const shared_ptr<Interaction>& contact = *ii;
			int id1 = contact->getId1();
			int id2 = contact->getId2();



//			if(!(((id1 == 17 && id2 == 13) || (id1 == 13 && id2 == 17))))
//				continue;

			if ( !( (*bodies)[id1]->getGroupMask() & (*bodies)[id2]->getGroupMask() & sdecGroupMask)  )
				continue; // skip other groups,

			std::cerr << __FILE__ << " " << id1 << " " << id2 << "\n";

			CohFrictMat* de1 			= YADE_CAST<CohFrictMat*>((*bodies)[id1]->physicalParameters.get());
			CohFrictMat* de2 			= YADE_CAST<CohFrictMat*>((*bodies)[id2]->physicalParameters.get());
//			ScGeom* currentContactGeometry		= YADE_CAST<ScGeom*>(contact->interactionGeometry.get());
			CohFrictPhys* currentContactPhysics = YADE_CAST<CohFrictPhys*> (contact->interactionPhysics.get());

			BssSnowGrain* b1 = dynamic_cast<BssSnowGrain*>((*bodies)[id1]->shape.get());
			BssSnowGrain* b2 = dynamic_cast<BssSnowGrain*>((*bodies)[id2]->shape.get());
			
			BshSnowGrain* B1 = dynamic_cast<BshSnowGrain*>((*bodies)[id1]->geometricalModel.get());
			BshSnowGrain* B2 = dynamic_cast<BshSnowGrain*>((*bodies)[id2]->geometricalModel.get());

			Vector3r F = currentContactPhysics->shearForce + currentContactPhysics->normalForce;

			//FIXME: moment is still not used...
			Vector3r M = currentContactPhysics->moment_twist + currentContactPhysics->moment_bending;

			if(b1)
			{
				Vector3r c1 = b1->m_copy.c_axis;
				c1 /= c1.Length();
				//FIXME: make sure that F is always in correct direction
				F *= -1.0;
				Vector3r F1 = de1->se3.orientation.Conjugate() * F;
				F1 = F1 - c1*F1.Dot(c1);
				b1->m_copy.has_deformed();
				B1->has_deformed();

				std::vector<Real> layer_totals;
				Real all_layers_total(0);
				layer_totals.resize(b1->m_copy.slices.size() , 0);

				BOOST_FOREACH(const depth_one& t,b1->depths[id2])
				{
					int i  = t.i;
					//int j  = t.j;

					layer_totals[i] += 1.0;//+= t.current_depth - t.original_depth;
					++all_layers_total;
					//b1->m_copy.slices[i][j] += F1/(creep_viscosity);
					//B1->       slices[i][j] += F1/(creep_viscosity);
				}

				for(size_t i = 0 ; i < b1->m_copy.slices.size() ; ++i )
				{
					if(layer_totals[i] > 0)
					{
						for(size_t j = 0 ; j < b1->m_copy.slices[0].size() ; ++j )
						{
							b1->m_copy.slices[i][j] += (F1*layer_totals[i])/(creep_viscosity*all_layers_total);
							B1->       slices[i][j] += (F1*layer_totals[i])/(creep_viscosity*all_layers_total);
						}
					}
				}

			//	for(size_t i=0;i < b1->m_copy.slices.size();++i)
			//	{
			//		for(size_t j=0 ; j < b1->m_copy.slices[i].size() ; ++j)
			//		{
			//			//if(b1->m_copy.slices[i][j][2]>0)
			//			{
			//				b1->m_copy.slices[i][j] += F1/(creep_viscosity);
			//				B1->       slices[i][j] += F1/(creep_viscosity);
			//			}
			//		}
			//	}
			}

			if(b2)
			{
				Vector3r c2 = b2->m_copy.c_axis;
				c2 /= c2.Length();
				F *= -1.0;
				Vector3r F2 = de2->se3.orientation.Conjugate() * F;
				F2 = F2 - c2*F2.Dot(c2);
				b2->m_copy.has_deformed();
				B2->has_deformed();

				std::vector<Real> layer_totals;
				Real all_layers_total(0);
				layer_totals.resize(b2->m_copy.slices.size() , 0);

				BOOST_FOREACH(const depth_one& t,b2->depths[id1])
				{
					int i  = t.i;
					//int j  = t.j;

					layer_totals[i] += 1.0;//+= t.current_depth - t.original_depth; 
					// FIXME: divide force between all layers proportionally to sum of (t.current_depth - t.original_depth) of all nodes in each layer 
					++all_layers_total;
					//b2->m_copy.slices[i][j] += F2/(creep_viscosity);
					//B2->       slices[i][j] += F2/(creep_viscosity);
				}

				for(size_t i = 0 ; i < b2->m_copy.slices.size() ; ++i )
				{
					if(layer_totals[i] > 0)
					{
						for(size_t j = 0 ; j < b2->m_copy.slices[0].size() ; ++j )
						{
							b2->m_copy.slices[i][j] += (F2*layer_totals[i])/(creep_viscosity*all_layers_total);
							B2->       slices[i][j] += (F2*layer_totals[i])/(creep_viscosity*all_layers_total);
						}
					}
				}


			//	for(size_t i=0;i < b2->m_copy.slices.size();++i)
			//	{
			//		for(size_t j=0 ; j < b2->m_copy.slices[i].size() ; ++j)
			//		{
			//			//if(b2->m_copy.slices[i][j][2]>0)
			//			{
			//				b2->m_copy.slices[i][j] += F2/(creep_viscosity);
			//				B2->       slices[i][j] += F2/(creep_viscosity);
			//			}
			//		}
			//	}
			}
		}
	}
}

YADE_PLUGIN((ElawSnowLayersDeformation));

YADE_REQUIRE_FEATURE(PHYSPAR);

