/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once

#include<iostream>
#include"Shape.hpp"
#include"Bound.hpp"
#include"State.hpp"
#include"Material.hpp"

#include"InteractionContainer.hpp"
#include"Interaction.hpp"

#include<yade/lib-base/Math.hpp>
#include<yade/lib-serialization/Serializable.hpp>
#include<yade/lib-multimethods/Indexable.hpp>

#define ISDYNAMIC_REDEFINED

class Scene;

class Body: public Serializable{
	public	:
		/*! \brief Numerical type for ::Body::id.
		 * \bug  Current code mixes signed and unsigned int, this might be a way to enforce consistence. */
		// typedef unsigned int id_t;
		//! symbolic constant for body that doesn't exist.
		static const body_id_t ID_NONE;
		//! get Body pointer given its id. 
		static const shared_ptr<Body>& byId(body_id_t _id,Scene* rb=NULL);
		static const shared_ptr<Body>& byId(body_id_t _id,shared_ptr<Scene> rb);

		
		//! Whether this Body is a Clump.
		//! @note The following is always true: \code (Body::isClump() XOR Body::isClumpMember() XOR Body::isStandalone()) \endcode
		bool isClump() const {return clumpId!=ID_NONE && id==clumpId;}
		//! Whether this Body is member of a Clump.
		bool isClumpMember() const {return clumpId!=ID_NONE && id!=clumpId;}
		//! Whether this body is standalone (neither Clump, nor member of a Clump)
		bool isStandalone() const {return clumpId==ID_NONE;}
		/*! Hook for clump to update position of members when user-forced reposition and redraw (through GUI) occurs.
		 * This is useful only in cases when engines that do that in every iteration are not active - i.e. when the simulation is paused.
		 * (otherwise, GLViewer would depend on Clump and therefore Clump would have to go to core...) */
		virtual void userForcedDisplacementRedrawHook(){return;}

		body_id_t getId() const {return id;};

		int getGroupMask() {return groupMask; };
		bool maskOk(int mask){return (mask==0 || (groupMask&mask));}

		// only BodyContainer can set the id of a body
		friend class BodyContainer;

		// Constructor/Destructor
		Body (body_id_t newId, int newGroup);

	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Body,Serializable,"A particle, basic element of simulation; interacts with other bodies.",
		((body_id_t,id,Body::ID_NONE,"[will be overridden]"))
		((int,groupMask,1,"Bitmask for determining interactions."))
		((bool,isDynamic,true,"Whether this body will be moved by forces."))

		((shared_ptr<Material>,material,,":yref:`Material` instance associated with this body."))
		((shared_ptr<State>,state,new State,"Physical :yref:`state<State>`."))
		((shared_ptr<Shape>,shape,,"Geometrical :yref:`Shape`."))
		((shared_ptr<Bound>,bound,,":yref:`Bound`, approximating volume for the purposes of collision detection."))

		((int,clumpId,Body::ID_NONE,"[will be overridden]")),
		/* ctor */,
		/* py */
		.def_readonly("id",&Body::id,"Unique id of this body") // overwrites automatic def_readwrite("id",...) earlier
		.def_readonly("clumpId",&Body::clumpId,"Id of clump this body makes part of; invalid number if not part of clump; see :yref:`Body::isStandalone`, :yref:`Body::isClump`, :yref:`Body::isClumpMember` properties. \n\n This property is not meant to be modified directly from Python, use :yref:`O.bodies.appendClumped<BodyContainer.appendClumped>` instead.")
		.def_readwrite("mat",&Body::material,"Shorthand for :yref:`Body::material`")
		.def_readwrite("dynamic",&Body::isDynamic,"Shorthand for :yref:`Body::isDynamic`")
		.def_readwrite("mask",&Body::groupMask,"Shorthand for :yref:`Body::groupMask`")
		.add_property("isStandalone",&Body::isStandalone,"True if this body is neither clump, nor clump member; false otherwise.")
		.add_property("isClumpMember",&Body::isClumpMember,"True if this body is clump member, false otherwise.")
		.add_property("isClump",&Body::isClump,"True if this body is clump itself, false otherwise.");
	);
};
REGISTER_SERIALIZABLE(Body);
