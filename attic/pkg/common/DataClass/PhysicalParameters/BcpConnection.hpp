/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#pragma once

#include <yade/pkg-dem/SimpleViscoelasticBodyParameters.hpp>


class BcpConnection : public SimpleViscoelasticBodyParameters
{
	public :
		// parameters (no parameters)

		// state
		unsigned int id1, id2;

		BcpConnection();
		virtual ~BcpConnection();

/// Serializable
	REGISTER_ATTRIBUTES(SimpleViscoelasticBodyParameters,(id1)(id2));
	REGISTER_CLASS_NAME(BcpConnection);
	REGISTER_BASE_CLASS_NAME(SimpleViscoelasticBodyParameters);

/// Indexable
	REGISTER_CLASS_INDEX(BcpConnection,SimpleViscoelasticBodyParameters);

};

REGISTER_SERIALIZABLE(BcpConnection);


