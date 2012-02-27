/*************************************************************************
*  Copyright (C) 2006 by luc Scholtes                                    *
*  luc.scholtes@hmg.inpg.fr                                              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#include <yade/pkg-common/Recorder.hpp>

#include <string>
#include <fstream>

class TriaxialCompressionEngine;

class CapillaryStressRecorder : public Recorder
{
	private :
		shared_ptr<TriaxialCompressionEngine> triaxialCompressionEngine;
		
	public :
		virtual void action();
		
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(CapillaryStressRecorder,Recorder,"Records information from capillary meniscii on samples submitted to triaxial compressions. -> New formalism needs to be tested!!!",,initRun=true;);
	DECLARE_LOGGER;

};
REGISTER_SERIALIZABLE(CapillaryStressRecorder);
