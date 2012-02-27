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


class StrainRecorder : public DataRecorder
{
	private :
		std::ofstream ofile; 

	public :
		std::string	 outputFile;
		unsigned int	 interval;
		std::vector<unsigned int> subscribedBodies;
		Real 		 initialLength;
		bool		 only_stress;

		StrainRecorder ();

		virtual void action();
		virtual bool isActivated(Scene*);

	protected :
		virtual void postProcessAttributes(bool deserializing);
	REGISTER_ATTRIBUTES(DataRecorder,(outputFile)(interval)(subscribedBodies)(initialLength)(only_stress));
	REGISTER_CLASS_NAME(StrainRecorder);
	REGISTER_BASE_CLASS_NAME(DataRecorder);
};
REGISTER_SERIALIZABLE(StrainRecorder);


