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

class MeasurePoisson : public DataRecorder
{
	private :
		std::ofstream ofile; 
		Real	last20[20];
		int	idx;
	

	public :
		std::string	 outputFile;
		unsigned int	 interval;
		Real 		 vertical,horizontal;
		int		 bottom,upper,left,right;
		
		MeasurePoisson ();

		virtual void action();
		virtual bool isActivated(Scene*);

	REGISTER_ATTRIBUTES(DataRecorder,(outputFile)(interval)(horizontal)(vertical)(bottom)(upper)(left)(right));
	protected :
		virtual void postProcessAttributes(bool deserializing);
	REGISTER_CLASS_NAME(MeasurePoisson);
	REGISTER_BASE_CLASS_NAME(DataRecorder);
};

REGISTER_SERIALIZABLE(MeasurePoisson);


