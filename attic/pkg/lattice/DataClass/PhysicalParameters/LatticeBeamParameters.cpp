/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LatticeBeamParameters.hpp"
#include<yade/core/Omega.hpp>

LatticeBeamParameters::LatticeBeamParameters() : PhysicalParameters()
{
	createIndex();
	count = 0.0;
	bendingRotation = Quaternionr(1.0,0.0,0.0,0.0);
	torsionalRotation = Quaternionr(1.0,0.0,0.0,0.0);
	torsionAngle = 0;
	torsionalStiffness = 0.6;
	lastIter_ = -1;
	shearing_strain = Vector3r(0,0,0);
}


LatticeBeamParameters::~LatticeBeamParameters()
{

}

void LatticeBeamParameters::postProcessAttributes(bool deserializing)
{
	if(deserializing)
	{

	//	torsionalStiffness /= 10.0;

//		static bool first=true;
//		if(first)
//		{
//			std::cerr << "using 0.2 for k.b. So k.l=" << longitudalStiffness 
////				<< " dt="<< Omega::instance().getTimeStep() 
////				<< " k.b=" << longitudalStiffness * Omega::instance().getTimeStep() 
//				<< " k.b=" << longitudalStiffness * 0.2 
//				<<"\n";
////				//<<"\nAnd criticalTensileStrain*10 and bendingStiffness/10 and longitudalStiffness/10\n\n", first=false;
//			 first=false;
//		}
//
//		bendingStiffness    = longitudalStiffness * 0.2;
//
//		bendingStiffness		/= 10.0;
//		longitudalStiffness		/= 10.0;
//		criticalTensileStrain		*= 10.0;


/*		if(criticalTensileStrain > 0.00015) // E.l
		{ // CEMENT MATRIX
			longitudalStiffness       = 1.0;	// k.l
			bendingStiffness          = 0.7;	// k.b
			criticalCompressiveStrain = 1.0;	// E.c
			criticalTensileStrain     = 0.0002;	// E.l
		}
		else if(criticalTensileStrain > 0.00006) // E.l
		{ // AGGREGATE
			longitudalStiffness       = 3.0;	// k.l
			bendingStiffness          = 2.1;	// k.b
			criticalCompressiveStrain = 1.0;	// E.c
			criticalTensileStrain     = 0.0001333;	// E.l
		}
		else
		{ // BOND / INTERFACE
			longitudalStiffness       = 1.0;	// k.l
			bendingStiffness          = 0.7;	// k.b
			criticalCompressiveStrain = 1.0;	// E.c
			criticalTensileStrain     = 0.0002;	// E.l
		}
*/		lastIter_ = -10;


// checking

#define chk(x,y) std::cout << #y << ",      " << x << " = " << y << "\n"
		if(criticalTensileStrain > 0.0015) // E.l
		{ // CEMENT MATRIX
			static int cement=0;	++cement;
			static bool d1=true;
			if(d1)
			{
				d1=false;
				std::cout << "\nCEMENT:\n";
				chk("\tk.l",longitudalStiffness);               // k.l
				chk("\t\tk.b",bendingStiffness);                // k.b
				chk("\tk.t",torsionalStiffness);              // k.t
				chk("E.c",criticalCompressiveStrain);           // E.c
				chk("\tE.l",criticalTensileStrain);             // E.l
			}
		}
		else if(criticalTensileStrain > 0.0006) // E.l
		{ // AGGREGATE
			static int aggregate=0; ++aggregate;
			static bool d2=true;
			if(d2)
			{
				d2=false;
				std::cout << "\nAGGREGATE:\n";
				chk("\tk.l",longitudalStiffness);               // k.l
				chk("\t\tk.b",bendingStiffness);                // k.b
				chk("\tk.t",torsionalStiffness);              // k.t
				chk("E.c",criticalCompressiveStrain);           // E.c
				chk("\tE.l",criticalTensileStrain);             // E.l
			}
		}
		else
		{ // BOND / INTERFACE
			static int bond=0; ++bond;
			static bool d3=true;
			if(d3)
			{
				d3=false;
				std::cout << "\nBOND:\n";
				chk("\tk.l",longitudalStiffness);               // k.l
				chk("\t\tk.b",bendingStiffness);                // k.b
				chk("\tk.t",torsionalStiffness);              // k.t
				chk("E.c",criticalCompressiveStrain);           // E.c
				chk("\tE.l",criticalTensileStrain);             // E.l
			}
		}
#undef chk

        }
}


Real LatticeBeamParameters::strain()
{
        long int it = Omega::instance().getCurrentIteration();
        if( lastIter_ == it )
                return strain_;
        else
        {
                lastIter_ = it;
                strain_ = (length - initialLength) / initialLength;
                return strain_;
        }
}

YADE_PLUGIN((LatticeBeamParameters));

YADE_REQUIRE_FEATURE(PHYSPAR);

