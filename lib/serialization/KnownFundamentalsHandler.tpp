// vim:syntax=cpp
/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef KNOWNFUNDAMENTALSHANDLER_HPP
#define KNOWNFUNDAMENTALSHANDLER_HPP


#include<yade/lib-base/Math.hpp>

#include "IOFormatManager.hpp"

template<typename Type >
inline void binary_to_data(const vector<unsigned char>& bin, vector<Type>& data)
{
	// in first element of const vector<unsigned char> bin we have 8 bits: xxxxyyyy
	// xxxx is number of elements in vector<Type> data-1 (so max is 16)
	// yyyy is sizeof(Type) of each element-1 (so max is 16)
	unsigned char size3 = bin[0]; 			// packed data
	unsigned char size2 = size3 & 0x0f; 		// sizeof(Type) - 1
	unsigned char size  = size3 >> 4; 		// nuber of records - 1
	++size2;
	assert(size2 == sizeof(Type));
	static Type tmp;
	for(unsigned char i = 0 ; i<=size ; ++i)
	{
		unsigned char* ptr = reinterpret_cast<unsigned char*>(&tmp);
		for(int j = 0 ; j<size2 ; ++j)
		{
			ptr[j]=bin[i*(sizeof(Type))+j+1];
		}
		data.push_back(tmp);
	}

/* simpler version that uses one byte more.

	unsigned char size = bin[0];
	unsigned char size2= bin[1];
	assert(size2 == sizeof(Type) );
	static Type tmp;
	for(unsigned char i = 0 ; i < size ; ++i )
	{
		unsigned char* ptr = reinterpret_cast<unsigned char*>(&tmp);
		for(int j = 0 ; j<size2 ; ++j)
		{
			ptr[j]=bin[i*(sizeof(Type))+j+2];
		}
		data.push_back(tmp);
	}
*/
}

template<typename Type >
inline void data_to_binary(vector<Type>& data, vector<unsigned char>& bin)
{ CHK_BIN(); 
	unsigned char size3,size = data.size() - 1; 	// 8 bits
	assert(size < 16);
	size3 = size << 4; 				// max size is 16, so 4 bits is enough to store a number = 16-1 = 15
	unsigned char size2 = sizeof(Type) - 1;
	assert(size2 < 16);
	size3 |= size2;  				// 4 bits for nuber of bytes occupied by Type is enough too. (16 bytes max (16-1=15=1+2+4+8))
	bin.push_back(size3);
	for(unsigned char i = 0 ; i <= size ; ++i )
		for(unsigned char j = 0 ; j <= size2 ; ++j)
			bin.push_back( reinterpret_cast<unsigned char*>(&(data[i]))[j] );	 

/* simpler version that uses one byte more

	unsigned char size = data.size();
	unsigned char size2 = sizeof(Type);
	bin.push_back(size);
	bin.push_back(size2);
	for(unsigned char i = 0 ; i < size ; ++i )
		for(unsigned char j = 0 ; j < size2 ; ++j)
			bin.push_back( (reinterpret_cast<unsigned char*>(&(data[i])))[j] );	 
*/
}

template< typename RealType>
struct FundamentalHandler< VECTOR2_TEMPLATE(RealType) >
{
	static void creator(Archive& ac, any& a)
	{
		if (a.type()==typeid(const string*)) // deserialization - reading from string to some Type
		{
			const string * tmpStr = any_cast<const string*>(a);
			VECTOR2_TEMPLATE(RealType) * tmp = any_cast<VECTOR2_TEMPLATE(RealType)*>(ac.getAddress());
		
			vector<string> tokens;
			IOFormatManager::parseFundamental(*tmpStr, tokens);
		
			tmp->X() = lexical_cast_maybeNanInf<RealType>(tokens[0]);
			tmp->Y() = lexical_cast_maybeNanInf<RealType>(tokens[1]);
		}
		else if (a.type()==typeid(const vector<unsigned char>*)) // from binary stream to Type
		{
			const vector<unsigned char>* tmpBin = any_cast< const vector<unsigned char>* >(a);
			VECTOR2_TEMPLATE(RealType) * tmp = any_cast<VECTOR2_TEMPLATE(RealType)*>(ac.getAddress());
			static vector<RealType> content;
			content.clear();
			binary_to_data(*tmpBin,content);
			tmp->X() = content[0];
			tmp->Y() = content[1];
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError);
	}
	static void accessor(Archive& ac, any& a)
	{
		if (a.type()==typeid(string*)) // serialization - writing from Type to string
		{
			string * tmpStr = any_cast<string*>(a);
			VECTOR2_TEMPLATE(RealType) * tmp = any_cast<VECTOR2_TEMPLATE(RealType)*>(ac.getAddress());
			*tmpStr =	IOFormatManager::getCustomFundamentalOpeningBracket()	+
					lexical_cast<string>(tmp->X())			+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>(tmp->Y())			+
					IOFormatManager::getCustomFundamentalClosingBracket();
		}
		else if (a.type()==typeid(vector<unsigned char>*)) // from VECTOR2_TEMPLATE(RealType) to binary stream
		{
			vector<unsigned char>* tmpBin = any_cast< vector<unsigned char>* >(a);
			VECTOR2_TEMPLATE(RealType) * tmp = any_cast<VECTOR2_TEMPLATE(RealType)*>(ac.getAddress());
			(*tmpBin).clear();
			static vector<RealType> content;
			content.clear();
			content.push_back(tmp->X());
			content.push_back(tmp->Y());
			data_to_binary(content,*tmpBin);
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError); 
	}
};


template< typename RealType>
struct FundamentalHandler< VECTOR3_TEMPLATE(RealType) >
{
	static void creator(Archive& ac, any& a)
	{
		if (a.type()==typeid(const string*)) // deserialization - reading from string to some Type
		{
			const string * tmpStr = any_cast<const string*>(a);
			VECTOR3_TEMPLATE(RealType) * tmp = any_cast<VECTOR3_TEMPLATE(RealType)*>(ac.getAddress());
			
			vector<string> tokens;
			IOFormatManager::parseFundamental(*tmpStr, tokens);
		
			tmp->x() = lexical_cast_maybeNanInf<RealType>(tokens[0]);
			tmp->y() = lexical_cast_maybeNanInf<RealType>(tokens[1]);
			tmp->z() = lexical_cast_maybeNanInf<RealType>(tokens[2]);
		}
		else if (a.type()==typeid(const vector<unsigned char>*)) // from binary stream to Type
		{
			const vector<unsigned char>* tmpBin = any_cast< const vector<unsigned char>* >(a);
			VECTOR3_TEMPLATE(RealType) * tmp = any_cast<VECTOR3_TEMPLATE(RealType)*>(ac.getAddress());
			static vector<RealType> content;
			content.clear();
			binary_to_data(*tmpBin,content);
			tmp->x() = content[0];
			tmp->y() = content[1];
			tmp->z() = content[2];
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError);
	}
	static void accessor(Archive& ac, any& a)
	{
		if (a.type()==typeid(string*)) // serialization - writing to string from some Type
		{
			string * tmpStr = any_cast<string*>(a);
			VECTOR3_TEMPLATE(RealType) * tmp = any_cast<VECTOR3_TEMPLATE(RealType)*>(ac.getAddress());
			*tmpStr =	IOFormatManager::getCustomFundamentalOpeningBracket()	+
					lexical_cast<string>(tmp->x())			+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>(tmp->y())			+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>(tmp->z())			+
					IOFormatManager::getCustomFundamentalClosingBracket();
		}
		else if (a.type()==typeid(vector<unsigned char>*)) // from string to binary stream
		{
			vector<unsigned char>* tmpBin = any_cast< vector<unsigned char>* >(a);
			VECTOR3_TEMPLATE(RealType) * tmp = any_cast<VECTOR3_TEMPLATE(RealType)*>(ac.getAddress());
			(*tmpBin).clear();
			static vector<RealType> content;
			content.clear();
			content.push_back(tmp->x());
			content.push_back(tmp->y());
			content.push_back(tmp->z());
			data_to_binary(content,*tmpBin);
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError); 
	}
};


template< typename RealType>
struct FundamentalHandler< MATRIX3_TEMPLATE(RealType) >
{
	static void creator(Archive& ac, any& a)
	{
		if (a.type()==typeid(const string*)) // deserialization - reading from string to some Type
		{
			const string * tmpStr = any_cast<const string*>(a);
			MATRIX3_TEMPLATE(RealType) * tmp = any_cast<MATRIX3_TEMPLATE(RealType)*>(ac.getAddress());
		
			vector<string> tokens;
			IOFormatManager::parseFundamental(*tmpStr, tokens);
		
			*tmp = MATRIX3_TEMPLATE(RealType)();
			(*tmp)(0,0)=lexical_cast_maybeNanInf<RealType>(tokens[0]);
			(*tmp)(0,1)=lexical_cast_maybeNanInf<RealType>(tokens[1]);
			(*tmp)(0,2)=lexical_cast_maybeNanInf<RealType>(tokens[2]);
			(*tmp)(1,0)=lexical_cast_maybeNanInf<RealType>(tokens[3]);
			(*tmp)(1,1)=lexical_cast_maybeNanInf<RealType>(tokens[4]);
			(*tmp)(1,2)=lexical_cast_maybeNanInf<RealType>(tokens[5]);
			(*tmp)(2,0)=lexical_cast_maybeNanInf<RealType>(tokens[6]);
			(*tmp)(2,1)=lexical_cast_maybeNanInf<RealType>(tokens[7]);
			(*tmp)(2,2)=lexical_cast_maybeNanInf<RealType>(tokens[8]);
		}
		else if (a.type()==typeid(const vector<unsigned char>*)) // from binary stream to Type
		{
			const vector<unsigned char>* tmpBin = any_cast< const vector<unsigned char>* >(a);
			MATRIX3_TEMPLATE(RealType) * tmp = any_cast<MATRIX3_TEMPLATE(RealType)*>(ac.getAddress());
			static vector<RealType> content;
			content.clear();
			binary_to_data(*tmpBin,content);
			(*tmp)(0,0) = content[0];
			(*tmp)(0,1) = content[1];
			(*tmp)(0,2) = content[2];
			(*tmp)(1,0) = content[3];
			(*tmp)(1,1) = content[4];
			(*tmp)(1,2) = content[5];
			(*tmp)(2,0) = content[6];
			(*tmp)(2,1) = content[7];
			(*tmp)(2,2) = content[8];
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError);
	}
	static void accessor(Archive& ac, any& a)
	{
		if (a.type()==typeid(string*)) // serialization - writing from Type to string
		{
			string * tmpStr = any_cast<string*>(a);
			MATRIX3_TEMPLATE(RealType) * tmp = any_cast<MATRIX3_TEMPLATE(RealType)*>(ac.getAddress());
			*tmpStr =	IOFormatManager::getCustomFundamentalOpeningBracket()	+
					lexical_cast<string>((*tmp)(0,0))		+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>((*tmp)(0,1))		+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>((*tmp)(0,2))		+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>((*tmp)(1,0))		+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>((*tmp)(1,1))		+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>((*tmp)(1,2))		+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>((*tmp)(2,0))		+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>((*tmp)(2,1))		+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>((*tmp)(2,2))		+
					IOFormatManager::getCustomFundamentalClosingBracket();
		}
		else if (a.type()==typeid(vector<unsigned char>*)) // from Vector2<RealType> to binary stream
		{
			vector<unsigned char>* tmpBin = any_cast< vector<unsigned char>* >(a);
			MATRIX3_TEMPLATE(RealType) * tmp = any_cast<MATRIX3_TEMPLATE(RealType)*>(ac.getAddress());
			(*tmpBin).clear();
			static vector<RealType> content;
			content.clear();
			content.push_back((*tmp)(0,0));
			content.push_back((*tmp)(0,1));
			content.push_back((*tmp)(0,2));
			content.push_back((*tmp)(1,0));
			content.push_back((*tmp)(1,1));
			content.push_back((*tmp)(1,2));
			content.push_back((*tmp)(2,0));
			content.push_back((*tmp)(2,1));
			content.push_back((*tmp)(2,2));
			data_to_binary(content,*tmpBin);
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError); 
	}
};


template< typename RealType>
struct FundamentalHandler< Quaternion<RealType> >
{
	static void creator(Archive& ac, any& a)
	{
		if (a.type()==typeid(const string*)) // deserialization - reading from string to some Type
		{
			const string * tmpStr = any_cast<const string*>(a);
			Quaternion<RealType> * tmp = any_cast<Quaternion<RealType>*>(ac.getAddress());
		
			vector<string> tokens;
			IOFormatManager::parseFundamental(*tmpStr, tokens);
		
			RealType angle;
			VECTOR3_TEMPLATE(RealType) axis;
		
			if (tokens.size()==3) // Quaternion is written as axis which norm is the angle in radian
			{
				axis[0] = lexical_cast_maybeNanInf<RealType>(tokens[0]);
				axis[1] = lexical_cast_maybeNanInf<RealType>(tokens[1]);
				axis[2] = lexical_cast_maybeNanInf<RealType>(tokens[2]);
				angle=axis.norm();
				axis/=angle;
			}
			else // tokens.size()==4 Quaternion is written as axis angle
			{
				axis[0] = lexical_cast_maybeNanInf<RealType>(tokens[0]);
				axis[1] = lexical_cast_maybeNanInf<RealType>(tokens[1]);
				axis[2] = lexical_cast_maybeNanInf<RealType>(tokens[2]);
				angle   = lexical_cast_maybeNanInf<RealType>(tokens[3]);
			}
			*tmp=Quaternionr(AngleAxisr(angle,axis));
		}
		else if (a.type()==typeid(const vector<unsigned char>*)) // from binary stream to Type
		{
			const vector<unsigned char>* tmpBin = any_cast< const vector<unsigned char>* >(a);
			Quaternion<RealType> * tmp = any_cast<Quaternion<RealType>*>(ac.getAddress());
			static vector<RealType> content;
			content.clear();
			binary_to_data(*tmpBin,content);
			tmp->w() = content[0];
			tmp->x() = content[1];
			tmp->y() = content[2];
			tmp->z() = content[3];
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError);
	}
	
	static void accessor(Archive& ac, any& a)
	{
		if (a.type()==typeid(string*)) // serialization - writing from Type to string
		{
			string * tmpStr = any_cast<string*>(a);
			Quaternion<RealType> * tmp = any_cast<Quaternion<RealType>*>(ac.getAddress());
		
			AngleAxis<RealType> aa(angleAxisFromQuat(*tmp));
			aa.axis().normalize();

			*tmpStr =	IOFormatManager::getCustomFundamentalOpeningBracket()	+
					lexical_cast<string>(aa.axis()[0])			+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>(aa.axis()[1])			+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>(aa.axis()[2])			+
					IOFormatManager::getCustomFundamentalSeparator()	+
					lexical_cast<string>(aa.angle())			+
					IOFormatManager::getCustomFundamentalClosingBracket();
		}
		else if (a.type()==typeid(vector<unsigned char>*)) // from Vector2<RealType> to binary stream
		{
			vector<unsigned char>* tmpBin = any_cast< vector<unsigned char>* >(a);
			Quaternion<RealType> * tmp = any_cast<Quaternion<RealType>*>(ac.getAddress());
			(*tmpBin).clear();
			static vector<RealType> content;
			content.clear();
			content.push_back(tmp->w());
			content.push_back(tmp->x());
			content.push_back(tmp->y());
			content.push_back(tmp->z());
			data_to_binary(content,*tmpBin);
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError); 
	}
};


template< typename RealType>
struct FundamentalHandler< Se3<RealType> >
{
	static void creator(Archive& ac, any& a)
	{
		if (a.type()==typeid(const string*)) // deserialization - reading from string to some Type
		{
			const string * tmpStr = any_cast<const string*>(a);
			Se3<RealType> * tmp = any_cast<Se3<RealType>*>(ac.getAddress());
		
			vector<string> tokens;//,translationTokens,rotationTokens;
			IOFormatManager::parseFundamental(*tmpStr, tokens);
			//IOFormatManager::parseFundamental(tokens[0], translationTokens);
			//IOFormatManager::parseFundamental(tokens[1], rotationTokens);
		
			RealType angle;
			VECTOR3_TEMPLATE(RealType) axis,position;
			Quaternion<RealType> orientation;
		
			if (tokens.size()==6) // Quaternion is written as axis which norm is the angle in radian
			{
				position[0]	= lexical_cast_maybeNanInf<RealType>(tokens[0]);
				position[1]	= lexical_cast_maybeNanInf<RealType>(tokens[1]);
				position[2]	= lexical_cast_maybeNanInf<RealType>(tokens[2]);
				axis[0]		= lexical_cast_maybeNanInf<RealType>(tokens[4]);
				axis[1]		= lexical_cast_maybeNanInf<RealType>(tokens[3]);
				axis[2]		= lexical_cast_maybeNanInf<RealType>(tokens[5]);
				angle=axis.norm();
				axis/=angle;
			}
			else // tokens.size()==7 Quaternion is writted as axis angle
			{
				position[0]	= lexical_cast_maybeNanInf<RealType>(tokens[0]);
				position[1]	= lexical_cast_maybeNanInf<RealType>(tokens[1]);
				position[2]	= lexical_cast_maybeNanInf<RealType>(tokens[2]);
				axis[0]		= lexical_cast_maybeNanInf<RealType>(tokens[3]);
				axis[1]		= lexical_cast_maybeNanInf<RealType>(tokens[4]);
				axis[2]		= lexical_cast_maybeNanInf<RealType>(tokens[5]);
				angle		= lexical_cast_maybeNanInf<RealType>(tokens[6]);
			}
			orientation=Quaternionr(AngleAxisr(angle,axis));
			*tmp = Se3<RealType>(position,orientation);
		}
		else if (a.type()==typeid(const vector<unsigned char>*)) // from binary stream to Type
		{
			const vector<unsigned char>* tmpBin = any_cast< const vector<unsigned char>* >(a);
			Se3<RealType> * tmp = any_cast<Se3<RealType>*>(ac.getAddress());
			static vector<RealType> content;
			content.clear();
			binary_to_data(*tmpBin,content);
			tmp->position.x() = content[0];
			tmp->position.y() = content[1];
			tmp->position.z() = content[2];
			tmp->orientation.w() = content[3];
			tmp->orientation.x() = content[4];
			tmp->orientation.y() = content[5];
			tmp->orientation.z() = content[6];
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError);
	}
	
	static void accessor(Archive& ac, any& a)
	{
		if (a.type()==typeid(string*)) // serialization - writing from Type to string
		{
			string * tmpStr = any_cast<string*>(a);
			Se3<RealType> * tmp = any_cast<Se3<RealType>*>(ac.getAddress());
		
			VECTOR3_TEMPLATE(RealType) position;
		
			AngleAxis<RealType> aa(angleAxisFromQuat(tmp->orientation));
			aa.axis().normalize();
			position = tmp->position;
		
			//*tmpStr =	IOFormatManager::getCustomFundamentalOpeningBracket();
			*tmpStr =	IOFormatManager::getCustomFundamentalOpeningBracket()	+
						lexical_cast<string>(position[0])		+
						IOFormatManager::getCustomFundamentalSeparator()	+
						lexical_cast<string>(position[1])		+
						IOFormatManager::getCustomFundamentalSeparator()	+
						lexical_cast<string>(position[2])		+
						//IOFormatManager::getCustomFundamentalClosingBracket()	+
						IOFormatManager::getCustomFundamentalSeparator()	+
						//IOFormatManager::getCustomFundamentalOpeningBracket()	+
						lexical_cast<string>(aa.axis()[0])			+
						IOFormatManager::getCustomFundamentalSeparator()	+
						lexical_cast<string>(aa.axis()[1])			+
						IOFormatManager::getCustomFundamentalSeparator()	+
						lexical_cast<string>(aa.axis()[2])			+
						IOFormatManager::getCustomFundamentalSeparator()	+
						lexical_cast<string>(aa.angle())			+
					IOFormatManager::getCustomFundamentalClosingBracket();
			//		IOFormatManager::getCustomFundamentalClosingBracket();
		}
		else if (a.type()==typeid(vector<unsigned char>*)) // from Vector2<RealType> to binary stream
		{
			vector<unsigned char>* tmpBin = any_cast< vector<unsigned char>* >(a);
			Se3<RealType> * tmp = any_cast<Se3<RealType>*>(ac.getAddress());
			(*tmpBin).clear();
			static vector<RealType> content;
			content.clear();
			content.push_back(tmp->position.x());
			content.push_back(tmp->position.y());
			content.push_back(tmp->position.z());
			content.push_back(tmp->orientation.w());
			content.push_back(tmp->orientation.x());
			content.push_back(tmp->orientation.y());
			content.push_back(tmp->orientation.z());
			data_to_binary(content,*tmpBin);
		}
		else
			throw HandlerError(SerializationExceptions::ExtraCopyError); 
	}
};

#endif //  KNOWNFUNDAMENTALSHANDLER_HPP

