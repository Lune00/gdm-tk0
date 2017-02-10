#ifndef _stringconv_hpp
#define _stringconv_hpp

#include <string>
#include "control.hpp"

namespace gdm
{	
	//
	bool stringToBool(string s)
	{
		if (s == "TRUE") return true;
		return false;
	}

	//
	string boolToString(bool b)
	{
		string s("FALSE");
		if (b) s = "TRUE";
		return s;
	}
	
	// 
	unsigned int stringToMode(string s)
	{
		if (s == "PRESSURE") return _PRESSURE;
		if (s == "VELOCITY") return _VELOCITY;
		if (s == "FORCE")    return _FORCE;
		return _VELOCITY;
	}

	//
	string mode2String(unsigned int which)
	{
		string mode("VELOCITY");
		if (which == _PRESSURE) mode = "PRESSURE";
		if (which == _VELOCITY) mode = "VELOCITY";
		if (which == _FORCE)    mode = "FORCE";
		return mode;
	}
}

#endif //_stringconv_hpp
