/*
 * normal.h
 *
 *  Created on: Jun 21, 2011
 *      Author: gpeng
 */

#ifndef NORMAL_H_
#define NORMAL_H_

#include <string>
#include <vector>
#include <limits>

//split string into several strings according to tok
//src: the string to be split
//tok: split scr by tok
//trim, NulStr: when there is nothing between two tok, if trim is true, nothing will be added to the result, if trim is false, NulStr will be added to the result
std::vector<std::string>  split(const std::string & src, const std::string & tok, bool trim=false, std::string NulStr="");

//split string into several strings according to any characters in toks
//src: the string to be split
//tok: split scr by tok
//trim, NulStr: when there is nothing between two tok, if trim is true, nothing will be added to the result, if trim is false, NulStr will be added to the result
std::vector<std::string>  split2(const std::string & src, const std::string & toks, bool trim=false, std::string NulStr="");

template<class Type>
int binSearch(const std::vector<Type> & src, Type value)
{
	if(src.size()==0)
	{
		return -1;
	}
	size_t st=0,ed=src.size();
	while(true)
	{
		size_t mid=(st+ed)/2;
		if(mid==st)
		{
			if(src[st]==value)
			{
				return int(st);
			}
			else if(src[ed]==value)
			{
				return int(ed);
			}
			else
			{
				//return std::numeric_limits<unsigned int>::quiet_NaN();
				return -1;
			}
		}

		if(src[mid]==value)
		{
			return int(mid);
		}
		else
		{
			if(src[mid]>value)
			{
				ed=mid;
			}
			else
			{
				st=mid;
			}
		}
	}
}

#endif /* NORMAL_H_ */
