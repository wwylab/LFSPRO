/*
 * normal.cpp
 *
 *  Created on: Jun 21, 2011
 *      Author: gpeng
 */

#include "normal.h"

using namespace std;

std::vector<std::string> split(const std::string & src, const std::string & tok, bool trim, std::string NulStr)
{
	std::vector <std::string> rlt;
	if(src.empty() || tok.empty())
	{
		return rlt;
	}


	size_t head=0, tail=0,length=0;
	while((tail=(src.find(tok,head))) != std::string::npos)
	{
		if((length=tail-head)>0)
		{
			rlt.push_back(src.substr(head,length));
		}
		else if(trim==false)
		{
			rlt.push_back(NulStr);
		}

		head=tail+1;
	}

	std::string endStr=src.substr(head);
	if(endStr.empty())
	{
		if(trim==false)
		{
			rlt.push_back(NulStr);
		}
	}
	else
	{
		rlt.push_back(endStr);
	}

	return rlt;
}


std::vector<std::string> split2(const std::string & src, const std::string & toks, bool trim, std::string NulStr)
{

	std::vector<std::string> rlt;
	if(src.empty() || toks.empty())
	{
		return  rlt;
	}


	size_t head=0, tail=0,length=0;
	while((tail=(src.find_first_of(toks,head))) != std::string::npos)
	{
		if((length=tail-head)>0)
		{
			rlt.push_back(src.substr(head,length));
		}
		else if(trim==false)
		{
			rlt.push_back(NulStr);
		}

		head=tail+1;
	}

	std::string endStr=src.substr(head);
	if(endStr.empty())
	{
		if(trim==false)
		{
			rlt.push_back(NulStr);
		}
	}
	else
	{
		rlt.push_back(endStr);
	}

	return rlt;
}

