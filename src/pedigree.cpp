/*
 * pedigree.cpp
 *
 *  Created on: Nov 6, 2012
 *      Author: gpeng1
 */

#include <iostream>

#include "pedigree.h"

using namespace std;

individual::individual(int id, int fid, int mid, int gender)
{
	m_id = id;
	m_gender = gender;
	m_fid = fid;
	m_mid = mid;
}

bool individual::set_id(int id)
{
	m_id = id;
	return true;
}

bool individual::set_gender(int gender)
{
	m_gender = gender;
	return true;
}

bool individual::set_fid(int fid)
{
	m_fid = fid;
	return true;
}

bool individual::set_mid(int mid)
{
	m_mid=mid;
	return true;
}


pedigree::pedigree(const vector<individual> & mem)
{
	m_mem = mem;

	m_nInd = m_mem.size();

	fulfillPed();

	setRelation();

	checkPed();

	//cout<<*this<<endl;
}

/*
bool pedigree::fulfillPed()
{
	size_t numInd=m_mem.size();

	int addMem=-1;
	for(size_t i=0;i<numInd;i++)
	{
		if(m_mem[i].get_fid()==0 && m_mem[i].get_mid()>0)
		{
//			bool findMiss=false;
//			for(size_t j=0;j<numInd;j++)
//			{
//				if(m_mem[i].get_mid()==m_mem[j].get_mid() && m_mem[j].get_fid()!=0)
//				{
//					m_mem[i].set_fid(m_mem[j].get_fid());
//					findMiss=true;
//					break;
//				}
//			}
//			if(!findMiss)
//			{
//				individual indTmp(addMem,0,0,1);
//				m_mem[i].set_fid(addMem);
//				m_mem.push_back(indTmp);
//				addMem--;
//			}

			individual indTmp(addMem,0,0,1);
			m_mem[i].set_fid(addMem);
			m_mem.push_back(indTmp);
			addMem--;
		}

		if(m_mem[i].get_mid()==0 && m_mem[i].get_fid()>0)
		{
//			bool findMiss=false;
//			for(size_t j=0;j<numInd;j++)
//			{
//				if(m_mem[i].get_fid()==m_mem[j].get_fid() && m_mem[j].get_mid()!=0)
//				{
//					m_mem[i].set_mid(m_mem[j].get_mid());
//					findMiss=true;
//					break;
//				}
//			}
//			if(!findMiss)
//			{
//				individual indTmp(addMem,0,0,2);
//				m_mem[i].set_mid(addMem);
//				m_mem.push_back(indTmp);
//				addMem--;
//			}
			individual indTmp(addMem,0,0,2);
			m_mem[i].set_mid(addMem);
			m_mem.push_back(indTmp);
			addMem--;
		}
	}

	m_nIndAll=m_mem.size();

	return true;
}
*/

bool pedigree::fulfillPed()
{
	size_t numInd=m_mem.size();

	vector<int> spouse(numInd,0);
	vector<int> idAll(numInd);
	for(size_t i=0;i<numInd;i++){
		idAll[i] = m_mem[i].get_id();
	}

	int addMem=-1;
	for(size_t i=0;i<numInd;i++){
		int fid = m_mem[i].get_fid();
		int mid = m_mem[i].get_mid();
		if(fid==0 && mid>0){
			size_t idIndex=-1;
			int idSP = 0;
			for(size_t j=0;j<numInd;j++){
				if(mid==idAll[j]){
					idIndex = j;
					if(spouse[j]!=0){
						idSP = spouse[j];
					}
					break;
				}
			}
			if(idSP!=0){
				m_mem[i].set_fid(idSP);
			}
			else{
				individual indTmp(addMem,0,0,1);
				m_mem[i].set_fid(addMem);
				m_mem.push_back(indTmp);
				spouse[idIndex] = addMem;
				addMem--;
			}
		}

		if(mid==0 && fid>0){
			size_t idIndex=-1;
			int idSP = 0;
			for(size_t j=0;j<numInd;j++){
				if(fid==idAll[j]){
					idIndex = j;
					if(spouse[j]!=0){
						idSP = spouse[j];
					}
					break;
				}
			}
			if(idSP!=0){
				m_mem[i].set_mid(idSP);
			}
			else{
				individual indTmp(addMem,0,0,2);
				m_mem[i].set_mid(addMem);
				m_mem.push_back(indTmp);
				spouse[idIndex] = addMem;
				addMem--;
			}
		}
	}

	m_nIndAll=m_mem.size();

	//cout<<"1"<<endl;

	return true;
}

bool pedigree::setRelation()
{
	child.clear();
	parent.clear();
	spouse.clear();

	child.resize(m_nIndAll);
	parent.resize(m_nIndAll);
	spouse.resize(m_nIndAll);

	for(unsigned int i=0;i<m_nIndAll;i++)
	{
		int fid=m_mem[i].get_fid();
		int mid=m_mem[i].get_mid();
		if(fid!=0 and mid!=0)
		{
			size_t indF=100000, indM=100000;
			for(size_t j=0;j<m_nIndAll;j++)
			{
				if(fid==m_mem[j].get_id())
				{
					indF=j;
				}
				if(mid==m_mem[j].get_id())
				{
					indM=j;
				}
			}
			if(indF==100000 || indM==100000)
			{
				cout<<"Cannot find parents' id for sample "<<m_mem[i].get_id()<<endl;
			}
			else
			{
				parent[i].push_back(indF);
				parent[i].push_back(indM);

				child[indF].push_back(i);
				child[indM].push_back(i);

				bool find=false;
				for(size_t j=0;j<spouse[indF].size();j++)
				{
					if(spouse[indF][j]==indM)
					{
						find=true;
						break;
					}
				}
				if(!find)
				{
					spouse[indF].push_back(indM);
				}

				find=false;
				for(size_t j=0;j<spouse[indM].size();j++)
				{
					if(spouse[indM][j]==indF)
					{
						find=true;
						break;
					}
				}
				if(!find)
				{
					spouse[indM].push_back(indF);
				}
			}
		}
	}
	return true;
}

bool pedigree::checkPed()
{
	for(unsigned int i=0;i<m_nIndAll;i++)
	{
		/*
		if(parent[i].size()==0 && child[i].size()==0)
		{
			cout<<"Sample "<<m_mem[i].get_id()<<" is not in the family."<<endl;
		}
		*/

		if(m_mem[i].get_gender()==1)
		{
			for(size_t j=0;j<child[i].size();j++)
			{
				if(parent[child[i][j]][0]!=i)
				{
					cout<<"Gender error: sample "<<m_mem[child[i][j]].get_id()<<"'s mother "<<m_mem[i].get_id()<<" is male."<<endl;
					return false;
				}
			}
		}
		else if(m_mem[i].get_gender()==2)
		{
			for(size_t j=0;j<child[i].size();j++)
			{
				if(parent[child[i][j]][1]!=i)
				{
					cout<<"Gender error: sample "<<m_mem[child[i][j]].get_id()<<"'s father "<<m_mem[i].get_id()<<" is female."<<endl;
					return false;
				}
			}
		}
		else
		{
			if(child[i].size()>0)
			{
				if(child[i].size()==1)
				{
					if(parent[child[i][0]][0]==i)
					{
						m_mem[i].set_gender(1);
					}
					else
					{
						m_mem[i].set_gender(2);
					}
				}
				else
				{
					if(parent[child[i][0]][0]==i)
					{
						for(size_t j=1;j<child[i].size();j++)
						{
							if(i != parent[child[i][j]][0])
							{
								cout<<"Gender error: sample "<<m_mem[child[i][0]].get_id()<<"'s father is sample "<<m_mem[child[i][j]].get_id()<<"'s mother."<<endl;
								return false;
							}
						}
					}
					else
					{
						for(size_t j=1;j<child[i].size();j++)
						{
							if(i != parent[child[i][j]][1])
							{
								cout<<"Gender error: sample "<<m_mem[child[i][0]].get_id()<<"'s mother is sample "<<m_mem[child[i][j]].get_id()<<"'s father."<<endl;
								return false;
							}
						}
					}

				}
			}
		}
	}
	return true;
}

std::ostream & operator<< (std::ostream & os, const pedigree & ped)
{
	for(size_t i=0;i<ped.get_nIndAll();i++)
	{
		os<<ped[i].get_id()<<'\t'<<ped[i].get_fid()<<'\t'<<ped[i].get_mid()<<'\t'<<ped[i].get_gender()<<endl;
	}
	return os;
}

std::ofstream & operator<< (std::ofstream & of, const pedigree & ped)
{
	for(size_t i=0;i<ped.get_nIndAll();i++)
	{
		of<<ped[i].get_id()<<'\t'<<ped[i].get_fid()<<'\t'<<ped[i].get_mid()<<'\t'<<ped[i].get_gender()<<endl;
	}
	return of;
}
