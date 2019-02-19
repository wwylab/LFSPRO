/*
 * Mendel.cpp
 *
 *  Created on: Nov 6, 2012
 *      Author: gpeng1
 */

#include "Mendel.h"

using namespace std;

Mendel::Mendel(const vector<individual> & mem, const vector<double> & genoProb, const vector<unsigned int> & numAllele, double mRate, bool store): pedigree(mem)
{
	m_genoProb=genoProb;
	m_pccp=PCCP(numAllele,mRate,store);
	m_error=!check_nGeno();
	m_chromSex=false;
}

Mendel::Mendel(const vector<individual> & mem, const vector<double> & genoProb, const vector<unsigned int> & numAllele, const vector<double> & mRate, bool store): pedigree(mem)
{
	m_genoProb=genoProb;
	m_pccp=PCCP(numAllele,mRate,store);
	m_error=!check_nGeno();
	m_chromSex=false;
}

Mendel::Mendel(const vector<individual> & mem, const vector<double> & genoProb, const vector<unsigned int> & numAllele, double mRate, const vector<int> & chromType, bool store): pedigree(mem)
{
	m_genoProb=genoProb;
	m_error=!check_nGeno();
	m_chromSex=false;
	for(size_t i=0;i<chromType.size();i++)
	{
		if(chromType[i]!=0)
		{
			m_chromSex=true;
		}
	}
	if(m_chromSex)
	{
		m_pccpMale=PCCP(numAllele,mRate,chromType,store);
		vector<int> chromTypeFemale=chromType;
		for(size_t i=0;i<chromType.size();i++)
		{
			if(chromTypeFemale[i]==1)
			{
				chromTypeFemale[i]=2;
			}
		}
		m_pccpFemale=PCCP(numAllele,mRate,chromTypeFemale,store);
	}
	else
	{
		m_pccp=PCCP(numAllele,mRate,chromType,store);
	}
}

Mendel::Mendel(const vector<individual> & mem, const vector<double> & genoProb, const vector<unsigned int> & numAllele, const vector<double> & mRate, const vector<int> & chromType, bool store): pedigree(mem)
{
	m_genoProb=genoProb;
	m_error=!check_nGeno();
	m_chromSex=false;
	for(size_t i=0;i<chromType.size();i++)
	{
		if(chromType[i]!=0)
		{
			m_chromSex=true;
		}
	}
	if(m_chromSex)
	{
		m_pccpMale=PCCP(numAllele,mRate,chromType,store);
		vector<int> chromTypeFemale=chromType;
		for(size_t i=0;i<chromType.size();i++)
		{
			if(chromTypeFemale[i]==1)
			{
				chromTypeFemale[i]=2;
			}
		}
		m_pccpFemale=PCCP(numAllele,mRate,chromTypeFemale,store);
	}
	else
	{
		m_pccp=PCCP(numAllele,mRate,chromType,store);
	}
}

bool Mendel::set_genoProb(const vector<double> & genoProb)
{
	if(genoProb.size()==m_pccp.get_nGenotype())
	{
		m_genoProb=genoProb;
		m_error=false;
		return true;
	}
	else
	{
		return false;
	}
}

bool Mendel::set_pccp(const PCCP & pccp)
{
	if(pccp.get_nGenotype()==m_genoProb.size())
	{
		m_pccp=pccp;
		m_error=false;
		return true;
	}
	else
	{
		return false;
	}
}

bool Mendel::set_par(const PCCP & pccp, const vector<double> & genoProb)
{
	if(pccp.get_nGenotype() != genoProb.size())
	{
		return false;
	}
	else
	{
		m_pccp=pccp;
		m_genoProb=genoProb;
		m_error=false;
		return true;
	}
}


dMatrix<double> Mendel::calPostProb(const dMatrix<double> & lk, calMethod method)
{
	if(m_error)
	{
		//cout<<"Error in parameters!"<<endl;
		return dMatrix<double> (1,1,-1);
	}

	if(lk.get_column() != m_genoProb.size())
	{
		//cout<<"Numbers of columns in the likelihood matrix doesn't match the number of genotype."<<endl;
		return dMatrix<double> (1,1,-1);
	}

	if(lk.get_row() != m_nInd)
	{
		//cout<<"Number of rows in the likelihood matrix doesn;t match the number of samples."<<endl;
		return dMatrix<double> (1,1,-1);
	}

	if(m_nInd==m_nIndAll)
	{
		m_lk=lk;
		if(method==Peeling)
		{
			return calPostProbPeeling();
		}
		else if(method==BN)
		{
			return calPostProbBN();
		}
		else
		{
			return calPostProbMCMC();
		}
	}
	else
	{
		dMatrix<double> lkNew(m_nIndAll,lk.get_column());
		for(unsigned int i=0;i<m_nInd;i++)
		{
			for(unsigned int j=0;j<lk.get_column();j++)
			{
				lkNew(i,j)=lk(i,j);
			}
		}
		for(unsigned int i=m_nInd;i<m_nIndAll;i++)
		{
			for(unsigned int j=0;j<lk.get_column();j++)
			{
				lkNew(i,j)=1;
			}
		}
		m_lk=lkNew;
		if(method==Peeling)
		{
			dMatrix<double> rltTmp=calPostProbPeeling();
			if(rltTmp.get_column()==1){
				return rltTmp;
			}
			dMatrix<double> rlt(lk.get_row(),lk.get_column());
			for(unsigned int i=0;i<lk.get_row();i++)
			{
				for(unsigned int j=0;j<lk.get_column();j++)
				{
					rlt(i,j)=rltTmp(i,j);
				}
			}
			return rlt;
		}
		else if(method==BN)
		{
			dMatrix<double> rltTmp=calPostProbBN();
			dMatrix<double> rlt(lk.get_row(),lk.get_column());
			for(unsigned int i=0;i<lk.get_row();i++)
			{
				for(unsigned int j=0;j<lk.get_column();j++)
				{
					rlt(i,j)=rltTmp(i,j);
				}
			}
			return rlt;
		}
		else
		{
			dMatrix<double> rltTmp=calPostProbMCMC();
			dMatrix<double> rlt(lk.get_row(),lk.get_column());
			for(unsigned int i=0;i<lk.get_row();i++)
			{
				for(unsigned int j=0;j<lk.get_column();j++)
				{
					rlt(i,j)=rltTmp(i,j);
				}
			}
			return rlt;
		}
	}
}

dMatrix<double> Mendel::calPostProb(const dMatrix<double> & lk, const vector<int> & sampleId, calMethod method)
{
	if(m_error)
	{
		//cout<<"Error in parameters!"<<endl;
		return dMatrix<double> (1,1,-1);
	}

	if(lk.get_column() != m_genoProb.size())
	{
		//cout<<"Numbers of columns in the likelihood matrix doesn't match the number of genotype."<<endl;
		return dMatrix<double> (1,1,-1);
	}

	if(lk.get_row() != m_nInd)
	{
		//cout<<"Number of rows in the likelihood matrix doesn;t match the number of samples."<<endl;
		return dMatrix<double> (1,1,-1);
	}

	vector<unsigned int> sampleIndex;
	for(size_t i=0;i<sampleId.size();i++)
	{
		bool find=false;
		for(size_t j=0;j<m_mem.size();j++)
		{
			if(m_mem[j].get_id()==sampleId[i])
			{
				sampleIndex.push_back(j);
				find = true;
				break;
			}
		}
		if(!find)
		{
			//cout<<"Warning: cannot find sample id: "<<sampleId[i]<<endl;
		}
	}

	if(sampleIndex.size()==0)
	{
		return dMatrix<double> (1,1,-1);
	}

	if(m_nInd==m_nIndAll)
	{
		m_lk=lk;
		if(method==Peeling)
		{
			return calPostProbPeeling(sampleIndex);
		}
		else if(method==BN)
		{
			dMatrix<double> rltTmp=calPostProbBN();
			dMatrix<double> rlt(sampleIndex.size(),lk.get_column());
			for(size_t i=0;i<sampleIndex.size();i++)
			{
				for(unsigned int j=0;j<lk.get_column();j++)
				{
					rlt(i,j)=rltTmp(sampleIndex[i],j);
				}
			}
			return rlt;
		}
		else
		{
			dMatrix<double> rltTmp = calPostProbMCMC();
			dMatrix<double> rlt(sampleIndex.size(),lk.get_column());
			for(size_t i=0;i<sampleIndex.size();i++)
			{
				for(unsigned int j=0;j<lk.get_column();j++)
				{
					rlt(i,j)=rltTmp(sampleIndex[i],j);
				}
			}
			return rlt;
		}
	}
	else
	{
		dMatrix<double> lkNew(m_nIndAll,lk.get_column());
		for(unsigned int i=0;i<m_nInd;i++)
		{
			for(unsigned int j=0;j<lk.get_column();j++)
			{
				lkNew(i,j)=lk(i,j);
			}
		}
		for(unsigned int i=m_nInd;i<m_nIndAll;i++)
		{
			for(unsigned int j=0;j<lk.get_column();j++)
			{
				lkNew(i,j)=1;
			}
		}
		m_lk=lkNew;
		if(method==Peeling)
		{
			return calPostProbPeeling(sampleIndex);
		}
		else if(method==BN)
		{
			dMatrix<double> rltTmp=calPostProbBN();
			dMatrix<double> rlt(sampleIndex.size(),lk.get_column());
			for(size_t i=0;i<sampleIndex.size();i++)
			{
				for(unsigned int j=0;j<lk.get_column();j++)
				{
					rlt(i,j)=rltTmp(sampleIndex[i],j);
				}
			}
			return rlt;
		}
		else
		{
			dMatrix<double> rltTmp=calPostProbMCMC();
			dMatrix<double> rlt(sampleIndex.size(),lk.get_column());
			for(size_t i=0;i<sampleIndex.size();i++)
			{
				for(unsigned int j=0;j<lk.get_column();j++)
				{
					rlt(i,j)=rltTmp(sampleIndex[i],j);
				}
			}
			return rlt;
		}
	}
}


bool Mendel::check_nGeno()
{
	if(m_genoProb.size()!=m_pccp.get_nGenotype())
	{
		//cout<<"Number of genotype in genotype probability ("<<m_genoProb.size()<<") and number of genotype calculated from number of alleles ("<<m_pccp.get_nGenotype()<<") is not match."<<endl;
		return false;
	}
	m_nGenotype=m_pccp.get_nGenotype();
	return true;
}

dMatrix<double> Mendel::calPostProbPeeling()
{
	dMatrix<double> antProb(m_nIndAll,m_nGenotype,-1);
	dMatrix<bool> checkLoop(m_nIndAll,m_nGenotype,false);
	m_loop=false;

	vector<dMatrix<double> > posProb;
	for(unsigned int i=0;i<m_nIndAll;i++)
	{
		dMatrix<double> dmTmp(m_nIndAll,m_nGenotype,-1);
		posProb.push_back(dmTmp);
	}

	//prepare anterior probability
	//set founder's anterior probability
	for(unsigned int i=0;i<m_nIndAll;i++)
	{
		if(parent[i].size()==0)
		{
			for(unsigned int j=0;j<m_nGenotype;j++)
			{
				antProb(i,j)=m_genoProb[j];
			}
		}
	}

	dMatrix<double> postProb(m_nIndAll,m_nGenotype);

	for(unsigned int i=0;i<m_nIndAll;i++)
	{
		for(unsigned int j=0;j<m_nGenotype;j++)
		{
			double postTmp=1;
			for(size_t k=0;k<spouse[i].size();k++)
			{
				postTmp=postTmp*calPosProb(i,spouse[i][k],j,antProb,posProb,checkLoop);
			}
			postTmp=postTmp*m_lk(i,j)*calAntProb(i,j,antProb,posProb,checkLoop);
			postProb(i,j)=postTmp;
		}

		double sTmp=postProb.sum_row(i);
		if(sTmp<=0)
		{
			if(m_loop){
				//cout<<"There are loops in the pedigree C!"<<endl;
				return dMatrix<double> (1,1,-2);
			}
			//cout<<"Error during peeling."<<endl;
			return dMatrix<double> (1,1,-1);
		}

		for(unsigned int j=0;j<m_nGenotype;j++)
		{
			postProb(i,j)=postProb(i,j)/sTmp;
		}
	}

	if(m_loop){
		//cout<<"There are loops in the pedigree C!"<<endl;
		return dMatrix<double> (1,1,-2);
	}

	return postProb;
}

dMatrix<double> Mendel::calPostProbPeeling(const vector<unsigned int> & sampleIndex)
{
	dMatrix<double> antProb(m_nIndAll,m_nGenotype,-1);
	dMatrix<bool> checkLoop(m_nIndAll,m_nGenotype,false);
	m_loop=false;

	vector<dMatrix<double> > posProb;
	for(unsigned int i=0;i<m_nIndAll;i++)
	{
		dMatrix<double> dmTmp(m_nIndAll,m_nGenotype,-1);
		posProb.push_back(dmTmp);
	}

	//prepare anterior probability
	//set founder's anterior probability
	for(unsigned int i=0;i<m_nIndAll;i++)
	{
		if(parent[i].size()==0)
		{
			for(unsigned int j=0;j<m_nGenotype;j++)
			{
				antProb(i,j)=m_genoProb[j];
			}
		}
	}

	dMatrix<double> postProb(sampleIndex.size(),m_nGenotype);
	for(unsigned int i=0;i<sampleIndex.size();i++)
	{
		for(unsigned int j=0;j<m_nGenotype;j++)
		{
			double postTmp=1;
			for(size_t k=0;k<spouse[sampleIndex[i]].size();k++)
			{
				postTmp=postTmp*calPosProb(sampleIndex[i],spouse[sampleIndex[i]][k],j,antProb,posProb,checkLoop);
			}
			postTmp=postTmp*m_lk(sampleIndex[i],j)*calAntProb(sampleIndex[i],j,antProb,posProb,checkLoop);
			postProb(i,j)=postTmp;
		}

		double sTmp=postProb.sum_row(i);
		if(sTmp<=0)
		{
			/*
			//cout<< postProb <<endl;
			for(size_t j=0;j<sampleIndex.size();j++)
			{
				//cout<<sampleIndex[j]<<'\t';
				for(size_t k=0;k<m_lk.get_column();k++)
				{
					//cout<<m_lk(sampleIndex[j],k)<<'\t';
				}
				//cout<<endl;
			}
			//cout<< *this <<endl;
			*/
			if(m_loop){
				//cout<<"There are loops in the pedigree C!"<<endl;
				return dMatrix<double> (1,1,-2);
			}
			//cout<<"Error during peeling."<<endl;
			return dMatrix<double> (1,1,-1);
		}

		for(unsigned int j=0;j<m_nGenotype;j++)
		{
			postProb(i,j)=postProb(i,j)/sTmp;
		}
	}

	if(m_loop){
		//cout<<"There are loops in the pedigree C!"<<endl;
		return dMatrix<double> (1,1,-2);
	}

	return postProb;
}

dMatrix<double> Mendel::calPostProbBN()
{
	dMatrix<double> postProb(m_nIndAll,m_nGenotype,0);

	//initial genotype
	vector<unsigned int> genotype(m_nIndAll,0);
	vector<double> indProb(m_nIndAll,0);

	while(true)
	{
		for(unsigned int i=0;i<m_nIndAll;i++)
		{
			if(parent[i].size()==0)
			{
				indProb[i]=m_genoProb[genotype[i]]*m_lk(i,genotype[i]);
			}
			else
			{
				if(m_chromSex)
				{
					if(m_mem[i].get_gender()==1)
					{
						indProb[i]=m_pccpMale(genotype[i],genotype[parent[i][0]],genotype[parent[i][1]])*m_lk(i,genotype[i]);
					}
					else
					{
						indProb[i]=m_pccpMale(genotype[i],genotype[parent[i][0]],genotype[parent[i][1]])*m_lk(i,genotype[i]);
					}
				}
				else
				{
					indProb[i]=m_pccp(genotype[i],genotype[parent[i][0]],genotype[parent[i][1]])*m_lk(i,genotype[i]);
				}
			}
		}

		double probAll=100000000;
		for(unsigned int i=0;i<m_nIndAll;i++)
		{
			probAll=probAll*indProb[i];
		}

		for(unsigned int i=0;i<m_nIndAll;i++)
		{
			postProb(i,genotype[i])=postProb(i,genotype[i])+probAll;
		}

		unsigned int index=0;
		while(index<m_nIndAll)
		{
			genotype[index]=genotype[index]+1;
			if(genotype[index]==m_nGenotype)
			{
				genotype[index]=0;
				index++;
			}
			else
			{
				break;
			}
		}

		if(index==m_nIndAll)
		{
			break;
		}
	}

	for(unsigned int i=0;i<m_nIndAll;i++)
	{
		double sTmp=postProb.sum_row(i);
		if(sTmp<=0)
		{
			//cout<<"Error during BN."<<endl;
			return dMatrix<double> (1,1,-1);
		}

		for(unsigned int j=0;j<m_nGenotype;j++)
		{
			postProb(i,j)=postProb(i,j)/sTmp;
		}
	}
	return postProb;
}

dMatrix<double> Mendel::calPostProbMCMC()
{
	dMatrix<double> postProb(m_nIndAll,m_nGenotype);
	return postProb;
}

double Mendel::calAntProb(unsigned int iInd, unsigned int iGeno, dMatrix<double> & antProb, vector<dMatrix<double> > & posProb, dMatrix<bool> & checkLoop)
{
	if(antProb(iInd,iGeno)>=0)
	{
		return antProb(iInd,iGeno);
	}

	if(checkLoop(iInd,iGeno)){
		m_loop = true;
		return 0;
	}

	checkLoop(iInd,iGeno) = true;

	unsigned int p1 = parent[iInd][0];
	unsigned int p2 = parent[iInd][1];

	vector<unsigned int> otherChild;
	for(size_t i=0;i<child[p1].size();i++)
	{
		for(size_t j=0;j<child[p2].size();j++)
		{
			if(child[p1][i]==child[p2][j] && child[p1][i]!=iInd)
			{
				otherChild.push_back(child[p1][i]);
				break;
			}
		}
	}

	vector<unsigned int> otherSpouse1;
	for(size_t i=0;i<spouse[p1].size();i++)
	{
		if(spouse[p1][i]!=p2)
		{
			otherSpouse1.push_back(spouse[p1][i]);
		}
	}

	vector<unsigned int> otherSpouse2;
	for(size_t i=0;i<spouse[p2].size();i++)
	{
		if(spouse[p2][i]!=p1)
		{
			otherSpouse2.push_back(spouse[p2][i]);
		}
	}

	double sm=0;
	for(unsigned int i=0;i<m_nGenotype;i++)
	{
		double sf=0;
		for(unsigned int j=0;j<m_nGenotype;j++)
		{
			double mc=1;
			for(size_t k=0;k<otherChild.size();k++)
			{
				double sc=0;
				for(unsigned int l=0;l<m_nGenotype;l++)
				{
					double mcs=1;
					for(size_t m=0;m<spouse[otherChild[k]].size();m++)
					{
						mcs=mcs*calPosProb(otherChild[k],spouse[otherChild[k]][m],l,antProb,posProb,checkLoop);
					}
					if(m_chromSex)
					{
						if(m_mem[otherChild[k]].get_gender()==1)
						{
							sc=sc+mcs*m_lk(otherChild[k],l)*m_pccpMale(l,i,j);
						}
						else
						{
							sc=sc+mcs*m_lk(otherChild[k],l)*m_pccpFemale(l,i,j);
						}
					}
					else
					{
						sc=sc+mcs*m_lk(otherChild[k],l)*m_pccp(l,i,j);

					}
				}
				mc=mc*sc;
			}

			double mf=1;
			for(size_t k=0;k<otherSpouse2.size();k++)
			{
				mf=mf*calPosProb(p2,otherSpouse2[k],j,antProb,posProb,checkLoop);
			}

			if(m_chromSex)
			{
				if(m_mem[iInd].get_gender()==1)
				{
					sf=sf+calAntProb(p2,j,antProb,posProb,checkLoop)*m_lk(p2,j)*mf*m_pccpMale(iGeno,i,j)*mc;
				}
				else
				{
					sf=sf+calAntProb(p2,j,antProb,posProb,checkLoop)*m_lk(p2,j)*mf*m_pccpFemale(iGeno,i,j)*mc;
				}
			}
			else
			{
				sf=sf+calAntProb(p2,j,antProb,posProb,checkLoop)*m_lk(p2,j)*mf*m_pccp(iGeno,i,j)*mc;
			}
		}
		double mm=1;
		for(size_t j=0;j<otherSpouse1.size();j++)
		{
			mm=mm*calPosProb(p1,otherSpouse1[j],i,antProb,posProb, checkLoop);
		}
		sm=sm+calAntProb(p1,i,antProb,posProb,checkLoop)*m_lk(p1,i)*mm*sf;
	}

	antProb(iInd,iGeno)=sm;
	return sm;
}

double Mendel::calPosProb(unsigned int iInd, unsigned int jInd, unsigned int iGeno, dMatrix<double> & antProb, vector<dMatrix<double> > & posProb, dMatrix<bool> & checkLoop)
{
	if(posProb[iInd](jInd,iGeno)>=0)
	{
		return posProb[iInd](jInd,iGeno);
	}

	//jInd's other spouse beside iInd
	vector<unsigned int> otherSpouse;
	for(size_t i=0;i<spouse[jInd].size();i++)
	{
		if(spouse[jInd][i]!=iInd)
		{
			otherSpouse.push_back(spouse[jInd][i]);
		}
	}

	//iInd and jInd's children
	vector<unsigned int > allChild;
	for(size_t i=0;i<child[iInd].size();i++)
	{
		for(size_t j=0;j<child[jInd].size();j++)
		{
			if(child[iInd][i]==child[jInd][j])
			{
				allChild.push_back(child[iInd][i]);
				break;
			}
		}
	}

	double sj=0;
	for(unsigned int i=0;i<m_nGenotype;i++)
	{
		double ms=1;
		for(size_t j=0;j<otherSpouse.size();j++)
		{
			ms=ms*calPosProb(jInd,otherSpouse[j],i,antProb,posProb,checkLoop);
		}
		double mc=1;
		for(size_t j=0;j<allChild.size();j++)
		{
			double sc=0;
			for(unsigned int k=0;k<m_nGenotype;k++)
			{
				double msc=1;
				for(size_t l=0;l<spouse[allChild[j]].size();l++)
				{
					msc=msc*calPosProb(allChild[j],spouse[allChild[j]][l],k,antProb,posProb,checkLoop);
				}
				if(m_chromSex)
				{
					if(m_mem[allChild[j]].get_gender()==1)
					{
						if(m_mem[iInd].get_gender()==1)
						{
							sc=sc+m_pccpMale(k,iGeno,i)*m_lk(allChild[j],k)*msc;
						}
						else
						{
							sc=sc+m_pccpMale(k,i,iGeno)*m_lk(allChild[j],k)*msc;
						}
					}
					else
					{
						if(m_mem[iInd].get_gender()==1)
						{
							sc=sc+m_pccpFemale(k,iGeno,i)*m_lk(allChild[j],k)*msc;
						}
						else
						{
							sc=sc+m_pccpFemale(k,i,iGeno)*m_lk(allChild[j],k)*msc;
						}
					}
				}
				else
				{
					sc=sc+m_pccp(k,i,iGeno)*m_lk(allChild[j],k)*msc;
				}
			}
			mc=mc*sc;
		}
		sj=sj+calAntProb(jInd,i,antProb,posProb,checkLoop)*m_lk(jInd,i)*ms*mc;
	}

	posProb[iInd](jInd,iGeno)=sj;
	return sj;
}
