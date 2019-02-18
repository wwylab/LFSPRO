/*
 * PCCP.cpp
 *
 *  Created on: Oct 16, 2012
 *      Author: gpeng1
 */

#include "PCCP.h"


using namespace std;

PCCP::PCCP()
{
	m_store=false;
	m_error=false;
	m_nGene=0;
	m_numAllele.clear();
	m_nGeno.clear();
	m_nGenotype=0;
	m_mut=false;
	m_mRate.clear();
	m_chromType.clear();
}

PCCP::PCCP(const vector<unsigned int> & numAllele, double mRate, bool store)
{
	m_store=store;
	m_error=false;
	m_nGene=numAllele.size();
	m_numAllele=numAllele;
	m_nGeno.resize(m_nGene);
	m_nGenotype=1;
	for(size_t i=0;i<numAllele.size();i++)
	{
		m_nGeno[i]=(numAllele[i]+1)*numAllele[i]/2;
		m_nGenotype=m_nGenotype*m_nGeno[i];
	}

	if(mRate==0)
	{
		m_mut=false;
	}
	else
	{
		m_mut=true;
	}

	m_mRate.clear();
	for(unsigned int i=0;i<m_nGene;i++)
	{
		m_mRate.push_back(mRate);
	}

	m_chromType=vector<int> (m_nGene,0);

	init();
}

PCCP::PCCP(const vector<unsigned int> & numAllele, const vector<double> & mRate, bool store)
{
	m_store=store;
	m_error=false;
	m_nGene=numAllele.size();
	m_numAllele=numAllele;
	m_nGeno.resize(m_nGene);
	m_nGenotype=1;
	for(size_t i=0;i<numAllele.size();i++)
	{
		m_nGeno[i]=(numAllele[i]+1)*numAllele[i]/2;
		m_nGenotype=m_nGenotype*m_nGeno[i];
	}

	m_mRate=mRate;
	m_mut=false;
	for(unsigned int i=0;i<m_nGene;i++)
	{
		if(m_mRate[i]!=0)
		{
			m_mut=true;
			break;
		}
	}

	m_chromType=vector<int> (m_nGene,0);

	init();
}

PCCP::PCCP(const vector<unsigned int> & numAllele, double mRate, const vector<int> & chromType, bool store)
{
	m_store=store;
	m_error=false;
	m_nGene=numAllele.size();
	m_numAllele=numAllele;
	m_nGeno.resize(m_nGene);
	m_nGenotype=1;
	for(size_t i=0;i<numAllele.size();i++)
	{
		m_nGeno[i]=(numAllele[i]+1)*numAllele[i]/2;
		m_nGenotype=m_nGenotype*m_nGeno[i];
	}

	if(mRate==0)
	{
		m_mut=false;
	}
	else
	{
		m_mut=true;
	}

	m_mRate.clear();
	for(unsigned int i=0;i<m_nGene;i++)
	{
		m_mRate.push_back(mRate);
	}

	m_chromType=chromType;

	init();
}

PCCP::PCCP(const vector<unsigned int> & numAllele, const vector<double> & mRate, const vector<int> & chromType, bool store)
{
	m_store=store;
	m_error=false;
	m_nGene=numAllele.size();
	m_numAllele=numAllele;
	m_nGeno.resize(m_nGene);
	m_nGenotype=1;
	for(size_t i=0;i<numAllele.size();i++)
	{
		m_nGeno[i]=(numAllele[i]+1)*numAllele[i]/2;
		m_nGenotype=m_nGenotype*m_nGeno[i];
	}

	m_mRate=mRate;
	m_mut=false;
	for(unsigned int i=0;i<m_nGene;i++)
	{
		if(m_mRate[i]!=0)
		{
			m_mut=true;
			break;
		}
	}

	m_chromType=chromType;

	init();
}


unsigned int PCCP::encode(const vector<unsigned int> & geneCode)
{
	if(geneCode.size()!=m_nGene)
	{
		return false;
	}
	unsigned int rlt=0;
	unsigned int base=1;
	for(unsigned int i=0; i< m_nGene;i++)
	{
		rlt=rlt+base*geneCode[i];
		base=base*m_nGeno[i];
	}
	return rlt;
}

bool PCCP::decode(unsigned int geno, vector<unsigned int > & geneCode)
{
	geneCode.resize(m_nGene,0);
	for(size_t i=0;i<m_nGene;i++)
	{
		unsigned int rm=geno%m_nGeno[i];
		geneCode[i]=rm;
		geno=(geno-rm)/m_nGeno[i];
	}
	return true;
}

bool PCCP::alleleCodeTable()
{
	a2gTable.clear();
	g2aTable.clear();
	for(size_t i=0;i<m_nGene;i++)
	{
		dMatrix<unsigned int> a2gTableTmp(m_numAllele[i],m_numAllele[i]);
		vector<allelePair> g2aTableTmp;

		unsigned int count=0;
		for(unsigned int j=0;j<m_numAllele[i];j++)
		{
			for(unsigned int k=j;k<m_numAllele[i];k++)
			{
				a2gTableTmp(j,k)=count;
				a2gTableTmp(k,j)=count;
				g2aTableTmp.push_back(allelePair(j,k));
				count++;
			}
		}
		a2gTable.push_back(a2gTableTmp);
		g2aTable.push_back(g2aTableTmp);
	}
	return true;
}

bool PCCP::calGenoDecode()
{
	genoDecodeTable=dMatrix<unsigned int> (m_nGenotype,m_nGene);
	for(unsigned int i=0;i<m_nGenotype;i++)
	{
		vector<unsigned int> genoCodeTmp;
		if(decode(i,genoCodeTmp))
		{
			for(unsigned int j=0;j<m_nGene;j++)
			{
				genoDecodeTable(i,j)=genoCodeTmp[j];
			}
		}
		else
		{
			return false;
		}
	}
	return true;
}

bool PCCP::calPCCPSAll()
{
	PCCPS.clear();
	for(unsigned int i=0;i<m_nGene;i++)
	{
		vector<dMatrix<double> > PCCPSTmp;
		if(calPCCPS(i,PCCPSTmp))
		{
			PCCPS.push_back(PCCPSTmp);
		}
		else
		{
			return false;
		}
	}
	return true;
}

bool PCCP::calPCCPS(int index, vector<dMatrix<double> > & PCCPSTmp)
{
	for(unsigned int i=0;i<m_nGeno[index];i++)
	{
		dMatrix<double> dmTmp(m_nGeno[index],m_nGeno[index],0);
		PCCPSTmp.push_back(dmTmp);
	}

	if(m_mRate[index]==0)
	{
		if(m_chromType[index]==0)
		{
			//haplotype
			for(unsigned int i=0;i<m_nGeno[index];i++)
			{
				for(unsigned int j=0;j<m_nGeno[index];j++)
				{
					unsigned int i0,i1,j0,j1;
					i0=g2aTable[index][i].first;
					i1=g2aTable[index][i].second;
					j0=g2aTable[index][j].first;
					j1=g2aTable[index][j].second;
					PCCPSTmp[a2gTable[index](i0,j0)](i,j)=PCCPSTmp[a2gTable[index](i0,j0)](i,j)+0.25;
					PCCPSTmp[a2gTable[index](i0,j1)](i,j)=PCCPSTmp[a2gTable[index](i0,j1)](i,j)+0.25;
					PCCPSTmp[a2gTable[index](i1,j0)](i,j)=PCCPSTmp[a2gTable[index](i1,j0)](i,j)+0.25;
					PCCPSTmp[a2gTable[index](i1,j1)](i,j)=PCCPSTmp[a2gTable[index](i1,j1)](i,j)+0.25;
				}
			}
		}
		else if(m_chromType[index]==1)
		{
			for(unsigned int i=0;i<m_numAllele[index];i++)
			{
				for(unsigned int j=0;j<m_nGeno[index];j++)
				{
					unsigned int j0,j1;
					//there are only one allele for male
					j0=g2aTable[index][j].first;
					j1=g2aTable[index][j].second;
					PCCPSTmp[a2gTable[index](j0,j0)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](j0,j0)](a2gTable[index](i,i),j)+0.5;
					PCCPSTmp[a2gTable[index](j1,j1)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](j1,j1)](a2gTable[index](i,i),j)+0.5;
				}
			}
		}
		else if(m_chromType[index]==2)
		{
			for(unsigned int i=0;i<m_numAllele[index];i++)
			{
				for(unsigned int j=0;j<m_nGeno[index];j++)
				{
					unsigned int j0,j1;
					j0=g2aTable[index][j].first;
					j1=g2aTable[index][j].second;
					PCCPSTmp[a2gTable[index](i,j0)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](i,j0)](a2gTable[index](i,i),j)+0.5;
					PCCPSTmp[a2gTable[index](i,j1)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](i,j1)](a2gTable[index](i,i),j)+0.5;
				}
			}
		}
		else if(m_chromType[index]==3)
		{
			for(unsigned int i=0;i<m_numAllele[index];i++)
			{
				for(unsigned int j=0;j<m_nGeno[index];j++)
				{
					PCCPSTmp[a2gTable[index](i,i)](a2gTable[index](i,i),j)=1;
				}
			}
		}
	}
	else
	{
		if(m_chromType[index]==0)
		{
			for(unsigned int i=0;i<m_nGeno[index];i++)
			{
				for(unsigned int j=0;j<m_nGeno[index];j++)
				{
					unsigned int i0,i1,j0,j1;
					i0=g2aTable[index][i].first;
					i1=g2aTable[index][i].second;
					j0=g2aTable[index][j].first;
					j1=g2aTable[index][j].second;

					vector<double> pi(m_numAllele[index],m_mRate[index]/(2*(m_numAllele[index]-1)));
					vector<double> pj(m_numAllele[index],m_mRate[index]/(2*(m_numAllele[index]-1)));

					//i0 and j0
					pi[i0]=(1-m_mRate[index])/2;
					pj[j0]=(1-m_mRate[index])/2;
					for(unsigned int k=0;k<m_numAllele[index];k++)
					{
						for(unsigned int l=0;l<m_numAllele[index];l++)
						{
							PCCPSTmp[a2gTable[index](k,l)](i,j)=PCCPSTmp[a2gTable[index](k,l)](i,j)+pi[k]*pj[l];
						}
					}

					//i0 and j1
					pj[j0]=m_mRate[index]/(2*(m_numAllele[index]-1));
					pj[j1]=(1-m_mRate[index])/2;
					for(unsigned int k=0;k<m_numAllele[index];k++)
					{
						for(unsigned int l=0;l<m_numAllele[index];l++)
						{
							PCCPSTmp[a2gTable[index](k,l)](i,j)=PCCPSTmp[a2gTable[index](k,l)](i,j)+pi[k]*pj[l];
						}
					}

					//i1 and j0
					pi[i0]=m_mRate[index]/(2*(m_numAllele[index]-1));
					pi[i1]=(1-m_mRate[index])/2;
					pj[j1]=m_mRate[index]/(2*(m_numAllele[index]-1));
					pj[j0]=(1-m_mRate[index])/2;
					for(unsigned int k=0;k<m_numAllele[index];k++)
					{
						for(unsigned int l=0;l<m_numAllele[index];l++)
						{
							PCCPSTmp[a2gTable[index](k,l)](i,j)=PCCPSTmp[a2gTable[index](k,l)](i,j)+pi[k]*pj[l];
						}
					}

					//i1 and j1
					pj[j0]=m_mRate[index]/(2*(m_numAllele[index]-1));
					pj[j1]=(1-m_mRate[index])/2;
					for(unsigned int k=0;k<m_numAllele[index];k++)
					{
						for(unsigned int l=0;l<m_numAllele[index];l++)
						{
							PCCPSTmp[a2gTable[index](k,l)](i,j)=PCCPSTmp[a2gTable[index](k,l)](i,j)+pi[k]*pj[l];
						}
					}
				}
			}
		}
		else if(m_chromType[index]==1)
		{
			for(unsigned int i=0;i<m_numAllele[index];i++)
			{
				for(unsigned int j=0;j<m_nGeno[index];j++)
				{
					unsigned int j0,j1;
					//there are only one allele for male from mother
					j0=g2aTable[index][j].first;
					j1=g2aTable[index][j].second;

					//j0
					for(unsigned int k=0;k<m_numAllele[index];k++)
					{
						if(k==j0)
						{
							PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)+(1-m_mRate[index])/2;
						}
						else
						{
							PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)+m_mRate[index]/(2*(m_numAllele[index]-1));;
						}
					}

					//j1
					for(unsigned int k=0;k<m_numAllele[index];k++)
					{
						if(k==j1)
						{
							PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)+(1-m_mRate[index])/2;
						}
						else
						{
							PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)+m_mRate[index]/(2*(m_numAllele[index]-1));;
						}
					}
				}
			}
		}
		else if(m_chromType[index]==2)
		{
			for(unsigned int i=0;i<m_numAllele[index];i++)
			{
				for(unsigned int j=0;j<m_nGeno[index];j++)
				{
					unsigned int j0,j1;
					j0=g2aTable[index][j].first;
					j1=g2aTable[index][j].second;

					vector<double> pi(m_numAllele[index],m_mRate[index]/(m_numAllele[index]-1));
					vector<double> pj(m_numAllele[index],m_mRate[index]/(2*(m_numAllele[index]-1)));

					//i and j0
					pi[i]=1-m_mRate[index];
					pj[j0]=(1-m_mRate[index])/2;
					for(unsigned int k=0;k<m_numAllele[index];k++)
					{
						for(unsigned int l=0;l<m_numAllele[index];l++)
						{
							PCCPSTmp[a2gTable[index](k,l)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](k,l)](a2gTable[index](i,i),j)+pi[k]*pj[l];
						}
					}

					//i and j1
					pj[j0]=m_mRate[index]/(2*(m_numAllele[index]-1));
					pj[j1]=(1-m_mRate[index])/2;
					for(unsigned int k=0;k<m_numAllele[index];k++)
					{
						for(unsigned int l=0;l<m_numAllele[index];l++)
						{
							PCCPSTmp[a2gTable[index](k,l)](a2gTable[index](i,i),j)=PCCPSTmp[a2gTable[index](k,l)](a2gTable[index](i,i),j)+pi[k]*pj[l];
						}
					}
				}
			}
		}
		else if(m_chromType[index]==3)
		{
			for(unsigned int i=0;i<m_numAllele[index];i++)
			{
				for(unsigned int j=0;j<m_nGeno[index];j++)
				{
					for(unsigned int k=0;k<m_numAllele[index];i++)
					{
						if(k==i)
						{
							PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)=1-m_mRate[index];
						}
						else
						{
							PCCPSTmp[a2gTable[index](k,k)](a2gTable[index](i,i),j)=m_mRate[index]/(m_numAllele[index]-1);
						}
					}
				}
			}
		}
	}
	return true;
}

bool PCCP::init()
{
	//preparation
	if(!alleleCodeTable())
	{
		m_error=true;
	}
	if(!calGenoDecode())
	{
		m_error=true;
	}
	if(!calPCCPSAll())
	{
		m_error=true;
	}

//	cout<<genoDecodeTable<<endl;
//	cout<<PCCPS[0][0]<<endl;
//	cout<<PCCPS[0][1]<<endl;
//	cout<<PCCPS[0][2]<<endl;
//
//	cout<<PCCPS[1][0]<<endl;
//	cout<<PCCPS[1][1]<<endl;
//	cout<<PCCPS[1][2]<<endl;

//	cout<<a2gTable[0]<<endl;
//	cout<<a2gTable[1]<<endl;

	if(m_store)
	{
		if(m_mut)
		{
			fullPCCP.clear();
			for(unsigned int i=0;i<m_nGenotype;i++)
			{
				dMatrix<double> dmTmp(m_nGenotype,m_nGenotype,0);
				for(unsigned int j=0;j<m_nGenotype;j++)
				{
					for(unsigned int k=0;k<m_nGenotype;k++)
					{
						double probTmp=1;
						for(unsigned int l=0;l<m_nGene;l++)
						{
							probTmp=probTmp*PCCPS[l][genoDecodeTable(i,l)](genoDecodeTable(j,l),genoDecodeTable(k,l));
						}
						dmTmp(j,k)=probTmp;
					}
				}
				fullPCCP.push_back(dmTmp);
			}
		}
		else
		{
			spPCCP.clear();
			for(unsigned int i=0;i<m_nGenotype;i++)
			{
				vector<map<int, double> > vmTmp;
				for(unsigned int j=0;j<m_nGenotype;j++)
				{
					map<int,double> mTmp;
					for(unsigned int k=0;k<m_nGenotype;k++)
					{
						double probTmp=1;
						for(unsigned int l=0;l<m_nGene;l++)
						{
							probTmp=probTmp*PCCPS[l][genoDecodeTable(k,l)](genoDecodeTable(i,l),genoDecodeTable(j,l));
							if(probTmp==0)
							{
								break;
							}
						}
						if(probTmp!=0)
						{
							mTmp[k]=probTmp;
						}
					}
					vmTmp.push_back(mTmp);
				}
				spPCCP.push_back(vmTmp);
			}
		}
		g2aTable.clear();
		a2gTable.clear();
		PCCPS.clear();
		genoDecodeTable=dMatrix<unsigned int> (1,1,0);
	}

	return true;
}

double PCCP::operator ()(unsigned int m, unsigned int n, unsigned int k) const
{
	if(m_error)
	{
		cout<<"Error in parameters!"<<endl;
		return -1;
	}
	if(m_store)
	{
		if(m_mut)
		{
			return fullPCCP[m](n,k);
		}
		else
		{
			map<int,double>::const_iterator it=spPCCP[n][k].find(m);
			if(it==spPCCP[n][k].end())
			{
				return 0;
			}
			else
			{
				return it->second;
			}
		}
	}
	else
	{
		double rlt=1;
		for(unsigned int i=0;i<m_nGene;i++)
		{
			rlt=rlt*PCCPS[i][genoDecodeTable(m,i)](genoDecodeTable(n,i),genoDecodeTable(k,i));
		}
		return rlt;
	}
}

PCCP & PCCP::operator =(const PCCP & pccp)
{
	if(this==&pccp)
	{
		return *this;
	}

	m_nGene=pccp.m_nGene;
	m_nGenotype=pccp.m_nGenotype;
	m_store=pccp.m_store;
	m_mut=pccp.m_mut;
	m_mRate=pccp.m_mRate;
	m_chromType=pccp.m_chromType;
	m_numAllele=pccp.m_numAllele;
	m_nGeno=pccp.m_nGeno;
	spPCCP=pccp.spPCCP;
	fullPCCP=pccp.fullPCCP;
	g2aTable=pccp.g2aTable;
	a2gTable=pccp.a2gTable;
	genoDecodeTable=pccp.genoDecodeTable;
	PCCPS=pccp.PCCPS;

	return *this;
}
