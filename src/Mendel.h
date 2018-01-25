/*
 * Mendel.h
 *
 *  Created on: Nov 6, 2012
 *      Author: gpeng1
 */

#ifndef MENDEL_H_
#define MENDEL_H_

#include "PCCP.h"
#include "pedigree.h"

enum calMethod {MCMC, BN, Peeling};

class Mendel: public pedigree
{
private:
	PCCP m_pccp;

	PCCP m_pccpMale;

	PCCP m_pccpFemale;

	//probability for each genotype
	//Pr(genotype)
	std::vector<double> m_genoProb;

	// likelihood
	dMatrix<double> m_lk;

	unsigned int m_nGenotype;

	// whether there is an error
	bool m_error;

	bool m_chromSex;

	bool m_loop;

	/////////////////////
	// check functions //
	/////////////////////

	// check whether the number of genotype in pccp matches the number of genotype in genoProb
	bool check_nGeno();

	/////////////////////////
	// calculate functions //
	/////////////////////////

	dMatrix<double> calPostProbPeeling();

	dMatrix<double> calPostProbPeeling(const std::vector<unsigned int> & sampleIndex);

	dMatrix<double> calPostProbBN();

	dMatrix<double> calPostProbMCMC();

	double calAntProb(unsigned int iInd, unsigned int iGeno, dMatrix<double> & antProb, std::vector<dMatrix<double> > & posProb, dMatrix<bool> & checkLoop);

	double calPosProb(unsigned int iInd, unsigned int jInd, unsigned int iGeno, dMatrix<double> & antProb, std::vector<dMatrix<double> > & posProb, dMatrix<bool> & checkLoop);

public:
	// constructor
	Mendel(const std::vector<individual> & mem, const std::vector<double> & genoProb, const std::vector<unsigned int> & numAllele, double mRate, bool store=true);
	Mendel(const std::vector<individual> & mem, const std::vector<double> & genoProb, const std::vector<unsigned int> & numAllele, const std::vector<double> & mRate, bool stroe=true);
	// for chromosome type
	// 0: autosome
	// 1: chromosome X
	Mendel(const std::vector<individual> & mem, const std::vector<double> & genoProb, const std::vector<unsigned int> & numAllele, double mRate, const std::vector<int> & chromType, bool store=true);
	Mendel(const std::vector<individual> & mem, const std::vector<double> & genoProb, const std::vector<unsigned int> & numAllele, const std::vector<double> & mRate, const std::vector<int> & chromType, bool store=true);

	///////////////////
	// get functions //
	///////////////////
	PCCP get_pccp() const { return m_pccp; }


	std::vector<double> get_genoProb() const { return m_genoProb; }

	unsigned int get_nGenotype() const{ return m_nGenotype;}

	dMatrix<double> get_lk() const { return m_lk; }

	///////////////////
	// set functions //
	///////////////////

	// set pccp
	// when pccp is not set properly, return false
	bool set_pccp(const PCCP & pccp);

	// set genoProb
	// when genoProb is not set properly, return false
	bool set_genoProb(const std::vector<double> & genoProb);

	//set pccp and genoProb
	bool set_par(const PCCP & pccp, const std::vector<double> & genoProb);


	//calculate the posterior probability
	dMatrix<double> calPostProb(const dMatrix<double> & lk, calMethod method=Peeling);

	dMatrix<double> calPostProb(const dMatrix<double> & lk, const std::vector<int> & sampleId, calMethod method=Peeling);
};


#endif /* MENDEL_H_ */
