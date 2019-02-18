/*
 * PCCP.h
 *
 * Probability of Child's genotype Conditional on Parents' genotype
 * Pr(Gc|Gf,Gm)
 *
 *  Created on: Oct 16, 2012
 *      Author: gpeng1
 */

#ifndef PCCP_H_
#define PCCP_H_

#include <vector>
#include <map>

#include "dMatrix.h"

typedef std::pair<int, double> idPair;
typedef std::pair<unsigned int, unsigned int> allelePair;

class PCCP
{
private:
	// number of genes
	unsigned int m_nGene;

	// number of all genotype
	unsigned int m_nGenotype;

	// indicate whether we need to store the PCCP
	bool m_store;

	// indicate whether considering mutaiton rate
	bool m_mut;

	// whether there is error
	bool m_error;

	// mutation rate for each gene
	std::vector<double> m_mRate;

	// chromosome type
	// 0: autosome
	// 1: chromosome X, child is male
	// 2: chromosome X, child is female
	// 3: chromosome Y, child is male
	std::vector<int> m_chromType;

	// number of allele for each gene
	std::vector<unsigned int> m_numAllele;

	// number of genotype for each gene
	std::vector<unsigned int> m_nGeno;

	// sparse matrix for PCCP, when mutation rate =0
	std::vector<std::vector <std::map<int, double> > > spPCCP;

	// matrix for PCCP, when mutation rate > 0
	std::vector<dMatrix<double> > fullPCCP;

	// code table for from gene to allele
	std::vector<std::vector<allelePair> > g2aTable;

	// code table for allele to gene
	std::vector<dMatrix<unsigned int> > a2gTable;

	// code table from whole code to code for each gene
	dMatrix<unsigned int> genoDecodeTable;

	// probability of child's genotype conditional on parents' genotype for each gene
	// there is no chromosome Y for female.
	// When calculate Pr(Gc | Gf, Gm) for chromosome Y, Gm could be set to any genotype
	std::vector<std::vector<dMatrix<double> > > PCCPS;


	/////////////////////
	//private functions//
	/////////////////////

	// encode the code for each gene to the whole code
	// input:
	// geneCode: code for each gene
	// output:
	// return: whole code
	unsigned int encode(const std::vector<unsigned int> & geneCode);

	// decode the whole code to the code for each gene
	// input:
	// gneo: whole code
	// output:
	// geneCode: code for each gene
	// return: true
	//
	// Code rule for multiple genes
	// The multiple gene code rule is based on the code of a single gene
	// Example:
	// Gene 1 (two alleles)
	// 0 0: 0; 0 1: 1; 1 1: 2;
	// Gene 2 (three alleles)
	// 0 0: 0; 0 1: 1; 0 2: 2; 1 1: 3; 1 2: 4; 2 2:5
	// The code for the whole genotype
	// 0 0 (the code for first gene is 0 and for second gene is 0): 0
	// 1 0 (the code for first gene is 1 and for second gene is 0): 1
	// 2 0: 2
	// 0 1: 3
	// 1 1: 4
	// 2 1: 5
	// 0 2: 7
	// ...
	// 2 5: 17
	bool decode(unsigned int geno, std::vector<unsigned int> & geneCode);

	// get encode and decode table for allele code to gene code and gene code to allele code
	// if there are n alleles in a gene, it is coded as following:
	// 0 0: 0
	// 0 1: 1
	// 0 2: 2
	// ...
	// 0 n-1: n
	// 1 1: n+1
	// 1 2: n+2
	// ...
	// 1 n-1: 2n-1
	// ...
	// n-1 n-1: n(n+1)/2
	// if it is chromosome X for male, it is coded the same way above.
	// 0 0 means he has allele 0; 1 1 means he has allele 1, etc.
	bool alleleCodeTable();

	// calculate PCCPS table for all genes
	bool calPCCPSAll();

	// calculate PCCPS table for one gene
	bool calPCCPS(int index, std::vector<dMatrix<double> > & PCCPSTmp);

	// calculate genoDecodeTable
	bool calGenoDecode();

	//initiate private members
	bool init();

public:
	// constructor
	PCCP();
	PCCP(const std::vector<unsigned int> & numAllele, double mRate, bool store=true);
	PCCP(const std::vector<unsigned int> & numAllele, const std::vector<double> & mRate, bool store=true);
	PCCP(const std::vector<unsigned int> & numAllele, double mRate, const std::vector<int> & chromType, bool store=true);
	PCCP(const std::vector<unsigned int> & numAllele, const std::vector<double> & mRate, const std::vector<int> & chromType, bool store=true);

	///////////////////
	// get functions //
	///////////////////
	unsigned int get_nGenotype() const {return m_nGenotype;}

	// Pr(Gc=m|Gf=n,Gm=k)
	double operator()(unsigned int m, unsigned int n, unsigned int k) const;

	PCCP & operator=(const PCCP & pccp);
};



#endif /* PCCP_H_ */
