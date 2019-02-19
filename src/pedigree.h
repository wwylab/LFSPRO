/*
 * pedigree.h
 *
 *  Created on: Nov 6, 2012
 *      Author: gpeng1
 */

#ifndef PEDIGREE_H_
#define PEDIGREE_H_

#include <vector>
#include <iostream>
#include <fstream>

class individual
{
private:
	// sample id
	int m_id;
	// gender
	// 1: male
	// 2: female
	// 0: unknown
	int m_gender;
	// father id
	int m_fid;
	// mother id
	int m_mid;

public:
	// constructor
	individual(int id, int fid, int mid, int gender);

	///////////////////
	// get functions //
	///////////////////
	int get_id() const {return m_id;}
	int get_gender() const {return m_gender;}
	int get_fid() const {return m_fid;}
	int get_mid() const {return m_mid;}

	///////////////////
	// set functions //
	///////////////////
	bool set_id(int id);
	bool set_gender(int gender);
	bool set_fid(int fid);
	bool set_mid(int mid);
};

class pedigree
{
protected:
	// sample id cannot be negative
	// 0: missing
	// members in a pedigree
	std::vector<individual> m_mem;

	// number of samples in original pedigree
	unsigned int m_nInd;

	// number of samples in after fulfill the parents
	unsigned int m_nIndAll;

	//one sample's child (only store the index)
	std::vector<std::vector<std::size_t> > child;

	//one sample's parent (only store the index)
	//first: father
	//second: mother
	std::vector<std::vector<std::size_t> > parent;

	//one sample's spouse
	std::vector<std::vector<std::size_t> > spouse;

private:
	//fulfill the pedigree
	bool fulfillPed();

	//set relationship in the pedigree
	bool setRelation();

	//check pedigree
	bool checkPed();

public:
	pedigree(const std::vector<individual> & mem);

	///////////////////
	// get functions //
	///////////////////
	unsigned int get_nInd() const { return m_nInd; }

	unsigned int get_nIndAll() const { return m_nIndAll; }

	std::vector<std::vector<std::size_t> > get_child() const { return child; }

	std::vector<std::size_t> get_child(unsigned int i) const {return child[i]; }

	std::vector<std::vector<std::size_t> > get_parent() const { return parent; }

	std::vector<std::size_t> get_parent(unsigned int i) const { return parent[i]; }

	std::vector<std::vector<std::size_t> > get_spouse() const { return spouse; }

	std::vector<std::size_t> get_spouse(unsigned int i) const { return spouse[i]; }

	individual operator[] (unsigned int i) const { return m_mem[i]; }


	friend std::ostream & operator<< (std::ostream & os, const pedigree & ped);

	friend std::ofstream & operator<< (std::ofstream & of, const pedigree & ped);
};

#endif /* PEDIGREE_H_ */
