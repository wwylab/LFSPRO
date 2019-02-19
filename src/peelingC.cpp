#include <R.h>

//#include "peelingC.h"

#include <cstdlib>
#include <cmath>

#include "Mendel.h"
#include "normal.h"

using namespace std;

extern "C"{
  void peelingC(double* genoprob, double* likelihood, int* id, int* gender, int* fid, int* mid, int* nSample, int* cid, int* ncid, int* nloci, int* nAllele, double* mrate, double* posteriorprob){
    
	int nGeno = 1;
	for(int i=0;i<*nloci;i++){
		nGeno = nGeno*(nAllele[i]*(nAllele[i]+1)/2);
	}

	/*
	for(int i=0;i<nGeno;i++){
		cout<<genoprob[i]<<'\t';
	}
	cout<<endl;
	*/


    vector<individual> mem;
    for(int i=0;i<*nSample;i++){
      individual indTmp(id[i], fid[i], mid[i],gender[i]);
      mem.push_back(indTmp);
    }
    
    vector<unsigned int> numAllele (nAllele,nAllele+*nloci);
    
    vector<int> numGeno(*nloci);
    for(int i=0;i<*nloci;i++){
    	numGeno[i] = numAllele[i]*(numAllele[i]+1)/2;
    }

    vector<double> mRate(mrate, mrate+*nloci);
    //for(int i=0;i<nloci;i++){
    //  mRate.push_back(mrate[i]);
    //}
    
    //cout<<mRate.size()<<endl;

    //vector<double> genoProb(genoprob,(genoprob+nGeno));
    //for(int i=0;i<3^nloci;i++){
    //  genoProb.push_back(genoprob[i])
    //}

    vector<vector<double> > genoProbSingle;
    int gpIndex = 0;
    for(int i=0;i<*nloci;i++){
    	vector<double> genoProbSingleTmp(numGeno[i]);
    	for(int j=0;j<numGeno[i];j++){
    		genoProbSingleTmp[j] = genoprob[gpIndex];
    		gpIndex++;
    	}
    	genoProbSingle.push_back(genoProbSingleTmp);
    }

    vector<double> genoProb(nGeno,1);
    //vector<int> genoCode(*nloci);
    for(int i=0;i<nGeno;i++){
    	int genoCodeAll = i;
    	for(int j=0;j<*nloci;j++){
    		int rm = genoCodeAll%numGeno[j];
    		//genoCode[j] = rm;
    		genoCodeAll = (genoCodeAll-rm)/numGeno[j];
    		genoProb[i] = genoProb[i]*genoProbSingle[j][rm];
    	}
    	//cout<<genoProb[i]<<endl;
    }

    Mendel mendel(mem,genoProb, numAllele, mRate);



    dMatrix<double> lk(*nSample,nGeno);
    for(int i=0;i<*nSample;i++){
      for(int j=0;j<nGeno;j++){
        lk(i,j)=likelihood[i*nGeno+j];
      }
    }
    
    //cout<<lk<<endl;

    vector<int> sampleID(cid, cid+*ncid);
    dMatrix<double> postProb = mendel.calPostProb(lk,sampleID);

    //cout<<postProb<<endl;

    //cout<<1<<endl;

    if(postProb.get_column()==1){
        for(int i=0;i<*ncid;i++){
          for(int j=0;j<nGeno;j++){
            posteriorprob[i*nGeno+j]=postProb(0,0);
          }
        }
    }
    else{
        for(int i=0;i<*ncid;i++){
          for(int j=0;j<nGeno;j++){
            posteriorprob[i*nGeno+j]=postProb(i,j);
          }
        }
    }
  }
}
