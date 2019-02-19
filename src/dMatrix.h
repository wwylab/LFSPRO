#ifndef DMATRIX_H_INCLUDED
#define DMATRIX_H_INCLUDED

//Created on March 5 2011
//Last updated in December 16, 2011
//Update:
//June 14
//add two friend function to output dMatrix (overload <<)
//December 16, 2011
//add transpose function tr()
//Gang Peng
//gpeng1@mdanderson.org
#include <iostream>
#include <fstream>

template <class Type>
class dMatrix
{
    private:
    //data
    Type* data;
    //row
    unsigned int m_row;
    //column
    unsigned int m_column;

    public:
    dMatrix();
    dMatrix(unsigned int row, unsigned int column, Type m=0);
    dMatrix(const dMatrix<Type> & dM);
    virtual ~dMatrix();

    unsigned int get_row() const { return m_row;}
    unsigned int get_column() const { return m_column;}


    bool setZero();

    dMatrix<Type> & operator=(const dMatrix<Type> & dM);
    const Type & operator()(unsigned int m, unsigned int n) const;
    Type & operator()(unsigned int m, unsigned int n);

    dMatrix<Type> tr() const;

    dMatrix<Type> operator * (const dMatrix<Type> & ma) const;

    //sumation of row m
    Type sum_row(unsigned int m);
    //sumation of column n
    Type sum_column(unsigned int n);
    //sumation of whole matrix
    Type sumAll();
    //output
    bool output() const;

    //output
    template <class T>
    friend std::ostream & operator<< (std::ostream & os, const dMatrix<T> & dM);

    template <class T>
    friend std::ofstream & operator<<(std::ofstream & of, const dMatrix<T> & dM);
};

template <class Type>
dMatrix<Type>::dMatrix()
{
    m_row=1;
    m_column=1;
    data=new Type[1];
    data[0]=0;
}

template <class Type>
dMatrix<Type>::dMatrix(unsigned int row, unsigned int column, Type m)
{
    m_row=row;
    m_column=column;
    data=new Type[row*column];
    for(unsigned int i=0;i<row;i++)
    {
        for(unsigned int j=0;j<column;j++)
        {
            (*this)(i,j)=m;
        }
    }
}

template <class Type>
dMatrix<Type>::dMatrix(const dMatrix<Type> & dM)
{
    m_row=dM.get_row();
    m_column=dM.get_column();
    data=new Type[m_row*m_column];
    for(unsigned int i=0;i<m_row;i++)
    {
        for(unsigned int j=0;j<m_column;j++)
        {
            (*this)(i,j)=dM(i,j);
        }
    }
}

template <class Type>
dMatrix<Type>::~dMatrix()
{
    delete[] data;
}

template <class Type>
dMatrix<Type> & dMatrix<Type>::operator=(const dMatrix<Type> & dM)
{
    if(this==&dM)
    {
        return *this;
    }

    delete[] data;
    m_row=dM.get_row();
    m_column=dM.get_column();
    data=new Type[m_row*m_column];
    for(unsigned int i=0;i<m_row;i++)
    {
        for(unsigned int j=0;j<m_column;j++)
        {
            (*this)(i,j)=dM(i,j);
        }
    }
    return *this;
}

template <class Type>
const Type & dMatrix<Type>::operator()(unsigned int m, unsigned int n) const
{
    return data[m*m_column+n];
}

template <class Type>
Type & dMatrix<Type>::operator()(unsigned int m, unsigned int n)
{
    return data[m*m_column+n];
}

template <class Type>
Type dMatrix<Type>::sum_row(unsigned int m)
{
    Type rlt=0;
    for(unsigned int i=0;i<m_column;i++)
    {
        rlt=rlt+(*this)(m,i);
    }
    return rlt;
}

template <class Type>
Type dMatrix<Type>::sum_column(unsigned int n)
{
    Type rlt=0;
    for(int i=0;i<m_row;i++)
    {
        rlt=rlt+(*this)(i,n);
    }
    return rlt;
}

template <class Type>
Type dMatrix<Type>::sumAll()
{
    Type rlt=0;
    for(int i=0;i<m_row;i++)
    {
        for(int j=0;j<m_column;j++)
        {
            rlt=rlt+(*this)(i,j);
        }
    }
    return rlt;
}

template <class Type>
bool dMatrix<Type>::output() const
{
    using namespace std;
    for(int i=0;i<m_row;i++)
    {
        for(int j=0;j<m_column;j++)
        {
            cout<<(*this)(i,j)<<'\t';
        }
        cout<<endl;
    }
    return true;
}

template <class T>
std::ostream & operator<<(std::ostream & os, const dMatrix<T> & dM)
{
	for(unsigned int i=0;i<dM.get_row();i++)
	{
		for(unsigned int j=0;j<dM.get_column();j++)
		{
			os<<dM(i,j)<<'\t';
		}
		os<<std::endl;
	}
	return os;
}

template <class T>
std::ofstream & operator<<(std::ofstream & of, const dMatrix<T> & dM)
{
	for(unsigned int i=0;i<dM.get_row();i++)
	{
		for(unsigned int j=0;j<dM.get_column();j++)
		{
			of<<dM(i,j)<<'\t';
		}
		of<<std::endl;
	}
	return of;
}

template <class Type>
bool dMatrix<Type>::setZero()
{
	for(int i=0;i<m_row;i++)
	{
		for(int j=0;j<m_column;j++)
		{
			(*this)(i,j)=0;
		}
	}

	return true;
}

template <class Type>
dMatrix<Type> dMatrix<Type>::tr() const
{
	dMatrix<Type> rlt(m_column,m_row);
	for(int i=0;i<m_row;i++)
	{
		for(int j=0;j<m_column;j++)
		{
			rlt(j,i)=(*this)(i,j);
		}
	}
	return rlt;
}

template <class Type>
dMatrix<Type> dMatrix<Type>::operator *(const dMatrix<Type> & ma) const
{
	dMatrix<Type> rlt(m_row,ma.get_column(),0);
	for(int i=0;i<m_row;i++)
	{
		for(int j=0;j<ma.get_column();j++)
		{
			for(int k=0;k<m_column;k++)
			{
				rlt(i,j)+=(*this)(i,k)*ma(k,j);
			}
		}
	}
	return rlt;
}
#endif // DMATRIX_H_INCLUDED
