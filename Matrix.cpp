#include "Matrix.hpp"
#include "mt.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>
#include <stdexcept>

void block::Matrix::init(int rows, int cols)
{
    this -> rows = rows;
    this -> cols = cols;
	this -> mat.resize(0);
    this -> mat.resize(rows,std::vector<double>(cols,0.0));
	this -> eigenDone = false;
	this -> eigenValue = 0.0;
	this -> eigenVector.resize(0);
}

int block::Matrix::checkInput(std::string input)
{
	int check = 0;
	size_t idx = 0;
	double val, tol = 1E-12;
	
 	try
	{
		val = std::stod(input,&idx);
	}
	catch(std::invalid_argument)
	{
		check = -1;
	}
	
	if(input.substr(idx,input.npos).size() != 0)
		check = -1;
	
	if(check != -1)
	{
		if(std::abs(val-std::floor(val+0.5)) > tol)
			check = -2;
		else if(val < 0 || std::abs(val) < tol)
			check = -3;
		else
			check = static_cast<int>(std::floor(val+0.5));
	}
	
	return check;
}

int block::Matrix::solicitUserInput()
{
	int check = 0;
	std::string input;
	while(check <= 0)
	{
		getline(std::cin,input);
		check = checkInput(input);
		switch(check)
		{
			case -1:
				std::cout << "You must enter a numeric value, specifically a postive integer." << std::endl;
				break;
			case -2:
				std::cout << "You must enter an integer, specifically a postive integer." << std::endl;
				break;
			case -3:
				std::cout << "You must enter a positive value, specifically a postive integer." << std::endl;
				break;
			default:
				std::cout << "Value <" << check << "> accepted." << std::endl;
				break;
		}
	}
	
	return check;
}

block::Matrix::Matrix()
{
	init(0,0);
}

block::Matrix::Matrix(int N)
{
    if (N <= 0)
    {
        std::cout << "The dimension of the matrix must be positive."
                  << std::endl;
        std::cout << "Enter a new dimension for the matrix to have:"
                  << std::endl;
        N = solicitUserInput();
    }
    
    init(N,N);
}

block::Matrix::Matrix(int N, int M)
{
    if (N <= 0)
    {
        std::cout << "The number of rows of the matrix must be positive."
        << std::endl;
        std::cout << "Enter the number of rows you want the matrix to have:"
        << std::endl;
		N = solicitUserInput();
    }
    
    if (M <= 0)
    {
        std::cout << "The number of columns of the matrix must be positive."
        << std::endl;
        std::cout << "Enter the number of columns you want the matrix to have:"
        << std::endl;
        M = solicitUserInput();
    }
    
    init(N,M);
}

block::Matrix::Matrix(std::string filename)
{
	setFromFile(filename);
}

block::Matrix::~Matrix() {}

const int block::Matrix::getRows() const
{
    return this -> rows;
}

const int block::Matrix::getCols() const
{
    return this -> cols;
}

const double block::Matrix::getEigenValue() const
{
	return this -> eigenValue;
}

const std::vector<double> block::Matrix::getEigenVector() const
{
	return this -> eigenVector;
}

const double block::Matrix::getVal(int n, int m) const
{
	return this -> mat.at(n-1).at(m-1);
}

double& block::Matrix::getVal(int n, int m)
{
   return mat.at(n-1).at(m-1);
}

void block::Matrix::setFromFile(std::string filename)
{
	int row, col;
	double val;
	std::ifstream inFile(filename);
	std::string line;
	char* lineC = new char[1024*1024];
	int numRows = 0;
	int numCols = 0;
	
	if(!inFile.is_open())
	{
		std::cout << "Problem opening inFile file." << std::endl;
		init(0,0);
	}
	else
	{
		getline(inFile,line);
		numRows++;
		strcpy(lineC,line.c_str());
		lineC = strtok(lineC," ");
		while(lineC != NULL)
		{
			lineC = strtok(NULL," ");
			numCols++;
		}
		delete[] lineC;
		while(getline(inFile,line))
			numRows++;
		
		init(numRows,numCols);
		
		inFile.clear();
		inFile.seekg(0);
		
		for (row = 1; row <= rows; row++)
		{
			for (col = 1; col <= cols; col++)
			{
				inFile >> val;
				getVal(row,col) = val;
			}
		}
		
		inFile.close();
	}
}

void block::Matrix::copyVals(const block::Matrix& matToCopy)
{
	for(int i=0; i<rows; i++)
		for(int j=0; j<cols; j++)
			mat.at(i).at(j) = matToCopy.mat.at(i).at(j);
}

std::vector<double>& block::Matrix::operator[](std::size_t idx)
{
	return mat.at(idx-1);
}

const std::vector<double> block::Matrix::operator[](std::size_t idx) const
{
	return mat.at(idx-1);
}

block::Matrix& block::Matrix::operator=(const block::Matrix& matToAssign)
{
	if (this != &matToAssign) { // self-assignment check expected
		if (matToAssign.rows != rows && matToAssign.cols != cols) {         // storage cannot be reused
			init(matToAssign.rows,matToAssign.cols);
			copyVals(matToAssign);
		} else {                          // storage can be reused
			copyVals(matToAssign);
		}
	}
	return *this;
}

block::Matrix& block::Matrix::operator+=(const block::Matrix& rhs)
{
	int row, col;
	int N1 = this -> rows;
	int M1 = this -> cols;
	int N2 = rhs.getRows();
	int M2 = rhs.getCols();
	
	if (N1 == N2 && M1 == M2)
		for (row = 1; row <= N1; row++)
			for (col = 1; col <= M1; col++)
				(*this).getVal(row,col) += rhs.getVal(row,col);
	else
		std::cout << "Matrices must have the same dimension to add them."
		<< std::endl;
	
	return *this;
}

block::Matrix& block::Matrix::operator-=(const block::Matrix& rhs)
{
	int row, col;
	int N1 = this -> rows;
	int M1 = this -> cols;
	int N2 = rhs.getRows();
	int M2 = rhs.getCols();
	
	if (N1 == N2 && M1 == M2)
		for (row = 1; row <= N1; row++)
			for (col = 1; col <= M1; col++)
				(*this).getVal(row,col) -= rhs.getVal(row,col);
	else
		std::cout << "Matrices must have the same dimension to subtract them."
		<< std::endl;
	
	return *this;
}

block::Matrix& block::Matrix::operator*=(const double rhs)
{
	int row, col;
	int N1 = this -> rows;
	int M1 = this -> cols;
	
	for (row = 1; row <= N1; row++)
		for (col = 1; col <= M1; col++)
			(*this).getVal(row,col) *= rhs;
	
	return *this;
}

block::Matrix& block::Matrix::operator*=(const block::Matrix& rhs)
{
	int row, col;
	int N1 = this -> rows;
	int M1 = this -> cols;
	int N2 = rhs.getRows();
	int M2 = rhs.getCols();
	block::Matrix temp;
	
	if (M1 == N2)
	{
		double prodsum;
		temp = *this;
		init(N1,M2);
		
		for (row = 1; row <= N1; row++)
			for (col = 1; col <= M2; col++)
			{
				prodsum = 0.0;
				for (int i = 1; i <= M1; i++)
					prodsum += temp.getVal(row,i) * rhs.getVal(i,col);
				(*this).getVal(row,col) = prodsum;
			}
	}
	else
	{
		std::cout << "The inner dimensions of the two matrices must match to"
		<< " multiply them." << std::endl;
	}
	
	return *this;
}

block::Matrix block::operator+(block::Matrix lhs, const block::Matrix& rhs)
{
	lhs += rhs;
	return lhs;
}

block::Matrix block::operator-(block::Matrix lhs, const block::Matrix& rhs)
{
	lhs -= rhs;
	return lhs;
}

block::Matrix block::operator*(block::Matrix lhs, const double rhs)
{
	lhs *= rhs;
	return lhs;
}

block::Matrix block::operator*(const double lhs, block::Matrix rhs)
{
	rhs *= lhs;
	return rhs;
}

block::Matrix block::operator*(block::Matrix lhs, const block::Matrix& rhs)
{
	lhs *= rhs;
	return lhs;
}

block::Matrix block::Matrix::T() const
{
	int i,j;
	block::Matrix t(cols,rows);
	for(i=1; i<=cols; i++)
		for(j=1; j<=rows; j++)
			t.getVal(i,j) = (*this).getVal(j,i);
	
	return t;
}

double block::Matrix::norm() const
{
	double norm = 0.0;
	
	if(cols == 1)
		for(int i=0; i<rows; i++)
			norm += mat.at(i).at(0) * mat.at(i).at(0);
	else if(rows == 1)
		for(int i=0; i<cols; i++)
			norm += mat.at(0).at(i) * mat.at(0).at(i);
	else
	{
		std::cout << "Cannot find the norm of a non-vector." << std::endl;
		return 0.0;
	}
	
	return std::sqrt(norm);
}

void block::Matrix::findLargestLambda(double shift)
{	
	// need block::Matrix for matrix operation
	// check if Matrix A is square, if not print error message and return
	if(rows != cols){
		return;}
	
	double tol = .0000001;
	int sumiteration = 0;
	Matrix V_0(rows,1);
	Matrix U_0(rows,1);
	Matrix V_1(rows,1);
	Matrix U_1(rows,1);
	Matrix test(rows,1);
	
	double N_0;
	double N_1;
	double N_new;
	double N_old;
	//Seeding V_old  dim M x 1 with random non zeros from Mersenne Twister
	// generating random real numbers between 0,1
	block::MT ran;
		auto rand_real = ran.rand_real(0.0,1.0);
		for (int i=1; i <= rows; i++)
		{
			V_0.getVal(i,1) = rand_real();
			//std::cout << "The ith element of V_old is " << i<< " " << V_0.getVal(i,1) << std::endl;
		}
	int minIterations = 25;
	double fracChange;
	N_0 = V_0.norm();			// No
	U_0 = 1/N_0 * V_0; 			// Uo
	do	
	{
		V_1 = (*this) * U_0;	// V1 
		N_1 = V_1.norm();		// N1
		fracChange = std::abs((N_1-N_0)/N_0);
		U_1 = 1/N_1 * V_1;		// U1
		sumiteration++;
		U_0 = U_1;				//over write U1 for loop
		N_0 = N_1;				//over write N1 for loop
	} 
	while(sumiteration < minIterations || fracChange > tol );
	std::cout << "The total sum of power iteration is :" << sumiteration << std::endl;
	// save eigenvector to eigenVector to vector object and print it out
	this -> eigenVector.resize(rows);
	std::cout << "The eigenvector is :" << std::endl;
	for (int i = 1; i <= rows; i++)
	{
		this -> eigenVector.at(i-1)=U_1.getVal(i,1);
		std::cout << eigenVector.at(i-1) << std::endl; //print out the vector
	}
	
	// Check if opposite signs present from test matrix to U_1
	//U_1 on Matrix A again to check signs
	test = (*this) * U_1; 
	if ( (test.getVal(1,1) < 0 )? (U_1.getVal(1,1) > 0): (U_1.getVal(1,1) < 0))
	{
		this -> eigenValue = -N_1 + shift;
		std::cout << "The eigenvalue of largest magnitude is " << eigenValue << std::endl;
	}
	else
	{
		this -> eigenValue = N_1 + shift;
		std::cout << "The eigenvalue of largest magnitude is " << eigenValue << std::endl;
	}
	eigenDone = true;

}

void block::Matrix::findLargestLambda()
{
	findLargestLambda(0.0);
}

block::Matrix block::Matrix::identity(int n)
{
	block::Matrix ident(n);
	for(int i=1; i<=n; i++)
		ident.getVal(i,i) = 1.0;
	return ident;
}

void block::Matrix::print() const
{
    int row, col;
    int N = this -> rows;
    int M = this -> cols;
    
    for (row = 1; row <= N; row++)
    {
        for (col = 1; col <= M; col++)
        {
            std::cout << getVal(row,col) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
