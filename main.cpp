#include "Matrix.hpp"
#include "mt.hpp"
#include <iostream>
#include <cstdlib>
#include <string>

using block::Matrix;

int main(int argc, char** argv)
{
	std::string filename;
	
		std::cout << "What is the full filename including "".txt"" of the Matrix to find its smallest eigenvalue and eigenvector?" << std::endl << std::endl;
		std::cin >> filename;
	
	Matrix mat1;
	mat1.setFromFile(filename);
	mat1.findLargestLambda();
	if (mat1.getEigenValue() < 0 )
	{
		std::cout << "You have found the smallest Eigenvalue and Eigenvector! Yay!" << std::endl << std::endl;
	}
	else{
		std::cout << "You have failed to find the smallest Eigenvalue, will rerun with new Matrix!" << std::endl << std::endl;
	Matrix Ident;
	Ident = Ident.identity(mat1.getRows());
	Matrix mat2 = mat1 - (mat1.getEigenValue() * Ident);
	mat2.findLargestLambda(mat1.getEigenValue());
	
	}
}