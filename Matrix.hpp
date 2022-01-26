#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <vector>

namespace block
{
    class Matrix
    {
        private:
        int rows; // number of rows in the Matrix
        int cols; // number of columns in the Matrix
        std::vector<std::vector<double> > mat; // internal data structure to store the Matrix values
        bool eigenDone; // flag denoting if the findLargestLambda() function has been run
        double eigenValue; // stores the eigenvalue found with findLargestLambda()
        std::vector<double> eigenVector; // stores the eigenvector associated with the above eigenvalue
        
        void init(int rows, int cols); // run every time a Matrix is created or changed
        int checkInput(std::string input); // validates a user's input and returns a code
        int solicitUserInput(); // continuously polls the user for input until it validates

        // helper function to copy Matrix values from a given Matrix to this Matrix
        void copyVals(const Matrix& matToCopy);
        
        public:
        // sets Matrix size and values based on a file; can be called independently or through the fourth constructor
        void setFromFile(std::string filename);
        
        // default constructor
        Matrix();
        
        // minimal non-default constructor
        Matrix(int N);
        
        // full constructor
        Matrix(int N, int M);
        
        // constructor that reads from file
        Matrix(std::string filename);
        
        // destructor; handles new Matrix from findLargestLambda
        ~Matrix();
        
        // gets
        const int getRows() const; // returns rows
        const int getCols() const; // returns cols
        const double getEigenValue() const; // returns eigenValue
        const std::vector<double> getEigenVector() const; // returns (*eigenVector)
        
        // get/set values
        const double getVal(int n, int m) const; // used for getting
        double& getVal(int n, int m); // used for setting
        
        // operator overloads
        std::vector<double>& operator[](std::size_t idx);
        const std::vector<double> operator[](std::size_t idx) const;
        Matrix& operator=(const Matrix& mat);
        Matrix& operator+=(const Matrix& rhs);
        Matrix& operator-=(const Matrix& rhs);
        Matrix& operator*=(const double rhs);
        Matrix& operator*=(const Matrix& rhs);
        friend Matrix operator+(Matrix lhs, const Matrix& rhs);
        friend Matrix operator-(Matrix lhs, const Matrix& rhs);
        friend Matrix operator*(Matrix lhs, const double rhs);
        friend Matrix operator*(const double lhs, Matrix rhs);
        friend Matrix operator*(Matrix lhs, const Matrix& rhs);
        
        // special matrix operations
        Matrix T() const; // returns the transpose of this Matrix
        double norm() const; // returns the norm of this Matrix (must be 1D)
        void findLargestLambda(double shift); // finds the eigenvalue of largest magnitude and its corresponding eigenvector; sets them to local attributes and changes the eigenDone flag to true when it runs successfully.
        void findLargestLambda(); // simply calls the function above with argument 0.0
        Matrix identity(int n);
        
        // class method to print out the matrix
        void print() const;

    };
}

#endif
