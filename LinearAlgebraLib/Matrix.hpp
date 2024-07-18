#ifndef LINEAR_ALGEBRA_LIBRARY_HPP
#define LINEAR_ALGEBRA_LIBRARY_HPP

#include <iostream>

namespace lalib { // starting namespace lalib



template <typename MatrixType = double>
class Matrix {
    public: 
        MatrixType **ptr;

    private:
        int n;
        int m;

    void Allocate(int nc,int mc);
    
    public: 

    // Matrix Constructor (size and verbal is optional)
    Matrix(int nc = 0, int mc = 0,bool verbal = false);

    // Matrix Constructor using write (verbal is optional)
    template <typename T, size_t m_size, size_t n_size>
    Matrix(T (&data)[m_size][n_size],bool verbal = false);

    // Returns N
    int GetN() const {return n;}
    // Returns M
    int GetM() const {return m;}

    // Alters the size of the matrix
    void ChangeSize(int nc,int mc);

    // Apply the data and size of an array to the matrix
    template <typename T, size_t m_size, size_t n_size>
    void Write(T (&data)[m_size][n_size]);

    // returns Kronecker Product
    Matrix<MatrixType> KroneckerProduct(const Matrix<MatrixType>& arr2) const;

    // returns Hadamard Product
    Matrix<MatrixType> HadamardProduct(const Matrix<MatrixType>& arr2) const;
    
    // reurns direct sum
    Matrix<MatrixType> DirectSum(const Matrix<MatrixType>& arr2) const;


    // Multiplies two matricies and returns the result
    Matrix<MatrixType> Multiply(const Matrix<MatrixType>& arr2) const;
    // (*) Multiplies two matricies and returns the result
    Matrix<MatrixType> operator* (const Matrix<MatrixType>& arr2) const;
    // (*=) Multiply two matricies and asigns it
    void operator*= (const Matrix<MatrixType>& arr2);


    // Adds two matricies and returns the result
    Matrix<MatrixType> Add(const Matrix<MatrixType> arr2) const;
    // (+) Adds two matricies and returns the result
    Matrix<MatrixType> operator+ (const Matrix<MatrixType>& arr2) const;
    // (+=) Adds two matricies and asigns it
    void operator+= (const Matrix<MatrixType>& arr2);
    // (-) Subtracts two matricies and returns the result
    Matrix<MatrixType> operator- (const Matrix<MatrixType>& arr2) const;
    // (-=) Subtracts two matricies and asigns it
    void operator-= (const Matrix<MatrixType>& arr2);


    // Multiplies matrix with a number and returns the result
    template< typename T = int>
    Matrix<MatrixType> Multiply(const T number) const;
    // (*) Multiplies matrix with a number and returns the result
    template< typename T = int>
    Matrix<MatrixType> operator* (const T number) const;
    // (*=) Multiplies matrix and asigns it
    template< typename T = int>
    void operator*= (const T number);

    // Divides matrix with a number and returns the result
    template< typename T = int>
    Matrix<MatrixType> divide(const T number) const;
    // (/) Divides matrix with a number and returns the result
    template< typename T = int>
    Matrix<MatrixType> operator/ (const T number) const;
    // (/=) Divides matrix with a number and and asigns it
    template< typename T = int>
    void operator/= (const T number);

    // Returns the power of a matrix
    Matrix<MatrixType> power(const int number) const;
    // (^) Returns the power of a matrix
    Matrix<MatrixType> operator^ (const int number);

    // Applies transpose
    void transpose();

    // returns the submatrix created by removing the rows and columns of the given index
    Matrix<MatrixType> Submatrix(const int iv,const int jv) const;

    // returns TRUE for equal FALSE for unequal
    bool equal(const Matrix<MatrixType> mat) const;
    // (==) returns TRUE for equal FALSE for unequal
    bool operator== (const Matrix<MatrixType> mat) const;
    // (!=) returns FALSE for equal TRUE for unequal
    bool operator!= (const Matrix<MatrixType> mat) const;


    // elementary row operation

    // Sawps two rows of a matrix
    void RowSwap(int j1, int j2);

    // Multiply row of a matrix with a number
    void RowScalarMult(int j1,MatrixType x);    

    // Multiply row of a matrix with a number then adds the result to another row
    void RowSum(int j1,int j2,MatrixType x);

    // Applies Gauss Elimination (Απαλοιφή Gauss) converting the matrix to row-echelon form / ref (κλιμακωτού πίνακα) 
    void Gauss();

    // Applies Gauss-Jordan Elimination (Απαλοιφή Gauss-Jordan) converting the matrix to reduced-row-echelon form / rref (κλιμακωτού πίνακα)
    void Gauss_Jordan();
};




template <typename MatrixType>
Matrix<MatrixType>::Matrix(int nc, int mc,bool verbal) { 

    Allocate(nc,mc);
    if (verbal) {
        if(n == 0 && m == 0) {
            std::cout << "Created zero size Matrix." << std::endl;
        }
        else {
            std::cout << "Created Matrix.\nm: " << m << ", n: " << n << std::endl;
        }
    }
};

template <typename MatrixType>
template <typename T, size_t m_size, size_t n_size>
Matrix<MatrixType>::Matrix(T (&data)[m_size][n_size],bool verbal) { 

    Allocate(n_size,m_size);
    Write(data);
    if (verbal) {
        if(n == 0 && m == 0) {
            std::cout << "Created zero size Matrix." << std::endl;
        }
        else {
            std::cout << "Created Matrix with Values.\nm: " << m << ", n: " << n << std::endl;
        }
    }
}


template <typename MatrixType>
template <typename T, size_t m_size, size_t n_size>
void Matrix<MatrixType>::Write(T (&data)[m_size][n_size]) {

    if(!(std::is_same<T,MatrixType>::value)) {
        std::cout << "Writing data failed. Wrong array type" << std::endl;
        return;
    }

    if(m_size != m || n_size != n) {
        std::cout << "Writing data failed. Wrong array size" << std::endl;
        return;
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            ptr[i][j] = data[i][j];
        }
    }
    return;
}



} // ending namespace lalib

#endif