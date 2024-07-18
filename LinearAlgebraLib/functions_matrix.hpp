#ifndef FUNCTIONS_MATRIX_HPP
#define FUNCTIONS_MATRIX_HPP

#include "Matrix.hpp"

namespace lalib { // starting namespace lalib

// returns the minor of a matrix
template <typename MatrixType>
MatrixType Minor(const Matrix<MatrixType> mat,const int iv,const int jv);

// returns the cofactor of a matrix
template <typename MatrixType>
MatrixType cofactor(const Matrix<MatrixType> mat,const int iv,const int jv);

// returns adjucate of a matrix
template <typename MatrixType>
Matrix<MatrixType> adj(const Matrix<MatrixType> mat);

// returns inverse of a matrix
template <typename MatrixType>
Matrix<MatrixType> inv(const Matrix<MatrixType> mat);

// Returns the transpose
template <typename MatrixType>
Matrix<MatrixType> transpose(const Matrix<MatrixType> mat);

// returns determinant of a matrix
template <typename MatrixType>
MatrixType det(const Matrix<MatrixType> mat);

// returns trace of a matrix
template <typename MatrixType>
MatrixType tr(const Matrix<MatrixType> mat);

// (*) Multiplies matrix with a number and returns the result
template<typename MatrixType , typename T >
Matrix<MatrixType> operator* (const T number,Matrix<MatrixType> mat);

// Prints a matrix in a nice formated way
template <typename T>
void printMatrix(Matrix<T> arr,int precisiond = 100);



} // ending namespace lalib

#endif