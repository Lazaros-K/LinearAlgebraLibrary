#include "Matrix.hpp"
#include <iostream>
#include <cmath>

namespace lalib { // starting namespace lalib

/*                                           */
/*    ===================================    */
/*       |                           |       */
/*       |  Matrix class functions   |       */
/*       |                           |       */
/*    ===================================    */
/*                                           */


template <typename MatrixType>
void Matrix<MatrixType>::Allocate(int nc,int mc) {

    n = nc;
    m = mc;
    ptr = (MatrixType**)calloc(n, sizeof(MatrixType));

    if(ptr == NULL) {
        std::cout << " Allocation FAILED" << std::endl;
        return;
    }

    for(int i = 0; i < n; i++) {
        ptr[i] = (MatrixType*)calloc(m, sizeof(MatrixType*));
        
        if(ptr[i] == NULL) {
            std::cout << " Allocation FAILED at column number " << i+1 << std::endl;
            return;
        }
    }

    return;
}


template <typename MatrixType>
void Matrix<MatrixType>::ChangeSize(int nc,int mc) {
    
    if(n == nc && m == mc) {
        std::cout << "Size changing operation skiped. Matrix is already the desired size." << std::endl;
        return;
    }


    MatrixType **new_ptr;

    new_ptr = (MatrixType**)calloc(nc , sizeof(MatrixType*));
    if(new_ptr == NULL) {
        std::cout << " Allocation FAILED" << std::endl;
        return;
    }        

    for(int i = 0; i < n; i++) {
        new_ptr[i] = (MatrixType*)realloc(ptr[i],mc * sizeof(MatrixType));

        if(new_ptr[i] == NULL) {
            std::cout << " RE-Allocation FAILED at column number " << i+1 << std::endl;
            return;
        }
        
    }

    for(int i = n; i < nc; i++) {
        new_ptr[i] = (MatrixType*)calloc(mc, sizeof(MatrixType));

        if(new_ptr[i] == NULL) {
            std::cout << " Allocation FAILED at column number " << i+1 << std::endl;
            return;
        }
        
    }

    
    ptr = new_ptr;
    

    n = nc;
    m = mc;

    return;
}


template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::KroneckerProduct(const Matrix<MatrixType>& arr2) const {

    const int Fm = m*arr2.m, Fn = n*arr2.n;

    Matrix<MatrixType> Farr(Fn,Fm,true);

    for (int j=0; j < m; j++) {
        for (int i=0; i < n; i++) {
            for (int k=0; k < arr2.m; k++) {
                for (int l=0; l < arr2.n; l++) {
                    Farr.ptr[i*arr2.n+l][j*arr2.m+k] = ptr[i][j] * arr2.ptr[l][k];
                }
            }
        }
    }
    return Farr;
} 


template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::HadamardProduct(const Matrix<MatrixType>& arr2) const {

    if(n != arr2.n || m != arr2.m) {
        std::cout << "Can't find Hadamard's Product. The matricies size is not equal" << std::endl;
    }

    Matrix<MatrixType> Farr(n,m,true);

    for (int j=0; j < m; j++) {
        for (int i=0; i < n; i++) {
            Farr.ptr[i][j] = ptr[i][j] * arr2.ptr[i][j];
        }
    }
    return Farr;
}


template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::DirectSum(const Matrix<MatrixType>& arr2) const {

    const int Fm = m+arr2.m, Fn = n+arr2.n;

    Matrix<MatrixType> Farr(Fn,Fm,true);

    for (int j=0; j < Fm; j++) {
        for (int i=0; i < Fn; i++) {
            if(i < n && j < m) {
                Farr.ptr[i][j] = ptr[i][j];
            }
            else if(i >= n && j >= m) {
                Farr.ptr[i][j] = arr2.ptr[i-n][j-m];
            }
            else {
                Farr.ptr[i][j] = 0;
            }
        }
    }
    return Farr;
}

template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::Multiply(const Matrix<MatrixType>& arr2) const {

    if(!(n == arr2.m)) {
        std::cout << "Can't multiply the Matricies. Wrong sizes." << std::endl;

        Matrix<MatrixType> Farr(0,0,true);
        return Farr;
    }

    const int Fm = m, Fn = arr2.n;

    Matrix<MatrixType> Farr(Fn,Fm,true);

    for (int j=0; j < Fm; j++) {
        for (int i=0; i < Fn; i++) {

            for (int k=0; k < arr2.m; k++) {
                //std::cout << Farr.ptr[i][j];
                Farr.ptr[i][j] += ptr[k][j] * arr2.ptr[i][k];
                //std::cout << " + (" << ptr[k][j] << "*" << arr2.ptr[i][k] << ") = " << Farr.ptr[i][j] << std::endl;
            }

        }
    }
    return Farr;
} 
template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::operator* (const Matrix<MatrixType>& arr2) const {
    return Multiply(arr2);
}
template <typename MatrixType>
void Matrix<MatrixType>::operator*= (const Matrix<MatrixType>& arr2) {
    *this = Multiply(arr2);
}


template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::Add(const Matrix<MatrixType> arr2) const {

    if(!(n == arr2.n && m == arr2.m)) {
        std::cout << "Can't add the Matrixs. Wrong sizes." << std::endl;

        Matrix<MatrixType> Farr(0,0);
        return Farr;
    }

    const int Fm = m, Fn = arr2.n;

    Matrix<MatrixType> Farr(Fn,Fm,true);

    for (int j=0; j < Fm; j++) {
        for (int i=0; i < Fn; i++) {

        Farr.ptr[i][j] = ptr[i][j] + arr2.ptr[i][j];

        }
    }
    return Farr;
}
template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::operator+ (const Matrix<MatrixType>& arr2) const {
    return Add(arr2);
}
template <typename MatrixType>
void Matrix<MatrixType>::operator+= (const Matrix<MatrixType>& arr2) {
    *this = Add(arr2);
}
template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::operator- (const Matrix<MatrixType>& arr2) const {
    return Add(arr2*-1);
}
template <typename MatrixType>
void Matrix<MatrixType>::operator-= (const Matrix<MatrixType>& arr2) {
    *this = Add(arr2*-1);
}



template <typename MatrixType>
template< typename T >
Matrix<MatrixType> Matrix<MatrixType>::Multiply(const T number) const {


    const int Fm = m, Fn = n;

    Matrix<MatrixType> Farr(Fn,Fm,true);

    for (int j=0; j < Fm; j++) {
        for (int i=0; i < Fn; i++) {

        Farr.ptr[i][j] = ptr[i][j] * number;

        }
    }
    return Farr;
}
template <typename MatrixType>
template< typename T >
Matrix<MatrixType> Matrix<MatrixType>::operator* (const T number) const {
    return Multiply(number);
}
template <typename MatrixType>
template< typename T >
void Matrix<MatrixType>::operator*= (const T number) {
    *this = Multiply(number);
}

template <typename MatrixType>
template< typename T >
Matrix<MatrixType> Matrix<MatrixType>::divide(const T number) const {
    return Multiply(MatrixType(1/number));
}
template <typename MatrixType>
template< typename T >
Matrix<MatrixType> Matrix<MatrixType>::operator/ (const T number) const {
    return divide(MatrixType(number));
}
template <typename MatrixType>
template< typename T >
void Matrix<MatrixType>::operator/= (const T number) {
    *this = divide(MatrixType(number));
}

template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::power(const int number) const {

    if(number == 1) {
        return *this;
    }

    if(!(n == m)) {
        std::cout << "Can't raise the Matrix to a number. Wrong sizes." << std::endl;

        Matrix<MatrixType> Farr(0,0,true);
        return Farr;
    }

    if(number == 0) {
        Matrix<MatrixType> Iarr(n,m,true);
        for(int i = 0; i<n; i++) {
            for(int j = 0; j<m; j++) {
                if(i == j) {
                    Iarr.ptr[i][j] = 1;
                }
                else {
                    Iarr.ptr[i][j] = 0;
                }
            }
        }
        return Iarr;
    }

    Matrix<MatrixType> Farr(n,m,true);

    Farr = Multiply(*this);
    for(int i; i<number-2; i++) {
        Farr = Farr.Multiply(*this);
    }
    return Farr;
}
template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::operator^ (const int number) {
    return power(number);
}

// template <typename MatrixType>
// void Matrix<MatrixType>::transpose() {
//     *this = transpose(*this);
//     return;
// }



template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::Submatrix(const int iv,const int jv) const {

    if(iv >= m || jv >= n || iv < 0 || jv < 0) {
        std::cout << "Can't find Submatrix of the Matrix. Wrong index." << std::endl;
        Matrix Submatrixzero;
        return Submatrixzero;
    }

    Matrix<MatrixType> Submatrix(n-1,m-1);

    int i2 = 0;
    int j2 = 0;
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            if(j != jv && i != iv) {
                Submatrix.ptr[i2][j2] = ptr[i][j]; 
            }
            if(j != jv) {
                j2 ++;
            }
        }
        if(i != iv) {
            i2 ++;
        }
        j2 = 0;
    }

    return Submatrix;
}


template <typename MatrixType>
bool Matrix<MatrixType>::equal(const Matrix<MatrixType> mat) const {


    const int Fm = m, Fn = n;
    int trueC = 0;

    for (int j=0; j < Fm; j++) {
        for (int i=0; i < Fn; i++) {
            if(ptr[i][j] == mat.ptr[i][j]) {
                trueC ++;
            }
        }
    }

    if(trueC == (Fm*Fn)) {
        return true;
    }
    return false;
}
template <typename MatrixType>
bool Matrix<MatrixType>::operator== (const Matrix<MatrixType> mat) const {
    return equal(mat);
}
template <typename MatrixType>
bool Matrix<MatrixType>::operator!= (const Matrix<MatrixType> mat) const {
    return !(equal(mat));
}


// elementary row operation

template <typename MatrixType>
void Matrix<MatrixType>::RowSwap(int j1, int j2) {
    if (j1 == j2 ) {
        return; 
    }
    else if ((j1 >= m && j1 < 0 ) || (j2 >= m && j2 < 0 )) {
        std::cout<< "Row number is not matching the matrix" <<std::endl;
        return; 
    }

    for(int i = 0; i < n; i++) {
        MatrixType x = ptr[i][j1];
        ptr[i][j1] = ptr[i][j2];
        ptr[i][j2] = x;
    }
}

template <typename MatrixType>
void Matrix<MatrixType>::RowScalarMult(int j1,MatrixType x) {
    
    for(int i = 0; i < n; i++) {
        if (ptr[i][j1] != 0) {
            ptr[i][j1] *= x;
        }
    }
}

template <typename MatrixType>
void Matrix<MatrixType>::RowSum(int j1,int j2,MatrixType x) {
    // j1 is the multiplied row
    // j2 is the mutated row (sum of j2 and j1*x)
    for(int i = 0; i < n; i++) {
        if (ptr[i][j1] != 0) {
            MatrixType Temp = ptr[i][j1] * x;
            ptr[i][j2] += Temp;
        }
    }
}

template <typename MatrixType>
void Matrix<MatrixType>::Gauss() {
    int k = 0; // number of top padding
    for(int i = 0; i < m; i++) {
        // find zero and bring it at the first row
        for(int j = 1+i; j < m; j++) {
            if(ptr[i][0+i] != 0) {
                break;
            }
            if(ptr[i][j] != 0) {
                RowSwap(0+i,j);
                break;
            }
            if(j == m-1) {
                i++;
                j = 1+k;
            }
        }
        
        // subtract to make a 0
        for(int j = 1+k; j < m; j++) {
            if (ptr[i][j] != 0) {
                RowSum(0+k,j,-(ptr[i][j]/ptr[i][0+k]));
            }
        }
        k++;
    }
}

template <typename MatrixType>
void Matrix<MatrixType>::Gauss_Jordan() {
    
    Gauss();
    
    
    for(int j = m-1; j >= 0; j--) {
        for(int i = 0; i < n; i++) {
            if(ptr[i][j] != 0) {
                RowScalarMult(j,1/ptr[i][j]);
                for(int k = 1; k <= j; k++) {
                    RowSum(j,j-k,-ptr[i][j-k]);
                }
                break;
            }
        }
    }

}



// ALL POSSIBLE CLASS TYPES

template class Matrix<float>;
template class Matrix<double>;
template class Matrix<int>;

template Matrix<double> Matrix<double>::Multiply<double>(double) const;
template Matrix<double> Matrix<double>::operator*<double>(double) const;
template Matrix<double> Matrix<double>::divide<double>(double) const;
template Matrix<double> Matrix<double>::operator/<double>(double) const;

template Matrix<int> Matrix<int>::Multiply<double>(double) const;
template Matrix<int> Matrix<int>::operator*<double>(double) const;
template Matrix<int> Matrix<int>::divide<double>(double) const;
template Matrix<int> Matrix<int>::operator/<double>(double) const;

template Matrix<float> Matrix<float>::Multiply<double>(double) const;
template Matrix<float> Matrix<float>::operator*<double>(double) const;
template Matrix<float> Matrix<float>::divide<double>(double) const;
template Matrix<float> Matrix<float>::operator/<double>(double) const;

// ----------------------------------------------------------------------

template Matrix<double> Matrix<double>::Multiply<int>(int) const;
template Matrix<double> Matrix<double>::operator*<int>(int) const;
template Matrix<double> Matrix<double>::divide<int>(int) const;
template Matrix<double> Matrix<double>::operator/<int>(int) const;

template Matrix<int> Matrix<int>::Multiply<int>(int) const;
template Matrix<int> Matrix<int>::operator*<int>(int) const;
template Matrix<int> Matrix<int>::divide<int>(int) const;
template Matrix<int> Matrix<int>::operator/<int>(int) const;

template Matrix<float> Matrix<float>::Multiply<int>(int) const;
template Matrix<float> Matrix<float>::operator*<int>(int) const;
template Matrix<float> Matrix<float>::divide<int>(int) const;
template Matrix<float> Matrix<float>::operator/<int>(int) const;

// ----------------------------------------------------------------------

template Matrix<double> Matrix<double>::Multiply<float>(float) const;
template Matrix<double> Matrix<double>::operator*<float>(float) const;
template Matrix<double> Matrix<double>::divide<float>(float) const;
template Matrix<double> Matrix<double>::operator/<float>(float) const;

template Matrix<int> Matrix<int>::Multiply<float>(float) const;
template Matrix<int> Matrix<int>::operator*<float>(float) const;
template Matrix<int> Matrix<int>::divide<float>(float) const;
template Matrix<int> Matrix<int>::operator/<float>(float) const;

template Matrix<float> Matrix<float>::Multiply<float>(float) const;
template Matrix<float> Matrix<float>::operator*<float>(float) const;
template Matrix<float> Matrix<float>::divide<float>(float) const;
template Matrix<float> Matrix<float>::operator/<float>(float) const;



} // ending namespace lalib