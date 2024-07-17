#ifndef LINEAR_ALGEBRA_LIBRARY_HPP
#define LINEAR_ALGEBRA_LIBRARY_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <chrono>
#include <thread>

#include "tools/Timer.hpp"

using namespace std;


template <typename MatrixType>
class Matrix {
    public: 
        MatrixType **ptr;

    private:
        int n;
        int m;

    void Allocate(int nc,int mc) {

        n = nc;
        m = mc;
        ptr = (MatrixType**)calloc(n, sizeof(MatrixType));

        if(ptr == NULL) {
            cout << " Allocation FAILED" << endl;
            return;
        }

        for(int i = 0; i < n; i++) {
            ptr[i] = (MatrixType*)calloc(m, sizeof(MatrixType*));
            
            if(ptr[i] == NULL) {
                cout << " Allocation FAILED at column number " << i+1 << endl;
                return;
            }
        }

        return;
    }
    
    public: 

    Matrix(int nc, int mc,bool verbal = false);

    template <typename T, size_t m_size, size_t n_size>
    Matrix(T (&data)[m_size][n_size],bool verbal = false);


    int GetN() const {return n;}
    int GetM() const {return m;}

    void ChangeSize(int nc,int mc);

    template <typename T, size_t m_size, size_t n_size>
    void Write(T (&data)[m_size][n_size]);

    Matrix<MatrixType> KroneckerProduct(const Matrix<MatrixType>& arr2) const;

    Matrix<MatrixType> HadamardProduct(const Matrix<MatrixType>& arr2) const;
    

    Matrix<MatrixType> directSum(const Matrix<MatrixType>& arr2) const;

    Matrix<MatrixType> multiply(const Matrix<MatrixType>& arr2) const;

    Matrix<MatrixType> operator* (const Matrix<MatrixType>& arr2) const;
    void operator*= (const Matrix<MatrixType>& arr2);


    Matrix<MatrixType> Add(const Matrix<MatrixType> arr2) const;
    Matrix<MatrixType> operator+ (const Matrix<MatrixType>& arr2) const;
    void operator+= (const Matrix<MatrixType>& arr2);
    Matrix<MatrixType> operator- (const Matrix<MatrixType>& arr2) const;
    void operator-= (const Matrix<MatrixType>& arr2);



    template< typename T >
    Matrix<MatrixType> multiply(const T number) const;
    template< typename T >
    Matrix<MatrixType> operator* (const T number) const;
    template< typename T >
    void operator*= (const T number);
    template< typename T >
    Matrix<MatrixType> divide(const T number) const;
    template< typename T >
    Matrix<MatrixType> operator/ (const T number) const;
    template< typename T >
    void operator/= (const T number);

    template< typename T >
    Matrix<MatrixType> power(T number) const;

    template< typename T >
    Matrix<MatrixType> operator^ (const T number);

    Matrix<MatrixType> gettranspose() const;
    void transpose();

    // returns determinant of a matrix
    MatrixType det() const;

    // returns trace of a matrix
    MatrixType tr() const;

    Matrix<MatrixType> Minor(const int iv,const int jv) const;

    MatrixType cofactor(const int iv,const int jv) const;

    // returns adjucate of a matrix
    Matrix<MatrixType> adj() const;

    // returns inverse of a matrix
    Matrix<MatrixType> inv() const;

    bool equal(Matrix<MatrixType> mat) const;
    bool operator== (Matrix<MatrixType> mat) const;
    bool operator!= (Matrix<MatrixType> mat) const;


    // elementary row operation

    void RowSwap(int j1, int j2);

    void RowScalarMult(int j1,MatrixType x);    

    void RowSum(int j1,int j2,MatrixType x);

    void Gauss();
    void Gauss_Jordan();
};












template<typename MatrixType , typename T >
Matrix<MatrixType> operator* (const T number,Matrix<MatrixType> mat);


template <typename T>
void printMatrix(Matrix<T> arr,int precisiond = 100);





























// from number
template<typename MatrixType , typename T >
Matrix<MatrixType> operator* (const T number,Matrix<MatrixType> mat) {
    return mat * number;
}

template <typename T>
void printMatrix(Matrix<T> arr,int precisiond) {


    if(arr.GetN() == 0 && arr.GetM()== 0) {
        cout << "||" << endl;
    }

    int maxlength = 0;
    int maxd = 0;
    int countd[arr.GetM()*arr.GetN()] = {0};
    int decimal_limit = numeric_limits<T>::max_digits10;

    if (decimal_limit > precisiond) {
        decimal_limit = precisiond;
    }

    for (int i = 0; i < arr.GetM(); i++)
    {
        for (int j = 0; j < arr.GetN(); j++)
        {
            
            // get desimal size
            T num = abs(arr.ptr[j][i]);

            int exponent = 1.0;
            T epsilon = numeric_limits<T>::epsilon() * num;
            T c = num - floor(num);

            

            while((c > epsilon && c < (1 - epsilon)) && countd[j+i*arr.GetN()] < decimal_limit)
            {
                exponent *= 10;

                c = num * exponent;
                c = c - floor(c);

                epsilon = numeric_limits<T>::epsilon() * exponent * num;

                countd[j+i*arr.GetN()]++;
            }

            
            string s = to_string((int)arr.ptr[j][i]);
            int sl = s.size();
            // get number size
            if((int)arr.ptr[j][i] == 0 && arr.ptr[j][i] < 0) {
                // add one more character
                sl ++;
            }

            if(countd[j+i*arr.GetN()] > 0) {
                
                sl += countd[j+i*arr.GetN()] + 1;
            }

            // max decimals digits
            if (countd[j+i*arr.GetN()] > maxd) {
                maxd = countd[j+i*arr.GetN()];
            }

            // max digits and decimal digits
            if (sl > maxlength) {
                maxlength = sl;
            }

        }
    }
    // cout << maxlength << endl;
    // cout << maxd << endl;


    string pavla;
    for (int i = 0; i < ((arr.GetN())*(maxlength)+(arr.GetN()-1)*3); i++)
    {
        pavla.append("-");
    }
    

    cout << "--" << pavla << "--" << endl;
    for (int i = 0; i < arr.GetM(); i++)
    {
        for (int j = 0; j < arr.GetN(); j++)
        {
            cout << "| " ;
            cout << setw(maxlength);
            cout << setprecision(countd[j+i*arr.GetN()]);
            cout << fixed;
            cout << T(arr.ptr[j][i]);
            cout << " ";
        }
        cout << "|" << endl;
        cout << "--" << pavla << "--" << endl;
    }
    cout << endl;
    return;
}


template <typename MatrixType>
Matrix<MatrixType>::Matrix(int nc, int mc,bool verbal) { // Constructor with parameters

    Allocate(nc,mc);
    if (verbal) {
        if(n == 0 && m == 0) {
            cout << "Created zero size Matrix." << endl;
        }
        else {
            cout << "Created Matrix.\nm: " << m << ", n: " << n << endl;
        }
    }
};

template <typename MatrixType>
template <typename T, size_t m_size, size_t n_size>
Matrix<MatrixType>::Matrix(T (&data)[m_size][n_size],bool verbal) { // Constructor with parameters and values

    Allocate(n_size,m_size);
    Write(data);
    if (verbal) {
        if(n == 0 && m == 0) {
            cout << "Created zero size Matrix." << endl;
        }
        else {
            cout << "Created Matrix with Values.\nm: " << m << ", n: " << n << endl;
        }
    }
}

template <typename MatrixType>
void Matrix<MatrixType>::ChangeSize(int nc,int mc) {
    
    if(n == nc && m == mc) {
        cout << "Size changing operation skiped. Matrix is already the desired size." << endl;
        return;
    }


    MatrixType **new_ptr;

    new_ptr = (MatrixType**)calloc(nc , sizeof(MatrixType*));
    if(new_ptr == NULL) {
        cout << " Allocation FAILED" << endl;
        return;
    }        

    for(int i = 0; i < n; i++) {
        new_ptr[i] = (MatrixType*)realloc(ptr[i],mc * sizeof(MatrixType));

        if(new_ptr[i] == NULL) {
            cout << " RE-Allocation FAILED at column number " << i+1 << endl;
            return;
        }
        
    }

    for(int i = n; i < nc; i++) {
        new_ptr[i] = (MatrixType*)calloc(mc, sizeof(MatrixType));

        if(new_ptr[i] == NULL) {
            cout << " Allocation FAILED at column number " << i+1 << endl;
            return;
        }
        
    }

    
    ptr = new_ptr;
    

    n = nc;
    m = mc;

    return;
}

template <typename MatrixType>
template <typename T, size_t m_size, size_t n_size>
void Matrix<MatrixType>::Write(T (&data)[m_size][n_size]) {

    if(!(is_same<T,MatrixType>::value)) {
        cout << "Writing data failed. Type of Matrix data is wrong" << endl;
        return;
    }

    if(m_size != m || n_size != n) {
        cout << "Writing data failed. Size of Matrix is wrong" << endl;
        return;
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            ptr[j][i] = data[i][j];
        }
    }
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
        cout << "Can't find Hadamard's Product. The matricies size is not equal" << endl;
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
Matrix<MatrixType> Matrix<MatrixType>::directSum(const Matrix<MatrixType>& arr2) const {

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
Matrix<MatrixType> Matrix<MatrixType>::multiply(const Matrix<MatrixType>& arr2) const {

    if(!(n == arr2.m)) {
        cout << "Can't multiply the Matricies. Wrong sizes." << endl;

        Matrix<MatrixType> Farr(0,0,true);
        return Farr;
    }

    const int Fm = m, Fn = arr2.n;

    Matrix<MatrixType> Farr(Fn,Fm,true);

    for (int j=0; j < Fm; j++) {
        for (int i=0; i < Fn; i++) {

            for (int k=0; k < arr2.m; k++) {
                //cout << Farr.ptr[i][j];
                Farr.ptr[i][j] += ptr[k][j] * arr2.ptr[i][k];
                //cout << " + (" << ptr[k][j] << "*" << arr2.ptr[i][k] << ") = " << Farr.ptr[i][j] << endl;
            }

        }
    }
    return Farr;
} 
template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::operator* (const Matrix<MatrixType>& arr2) const {
    return multiply(arr2);
}
template <typename MatrixType>
void Matrix<MatrixType>::operator*= (const Matrix<MatrixType>& arr2) {
    *this = multiply(arr2);
}


template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::Add(const Matrix<MatrixType> arr2) const {

    if(!(n == arr2.n && m == arr2.m)) {
        cout << "Can't add the Matrixs. Wrong sizes." << endl;

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
    return Add(-1*arr2);
}
template <typename MatrixType>
void Matrix<MatrixType>::operator-= (const Matrix<MatrixType>& arr2) {
    *this = Add(-1*arr2);
}



template <typename MatrixType>
template< typename T >
Matrix<MatrixType> Matrix<MatrixType>::multiply(const T number) const {


    const int Fm = m, Fn = n;

    Matrix<MatrixType> Farr(Fn,Fm,true);

    for (int j=0; j < Fm; j++) {
        for (int i=0; i < Fn; i++) {

        Farr.ptr[i][j] = ptr[i][j] * (MatrixType)number;

        }
    }
    return Farr;
}
template <typename MatrixType>
template< typename T >
Matrix<MatrixType> Matrix<MatrixType>::operator* (const T number) const {
    return multiply(number);
}
template <typename MatrixType>
template< typename T >
void Matrix<MatrixType>::operator*= (const T number) {
    *this = multiply(number);
}
template <typename MatrixType>
template< typename T >
Matrix<MatrixType> Matrix<MatrixType>::divide(const T number) const {
    return multiply(MatrixType(1/number));
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
template< typename T >
Matrix<MatrixType> Matrix<MatrixType>::power(T number) const {

    number = (int)number;

    if(number == 1) {
        return *this;
    }

    if(!(n == m)) {
        cout << "Can't raise the Matrix to a number. Wrong sizes." << endl;

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

    Farr = multiply(*this);
    for(int i; i<number-2; i++) {
        Farr = Farr.multiply(*this);
    }
    return Farr;
}
template <typename MatrixType>
template< typename T >
Matrix<MatrixType> Matrix<MatrixType>::operator^ (const T number) {
    return power(number);
}

template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::gettranspose() const {
    const Matrix<MatrixType> oldArr = *this;
    Matrix<MatrixType> Farr(m,n,true);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            Farr.ptr[j][i] = oldArr.ptr[i][j];
        }
    }
    return Farr;
}
template <typename MatrixType>
void Matrix<MatrixType>::transpose() {
    *this = gettranspose();
    return;
}

// returns determinant of a matrix
template <typename MatrixType>
MatrixType Matrix<MatrixType>::det() const {
    
    MatrixType det;

    if(!(n == m)) {
        cout << "Can't find supplement the Matrix to a number. Wrong sizes." << endl;

        Matrix<MatrixType> Farr(0,0,true);
        return det;
    }


    if (n == 1) {
        return ptr[0][0];
    }
    else if(n == 2) {
        
        det = ptr[0][0]*ptr[1][1] - ptr[1][0]*ptr[0][1];
        
        return det;
    }

    det = 0;
    for (int i = 0; i < n; i++) {
        det += ptr[i][0] * cofactor(i,0);
    }
    
    return det;

}

// returns trace of a matrix
template <typename MatrixType>
MatrixType Matrix<MatrixType>::tr() const {
    
    MatrixType tr;

    if(!(n == m)) {
        cout << "Can't find supplement the Matrix to a number. Wrong sizes." << endl;

        Matrix<MatrixType> Farr(0,0,true);
        return tr;
    }

    tr = 0;
    for (int i = 0; i < n; i++) {
        tr += ptr[i][i];
    }
    
    return tr;

}

template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::Minor(const int iv,const int jv) const {
    
    Matrix<MatrixType> Minor(n-1,m-1,true);

    int i2 = 0;
    int j2 = 0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(j != jv && i != iv) {
                Minor.ptr[i2][j2] = ptr[i][j];
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

    return Minor;
}

template <typename MatrixType>
MatrixType Matrix<MatrixType>::cofactor(const int iv,const int jv) const {
    
    MatrixType cofactor;

    if(!(n == m)) {
        cout << "Can't find cofactor the Matrix to a number. Wrong sizes." << endl;

        return cofactor;
    }

    cofactor = pow(-1,iv+jv) * Minor(iv,jv).det();
    return cofactor;
}

// returns adjucate of a matrix
template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::adj() const {
    
    if(!(n == m)) {
        cout << "Can't find adjucate of the Matrix. Wrong sizes." << endl;

        Matrix<MatrixType> Farr(0,0,true);
        return Farr;
    }

    Matrix<MatrixType> Farr(n,m,true);

    for(int j = 0; j < m; j++) {
        for(int i = 0; i < n; i++) {
            Farr.ptr[j][i] = cofactor(i,j);
        }
    }
    return Farr;
}

// returns inverse of a matrix
template <typename MatrixType>
Matrix<MatrixType> Matrix<MatrixType>::inv() const {
    
    if(!(n == m)) {
        cout << "Can't find inverse of the Matrix. Wrong sizes." << endl;

        Matrix<MatrixType> Farr(0,0,true);
        return Farr;
    }

    MatrixType det1 = det();
    if(det1 == 0) {
        cout << "Can't find inverse of the Matrix. Determinant is 0." << endl;

        Matrix<MatrixType> Farr(0,0,true);
        return Farr;
    }

    Matrix<MatrixType> Farr(n,m,true);

    Farr = adj().divide(det1);

    return Farr;
}

template <typename MatrixType>
bool Matrix<MatrixType>::equal(Matrix<MatrixType> mat) const {


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
bool Matrix<MatrixType>::operator== (Matrix<MatrixType> mat) const {
    return equal(mat);
}
template <typename MatrixType>
bool Matrix<MatrixType>::operator!= (Matrix<MatrixType> mat) const {
    return !(equal(mat));
}


// elementary row operation

template <typename MatrixType>
void Matrix<MatrixType>::RowSwap(int j1, int j2) {
    if (j1 == j2 ) {
        return; 
    }
    else if ((j1 >= m && j1 < 0 ) || (j2 >= m && j2 < 0 )) {
        cout<< "Row number is not matching the matrix" <<endl;
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



#endif