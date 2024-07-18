#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include "Matrix.hpp"
#include "functions_matrix.hpp"

namespace lalib { // starting namespace lalib


template <typename MatrixType>
MatrixType Minor(const Matrix<MatrixType> mat,const int iv,const int jv) {

    if(!(mat.GetN() == mat.GetM())) {
        std::cout << "Can't find Minor the Matrix to a number. Wrong sizes." << std::endl;

        MatrixType Minor;
        return Minor;
    }

    return det(mat.Submatrix(iv,jv));
}

template <typename MatrixType>
MatrixType cofactor(const Matrix<MatrixType> mat,const int iv,const int jv) {

    if(!(mat.GetN() == mat.GetM())) {
        std::cout << "Can't find cofactor the Matrix to a number. Wrong sizes." << std::endl;

        return 0;
    }

    MatrixType cofactor = Minor(mat,iv,jv);

    if(cofactor == 0) {
        return cofactor;
    }

    return cofactor * (pow(-1,iv+jv) );
}

template <typename MatrixType>
Matrix<MatrixType> adj(const Matrix<MatrixType> mat) {
    
    if(!(mat.GetN() == mat.GetM())) {
        std::cout << "Can't find adjucate of the Matrix. Wrong sizes." << std::endl;

        Matrix<MatrixType> Farr(0,0,true);
        return Farr;
    }

    Matrix<MatrixType> Farr(mat.GetN(),mat.GetM(),true);

    for(int j = 0; j < mat.GetM(); j++) {
        for(int i = 0; i < mat.GetN(); i++) {
            Farr.ptr[j][i] = cofactor(mat,i,j);
        }
    }
    return Farr;
}

// returns inverse of a matrix
template <typename MatrixType>
Matrix<MatrixType> inv(const Matrix<MatrixType> mat) {
    
    if(!(mat.GetN() == mat.GetM())) {
        std::cout << "Can't find inverse of the Matrix. Wrong sizes." << std::endl;

        Matrix<MatrixType> Farr;
        return Farr;
    }

    MatrixType det1 = det(mat);
    if(det1 == 0) {
        std::cout << "Can't find inverse of the Matrix. Determinant is 0." << std::endl;

        Matrix<MatrixType> Farr;
        return Farr;
    }

    Matrix<MatrixType> Farr(mat.GetN(),mat.GetM(),true);

    Farr = adj(mat).divide(det1);

    return Farr;
}

template <typename MatrixType>
Matrix<MatrixType> transpose(const Matrix<MatrixType> mat) {
    const Matrix<MatrixType> oldArr = mat;
    Matrix<MatrixType> Farr(mat.GetM(),mat.GetN(),true);

    for(int i = 0; i < mat.GetN(); i++) {
        for(int j = 0; j < mat.GetM(); j++) {
            Farr.ptr[j][i] = oldArr.ptr[i][j];
        }
    }
    return Farr;
}

// returns determinant of a matrix
template <typename MatrixType>
MatrixType det(const Matrix<MatrixType> mat) {
    
    MatrixType det;

    if(!(mat.GetN() == mat.GetM())) {
        std::cout << "Can't find supplement the Matrix to a number. Wrong sizes." << std::endl;

        Matrix<MatrixType> Farr;
        return det;
    }


    if (mat.GetN() == 1) {
        return mat.ptr[0][0];
    }
    else if(mat.GetN() == 2) {
        
        det = mat.ptr[0][0]*mat.ptr[1][1] - mat.ptr[1][0]*mat.ptr[0][1];
        
        return det;
    }

    det = 0;
    for (int i = 0; i < mat.GetN(); i++) {
        det += mat.ptr[i][0] * cofactor(mat,i,0);
    }
    
    return det;

}

template <typename MatrixType>
MatrixType tr(const Matrix<MatrixType> mat) {
    
    MatrixType tr;

    if(!(mat.GetN() == mat.GetM())) {
        std::cout << "Can't find supplement the Matrix to a number. Wrong sizes." << std::endl;

        Matrix<MatrixType> Farr(0,0,true);
        return tr;
    }

    tr = 0;
    for (int i = 0; i < mat.GetN(); i++) {
        tr += mat.ptr[i][i];
    }
    
    return tr;

}


template<typename MatrixType , typename T >
Matrix<MatrixType> operator* (const T number,Matrix<MatrixType> mat) {
    return mat * number;
}

template <typename MatrixType>
void printMatrix(Matrix<MatrixType> arr,int precisiond) {

    using namespace std;

    if(arr.GetN() == 0 && arr.GetM()== 0) {
        std::cout << "||" << std::endl;
    }

    int maxlength = 0;
    int maxd = 0;
    int countd[arr.GetM()*arr.GetN()] = {0};
    int decimal_limit = numeric_limits<MatrixType>::max_digits10;

    if (decimal_limit > precisiond) {
        decimal_limit = precisiond;
    }

    for (int i = 0; i < arr.GetM(); i++)
    {
        for (int j = 0; j < arr.GetN(); j++)
        {
            
            // get desimal size
            MatrixType num = abs(arr.ptr[j][i]);

            int exponent = 1.0;
            MatrixType epsilon = numeric_limits<MatrixType>::epsilon() * num;
            MatrixType c = num - floor(num);

            

            while((c > epsilon && c < (1 - epsilon)) && countd[j+i*arr.GetN()] < decimal_limit)
            {
                exponent *= 10;

                c = num * exponent;
                c = c - floor(c);

                epsilon = numeric_limits<MatrixType>::epsilon() * exponent * num;

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
    // std::cout << maxlength << std::endl;
    // std::cout << maxd << std::endl;


    string pavla;
    for (int i = 0; i < ((arr.GetN())*(maxlength)+(arr.GetN()-1)*3); i++)
    {
        pavla.append("-");
    }
    

    std::cout << "--" << pavla << "--" << std::endl;
    for (int i = 0; i < arr.GetM(); i++)
    {
        for (int j = 0; j < arr.GetN(); j++)
        {
            std::cout << "| " ;
            std::cout << setw(maxlength);
            std::cout << setprecision(countd[j+i*arr.GetN()]);
            std::cout << fixed;
            std::cout << MatrixType(arr.ptr[i][j]);
            std::cout << " ";
        }
        std::cout << "|" << std::endl;
        std::cout << "--" << pavla << "--" << std::endl;
    }
    std::cout << std::endl;
    return;
}

template double Minor<double>(Matrix<double>,int,int);
template double cofactor<double>(Matrix<double>,int,int);
template Matrix<double> adj<double>(Matrix<double>);
template Matrix<double> inv<double>(Matrix<double>);
template Matrix<double> transpose<double>(Matrix<double>);
template double det<double>(Matrix<double>);
template double tr<double>(Matrix<double>);
template void printMatrix<double>(Matrix<double>,int);

template int Minor<int>(Matrix<int>,int,int);
template int cofactor<int>(Matrix<int>,int,int);
template Matrix<int> adj<int>(Matrix<int>);
template Matrix<int> inv<int>(Matrix<int>);
template Matrix<int> transpose<int>(Matrix<int>);
template int det<int>(Matrix<int>);
template int tr<int>(Matrix<int>);
template void printMatrix<int>(Matrix<int>,int);

template float Minor<float>(Matrix<float>,int,int);
template float cofactor<float>(Matrix<float>,int,int);
template Matrix<float> adj<float>(Matrix<float>);
template Matrix<float> inv<float>(Matrix<float>);
template Matrix<float> transpose<float>(Matrix<float>);
template float det<float>(Matrix<float>);
template float tr<float>(Matrix<float>);
template void printMatrix<float>(Matrix<float>,int);



template Matrix<double> operator*<double,double>(double,Matrix<double>);
template Matrix<int> operator*<int,double>(double,Matrix<int>);
template Matrix<float> operator*<float,double>(double,Matrix<float>);


template Matrix<double> operator*<double,int>(int,Matrix<double>);
template Matrix<int> operator*<int,int>(int,Matrix<int>);
template Matrix<float> operator*<float,int>(int,Matrix<float>);


template Matrix<double> operator*<double,float>(float,Matrix<double>);
template Matrix<int> operator*<int,float>(float,Matrix<int>);
template Matrix<float> operator*<float,float>(float,Matrix<float>);



} // ending namespace lalib