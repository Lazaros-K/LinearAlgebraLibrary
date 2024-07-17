#include <iostream>
#include <thread>

#include "LinearAlgebraLibrary.hpp"
#include "tools/Timer.hpp"

using namespace std;


int main() {
    cout << "\n\n\n\n\n\n\n\n\n\n";
    
    double x[4][6] = 
        {
            {0,2,-1,1,0,1},
            {0,8,9,1,1,1},
            {1,-1,2,2,0,1},
            {2,3,-1,1,1,1}
        };

    double y[2][2] = 
        {
            {-2,2},
            {1,1},
        };

    double z[2][2] = 
        {
            {1,6},
            {-2,2}
        };
    
    Matrix<double> arr1(x,true);
    Matrix<double> arr2(y,true);
    Matrix<double> arr3(z,true);

    /*
    arr1.RowSum(0,1,-2);

    arr1.RowSum(0,2,-1);

    arr1.RowSum(1,2,-(6.0/15.0));
    

    arr1.RowScalarMult(2,1.0/0.6);
    arr1.RowSum(2,1,-1);

    arr1.RowScalarMult(1,-1.0/15);
    arr1.RowSum(1,0,-8);
    */

    Timer timer;
    arr1.Gauss_Jordan();
    timer.stopPrint();

    printMatrix(arr1,2);

    

    //arr2.ChangeSize(2,2);
    //arr1.Write(x);

    // arr2.reverse();
    //printMatrix(arr1^2);
    // cout << arr1.tr() <<endl;

    //this_thread::sleep_for(chrono::seconds(10));








    /*
    ---- Theory-3 ----

    bool ans = (arr3.KroneckerProduct(arr1.directSum(arr2)) == arr3.KroneckerProduct(arr1).directSum(arr3.KroneckerProduct(arr2)));
    bool ans2 = ((arr1.directSum(arr2)).KroneckerProduct(arr3) == arr1.KroneckerProduct(arr3).directSum(arr2.KroneckerProduct(arr3)));
    cout << ans2 << endl;
    cout << ans1 << endl;
    printMatrix((arr1.directSum(arr2)).KroneckerProduct(arr3));
    printMatrix(arr1.KroneckerProduct(arr3).directSum(arr2.KroneckerProduct(arr3)));
    */
   system("pause");
    
    return 0;
}