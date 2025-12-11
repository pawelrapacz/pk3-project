#include "Matrix.h"

int main() {
    CRSMatrix<double, 4> matrix1({
        { 0, 0, 0, 0 },
        { 5, 8, 0, 0 },
        { 0, 0, 3, 0 },
        { 0, 6, 0, 0 }
    });
    matrix1.print();

    CRSMatrix<double, 4, 6> matrix2({
        { 10, 20,  0,  0,  0,  0 },
        {  0, 30,  0,  4,  0,  0 },
        {  0,  0, 50, 60, 70,  0 },
        {  0,  0,  0,  0,  0, 80 }
    });
    matrix2.print();
    return 0;
}