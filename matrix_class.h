#ifndef MATRIX_CLASS_H
#define MATRIX_CLASS_H


typedef struct matrix_multiplication {
    int num; // thread index
    int t_quantity; // quantity of threads
    int n; // matrix size
    double *matrix;
    double *reverse;
    double *result;
} m_m;


#endif //MATRIX_CLASS_H
