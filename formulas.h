#ifndef QR_FORMULAS_H
#define QR_FORMULAS_H

#include <stdio.h>
#include <stdlib.h>


double formula(int i, int j, int size) {
    if (i == 0 && j == 0)
        return -1;
    if (i == size - 1 && j == size - 1)
        return -(size - 1.0) / size;
    if (i == j)
        return -2.0;
    if (i == j - 1)
        return 1.0;
    if (i == j + 1)
        return 1.0;
    return 0;
}


void print(double *matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            printf("%lf ", matrix[i * size + j]);
        printf("\n");
    }
    printf("\n");
}

#endif //QR_FORMULAS_H
