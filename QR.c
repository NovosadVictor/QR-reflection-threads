#include <malloc.h>
#include <time.h>
#include <pthread.h>

#include "qr_functions.h"
#include "formulas.h"



int main(int argc, char **argv) {
	if (argc != 4) {
		printf("You should choose mode(1 - formula, 2 - file)\n");
        printf("And write matrix size(for mode 1) of filename(for mode 2)\n");
        printf("And also write quantity of threads\n");
		return -1;
	}

    int k = atoi(argv[3]);
    pthread_t *threads = (pthread_t *) malloc(k * sizeof(pthread_t));


    if (atoi(argv[1]) == 2) {
        FILE *fi = fopen(argv[2], "r");
        if (fi == NULL) {
            printf("incorrect filename\n");
            return -1;
        }

        int n;
        if (fscanf(fi, "%d", &n) != 1) {
            printf("Error in matrix size\n");
            return -1;
        }

        m_m *classes = (m_m *) malloc(k * sizeof(m_m));
        double *matrix = (double *) malloc(n * n * sizeof(double));
        double *matrix_for_norma = (double *) malloc(n * n * sizeof(double));
        double *result = (double *) malloc(n * n * sizeof(double));
        double *for_norma = (double *) malloc(n * n * sizeof(double));
        double *d = (double *) malloc(n * sizeof(double));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double elem;
                if (fscanf(fi, "%lf", &elem) != 1) {
                    printf("Incorrect matrix in file\n");
                    return -1;
                }
                matrix[i * n + j] = elem;
                matrix_for_norma[i * n + j] = matrix[i * n + j];
                for_norma[i * n + j] = 0.0;
            }

        char c;
        if (fscanf(fi, "%c", &c) > 0) {
            printf("NOT EOF\n");
            return -1;
        }

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    result[i * n + j] = 1.0;
                else
                    result[i * n + j] = 0.0;
            }

//        printf("Input matrix\n");
//        print(matrix, n);

        if (QR_decomposition(n, matrix, result, d) == 1) {
            printf("\nSingular matrix\n");
            return -1;
        }

        for (int s = 0; s < k; s++) {
            classes[s].n = n;
            classes[s].matrix = matrix_for_norma;
            classes[s].reverse = result;
            classes[s].result = for_norma;
            classes[s].t_quantity = k;
            classes[s].num = s + 1;
        }

        for (int s = 0; s < k; s++)
            if (pthread_create(&threads[s], NULL, norma, (void *)&classes[s]) != 0) {
                printf("\nError, cant create %d thread\n", s);
                return -1;
            }

        for (int s = 0; s < k; s++)
            if (pthread_join(threads[s], NULL) != 0) {
                printf("\nError, cant join %d thread\n", s);
                return -1;
            }


 //       print(result, n);
 //       print(matrix_for_norma, n);
 //       print(for_norma, n);

        time_t start = clock();

        printf("Norma = %lf\n", find_norma(n, for_norma));

        time_t end = clock();
//        printf("Resulting matrix\n");
//        print(result, n);

        printf("execution time: %f sec\n", (double) (end - start) / CLOCKS_PER_SEC);
        free(matrix);
        free(result);
        free(d);
        free(matrix_for_norma);
        free(for_norma);
        free(classes);
        fclose(fi);
    }
    if (atoi(argv[1]) == 1) {
        int n = atoi(argv[2]);

        m_m *classes = (m_m *) malloc(k * sizeof(m_m));
        double *matrix = (double *) malloc(n * n * sizeof(double));
        double *matrix_for_norma = (double *) malloc(n * n * sizeof(double));
        double *result = (double *) malloc(n * n * sizeof(double));
        double *for_norma = (double *) malloc(n * n * sizeof(double));
        double *d = (double *) malloc(n * sizeof(double));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                matrix[i * n + j] = formula(i, j, n);
                matrix_for_norma[i * n + j] = matrix[i * n + j];
                for_norma[i * n + j] = 0.0;
            }


        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    result[i * n + j] = 1.0;
                else
                    result[i * n + j] = 0.0;
            }

//        printf("Input matrix\n");
//        print(matrix, n);

        if (QR_decomposition(n, matrix, result, d) == 1) {
                printf("\nSingular matrix\n");
                return -1;
            }

        for (int s = 0; s < k; s++) {
            classes[s].n = n;
            classes[s].matrix = matrix_for_norma;
            classes[s].reverse = result;
            classes[s].result = for_norma;
            classes[s].t_quantity = k;
            classes[s].num = s + 1;
        }

        for (int s = 0; s < k; s++)
            if (pthread_create(&threads[s], NULL, norma, (void *)&classes[s]) != 0) {
                printf("\nError, cant create %d thread\n", s);
                return -1;
            }

        for (int s = 0; s < k; s++)
            if (pthread_join(threads[s], NULL) != 0) {
                printf("\nError, cant join %d thread\n", s);
                return -1;
            }

 //       print(result, n);
 //       print(matrix_for_norma, n);
 //       print(for_norma, n);

        time_t start = clock();

        printf("Norma = %lf\n", find_norma(n, for_norma));

        time_t end = clock();
//        printf("Resulting matrix\n");
//        print(result, n);

        printf("execution time: %f sec\n", (double) (end - start) / CLOCKS_PER_SEC);
        free(matrix);
        free(result);
        free(d);
        free(matrix_for_norma);
        free(for_norma);
        free(classes);
    }
	return 0;
}

