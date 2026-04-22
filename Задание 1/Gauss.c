/*Решение СЛАУ Ax=b методом Гаусса и его анализ*/

#include <stdio.h>
#include <math.h>
#include <locale.h>
#include <stdlib.h>
#include <windows.h>

int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_NUMERIC, "C");

    int n;
    
    FILE *fp = fopen("matrixA.txt", "r");
    
    // Чтение размера матрицы
    if (fscanf(fp, "%d", &n) != 1) {
        printf("Ошибка чтения размера матрицы\n");
        fclose(fp);
        return 1;
    }

    // A --- Расширенная матрица коэффициентов (n x (n+1))
    double **A = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *)malloc((n + 1) * sizeof(double));
    }

    // X - массив для решения (вектор неизвестных)
    double *X = (double *)malloc(n * sizeof(double));

    // Чтение расширенной матрицы A|b из файла
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            if (fscanf(fp, "%lf", &A[i][j]) != 1) {
                printf("Ошибка чтения элемента матрицы A[%d][%d]\n", i, j);
                fclose(fp);
                return 1;
            }
        }
    }
    
    // Закрытие файла
    fclose(fp);

    printf("Матрица A|b:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            printf("%.2f ", A[i][j]);
        }
        printf("\n");
    }

    // Прямой ход метода Гаусса (приведение к верхнетреугольному виду)
    double m;
    for (int k = 1; k < n; k++) {  // k от 1 до n-1
        for (int j = k; j < n; j++) {  // j от k до n-1
            if (A[k][k] == 0) {
                printf("Ошибка: нулевой элемент на диагонали (k=%d)\n", k);
                return 1;
            }
            m = A[j][k-1] / A[k-1][k-1];
            for (int i = 0; i < n + 1; i++) {  // i от 0 до n
                A[j][i] -= m * A[k-1][i];
            }
        }
    }

    printf("Матрица после прямого хода:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            printf("%.2f ", A[i][j]);
        }
        printf("\n");
    }

    // Обратный ход: вычисление X
    for (int i = n - 1; i >= 0; i--) {
        if (A[i][i] == 0) {
            printf("Ошибка: нулевой элемент на диагонали при обратном ходе (i=%d)\n", i);
            return 1;
        }
        X[i] = A[i][n] / A[i][i];
        for (int c = n - 1; c > i; c--) { 
            X[i] -= A[i][c] * X[c] / A[i][i];
        }
    }

    // Вывод решения
    printf("Решение системы:\n");
    for (int i = 0; i < n; i++) {
        printf("X[%d] = %.6f\n", i, X[i]);
    }

    // Освобождение памяти
    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(X);

    return 0;
}




