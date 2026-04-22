/*Решение СЛАУ Ax=b методом Гаусса и его анализ*/

#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>


double* gauss(double **A, int n) {
    // Прямой ход метода Гаусса (приведение к верхнетреугольному виду)
    double m;
    for (int k = 1; k < n; k++) {  // k от 1 до n-1
        for (int j = k; j < n; j++) {  // j от k до n-1
            if (A[k][k] == 0) {
                printf("Ошибка: нулевой элемент на диагонали (k=%d)\n", k);
                break;
            }
            m = A[j][k-1] / A[k-1][k-1];
            for (int i = 0; i < n + 1; i++) {  // i от 0 до n
                A[j][i] -= m * A[k-1][i];
            }
        }
    }
    if (n < 10) {
        printf("Матрица после прямого хода:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n + 1; j++) {
                printf("%.2f ", A[i][j]);
            }
            printf("\n");
        }
    }

    // Обратный ход: вычисление X
    double *X = (double *)malloc(n * sizeof(double));

    for (int i = n - 1; i >= 0; i--) {
        if (A[i][i] == 0) {
            printf("Ошибка: нулевой элемент на диагонали при обратном ходе (i=%d)\n", i);
            break;
        }
        X[i] = A[i][n] / A[i][i];
        for (int c = n - 1; c > i; c--) { 
            X[i] -= A[i][c] * X[c] / A[i][i];
        }
    }

    return X;
}

double* matrix_x_vector_column(double **A, double *X, int n) {
    double *B = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        B[i] = 0.0;
        for (int j = 0; j < n; j++) {
            B[i] += A[i][j] * X[j];
        }
    }
    return B;
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_NUMERIC, "C");

    double **A;
    double *B;
    double *X;
    double *X_target;
    int n;

    int is_use_file = 0; // 0 - не использовать файл, 1 - использовать файл
    printf("Введите 0 для генерации матрицы, 1 для чтения из файла: ");
    scanf("%d", &is_use_file);

    if (is_use_file == 1) { // ИСПОЛЬЗУЕТСЯ ФАЙЛ
        FILE *fp = fopen("matrixA.txt", "r");
        // Чтение размера матрицы
        if (fscanf(fp, "%d", &n) != 1) {
            printf("Ошибка чтения размера матрицы\n");
            fclose(fp);
            return 1;
        }
        // A --- Расширенная матрица коэффициентов (n x (n+1))
        A = (double **)malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++) {
            A[i] = (double *)malloc((n + 1) * sizeof(double));
        }
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

        // Решение 
        printf("\n=================\nМатрица A|b:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n + 1; j++) {
                printf("%.2f ", A[i][j]);
            }
            printf("\n");
        }

        // Вызов функции Гаусса
        clock_t start = clock();
        X = gauss(A, n);
        clock_t end = clock();
        double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        printf("Время выполнения: %.6f секунд\n", time_spent);
        printf("Процессорные такты: %ld\n", (long)(end - start));

        // Вывод решения
        printf("Решение системы:\n");
        for (int i = 0; i < n; i++) {
            printf("X[%d] = %.6f\n", i, X[i]);
        }
        printf("=================\n");

        // Освобождение памяти
        for (int i = 0; i < n; i++) {
            free(A[i]);
        }
        free(A);
        free(X);
    }
    else { // НЕ ИСПОЛЬЗУЕТСЯ ФАЙЛ
        int nn = 100;
        //ввод желаемого размера матрицы через консоль
        printf("Введите размер матрицы n: ");
        scanf("%d", &nn);

        for (n = nn; n < 1500; n += 100)
        {
            // A --- Расширенная матрица коэффициентов (n x (n+1))
            A = (double **)malloc(n * sizeof(double *));
            for (int i = 0; i < n; i++) {
                A[i] = (double *)malloc((n + 1) * sizeof(double));
            }
            //заполнение матрицы A и B и их объединение в РМК
            X_target = (double *)malloc(n * sizeof(double));

            int a = 1;
            for (int i = 0; i < n; i++) {
                X_target[i] = a;

                for (int j = 0; j < n; j++) {
                    if (i == j)
                        A[i][j] = a; 
                    else
                        A[i][j] = 0.1 * a; 
                }
                a +=1;
            }

            B = matrix_x_vector_column(A, X_target, n);
            
            for (int i = 0; i < n; i++) {
                A[i][n] = B[i];
            }

            if (n < 10) {
                // Решение 
                printf("\n=================\nМатрица A|b:\n");
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n + 1; j++) {
                        printf("%.2f ", A[i][j]);
                    }
                    printf("\n");
                }
            }

            // Вызов функции Гаусса
            clock_t start = clock();
            X = gauss(A, n);
            clock_t end = clock();
            double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
            printf("Время выполнения: %.6f секунд\n", time_spent);
            printf("Процессорные такты: %ld\n", (long)(end - start));


            // Вывод решения
            printf("Решение системы:\n");
            if (n < 10)
                for (int i = 0; i < n; i++)
                    printf("X[%d] = %.6f\n", i, X[i]);
            //else for (int i = 0; i < n; i++) printf("%.2f ", X[i]);
            
            printf("=================\n");    

            // Освобождение памяти
            for (int i = 0; i < n; i++) {
                free(A[i]);
            }
            free(A);
            free(X);
            free(B);
            free(X_target);
        }     
    }

    return 0;
}
