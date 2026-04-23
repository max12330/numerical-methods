/*Решение СЛАУ Ax=b методом Гаусса и его анализ*/

#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <windows.h>
#include <math.h>
#include <time.h>


double* gauss(double **A, int n) {
    // Прямой ход с выбором главного элемента по столбцу
    for (int k = 0; k < n - 1; k++) {
        // 1. Поиск главного элемента в столбце k
        int max_row = k;
        double max_val = fabs(A[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > max_val) {
                max_val = fabs(A[i][k]);
                max_row = i;
            }
        }
        // Проверка на ноль
        if (max_val < 1e-12) {
            printf("Ошибка: столбец %d почти нулевой (главный элемент = %g)\n", k, max_val);
            return NULL;
        }
        // 2. Перестановка строк
        if (max_row != k) {
            double *tmp = A[k];
            A[k] = A[max_row];
            A[max_row] = tmp;
        }
        // 3. Исключение под диагональю
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n + 1; j++) {
                A[i][j] -= factor * A[k][j];
            }
        }
    }

    // Вывод матрицы после прямого хода
    if (n < 10) {
        printf("Матрица после прямого хода (верхнетреугольная):\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n + 1; j++) {
                printf("%.3f ", A[i][j]);
            }
            printf("\n");
        }
    }

    // Обратный ход с проверкой диагонали
    double *X = (double *)malloc(n * sizeof(double));
    for (int i = n - 1; i >= 0; i--) {
        if (fabs(A[i][i]) < 1e-12) {
            printf("Ошибка: нулевой диагональный элемент на шаге %d\n", i);
            free(X);
            return NULL;
        }
        X[i] = A[i][n];
        for (int j = i + 1; j < n; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
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
    double **A_for_check;
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
        //копируем матрицу A для проверки
        A_for_check = (double **)malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++) {
            A_for_check[i] = (double *)malloc((n + 1) * sizeof(double));
            for (int j = 0; j < n + 1; j++) {
                A_for_check[i][j] = A[i][j];
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

        double *check = (double *)malloc(n * sizeof(double));
        check = matrix_x_vector_column(A_for_check, X, n);
        for (int i = 0; i < n; i++) {
            check[i] -= A_for_check[i][n]; // вычитаем b[i]
        }
        for (int i = 0; i < n; i++) {
            printf("check[%d] = %.6f\n", i, check[i]);
        }

        // Освобождение памяти
        for (int i = 0; i < n; i++) {
            free(A[i]);
            free(A_for_check[i]);
        }
        free(A);
        free(X);
    }
    else { // НЕ ИСПОЛЬЗУЕТСЯ ФАЙЛ
        int nn = 100;
        //ввод желаемого размера матрицы через консоль
        printf("Введите размер матрицы n: ");
        scanf("%d", &nn);

        for (n = nn; n < 100; n += 100)
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
