#include <stdio.h>
#include <math.h>
#include <locale.h>
#include <windows.h>
#include <stdlib.h>
#include <string.h>

// ======================== ОСНОВНЫЕ ФУНКЦИИ МАТРИЧНЫХ ОПЕРАЦИЙ ========================

// Выделение памяти под квадратную матрицу
double** alloc_matrix(int n) {
    double** A = (double**)malloc(n * sizeof(double*));
    if (A == NULL) return NULL;
    for (int i = 0; i < n; i++) {
        A[i] = (double*)malloc(n * sizeof(double));
        if (A[i] == NULL) {
            for (int j = 0; j < i; j++) free(A[j]);
            free(A);
            return NULL;
        }
    }
    return A;
}

// Освобождение памяти матрицы
void free_matrix(double** A, int n) {
    if (A == NULL) return;
    for (int i = 0; i < n; i++) free(A[i]);
    free(A);
}

// Копирование матрицы
void copy_matrix(double** dest, double** src, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

// Умножение матриц: C = A * B
void matrix_multiply(double** A, double** B, double** C, int n) {
    double** temp = alloc_matrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            temp[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                temp[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    copy_matrix(C, temp, n);
    free_matrix(temp, n);
}

// Норма для проверки сходимости (максимум внедиагональных элементов)
double off_diagonal_norm(double** A, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) sum += A[i][j] * A[i][j];
        }
    }
    return sqrt(sum);
}

// Печать матрицы
void print_matrix(double** A, int n, const char* title) {
    printf("\n%s:\n", title);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%12.6f ", A[i][j]);
        }
        printf("\n");
    }
}

// Печать вектора
void print_vector(double* v, int n, const char* title) {
    printf("\n%s:\n", title);
    for (int i = 0; i < n; i++) {
        printf("λ%d = %12.6f\n", i + 1, v[i]);
    }
}

// ===== ПРЕОБРАЗОВАНИЯ ХАУСХОЛДЕРА 

// Приведение матрицы к верхней форме Хессенберга
// Используются преобразования Хаусхолдера
void hessenberg_reduction(double** A, int n) {
    for (int k = 0; k < n - 2; k++) {
        // Находим вектор Хаусхолдера для поддиагональной части столбца k
        double norm = 0.0;
        for (int i = k + 1; i < n; i++) {
            norm += A[i][k] * A[i][k];
        }
        norm = sqrt(norm);
        if (norm == 0.0) continue;
        
        // Вычисляем alpha = -sign(A[k+1][k]) * norm
        double alpha = (A[k + 1][k] > 0 ? -norm : norm);
        double beta = alpha * alpha - alpha * A[k + 1][k];
        beta = 1.0 / beta;
        
        // Формируем вектор u
        double* u = (double*)malloc(n * sizeof(double));
        memset(u, 0, n * sizeof(double));
        u[k + 1] = A[k + 1][k] - alpha;
        for (int i = k + 2; i < n; i++) {
            u[i] = A[i][k];
        }
        
        // Применяем преобразование слева: A = (I - beta*u*u^T) * A
        // Сначала вычисляем w = beta * A^T * u
        double* w = (double*)malloc(n * sizeof(double));
        memset(w, 0, n * sizeof(double));
        for (int j = k; j < n; j++) {
            double sum = 0.0;
            for (int i = k + 1; i < n; i++) {
                sum += A[i][j] * u[i];
            }
            w[j] = beta * sum;
        }
        
        // Обновляем A: A = A - u * w^T
        for (int i = k + 1; i < n; i++) {
            for (int j = k; j < n; j++) {
                A[i][j] -= u[i] * w[j];
            }
        }
        
        // Применяем преобразование справа: A = A * (I - beta*u*u^T)
        // Вычисляем w = beta * A * u
        memset(w, 0, n * sizeof(double));
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = k + 1; j < n; j++) {
                sum += A[i][j] * u[j];
            }
            w[i] = beta * sum;
        }
        
        // Обновляем A: A = A - w * u^T
        for (int i = 0; i < n; i++) {
            for (int j = k + 1; j < n; j++) {
                A[i][j] -= w[i] * u[j];
            }
        }
        
        free(u);
        free(w);
    }
}

// ======================== QR-АЛГОРИТМ ========================

// QR-разложение с помощью вращений Гивенса
// Матрица A на входе (разрушается), на выходе Q и R
void qr_decomposition(double** A, double** Q, double** R, int n) {
    // Копируем A в R
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            R[i][j] = A[i][j];
        }
    }
    // Q = единичная матрица
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Q[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // Обнуляем поддиагональные элементы с помощью вращений Гивенса
    for (int j = 0; j < n - 1; j++) {
        for (int i = n - 1; i > j; i--) {
            if (R[i][j] != 0.0) {
                double a = R[i - 1][j];
                double b = R[i][j];
                double r = sqrt(a * a + b * b);
                double c = a / r;
                double s = -b / r;
                
                // Применяем вращение к строкам R
                for (int k = j; k < n; k++) {
                    double temp = c * R[i - 1][k] - s * R[i][k];
                    R[i][k] = s * R[i - 1][k] + c * R[i][k];
                    R[i - 1][k] = temp;
                }
                // Применяем вращение к столбцам Q
                for (int k = 0; k < n; k++) {
                    double temp = c * Q[k][i - 1] - s * Q[k][i];
                    Q[k][i] = s * Q[k][i - 1] + c * Q[k][i];
                    Q[k][i - 1] = temp;
                }
            }
        }
    }
}

// QR-алгоритм для нахождения собственных значений и векторов
// Возвращает 1 в случае успеха, 0 при превышении максимального числа итераций
int qr_algorithm(double** A, double* eigenvalues, double** eigenvectors, int n, double eps, int max_iter) {
    // Создаем рабочие матрицы
    double** H = alloc_matrix(n);
    double** Q = alloc_matrix(n);
    double** R = alloc_matrix(n);
    double** new_H = alloc_matrix(n);
    double** V = alloc_matrix(n);  // накопление собственных векторов
    
    if (!H || !Q || !R || !new_H || !V) {
        printf("Ошибка выделения памяти.\n");
        return 0;
    }
    
    // Копируем исходную матрицу
    copy_matrix(H, A, n);
    
    // 1. Приводим к форме Хессенберга для ускорения сходимости
    hessenberg_reduction(H, n);
    
    // 2. Инициализируем матрицу накопления векторов (единичная)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            V[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    int iter = 0;
    double norm = off_diagonal_norm(H, n);
    
    // 3. Основной итерационный процесс
    while (norm > eps && iter < max_iter) {
        // QR-разложение текущей матрицы H
        qr_decomposition(H, Q, R, n);
        
        // H_new = R * Q
        matrix_multiply(R, Q, new_H, n);
        
        // Обновляем накопленные векторы: V = V * Q
        double** temp = alloc_matrix(n);
        matrix_multiply(V, Q, temp, n);
        copy_matrix(V, temp, n);
        free_matrix(temp, n);
        
        // Копируем new_H обратно в H
        copy_matrix(H, new_H, n);
        
        norm = off_diagonal_norm(H, n);
        iter++;
        
        // Для отладки (можно закомментировать)
        // printf("Итерация %d: норма внедиагональных элементов = %e\n", iter, norm);
    }
    
    if (iter >= max_iter) {
        printf("Предупреждение: достигнуто максимальное число итераций (%d).\n", max_iter);
        printf("Текущая норма внедиагональных элементов: %e\n", norm);
        free_matrix(H, n);
        free_matrix(Q, n);
        free_matrix(R, n);
        free_matrix(new_H, n);
        free_matrix(V, n);
        return 0;
    }
    
    // 4. Извлекаем собственные значения с диагонали
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = H[i][i];
    }
    
    // 5. Сортируем собственные значения и соответствующие векторы по убыванию
    // Простой пузырьковый метод для демонстрации
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (eigenvalues[i] < eigenvalues[j]) {
                double temp = eigenvalues[i];
                eigenvalues[i] = eigenvalues[j];
                eigenvalues[j] = temp;
                
                for (int k = 0; k < n; k++) {
                    double temp_vec = V[k][i];
                    V[k][i] = V[k][j];
                    V[k][j] = temp_vec;
                }
            }
        }
    }
    
    // 6. Нормируем собственные векторы
    for (int j = 0; j < n; j++) {
        double norm_vec = 0.0;
        for (int i = 0; i < n; i++) {
            norm_vec += V[i][j] * V[i][j];
        }
        norm_vec = sqrt(norm_vec);
        if (norm_vec > 1e-12) {
            for (int i = 0; i < n; i++) {
                V[i][j] /= norm_vec;
            }
        }
    }
    
    // Копируем результат в выходную матрицу eigenvectors
    copy_matrix(eigenvectors, V, n);
    
    free_matrix(H, n);
    free_matrix(Q, n);
    free_matrix(R, n);
    free_matrix(new_H, n);
    free_matrix(V, n);
    return 1;
}

// ======================== ВВОД МАТРИЦЫ ========================

// Выделение памяти под квадратную матрицу с проверкой на пустой файл
double** matrix_from_file(int* n) {
    const char* filename = "matrix.txt";
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Ошибка: не удалось открыть файл '%s'\n", filename);
        return NULL;
    }
    
    if (fscanf(fp, "%d", n) != 1 || *n <= 0) {
        printf("Ошибка: в файле '%s' не указан корректный размер матрицы.\n", filename);
        fclose(fp);
        return NULL;
    }
    
    double** A = alloc_matrix(*n);
    if (A == NULL) {
        printf("Ошибка выделения памяти для матрицы размером %d.\n", *n);
        fclose(fp);
        return NULL;
    }
    
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
            if (fscanf(fp, "%lf", &A[i][j]) != 1) {
                printf("Ошибка: недостаточно данных в файле для матрицы %dx%d.\n", *n, *n);
                fclose(fp);
                free_matrix(A, *n);
                return NULL;
            }
        }
    }
    fclose(fp);
    return A;
}

// Ввод матрицы с клавиатуры
double** matrix_from_keyboard(int* n) {
    printf("Введите размерность квадратной матрицы n: ");
    scanf("%d", n);
    if (*n <= 0) {
        printf("Размерность должна быть положительным числом.\n");
        return NULL;
    }
    
    double** A = alloc_matrix(*n);
    if (A == NULL) {
        printf("Ошибка выделения памяти для матрицы размером %d.\n", *n);
        return NULL;
    }
    
    printf("Введите элементы матрицы построчно:\n");
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
            scanf("%lf", &A[i][j]);
        }
    }
    return A;
}

// ======================== ОСНОВНАЯ ФУНКЦИЯ ========================

int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_NUMERIC, "C");

    int n, choice, input_choice;
    double eps;
    int max_iter;
    double** A = NULL;
    double* eigenvalues = NULL;
    double** eigenvectors = NULL;
    
    printf("===== QR-АЛГОРИТМ ДЛЯ НАХОЖДЕНИЯ СОБСТВЕННЫХ ЗНАЧЕНИЙ И ВЕКТОРОВ =====\n\n");
    
    // Ввод матрицы
    printf("Выберите способ ввода матрицы:\n");
    printf("1 - из файла 'matrix.txt'\n");
    printf("2 - ввод с клавиатуры\n");
    printf("Ваш выбор: ");
    scanf("%d", &input_choice);
    
    if (input_choice == 1) {
        A = matrix_from_file(&n);
        if (A == NULL) return 1;
    } else if (input_choice == 2) {
        A = matrix_from_keyboard(&n);
        if (A == NULL) return 1;
    } else {
        printf("Неверный выбор.\n");
        return 1;
    }
    
    printf("\nИсходная матрица (%dx%d):\n", n, n);
    print_matrix(A, n, "Исходная матрица");
    
    // Ввод параметров точности
    printf("\nВведите точность (eps > 0): ");
    scanf("%lf", &eps);
    printf("Введите максимальное число итераций: ");
    scanf("%d", &max_iter);
    
    // Выделение памяти под результаты
    eigenvalues = (double*)malloc(n * sizeof(double));
    eigenvectors = alloc_matrix(n);
    if (eigenvalues == NULL || eigenvectors == NULL) {
        printf("Ошибка выделения памяти.\n");
        return 1;
    }
    
    // Выполнение QR-алгоритма
    int success = qr_algorithm(A, eigenvalues, eigenvectors, n, eps, max_iter);
    
    if (success) {
        printf("\n===== РЕЗУЛЬТАТЫ =====\n");
        print_vector(eigenvalues, n, "Собственные значения (в порядке убывания)");
        print_matrix(eigenvectors, n, "Собственные векторы (по столбцам)");
    } else {
        printf("\nQR-алгоритм не сошёлся за заданное число итераций.\n");
    }
    
    // Освобождение памяти
    free_matrix(A, n);
    free_matrix(eigenvectors, n);
    free(eigenvalues);
    
    return 0;
}