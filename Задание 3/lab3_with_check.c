#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <windows.h>

//  ОСНОВНЫЕ ФУНКЦИИ 

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

void free_matrix(double** A, int n) {
    if (A == NULL) return;
    for (int i = 0; i < n; i++) free(A[i]);
    free(A);
}

void copy_matrix(double** dest, double** src, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            dest[i][j] = src[i][j];
}

void matrix_multiply(double** A, double** B, 
                double** C, int n) {
    double** temp = alloc_matrix(n);
    if (temp == NULL) return;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            temp[i][j] = 0.0;
            for (int k = 0; k < n; k++)
                temp[i][j] += A[i][k] * B[k][j];
        }
    }
    copy_matrix(C, temp, n);
    free_matrix(temp, n);
}

double off_diagonal_norm(double** A, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j) sum += A[i][j] * A[i][j];
    return sqrt(sum);
}

void print_matrix(double** A, int n, const char* title) {
    printf("\n%s:\n", title);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%12.6f ", A[i][j]);
        printf("\n");
    }
}

void print_eigenvalues(double* eigenvals, int n) {
    printf("\nСобственные значения (в порядке убывания):\n");
    for (int i = 0; i < n; i++)
        printf("lambda%d = %12.6f\n", i+1, eigenvals[i]);
}

void print_eigenvectors(double** V, int n) {
    printf("\nСобственные векторы (по столбцам):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%12.6f ", V[i][j]);
        printf("\n");
    }
}

// =============== QR-РАЗЛОЖЕНИЕ (ВРАЩЕНИЯ ГИВЕНСА) 
void qr_decomposition(double** A, double** Q, 
                double** R, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            R[i][j] = A[i][j];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Q[i][j] = (i == j) ? 1.0 : 0.0;
    
    for (int j = 0; j < n-1; j++) {
        for (int i = n-1; i > j; i--) {
            if (fabs(R[i][j]) > 1e-14) {
                double a = R[i-1][j];
                double b = R[i][j];
                double r = sqrt(a*a + b*b);
                double c = a / r;
                double s = -b / r;
                for (int k = j; k < n; k++) {
                    double temp = c * R[i-1][k] - s * R[i][k];
                    R[i][k] = s * R[i-1][k] + c * R[i][k];
                    R[i-1][k] = temp;
                }
                for (int k = 0; k < n; k++) {
                    double temp = c * Q[k][i-1] - s * Q[k][i];
                    Q[k][i] = s * Q[k][i-1] + c * Q[k][i];
                    Q[k][i-1] = temp;
                }
            }
        }
    }
}

// ============ QR-АЛГОРИТМ 
int qr_algorithm(double** A, double* eigenvalues, 
        double** eigenvectors, int n, double eps, 
        int max_iter, int* iterations_used) 
    {
    double** H = alloc_matrix(n);
    double** Q = alloc_matrix(n);
    double** R = alloc_matrix(n);
    double** newH = alloc_matrix(n);
    double** V = alloc_matrix(n);
    if (!H || !Q || !R || !newH || !V) {
        printf("Ошибка выделения памяти.\n");
        return 0;
    }
    
    copy_matrix(H, A, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            V[i][j] = (i == j) ? 1.0 : 0.0;
    
    int iter = 0;
    double norm = off_diagonal_norm(H, n);
    
    while (norm > eps && iter < max_iter) {
        qr_decomposition(H, Q, R, n);
        matrix_multiply(R, Q, newH, n);
        double** temp = alloc_matrix(n);
        matrix_multiply(V, Q, temp, n);
        copy_matrix(V, temp, n);
        free_matrix(temp, n);
        copy_matrix(H, newH, n);
        norm = off_diagonal_norm(H, n);
        iter++;
    }
    
    *iterations_used = iter;
    
    if (iter >= max_iter) {
        printf("Внимание: достигнуто максимальное число итераций (%d)\n", 
            max_iter);
        printf("Текущая норма внедиагональных элементов: %e\n", norm);
        free_matrix(H, n); free_matrix(Q, n); free_matrix(R, n);
        free_matrix(newH, n); free_matrix(V, n);
        return 0;
    }
    
    for (int i = 0; i < n; i++)
        eigenvalues[i] = H[i][i];
    
    // Сортировка по убыванию
    for (int i = 0; i < n-1; i++) {
        for (int j = i+1; j < n; j++) {
            if (eigenvalues[i] < eigenvalues[j]) {
                double tmp = eigenvalues[i];
                eigenvalues[i] = eigenvalues[j];
                eigenvalues[j] = tmp;
                for (int k = 0; k < n; k++) {
                    double tmp_vec = V[k][i];
                    V[k][i] = V[k][j];
                    V[k][j] = tmp_vec;
                }
            }
        }
    }
    
    // Нормировка векторов
    for (int j = 0; j < n; j++) {
        double norm_vec = 0.0;
        for (int i = 0; i < n; i++)
            norm_vec += V[i][j] * V[i][j];
        norm_vec = sqrt(norm_vec);
        if (norm_vec > 1e-12)
            for (int i = 0; i < n; i++)
                V[i][j] /= norm_vec;
    }
    
    copy_matrix(eigenvectors, V, n);
    
    free_matrix(H, n); free_matrix(Q, n); free_matrix(R, n);
    free_matrix(newH, n); free_matrix(V, n);
    return 1;
}

// ======== ВВОД МАТРИЦЫ 
double** read_matrix_from_file(int* n) {
    const char* fname = "matrix.txt";
    FILE* fp = fopen(fname, "r");
    if (!fp) {
        printf("Не удалось открыть файл %s\n", fname);
        return NULL;
    }
    if (fscanf(fp, "%d", n) != 1 || *n <= 0) {
        printf("Ошибка чтения размера матрицы.\n");
        fclose(fp);
        return NULL;
    }
    double** A = alloc_matrix(*n);
    if (!A) return NULL;
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
            if (fscanf(fp, "%lf", &A[i][j]) != 1) {
                printf("Ошибка чтения элемента матрицы.\n");
                fclose(fp);
                free_matrix(A, *n);
                return NULL;
            }
        }
    }
    fclose(fp);
    return A;
}

double** read_matrix_from_keyboard(int* n) {
    printf("Введите размерность квадратной матрицы n: ");
    scanf("%d", n);
    if (*n <= 0) {
        printf("Размерность должна быть положительной.\n");
        return NULL;
    }
    double** A = alloc_matrix(*n);
    if (!A) return NULL;
    printf("Введите элементы матрицы построчно:\n");
    for (int i = 0; i < *n; i++)
        for (int j = 0; j < *n; j++)
            scanf("%lf", &A[i][j]);
    return A;
}

// ======== ПРОВЕРКА A*v - lambda*v = 0 
void check_residuals(double** A, double* eigenvalues, 
    double** eigenvectors, int n) {
    double max_resid = 0.0;
    printf("\n=== ПРОВЕРКА A*v - lambda*v ===\n");
    for (int j = 0; j < n; j++) {
        double* Av = (double*)malloc(n * sizeof(double));
        double* resid = (double*)malloc(n * sizeof(double));
        // Вычисляем A * v_j
        for (int i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (int k = 0; k < n; k++)
                Av[i] += A[i][k] * eigenvectors[k][j];
            resid[i] = Av[i] - eigenvalues[j] * eigenvectors[i][j];
        }
        double norm_resid = 0.0;
        for (int i = 0; i < n; i++)
            norm_resid += resid[i] * resid[i];
        norm_resid = sqrt(norm_resid);
        free(Av);
        free(resid);
        printf("Невязка для lambda%d = %e\n", j+1, norm_resid);
        if (norm_resid > max_resid) max_resid = norm_resid;
    }
    printf("Максимальная невязка: %e\n", max_resid);
}

//  ОСНОВНАЯ ФУНКЦИЯ 
int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_NUMERIC, "C");
    int n, choice;
    double** A = NULL;
    
    printf("=== QR-АЛГОРИТМ ДЛЯ СИММЕТРИЧНОЙ МАТРИЦЫ ===\n");
    printf("Ввод матрицы:\n1 - из файла matrix.txt");
    printf("2 - с клавиатуры\nВаш выбор: ");
    scanf("%d", &choice);
    
    if (choice == 1) {
        A = read_matrix_from_file(&n);
    } else if (choice == 2) {
        A = read_matrix_from_keyboard(&n);
    } else {
        printf("Неверный выбор.\n");
        return 1;
    }
    
    if (A == NULL) {
        printf("Ошибка при чтении матрицы.\n");
        return 1;
    }
    
    printf("\nВведённая матрица:\n");
    print_matrix(A, n, "Исходная матрица");
    
    // Проверка симметричности (предупреждение)
    int symmetric = 1;
    for (int i = 0; i < n && symmetric; i++)
        for (int j = i+1; j < n; j++)
            if (fabs(A[i][j] - A[j][i]) > 1e-8) 
                { symmetric = 0; break; }
    if (!symmetric) {
        printf("\nПредупреждение: матрица не симметрична!");
        printf(" Алгоритм может работать некорректно.\n");
    }
    
    double eps;
    int max_iter;
    printf("\nВведите точность (eps > 0): ");
    scanf("%lf", &eps);
    printf("Введите максимальное число итераций: ");
    scanf("%d", &max_iter);
    
    double* eigenvalues = (double*)malloc(n * sizeof(double));
    double** eigenvectors = alloc_matrix(n);
    if (eigenvalues == NULL || eigenvectors == NULL) {
        printf("Ошибка выделения памяти.\n");
        free_matrix(A, n);
        return 1;
    }
    
    int iter_used = 0;
    int success = qr_algorithm(A, eigenvalues, eigenvectors, n, eps, max_iter, &iter_used);
    
    if (success) {
        printf("\n=== РЕЗУЛЬТАТЫ ===\n");
        printf("Количество итераций QR-алгоритма: %d\n", iter_used);
        print_eigenvalues(eigenvalues, n);
        print_eigenvectors(eigenvectors, n);
        
        check_residuals(A, eigenvalues, eigenvectors, n);
    } else {
        printf("\nQR-алгоритм не сошёлся за %d итераций.\n", max_iter);
    }
    
    free_matrix(A, n);
    free_matrix(eigenvectors, n);
    free(eigenvalues);
    
    return 0;
}