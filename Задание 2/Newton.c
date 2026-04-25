#include <stdio.h>
#include <math.h>
#include <locale.h>
#include <stdlib.h>
#include <windows.h>

/* Определение системы нелинейных уравнений:
   f1(x, y) = x^2 + y^2 - 4 = 0
   f2(x, y) = x*y - 1 = 0
*/

double f1(double x, double y) {
    return x * x + y * y - 4.0;
}

double f2(double x, double y) {
    return x * y - 1.0;
}

double df1dx(double x, double y) {
    return 2.0 * x;
}

double df1dy(double x, double y) {
    return 2.0 * y;
}

double df2dx(double x, double y) {
    return y;
}

double df2dy(double x, double y) {
    return x;
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_NUMERIC, "C");

    double x0, y0, eps;
    int max_iter;

    printf("Решение системы нелинейных уравнений методом Ньютона:\n");
    printf("  f1(x,y) = x^2 + y^2 - 4 = 0\n");
    printf("  f2(x,y) = x*y - 1 = 0\n\n");

    // Ввод начального приближения
    printf("Введите начальное приближение (x0 y0): ");
    if (scanf("%lf %lf", &x0, &y0) != 2) {
        printf("Ошибка: некорректный ввод начального приближения.\n");
        return 1;
    }

    // Ввод точности
    printf("Введите точность (eps > 0): ");
    if (scanf("%lf", &eps) != 1 || eps <= 0.0) {
        printf("Ошибка: точность должна быть положительным числом.\n");
        return 1;
    }

    // Ввод максимального числа итераций
    printf("Введите максимальное число итераций: ");
    if (scanf("%d", &max_iter) != 1 || max_iter <= 0) {
        printf("Ошибка: максимальное число итераций должно быть положительным целым числом.\n");
        return 1;
    }

    printf("\nНачальное приближение: (%.6lf, %.6lf)\n", x0, y0);
    printf("Точность: %e\n", eps);
    printf("Максимум итераций: %d\n\n", max_iter);

    double x = x0, y = y0;
    int iter;
    // норма разности между двумя последними приближениями
    double norm_delta_last = 0.0;
    // норма невязки в последней точке
    double norm_F_last = 0.0;

    for (iter = 0; iter < max_iter; ++iter) {
        // Вычисление значений функций и матрицы Якоби
        double f1_val = f1(x, y);
        double f2_val = f2(x, y);
        double J11 = df1dx(x, y);
        double J12 = df1dy(x, y);
        double J21 = df2dx(x, y);
        double J22 = df2dy(x, y);

        double det = J11 * J22 - J12 * J21;
        if (fabs(det) < 1e-12) {
            printf("Ошибка: якобиан вырожден (det = %e) на итерации %d\n", det, iter);
            break;
        }

        // Решение системы J * delta = -F
        double dx = (-J22 * f1_val + J12 * f2_val) / det;
        double dy = (J21 * f1_val - J11 * f2_val) / det;

        double x_new = x + dx;
        double y_new = y + dy;

        // Норма разности между приближениями
        norm_delta_last = sqrt(dx * dx + dy * dy);

        // Вычисление невязки в новом приближении
        double f1_new = f1(x_new, y_new);
        double f2_new = f2(x_new, y_new);
        norm_F_last = sqrt(f1_new * f1_new + f2_new * f2_new);

        // Печать информации об итерации
        printf("Итерация %2d: x = %.8lf, y = %.8lf, |Δ| = %e, |F| = %e\n",
               iter + 1, x_new, y_new, norm_delta_last, norm_F_last);

        // Обновление переменных
        x = x_new;
        y = y_new;

        // Проверка критериев остановки
        if (norm_delta_last < eps || norm_F_last < eps) {
            printf("\nРешение найдено.\n");
            break;
        }
    }

    // Вывод результатов
    printf("\n=================================\n");
    printf("Приближённое решение:\n");
    printf("x = %.10lf\n", x);
    printf("y = %.10lf\n", y);
    printf("Количество итераций: %d\n", iter + 1);
    printf("Норма разности между двумя последними приближениями: %e\n", norm_delta_last);
    printf("Погрешность (норма невязки |F|): %e\n", norm_F_last);

    if (iter == max_iter && norm_delta_last >= eps && norm_F_last >= eps) {
        printf("\nПредупреждение: достигнуто максимальное число итераций без выполнения критерия остановки.\n");
    }

    return 0;
}