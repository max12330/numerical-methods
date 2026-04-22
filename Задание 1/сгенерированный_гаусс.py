import numpy as np

def gauss_elimination(A, b, n):
    """Метод Гаусса"""
    Ab = np.hstack((A.astype(float), b.reshape(-1, 1)))
    
    # Выводим исходную матрицу только если n <= 10
    if n <= 10:
        print("Исходная расширенная матрица [A|b]:")
        for row in Ab:
            formatted_row = [f"{x:.2f}" for x in row]
            print("  ".join(formatted_row))
    
    # Прямой ход
    for i in range(n):
        # Выбор ведущего элемента
        max_row = i
        for k in range(i + 1, n):
            if abs(Ab[k, i]) > abs(Ab[max_row, i]):
                max_row = k
        
        if max_row != i:
            Ab[[i, max_row]] = Ab[[max_row, i]]
            
            if n <= 4:
                print(
                f"\nШаг {i+1}: Перестановка строк {i+1} и {max_row+1}")
                for row in Ab:
                    formatted_row = [f"{x:.2f}" for x in row]
                    print("  ".join(formatted_row))
        
        # Исключение
        for k in range(i + 1, n):
            factor = Ab[k, i] / Ab[i, i]
            Ab[k, i:] -= factor * Ab[i, i:]
        
        if n <= 4:
            print(f"\nШаг {i+1}: После исключения переменной x{i+1}")
            for row in Ab:
                formatted_row = [f"{x:.2f}" for x in row]
                print("  ".join(formatted_row))
    
    # Выводим матрицу после прямого хода только для n от 5 до 10
    if 5 <= n <= 10:
        print("\nМатрица после прямого хода:")
        for row in Ab:
            formatted_row = [f"{x:.2f}" for x in row]
            print("  ".join(formatted_row))
    
    # Обратный ход
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = Ab[i, -1]
        for j in range(i + 1, n):
            x[i] -= Ab[i, j] * x[j]
        x[i] /= Ab[i, i]
    
    return x

def main():
    n = int(input("Введите размерность матрицы n: "))
    
    # Создание матрицы A
    A = np.zeros((n, n))
    for i in range(n):
        diag_value = 34 + i
        row_value = diag_value * 0.01
        for j in range(n):
            if i == j:
                A[i, j] = diag_value
            else:
                A[i, j] = row_value
    
    # Создание вектора b
    x_true = np.array([34 + i for i in range(n)])
    b = A @ x_true
    
    # Решение
    solution = gauss_elimination(A.copy(), b.copy(), n)
    
    # Вывод результатов
    if n > 10:
        # Только значения x в строку
        print(" ".join([f"{val:.6f}" for val in solution]))
    else:
        print("\nВектор решения x:")
        for i, val in enumerate(solution):
            print(f"x[{i+1}] = {val:.6f}")

if __name__ == "__main__":
    main()