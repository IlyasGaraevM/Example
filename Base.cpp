#include <iostream>
#include <cstdlib>
#include <math.h>
#include <random>

//*****Задаваемые переменные*****//
#define N 20				// Количество точек
#define eps 0.00005			// Погрешность фи 0.005%
#define shag 1				// Шаг градиентного спуска
#define sigma_fi 1			// Сигма Фи
#define iterNumber 100		// Количество общих итераций
#define eps_lambda 0.001	// Ошибка в расчете собственных чисел

using namespace std;
default_random_engine generator;
normal_distribution<double> distribution(0, eps*1000);

class coord
{
private:
	const double PI = 3.1415926535897932384626433832795;
// Теоретические значения, которые генерируют фи	
	double Vx_0 = -2;						// Проекция скорости на Ох
	double Vy_0 = 3;						// Проекция скорости на Оу
	double x_0 = 40;						// Начальная координата по Ох
	double y_0 = 4;							// Начальная координата по Оу
// Используемые внутренние переменные
	double Oshibka = 100;					// Рассматриваемая ошибка в методе переборки				
	double vectorError[4] = { 0,0,0,0 };	// Вектор ошибок
	double matrixError[4][4];				// Матрица ошибок

public:
	
	double alpha = abs(atan(Vy_0 / Vx_0));	// Угол наклона скорости

	void Result_r_V(double x0, double y0, double Vx0, double Vy0)	// Функция выводящая результат
	{
		cout << "x_teor = " << x_0 << "  x_exp = " << x0 << endl;
		cout << "y_teor = " << y_0 << "  y_exp = " << y0 << endl;
		cout << "Vx_teor = " << Vx_0 << "  Vx_exp = " << Vx0 << endl;
		cout << "Vy_teor = " << Vy_0 << "  Vy_exp = " << Vy0 << endl;
	}

	double Fi_ishodnoe(int t)			// Функция задающая Фи по начальным значениям
	{
		double x_i = x_0 + Vx_0 * t;	// Расчет х в момент t
		double y_i = y_0 + Vy_0 * t;	// Расчет y в момент t
		return atan(y_i / x_i);			// Возвращает фи как arctg(y/x)
	}

	double Fi_s_pogresh(int t)			// Углы фи с добавлением гауссовой погрешности (0, eps)
	{
		return Fi_ishodnoe(t)*(1 + distribution(generator));
	}
//*********************************//
//******Блок метода переборки******//
//*********************************//
	double Proverka(double x0, double y0, double Vx, double Vy, double* Fi)		//Проверяем на минимальную ошибку
	{
		double Sredn_oshibka = 0;
		double iteration;				// Значение на данной итерации
		for (int t = 0; t < N; t++)		// Проводим итерации
		{
			iteration = (Fi[t] - atan((y0 + Vy * t) / (x0 + Vx * t + 0.0000001))) / sigma_fi;
			Sredn_oshibka = Sredn_oshibka + iteration * iteration / 2;			// Суммарная ошибка
		}
		if (Oshibka > Sredn_oshibka / N)										// Выделение случая с наименьшей ошибкой
		{
			Oshibka = Sredn_oshibka / N;										// В случае ошибки в программе
			cout << "oshibka = " << Oshibka << endl;
		}
		else
			return -1;
	} 
		
	void Pereborka_r_V(double* Fi)												// Функция нахождения координат методом переборки
	{
		double x0, y0, Vx0, Vy0;
		for (double x = 1; x < 50; x = x + 0.2)
			for (double y = 1; y < 50; y = y + 0.2)
				for (double Vx = -5; Vx < 10; Vx = Vx + 0.2)
					for (double Vy = -5; Vy < 10; Vy = Vy + 0.2)
						if (Proverka(x, y, Vx, Vy, Fi) != -1)					// Проверяем на минимальную ошибку
						{
							x0 = x;
							y0 = y; 
							Vx0 = Vx;
							Vy0 = Vy;
						}
		Result_r_V(x0,y0,Vx0,Vy0);												// Выводим результат
	}

//************************************//
//******Блок градиентного спуска******//
//************************************//
	double func(double x, double y, double Vx, double Vy, double* Fi)			// Функция разницы Фи
	{
		double sum = 0;
		for (int t = 0; t < N; t++)
			sum = sum + pow((Fi[t] - atan((y + Vy*t) / (x + Vx*t))), 2);
		return sum;
	}
//*****Производные по переменным*****//
	double Gradient_x(double x, double y, double Vx, double Vy, double* Fi)		// Производная по x
	{
		double sum = 0;
		for (int t = 0; t < N; t++)
			sum = sum + 2 * (Fi[t] - atan((y + Vy*t) / (x + Vx*t)))*(y + Vy*t) / ((x + Vx*t)*(x + Vx*t) + (y + Vy*t)*(y + Vy*t));
		return sum;
	}
	
	double Gradient_y(double x, double y, double Vx, double Vy, double* Fi)		// Производная по y
	{
		double sum = 0;
		for (int t = 0; t < N; t++)
			sum = sum - 2 * (Fi[t] - atan((y + Vy*t) / (x + Vx*t)))*(x + Vx*t) / ((x + Vx*t)*(x + Vx*t) + (y + Vy*t)*(y + Vy*t));
		return sum;
	}
	
	double Gradient_Vx(double x, double y, double Vx, double Vy, double* Fi)	// Производная по Vx
	{
		double sum = 0;
		for (int t = 0; t < N; t++)
			sum = sum + 2 * (Fi[t] - atan((y + Vy*t) / (x + Vx*t)))*t*(y + Vy*t) / ((x + Vx*t)*(x + Vx*t) + (y + Vy*t)*(y + Vy*t));
		return sum;
	}
	
	double Gradient_Vy(double x, double y, double Vx, double Vy, double* Fi)	// Производная по Vy
	{
		double sum = 0;
		for (int t = 0; t < N; t++)
			sum = sum - 2 * (Fi[t] - atan((y + Vy*t) / (x + Vx*t)))*t*(x + Vx*t) / ((x + Vx*t)*(x + Vx*t) + (y + Vy*t)*(y + Vy*t));
		return sum;
	}
//*************************//

	int Gradient_spusk(double* Fi)												// Основная программа расчета
	{
		int iter = 0;															// Количество итераций
//***Начальные координаты для градиентного спуска***//
		double x = 40;
		double y = 4;
		double Vx = -2;
		double Vy = 3;
//**************************//
		double func_0 = - 100;													// Предыдущее значение функции
		while ((abs(func(x, y, Vx, Vy, Fi) - func_0) >= eps) && (iter < 1000))	// Цикл для градиентного спуска
			{
				func_0 = func(x, y, Vx, Vy, Fi);
				x = x - shag*Gradient_x(x, y, Vx, Vy, Fi);						// Покоординатное изменение в градиентном спуске
				y = y - shag*Gradient_y(x, y, Vx, Vy, Fi);
				Vx = Vx - shag*Gradient_Vx(x, y, Vx, Vy, Fi);
				Vy = Vy - shag*Gradient_Vy(x, y, Vx, Vy, Fi);

				iter++;															// Подсчет количества итераций
			}
		Result_r_V(x, y, Vx, Vy);												// Полученные результаты
		Error_Delta(x - x_0, y - y_0, Vx - Vx_0, Vy - Vy_0);					// Заполняем данные об ошибке
		return iter;															// Количество итераций
	}

//***********************************//
//****Блок расчета матрицы ошибок****//
//***********************************//
	void Error_Delta(double delta_x, double delta_y, double delta_Vx, double delta_Vy)	// Функция записи ошибки после
	{																					// каждой генерации углов
		vectorError[0] = vectorError[0] + delta_x;								// Записываем ошибку как 4-вектор
		vectorError[1] = vectorError[1] + delta_y;
		vectorError[2] = vectorError[2] + delta_Vx;
		vectorError[3] = vectorError[3] + delta_Vy;
	}

	void Error_Matrix()															// Составляем матрицу ошибок
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				matrixError[i][j] = (vectorError[i] * vectorError[j]) / (iterNumber*iterNumber);
				cout << matrixError[i][j] << ' ';								// Выводим матрицу на экран
			}
			cout << endl;
		}
	}

	void Error_Lambda()															// Расчет собственных значений
	{
// Используемые величины //
		double Stroke0[4] = { 1,0,0,0 };
		double Stroke[4];
		double summ = 0;
		double func0, func1, delta;
//***********************//
		double Lambda[4];														// Значения собственных значений

		do																		// Расчитываем через метож наименьших
		{																		// квадратов собственные значения
			for (int i = 0; i < 4; i++)
				summ = summ + Stroke0[i] * Stroke0[i];
			func0 = sqrt(summ);
			for (int i = 0;	i < 4;	i++)
				Lambda[i] = Stroke0[i] / func0;
			for (int i = 0; i < 4; i++)
			{
				Stroke[i] = 0;
				for (int j = 0; j < 4; j++)
					Stroke[i] = Stroke[i] + matrixError[i][j] * Lambda[j];
			}
			summ = 0;
			for (int i = 0; i < 4; i++)
				summ = summ + Stroke[i] * Stroke[i];
			func1 = sqrt(summ);
			delta = abs(func1 - func0);
			for (int i = 0; i < 4; i++)
				Stroke0[i] = Stroke[i];
			summ = 0;
		} while (delta > eps_lambda);											// Задаем точность измерений
		cout << "func = " << func1 << endl;
		for (int i = 0; i < 4; i++)
			cout << "Lambda = " << Lambda[i] << endl;							// Вывод собственных чисел
	}

/*	double Norm(double x, double y, double Vx, double Vy, double* Fi)				// Выдает число для нормировки шага градиента
	{
		return sqrt(pow(Gradient_x(x, y, Vx, Vy, Fi), 2) + pow(Gradient_y(x, y, Vx, Vy, Fi), 2) + pow(Gradient_Vx(x, y, Vx, Vy, Fi), 2) + pow(Gradient_Vy(x, y, Vx, Vy, Fi), 2));
	}
	double Sigma(double* Fi)
	{
		int flag;									// Позиция наиболее близкого угла
		double max_delta = 0;						// Максимальная разница между углами фи
		double Fi_perpendicular;					// Расчетное значение фи перпендикулярного
		for (int i = 0; i < N - 1; i++)
		{
			delta[i] = Fi[i + 1] - Fi[i];			// Заполняем массив из разниц между углами
			if (delta[i] > max_delta)
			{
				max_delta = delta[i];				// Фиксируем максимальную разницу
				flag = i;							// Номер угла с максимальной разницей
			//	cout << "maxdelta =" << max_delta << endl;
			}
		}
		if (flag <= 0)
			flag = 1;
	//	cout << "Fi_flag=" << Fi[flag] << endl;
		// Случай когда движемся в большую сторону
		if (delta[flag + 1] > delta[flag - 1]) 
			Fi_perpendicular = Fi[flag] + (delta[flag + 1] - delta[flag]) / 2;
		// Случай когда движемся в меньшую сторону
		else
			Fi_perpendicular = Fi[flag] + (delta[flag] - delta[flag - 1]) / 2;
		cout << "Fi perpendicular = " << Fi_perpendicular << endl;
	//	cout << "alpha =" << PI /2 - alpha << endl;
		return abs((PI / 2 - alpha) - Fi_perpendicular);							
	}*/
};

int main()
{
	int waiting;
	double Fi_teor[N];									// Фи без погрешности
	double Fi_exp[N];									// Фи с погрешностью
	coord raschet;

	for (int i = 0; i < N; i++)							// Заполняем массив теоретическими Фи
	{
		Fi_teor[i] = raschet.Fi_ishodnoe(i);
//		cout << "Fi(" << i << ")=" << Fi_teor[i] << endl;
	}
	cout << "New Fi:" << endl;
	
	for (int j = 0; j < 100; j++)										// Генерируем разные ситуации
	{
		for (int i = 0; i < N; i++)										// Считаем новые фи
		{
			Fi_exp[i] = raschet.Fi_s_pogresh(i);						// Через функцию разброса
//			cout << "Fi(" << i << ")=" << Fi_exp[i] << endl;
		}
		cout << "iter = " << raschet.Gradient_spusk(Fi_exp) << endl;	// Считаем по градиентному спуску
		raschet.Error_Matrix();											// Создаем матрицу ошибок
		raschet.Error_Lambda();											// Расчет собственных чисел
	}
	cin >> waiting;
}
