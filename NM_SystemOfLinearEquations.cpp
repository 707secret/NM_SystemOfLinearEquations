#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

//7 вариант - N = 700, l = 16
const int N = 700;
const int l = 16;
const double eps = 0.000001;

using namespace std;


void printMatrix(vector<vector<double>> vec)
{
	for (int i = 0; i < vec.size(); i++)
	{
		for (int j = 0; j < vec[i].size(); j++)
		{
			cout << fixed << showpoint;
			cout << setprecision(4) << " " << vec[i][j];
		}
		cout << endl;
	}
	cout << endl;
}

void printVector(vector<double> vec)
{
	for (int i = 0; i < vec.size(); i++)
	{
		cout << " " << vec[i];
	}
	cout << endl;
}

vector<double> MatrixVectorMultiply(int N, vector<vector<double>> Matrix, vector<double> Vector)
{
	vector<double> VectorSolution(N, 0);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			VectorSolution[i] += Matrix[i][j] * Vector[j];
		}
	}
	return VectorSolution;
}

vector<double> MatrixLentaVectorMultiply(int N, int lw, vector<vector<double>> MatrixLenta, vector<double> Vector)
{
	vector<double> VectorSolution(N, 0);
	for (int i = 0; i < N; i++)
	{

		for (int j = 0; j < N; j++)
		{
			int ind = j - i + lw;
			if ((ind >= 0) && (ind < 2 * lw + 1))
			{
				VectorSolution[i] += MatrixLenta[i][ind] * Vector[j];
			}
		}
	}
	return VectorSolution;
}

double VectorVectorMult(int N, vector<double> Vector1, vector<double> Vector2)
{
	double result = 0;
	for (int i = 0; i < N; i++)
	{
		result += Vector1[i] * Vector2[i];
	}
	return result;
}

vector<vector<double>> MatrixMatrixmult(int N, vector<vector<double>> M1, vector<vector<double>> M2)
{
	vector<vector<double>> Mres(N, vector<double>(N, 0));
	double num;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			num = 0;
			for (int k = 0; k < N; k++)
			{
				num += M1[i][k] * M2[k][j];
			}
			Mres[i][j] = num;
		}
	}
	return Mres;
}

vector<vector<double>> MatrixTranspose(int N, vector<vector<double>> M)
{
	double num;
	double t;
	for (int i = 0; i < N; ++i)
	{
		for (int j = i; j < N; ++j)
		{
			t = M[i][j];
			M[i][j] = M[j][i];
			M[j][i] = t;
		}
	}

	return M;
}

vector<vector<double>> MatrixtoLenta(int N, int l, vector<vector<double>> Matrix, vector<vector<double>> MatrixLenta)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (Matrix[i][j] != 0)
				MatrixLenta[i][j - i + l] = Matrix[i][j];
		}
	}
	return MatrixLenta;
}

//Метод Якоби!
vector<double> JacobiMethod(int N, int lw, vector<vector<double>> A, vector<double> b)
{
	vector<double> xprev(N, 0); 
	vector<double> x(N, 0);
	double diff;
	int iter = 0;
	do
	{
		iter++;
		for (int i = 0; i < N; i++)
		{
			x[i] = b[i];
			for (int j = 0; j < N; j++)
			{
				int ind = j - i + lw;
				if ((j != i) && (ind >= 0) && (ind < 2 * lw + 1))
					x[i] -= A[i][ind] * xprev[j];
			}
			x[i] /= A[i][lw];
		}
		diff = fabs(x[0] - xprev[0]);
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - xprev[i]) > diff)
				diff = fabs(x[i] - xprev[i]);
			xprev[i] = x[i];
		}
	} while (diff > eps);
	cout << "Iterations: " << iter << endl;
	return x;
}

//Вычисление нормы
void Norm(int N, vector<double> solTrue, vector<double> solOur)
{
	double max = fabs(solOur[0] - solTrue[0]);
	for (int i = 0; i < N; i++)
	{
		if (fabs(solOur[i] - solTrue[i]) > max)
			max = fabs(solOur[i] - solTrue[i]);
	}
	cout << "Norm = " << max << endl;
}

vector<double> SORLenta(int N, int l, vector<vector<double>> A, vector<double> b, double relaxation)
{
	vector<double> xprev(N, 0); // Начальнoe приближение
	vector<double> x(N, 0); // Текущее приближение
	double diff; // Разница между текущим и предыдущим приближением, итерации вплоть до разницы меньше eps
	int iter = 0;

	do
	{
		iter++;
		for (int i = 0; i < N; i++)
		{
			x[i] = b[i];
			for (int j = 0; j < i; j++)
			{
				int ind = j - i + l;
				if ((j != i) && (ind >= 0) && (ind < 2 * l + 1))
					x[i] -= A[i][ind] * x[j];
			}
			for (int j = i + 1; j < N; j++)
			{
				int ind = j - i + l;
				if ((j != i) && (ind >= 0) && (ind < 2 * l + 1))
					x[i] -= A[i][ind] * xprev[j];
			}
			x[i] /= A[i][l];
			x[i] *= relaxation;
			x[i] += (1 - relaxation) * xprev[i];
		}
		diff = fabs(x[0] - xprev[0]);
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - xprev[i]) > diff)
				diff = fabs(x[i] - xprev[i]);
			xprev[i] = x[i];
		}
	} while (diff > eps);
	cout << " Iterations: " << iter << endl;
	return x;
}


int main(void)
{
	srand(time(NULL));
	vector<vector<double>> A1(N, vector<double>(N, 0)); //ленточные матрицы
	vector<vector<double>> A2(N, vector<double>(N, 0));
	vector<vector<double>> A3(N, vector<double>(N, 0));
	vector<vector<double>> AT1(N, vector<double>(N, 0)); //транспонированные матрицы
	vector<vector<double>> AT2(N, vector<double>(N, 0));
	vector<vector<double>> AT3(N, vector<double>(N, 0));
	vector<vector<double>> M1(N, vector<double>(N, 0)); //матрицы АТ * А
	vector<vector<double>> M2(N, vector<double>(N, 0));
	vector<vector<double>> M3(N, vector<double>(N, 0));

	double q1 = 1.1;
	double q2 = 2;
	double q3 = 10;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if ((l + 1) - abs(i - j) > 0) {
				//генерация внедиагональных элементов ленточной матрицы А
				A1[i][j] = (double)(rand()) / RAND_MAX * (1 - (-1)) + (-1);
			}
			A2[i][j] = A1[i][j];
			A3[i][j] = A1[i][j];
		}
	}

	//создаем диагональное преобладание
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A1[i][i] += fabs(A1[i][j]);
			A2[i][i] += fabs(A2[i][j]);
			A3[i][i] += fabs(A3[i][j]);
		}
		A1[i][i] *= q1;
		A2[i][i] *= q2;
		A3[i][i] *= q3;
	}

	//Вывод матриц, тут видно, что они отличаются только диагоналями
	//printMatrix(A1);
	//printMatrix(A2);
	//printMatrix(A3);


	vector<double> sol(N, 0); //cоздание и генерация вектора точного решения СЛАУ
	for (int i = 0; i < N; i++)
	{
		sol[i] = (double)(rand()) / RAND_MAX * (1 - (-1)) + (-1);
	}

	//cout << "Solution" << endl;
	//printVector1dDouble(sol);

	//вычисление векторов правой части системы
	vector<double> rhs1 = MatrixVectorMultiply(N, A1, sol);
	vector<double> rhs2 = MatrixVectorMultiply(N, A2, sol);
	vector<double> rhs3 = MatrixVectorMultiply(N, A3, sol);

	//транспонированные матрицы
	AT1 = MatrixTranspose(N, A1);
	AT2 = MatrixTranspose(N, A2);
	AT3 = MatrixTranspose(N, A3);

	//AT*A : 
	M1 = MatrixMatrixmult(N, AT1, A1);
	M2 = MatrixMatrixmult(N, AT2, A2);
	M3 = MatrixMatrixmult(N, AT3, A3);

	//b*:
	vector<double> b1 = MatrixVectorMultiply(N, AT1, rhs1);
	vector<double> b2 = MatrixVectorMultiply(N, AT2, rhs2);
	vector<double> b3 = MatrixVectorMultiply(N, AT3, rhs3);

	vector<double> x1(N, 0);
	vector<double> x2(N, 0);
	vector<double> x3(N, 0);
	cout << endl;

	// Создаем ленточные матрицы c шириной ленты 4l+1
	vector<vector<double>> Lenta1(N, vector<double>(4 * l + 1, 0));
	vector<vector<double>> Lenta2(N, vector<double>(4 * l + 1, 0));
	vector<vector<double>> Lenta3(N, vector<double>(4 * l + 1, 0));

	//заполняем ленточные матрицы, на этом у нас 1 задание сделано
	Lenta1 = MatrixtoLenta(N, 2 * l , M1, Lenta1);
	Lenta2 = MatrixtoLenta(N, 2 * l , M2, Lenta2);
	Lenta3 = MatrixtoLenta(N, 2 * l , M3, Lenta3);
	//printMatrix(Lenta1);
	cout << "__________________" << endl;
	cout << endl;


	cout << "2. Jacobi method iterations and solutions: " << endl;
	cout << " -------------------------- " << endl;
	cout << "For q1 = " << q1 << " ";
	x1 = JacobiMethod(N, 2 * l , Lenta1, b1);
	Norm(N, x1, sol);

	cout << "For q2 = " << q2 << " ";
	x2 = JacobiMethod(N, 2 * l , Lenta2, b2);
	Norm(N, x2, sol);

	cout << "For q3 = " << q3 << " ";
	x3 = JacobiMethod(N, 2 * l , Lenta3, b3);
	Norm(N, x3, sol);
	cout << endl;

	cout << "__________________" << endl;
	cout << endl;


	cout << "3. SOR method iterations and solutions: " << endl;
	double w1 = 0.2; //наши значения w от 0 до 2
	double w2 = 1;
	double w3 = 1.9;

	cout << " relax = " << w1 << endl;
	cout << " For q1 = " << q1 << " "; x1 = SORLenta(N, 2 * l, Lenta1, b1, w1);
	Norm(N, x1, sol);
	cout << " For q2 = " << q2 << " "; x2 = SORLenta(N, 2 * l, Lenta2, b2, w1);
	Norm(N, x2, sol);
	cout << " For q3 = " << q3 << " "; x3 = SORLenta(N, 2 * l, Lenta3, b3, w1);
	Norm(N, x3, sol);

	cout << endl;
	cout << " relax = " << w2 << endl;
	cout << " For q1 = " << q1 << " "; x1 = SORLenta(N, 2 * l, Lenta1, b1, w2);
	Norm(N, x1, sol);
	cout << " For q2 = " << q2 << " "; x2 = SORLenta(N, 2 * l, Lenta2, b2, w2);
	Norm(N, x2, sol);
	cout << " For q3 = " << q3 << " "; x3 = SORLenta(N, 2 * l, Lenta3, b3, w2);
	Norm(N, x3, sol);
	cout << endl;

	cout << " relax = " << w3 << endl;
	cout << " For q1 = " << q1 << " "; x1 = SORLenta(N, 2 * l, Lenta1, b1, w3);
	Norm(N, x1, sol);
	cout << " For q2 = " << q2 << " "; x2 = SORLenta(N, 2 * l, Lenta2, b2, w3);
	Norm(N, x2, sol);
	cout << " For q3 = " << q3 << " "; x3 = SORLenta(N, 2 * l, Lenta3, b3, w3);
	Norm(N, x3, sol);
	/*printVector(x1);
	printVector(x2);
	printVector(x3);*/
	cout << endl;


	return 0;
}
