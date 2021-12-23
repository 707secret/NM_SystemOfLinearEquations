#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

//7 вариант - N = 700, l = 16
const int N = 10;
const int l = 5;
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
vector<double> JacobiLenta(int N, int lw, vector<vector<double>> A, vector<double> b)
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
void norma(int N, vector<double> solTrue, vector<double> solOur)
{
	double max = fabs(solOur[0] - solTrue[0]);
	for (int i = 0; i < N; i++)
	{
		if (fabs(solOur[i] - solTrue[i]) > max)
			max = fabs(solOur[i] - solTrue[i]);
	}
	cout << "Norm = " << max << endl;
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
	Lenta1 = MatrixtoLenta(N, 2 * l + 1, M1, Lenta1);
	Lenta2 = MatrixtoLenta(N, 2 * l + 1, M2, Lenta2);
	Lenta3 = MatrixtoLenta(N, 2 * l + 1, M3, Lenta3);
	//printMatrix(Lenta1);
	cout << "__________________" << endl;
	cout << endl;


	cout << "2. Jacobi method iterations and solutions: " << endl;
	cout << " -------------------------- " << endl;
	cout << "For q1 = " << q1 << " ";
	x1 = JacobiLenta(N, 2 * l + 1, Lenta1, b1);
	norma(N, x1, sol);

	cout << "For q2 = " << q2 << " ";
	x2 = JacobiLenta(N, 2 * l + 1, Lenta2, b2);
	norma(N, x2, sol);

	cout << "For q3 = " << q3 << " ";
	x3 = JacobiLenta(N, 2 * l + 1, Lenta3, b3);
	norma(N, x3, sol);
	cout << endl;

	return 0;
}
