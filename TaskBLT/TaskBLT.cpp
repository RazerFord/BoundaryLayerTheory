#include <iostream>
#include <fstream>

using namespace std;

void Tom(double* A, double* B, double* C, double* F, int n, int nXStart, double* Fun, int nAll);
double** getArray(int nY, int nX);
bool deleteArray(double**& arr, int nY);
void borderConditions(double** arr, int nY, int nX, double hy, double hx, double xBarrierLenght, double yDistancesToBarrierDown, double yDistancesToBarrierUp, double xDistancesToBarrier, double** speedX, double** speedY);

int main()
{
	setlocale(LC_ALL, "Russian");

	/** Объект **/
	//размеры образца по x, y
	double xLenght = 3.0, yLenght = 3.0;
	//количество областей по x
	const int n = 3;
	//количество областей по y
	const int m = 3;
	//количество узлов по x
	const int nX = 10;
	//количество узлов по y
	const int nY = 12;

	double** psi = getArray(nY + 1, nX + 1);
	double** swirl = getArray(nY + 1, nX + 1);
	double** speedX = getArray(nY + 1, nX + 1);
	double** speedY = getArray(nY + 1, nX + 1);

	/** Препятствие **/
	//размеры препятствия по x, y
	double xBarrierLenght = 1.0, yBarrierLenght = 1.0;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ X
	double xDistancesToBarrier = 0.0;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ Y
	double yDistancesToBarrierDown = 1.0;
	double yDistancesToBarrierUp = yBarrierLenght + yDistancesToBarrierDown;

	//шаг по времени
	const double dt = 0.01;
	//шаг по x
	double hx = xLenght / nX;
	//шаг по y
	double hy = yLenght / nY;

	borderConditions(psi, nY + 1, nX + 1, hy, hx, xBarrierLenght, yDistancesToBarrierDown, yDistancesToBarrierUp, xDistancesToBarrier, speedX, speedY);
	for (int i = nY; i >= 0; i--) {
		for (int j = 0; j <= nX; j++) {
			cout << psi[i][j] << "   \t";
		}
		cout << endl;
	}

	int i1, i2, j1;

	i1 = static_cast<int>(yDistancesToBarrierDown / hy);

	i2 = static_cast<int>((yDistancesToBarrierUp) / hy) + 1;

	j1 = static_cast<int>((xDistancesToBarrier + xBarrierLenght) / hx);
	cout << i1 << "   " << i2 << "   " << j1 << endl;
	cin.get();
	ofstream out("D://datka.dat");
	double timer = 0.0;
	int k = 0;
	while (timer < 11) {
		timer = dt * k++;
		cout << timer << endl;

		int max = std::max(nX, nY);
		double* A = new double[max + 1];
		double* B = new double[max + 1];
		double* C = new double[max + 1];
		double* F = new double[max + 1];

		//прогонка по X 1 область
		for (int i = 1; i < i1; i++) {
			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = psi[i][0];

			for (int j = 1; j < nX; j++) {
				A[j] = dt / pow(hx, 2);
				B[j] = -1.0 - 2.0 * dt / pow(hx, 2);
				C[j] = dt / pow(hx, 2);
				F[j] = -psi[i][j] - dt * swirl[i][j];
			}

			A[nX] = -1.0;
			B[nX] = 1.0;
			C[nX] = 0.0;
			F[nX] = 0.0;

			Tom(A, B, C, F, nX, 0, psi[i], nX);
		}

		//прогонка по X 2 область
		for (int i = i1; i < i2; i++) {
			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = psi[i][j1];

			for (int j = 1; j < nX - j1; j++) {
				A[j] = dt / pow(hx, 2);
				B[j] = -1.0 - 2.0 * dt / pow(hx, 2);
				C[j] = dt / pow(hx, 2);
				F[j] = -psi[i][j + j1] - dt * swirl[i][j + j1];
			}

			A[nX - j1] = -1.0;
			B[nX - j1] = 1.0;
			C[nX - j1] = 0.0;
			F[nX - j1] = 0.0;

			Tom(A, B, C, F, nX - j1, j1, psi[i], nX);
		}

		//прогонка по X 3 область
		for (int i = i2; i < nY; i++) {
			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = psi[i][0];

			for (int j = 1; j < nX; j++) {
				A[j] = dt / pow(hx, 2);
				B[j] = -1.0 - 2.0 * dt / pow(hx, 2);
				C[j] = dt / pow(hx, 2);
				F[j] = -psi[i][j] - dt * swirl[i][j];
			}

			A[nX] = -1.0;
			B[nX] = 1.0;
			C[nX] = 0.0;
			F[nX] = 0.0;

			Tom(A, B, C, F, nX, 0, psi[i], nX);
		}

		//прогонка по Y
		for (int j = 1; j < nX; j++) {

			int nYStart, nYEnd, size;
			double* A = new double[nY + 1]{ 0 };
			double* B = new double[nY + 1]{ 0 };
			double* C = new double[nY + 1]{ 0 };
			double* F = new double[nY + 1]{ 0 };

			double* psiN = new double[nY + 1]{ 0 };

			for (int k = 0; k <= nY; k++) {
				psiN[k] = psi[k][j];
			}

			if (j * hx <= xBarrierLenght + xDistancesToBarrier) {
				nYStart = 0;
				nYEnd = static_cast<int>(yDistancesToBarrierDown / hy);
				size = nYEnd - nYStart;

				A[0] = 0.0;
				B[0] = 1.0;
				C[0] = 0.0;
				F[0] = psi[0][j];

				for (int i = 1; i < size; i++) {
					A[i] = dt / pow(hy, 2);
					B[i] = -1.0 - 2.0 * dt / pow(hy, 2);
					C[i] = dt / pow(hy, 2);
					F[i] = -psi[i + nYStart][j];
				}

				A[size] = 0.0;
				B[size] = 1.0;
				C[size] = 0.0;
				F[size] = psi[size][j];

				Tom(A, B, C, F, size, nYStart, psiN, nYEnd);

				for (int k = nYStart; k <= nYEnd; k++) {
					psi[k][j] = psiN[k];
				}

				nYStart = static_cast<int>(yDistancesToBarrierUp / hy);
				nYEnd = nY;
				size = nYEnd - nYStart;

				A[0] = 0.0;
				B[0] = 1.0;
				C[0] = 0.0;
				F[0] = psi[nYStart][j];

				for (int i = 1; i < size; i++) {
					A[i] = dt / pow(hy, 2);
					B[i] = -1.0 - 2.0 * dt / pow(hy, 2);
					C[i] = dt / pow(hy, 2);
					F[i] = -psi[i + nYStart][j];
				}

				A[size] = 0.0;
				B[size] = 1.0;
				C[size] = 0.0;
				F[size] = psi[nYEnd][j];

				Tom(A, B, C, F, size, nYStart, psiN, nYEnd);

				for (int k = nYStart; k <= nYEnd; k++) {
					psi[k][j] = psiN[k];
				}
			}
			if (j * hx > xBarrierLenght + xDistancesToBarrier) {
				nYStart = 0;
				nYEnd = nY;
				size = nYEnd - nYStart;

				A[0] = 0.0;
				B[0] = 1.0;
				C[0] = 0.0;
				F[0] = psi[nYStart][j];

				for (int i = 1; i < size; i++) {
					A[i] = dt / pow(hy, 2);
					B[i] = -1.0 - 2.0 * dt / pow(hy, 2);
					C[i] = dt / pow(hy, 2);
					F[i] = -psi[i + nYStart][j];
				}

				A[size] = 0.0;
				B[size] = 1.0;
				C[size] = 0.0;
				F[size] = psi[nYEnd][j];

				Tom(A, B, C, F, size, nYStart, psiN, nYEnd);

				for (int k = nYStart; k <= nYEnd; k++) {
					psi[k][j] = psiN[k];
				}
			}

			delete[] A;
			delete[] B;
			delete[] C;
			delete[] F;
			delete[] psiN;
		}

		delete[] A;
		delete[] B;
		delete[] C;
		delete[] F;

		for (int i = 1; i < nY; i++)
		{
			psi[i][nX] = psi[i][nX - 1];
		}

		//Короче тут вычисляется скорость
		for (int i = 1; i < nY; i++) {
			for (int j = 1; j < nX; j++) {
				speedX[i][j] = (psi[i + 1][j] - psi[i][j]) / hy;
				speedY[i][j] = -(psi[i][j + 1] - psi[i][j]) / hx;
			}
			speedX[i][nX] = speedX[i][nX - 1];
			speedY[i][nX] = speedY[i][nX - 1];
		}

		//Граничные условия для функции вихря
		//for (int i = 0; i <= nY; i++) {
		//	for (int j = 0; j <= nX; j++) {

		//	}
		//}

		//Прогонка для вихрей по X
		//for (int i = 1; i < nY; i++) {

		//	int nXStart;

		//	if (i * hy < yDistancesToBarrierDown) {
		//		nXStart = static_cast<int>(xDistancesToBarrier / hx);
		//	}
		//	if (i * hy >= yDistancesToBarrierDown && i * hy <= yDistancesToBarrierUp) {
		//		nXStart = static_cast<int>((xDistancesToBarrier + xBarrierLenght) / hx);
		//	}
		//	if (i * hy > yDistancesToBarrierUp) {
		//		nXStart = static_cast<int>(xDistancesToBarrier / hx);
		//	}

		//	int size = nX - nXStart;

		//	double* A = new double[size + 1];
		//	double* B = new double[size + 1];
		//	double* C = new double[size + 1];
		//	double* F = new double[size + 1];

		//	A[0] = 0.0;
		//	B[0] = 1.0;
		//	C[0] = 0.0;
		//	F[0] = psi[i][nXStart];

		//	for (int j = 1; j < size; j++) {
		//		A[j] = dt / pow(hx, 2);
		//		B[j] = -1.0 - 2.0 * dt / pow(hx, 2);
		//		C[j] = dt / pow(hx, 2);
		//		F[j] = -psi[i][j + nXStart] - dt * swirl[i][j + nXStart];
		//	}

		//	A[size] = -1.0;
		//	B[size] = 1.0;
		//	C[size] = 0.0;
		//	F[size] = 0.0;

		//	Tom(A, B, C, F, size, nXStart, psi[i], nX);

		//	delete[] A;
		//	delete[] B;
		//	delete[] C;
		//	delete[] F;
		//}
	}
	cout << "finish \n";

	for (int i = nY; i >= 0; i--) {
		for (int j = 0; j <= nX; j++) {
			out << j * hx << " " << i * hy << " " << psi[i][j] << endl;
		}
	}
	deleteArray(psi, nY + 1);
	out.close();

	return 0;
}


void Tom(double* A, double* B, double* C, double* F, int n, int nXStart, double* Fun, int nAll)
{
	double* ai = new double[n] {0};
	double* bi = new double[n] {0};
	ai[0] = -C[0] / B[0];
	bi[0] = F[0] / B[0];

	for (int i = 1; i < n; i++)
	{
		ai[i] = -C[i] / (A[i] * ai[i - 1] + B[i]);
		bi[i] = (F[i] - A[i] * bi[i - 1]) / (A[i] * ai[i - 1] + B[i]);
	}

	Fun[nAll] = (F[n] - A[n] * bi[n - 1]) / (A[n] * ai[n - 1] + B[n]);

	for (int i = n - 1; i >= 0; i--) {
		Fun[nXStart + i] = ai[i] * Fun[nXStart + i + 1] + bi[i];
	}
	delete[] ai;
	delete[] bi;
}

double** getArray(int nY, int nX)
{
	double** arr = new double* [nY];
	for (int i = 0; i < nY; i++) {
		arr[i] = new double[nX] { 0 };
	}
	return arr;
}

bool deleteArray(double**& arr, int nY)
{
	try {
		for (int i = 0; i < nY; i++) {
			delete[] arr[i];
		}
		delete[] arr;
		return true;
	}
	catch (exception) {
		return false;
	}
}

void borderConditions(double** arr, int nY, int nX, double hy, double hx, double xBarrierLenght, double yDistancesToBarrierDown, double yDistancesToBarrierUp, double xDistancesToBarrier, double** speedX, double** speedY)
{
	/** Граничные условия **/
	//на левой границе
	double Ux = 1.0, Uy = 0.0;

	//граничные условия на нижней границе
	for (int j = 0; j < nX; j++) {
		arr[0][j] = 0.0;
	}

	//граничные условия на левой границе
	for (int i = 0; i < nY; i++) {
		int nXStart;

		if (i * hy < yDistancesToBarrierDown) {
			nXStart = static_cast<int>(xDistancesToBarrier / hx);

			arr[i][nXStart] = i * hy * Ux;
		}
		if (i * hy >= yDistancesToBarrierDown && i * hy <= yDistancesToBarrierUp) {
			nXStart = static_cast<int>((xDistancesToBarrier + xBarrierLenght) / hx);
			arr[i][nXStart] = yDistancesToBarrierDown;
			if (abs(i * hy - yDistancesToBarrierDown) < 0.0000001 || abs(i * hy - yDistancesToBarrierUp) < 0.0000001) {
				for (int j = 0; j <= nXStart; j++) {
					arr[i][j] = yDistancesToBarrierDown;
				}
			}
		}
		if (i * hy > yDistancesToBarrierUp) {
			nXStart = static_cast<int>(xDistancesToBarrier / hx);

			arr[i][nXStart] = Ux * (i * hy - yDistancesToBarrierUp) + yDistancesToBarrierDown;
		}
	}

	//граничные условия на верхней границе
	for (int j = 1; j < nX; j++) {
		arr[nY - 1][j] = arr[nY - 1][0];
	}

	//граничные условия на правой границе
	for (int i = nY - 2; i > 1; i--) {
		arr[i][nX - 1] = 0.0;
	}

	//граничные условия для скорости на нижней границе
	for (int j = 0; j < nX; j++) {
		speedX[0][j] = 0.0;
		speedY[0][j] = 0.0;
	}
	//граничные условия для скорости на левой границе
	for (int i = 0; i < nY; i++) {
		speedX[i][0] = Ux;
		speedY[i][0] = 0.0;
	}
	//граничные условия для скорости на верхней границе
	for (int j = 0; j < nX; j++) {
		speedX[nY - 1][j] = 0.0;
		speedY[nY - 1][j] = 0.0;
	}
	//граничные условия для скорости на правой
	for (int i = 0; i < nY; i++) {
		speedX[i][0] = 0.0;
		speedY[i][0] = 0.0;
	}
}