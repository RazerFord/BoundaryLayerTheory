#include <iostream>
#include <fstream>
#include <iomanip>

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::fixed;
using std::setw;
using std::setprecision;

void Tom(double* A, double* B, double* C, double* F, int n, int nXStart, double* Fun, int nAll);
double** getArray(int nY, int nX);
bool deleteArray(double**& arr, int nY);
void borderConditionsPsi(double** arr, int nY, int nX, double hy, double hx, double xBarrierLenght, double yDistancesToBarrierDown, double yDistancesToBarrierUp, double xDistancesToBarrier, double** speedX, double** speedY);

int main()
{
	setlocale(LC_ALL, "Russian");

	//число Рейнольдса
	double Re = 1.0;

	/** Объект **/
	//размеры образца по x, y
	double xLenght = 3.0, yLenght = 3.0;
	//количество областей по x
	const int n = 3;
	//количество областей по y
	const int m = 3;
	//количество узлов по x
	int nX = 100;
	//количество узлов по y
	int nY = 100;

	/** Препятствие **/
	//размеры препятствия по x, y
	double xBarrierLenght = 1.0, yBarrierLenght = 1.0;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ X
	double xDistancesToBarrier = 0.0;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ Y
	double yDistancesToBarrierDown = 1.0;
	double yDistancesToBarrierUp = yBarrierLenght + yDistancesToBarrierDown;

	double hx = xLenght / nX;
	int testX = static_cast<int>(xBarrierLenght / hx);
	while (true) {
		if (!(abs((testX * hx) - xBarrierLenght) > 0.0000001))
			break;
		nX++;
		hx = xLenght / nX;
		testX = static_cast<int>(xBarrierLenght / hx);
	}
	double hy = yLenght / nY;
	int testY = static_cast<int>(yBarrierLenght / hy);
	while (true) {
		if (!(abs(testY * hy - yBarrierLenght) > 0.0000001))
			break;
		nY++;
		hy = yLenght / nY;
		testY = static_cast<int>(xBarrierLenght / hy);
	}

	cout << "Количество точек по X: " << nX << endl;
	cout << "Количество точек по Y: " << nY << endl;

	double** psi = getArray(nY + 1, nX + 1);
	double** omega = getArray(nY + 1, nX + 1);
	double** speedX = getArray(nY + 1, nX + 1);
	double** speedY = getArray(nY + 1, nX + 1);

	//шаг по времени
	const double dt = 0.005;
	//шаг по x
	hx = xLenght / nX;
	//шаг по y
	hy = yLenght / nY;
	cin.get();
	borderConditionsPsi(psi, nY + 1, nX + 1, hy, hx, xBarrierLenght, yDistancesToBarrierDown, yDistancesToBarrierUp, xDistancesToBarrier, speedX, speedY);
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
	cout << j1 * hx << endl;
	cout << i1 << "   " << i2 << "   " << j1 << endl;
	cin.get();
	ofstream out("D://datka.dat");
	double timer = 0.0;
	int k = 0;
	int max = std::max(nX, nY);
	double* A = new double[max + 1];
	double* B = new double[max + 1];
	double* C = new double[max + 1];
	double* F = new double[max + 1];
	while (timer < 10) {
		timer = dt * k++;
		cout << timer << endl;


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
				F[j] = -psi[i][j] + dt * omega[i][j];
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

			for (int j = j1 + 1; j < nX; j++) {
				A[j - j1] = dt / pow(hx, 2);
				B[j - j1] = -1.0 - 2.0 * dt / pow(hx, 2);
				C[j - j1] = dt / pow(hx, 2);
				F[j - j1] = -psi[i][j] + dt * omega[i][j];
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
				F[j] = -psi[i][j] + dt * omega[i][j];
			}

			A[nX] = -1.0;
			B[nX] = 1.0;
			C[nX] = 0.0;
			F[nX] = 0.0;

			Tom(A, B, C, F, nX, 0, psi[i], nX);
		}

		double* e = new double[nY + 1]{ 0 };

		//прогонка по Y 1 и 2 область
		for (int j = 1; j <= j1; j++) {


			for (int k = 0; k <= nY; k++) {
				e[k] = psi[k][j];
			}

			if (j * hx <= xBarrierLenght + xDistancesToBarrier) {

				A[0] = 0.0;
				B[0] = 1.0;
				C[0] = 0.0;
				F[0] = psi[0][j];

				for (int i = 1; i < i1; i++) {
					A[i] = dt / pow(hy, 2);
					B[i] = -1.0 - 2.0 * dt / pow(hy, 2);
					C[i] = dt / pow(hy, 2);
					F[i] = -psi[i][j];
				}

				A[i1] = 0.0;
				B[i1] = 1.0;
				C[i1] = 0.0;
				F[i1] = psi[i1][j];

				Tom(A, B, C, F, i1, 0, e, i1);

				A[0] = 0.0;
				B[0] = 1.0;
				C[0] = 0.0;
				F[0] = psi[i2 - 1][j];

				for (int i = i2; i < nY; i++) {
					A[i - i2 + 1] = dt / pow(hy, 2);
					B[i - i2 + 1] = -1.0 - 2.0 * dt / pow(hy, 2);
					C[i - i2 + 1] = dt / pow(hy, 2);
					F[i - i2 + 1] = -psi[i][j];
				}

				A[nY - i2 + 1] = 0.0;
				B[nY - i2 + 1] = 1.0;
				C[nY - i2 + 1] = 0.0;
				F[nY - i2 + 1] = psi[nY][j];

				Tom(A, B, C, F, nY - i2 + 1, i2 - 1, e, nY);

				for (int k = 0; k <= nY; k++) {
					psi[k][j] = e[k];
				}
			}
		}

		//прогонка по Y 3 область
		for (int j = j1 + 1; j < nX; j++) {

			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = psi[0][j];

			for (int i = 1; i < nY; i++) {
				A[i] = dt / pow(hy, 2);
				B[i] = -1.0 - 2.0 * dt / pow(hy, 2);
				C[i] = dt / pow(hy, 2);
				F[i] = -psi[i][j];
			}

			A[nY] = 0.0;
			B[nY] = 1.0;
			C[nY] = 0.0;
			F[nY] = psi[nY][j];

			Tom(A, B, C, F, nY, 0, e, nY);

			for (int k = 0; k <= nY; k++) {
				psi[k][j] = e[k];
			}
		}

		for (int i = 1; i < nY; i++)
		{
			psi[i][nX] = psi[i][nX - 1];
		}

		//Вычисления скорости по X
		for (int i = 1; i < nY; i++) {
			for (int j = 1; j < nX; j++) {
				speedX[i][j] = (psi[i + 1][j] - psi[i - 1][j]) / hy / 2;
			}
			speedX[i][nX] = speedX[i][nX - 1];
		}

		//for (int i = i1; i < i2; i++) {
		//	for (int j = j1 + 1; j < nX; j++) {
		//		speedX[i][j] = (psi[i + 1][j] - psi[i][j]) / hy;
		//	}
		//	speedX[i][nX] = speedX[i][nX - 1];
		//}

		//for (int i = i2; i < nY; i++) {
		//	for (int j = 1; j < nX; j++) {
		//		speedX[i][j] = (psi[i + 1][j] - psi[i][j]) / hy;
		//	}
		//	speedX[i][nX] = speedX[i][nX - 1];
		//}

		//Вычисления скорости по Y
		for (int i = 1; i < nY; i++) {
			for (int j = 1; j < nX; j++) {
				speedY[i][j] = -(psi[i][j + 1] - psi[i][j - 1]) / hx / 2;
			}
			speedY[i][nX] = speedY[i][nX - 1];
		}

		//for (int i = i1; i < i2; i++) {
		//	for (int j = j1; j < nX; j++) {
		//		speedY[i][j] = -(psi[i][j + 1] - psi[i][j]) / hx;
		//	}
		//	speedY[i][nX] = speedY[i][nX - 1];
		//}

		//for (int i = i2 - 1; i < nY; i++) {
		//	for (int j = 0; j < nX; j++) {
		//		speedY[i][j] = -(psi[i][j + 1] - psi[i][j]) / hx;
		//	}
		//	speedY[i][nX] = speedY[i][nX - 1];
		//	}

		//прогонка по X 1 область для вихрей
		for (int i = 1; i < i1; i++) {
			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = 0.0;

			for (int j = 1; j < nX; j++) {
				A[j] = -dt * (abs(speedX[i][j]) + speedX[i][j]) / (2.0 * hx) - dt / (Re * pow(hx, 2));
				B[j] = 1.0 + dt * abs(speedX[i][j]) / hx + 2.0 * dt / (Re * pow(hx, 2));
				C[j] = -dt * (abs(speedX[i][j]) - speedX[i][j]) / (2.0 * hx) - dt / (Re * pow(hx, 2));
				F[j] = omega[i][j];

				/*A[j] = -dt * speedX[i][j] / (2.0 * hx) - dt / (Re * pow(hx, 2));
				B[j] = 1.0 + 2.0 * dt / (Re * pow(hx, 2));
				C[j] = dt * speedX[i][j] / (2.0 * hx) - dt / (Re * pow(hx, 2));
				F[j] = omega[i][j];*/
			}

			A[nX] = -1.0;
			B[nX] = 1.0;
			C[nX] = 0.0;
			F[nX] = 0.0;

			Tom(A, B, C, F, nX, 0, omega[i], nX);
		}

		//прогонка по X 2 область для вихрей
		for (int i = i1; i < i2; i++) {
			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = 2.0 * (psi[i][j1 + 1] - psi[i][j1]) / pow(hx, 2);

			for (int j = j1 + 1; j < nX; j++) {
				A[j - j1] = -dt * (abs(speedX[i][j]) + speedX[i][j]) / (2.0 * hx) - dt / (Re * pow(hx, 2));
				B[j - j1] = 1.0 + dt * abs(speedX[i][j]) / hx + 2.0 * dt / (Re * pow(hx, 2));
				C[j - j1] = -dt * (abs(speedX[i][j]) - speedX[i][j]) / (2.0 * hx) - dt / (Re * pow(hx, 2));
				F[j - j1] = omega[i][j];

				/*A[j - j1] = -dt * speedX[i][j] / (2.0 * hx) - dt / (Re * pow(hx, 2));
				B[j - j1] = 1.0 + 2.0 * dt / (Re * pow(hx, 2));
				C[j - j1] = dt * speedX[i][j] / (2.0 * hx) - dt / (Re * pow(hx, 2));
				F[j - j1] = omega[i][j];*/
			}

			A[nX - j1] = 1.0;
			B[nX - j1] = 1.0;
			C[nX - j1] = 0.0;
			F[nX - j1] = 0.0;

			Tom(A, B, C, F, nX - j1, j1, omega[i], nX);
		}

		//прогонка по X 3 область для вихрей
		for (int i = i2; i < nY; i++) {
			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = 0.0;

			for (int j = 1; j < nX; j++) {
				A[j] = -dt * (abs(speedX[i][j]) + speedX[i][j]) / (2.0 * hx) - dt / (Re * pow(hx, 2));
				B[j] = 1.0 + dt * abs(speedX[i][j]) / hx + 2.0 * dt / (Re * pow(hx, 2));
				C[j] = -dt * (abs(speedX[i][j]) - speedX[i][j]) / (2.0 * hx) - dt / (Re * pow(hx, 2));
				F[j] = omega[i][j];

				/*A[j] = -dt * speedX[i][j] / (2.0 * hx) - dt / (Re * pow(hx, 2));
				B[j] = 1.0 + 2.0 * dt / (Re * pow(hx, 2));
				C[j] = dt * speedX[i][j] / (2.0 * hx) - dt / (Re * pow(hx, 2));
				F[j] = omega[i][j];*/
			}

			A[nX] = 1.0;
			B[nX] = -1.0;
			C[nX] = 0.0;
			F[nX] = 0.0;

			Tom(A, B, C, F, nX, 0, omega[i], nX);
		}

		//прогонка по Y 1 и 2 область для вихрей
		for (int j = 1; j <= j1; j++) {


			for (int k = 0; k <= nY; k++) {
				e[k] = omega[k][j];
			}

			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = 2.0 * (psi[1][j] - psi[0][j]) / pow(hy, 2);

			for (int i = 1; i < i1; i++) {
				A[i] = -dt * (abs(speedY[i][j]) + speedY[i][j]) / (2.0 * hy) - dt / (Re * pow(hy, 2));
				B[i] = 1.0 + dt * abs(speedY[i][j]) / hy + 2.0 * dt / (Re * pow(hy, 2));
				C[i] = -dt * (abs(speedY[i][j]) - speedY[i][j]) / (2.0 * hy) - dt / (Re * pow(hy, 2));
				F[i] = omega[i][j];

				/*A[i] = -dt * speedY[i][j] / (2.0 * hy) - dt / (Re * pow(hy, 2));
				B[i] = 1.0  + 2.0 * dt / (Re * pow(hy, 2));
				C[i] = dt * speedY[i][j] / (2.0 * hy) - dt / (Re * pow(hy, 2));
				F[i] = omega[i][j];*/
			}

			A[i1] = 0.0;
			B[i1] = 1.0;
			C[i1] = 0.0;
			F[i1] = -2.0 * (psi[i1][j] - psi[i1 - 1][j]) / pow(hy, 2);

			Tom(A, B, C, F, i1, 0, e, i1);

			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = 2.0 * (psi[i2][j] - psi[i2 - 1][j]) / pow(hy, 2);

			for (int i = i2; i < nY; i++) {
				A[i - i2 + 1] = -dt * (abs(speedY[i][j]) + speedY[i][j]) / (2.0 * hy) - dt / (Re * pow(hy, 2));
				B[i - i2 + 1] = 1.0 + dt * abs(speedY[i][j]) / hy + 2.0 * dt / (Re * pow(hy, 2));
				C[i - i2 + 1] = -dt * (abs(speedY[i][j]) - speedY[i][j]) / (2.0 * hy) - dt / (Re * pow(hy, 2));
				F[i - i2 + 1] = omega[i][j];

			/*	A[i - i2 + 1] = -dt * speedY[i][j] / (2.0 * hy) - dt / (Re * pow(hy, 2));
				B[i - i2 + 1] = 1.0 + 2.0 * dt / (Re * pow(hy, 2));
				C[i - i2 + 1] = dt * speedY[i][j] / (2.0 * hy) - dt / (Re * pow(hy, 2));
				F[i - i2 + 1] = omega[i][j];*/
			}

			A[nY - i2 + 1] = 0.0;
			B[nY - i2 + 1] = 1.0;
			C[nY - i2 + 1] = 0.0;
			F[nY - i2 + 1] = -2.0 * (psi[nY][j] - psi[nY - 1][j]) / pow(hy, 2);

			Tom(A, B, C, F, nY - i2 + 1, i2 - 1, e, nY);

			for (int k = 0; k <= nY; k++) {
				omega[k][j] = e[k];
			}
		}

		//прогонка по Y 3 область
		for (int j = j1 + 1; j < nX; j++) {

			A[0] = 0.0;
			B[0] = 1.0;
			C[0] = 0.0;
			F[0] = 2.0 * (psi[1][j] - psi[0][j]) / pow(hy, 2);

			for (int i = 1; i < nY; i++) {
				A[i] = -dt * (abs(speedY[i][j]) + speedY[i][j]) / (2.0 * hy) - dt / (Re * pow(hy, 2));
				B[i] = 1.0 + dt * abs(speedY[i][j]) / hy + 2.0 * dt / (Re * pow(hy, 2));
				C[i] = -dt * (abs(speedY[i][j]) - speedY[i][j]) / (2.0 * hy) - dt / (Re * pow(hy, 2));
				F[i] = omega[i][j];

				/*A[i] = -dt * speedY[i][j] / (2.0 * hy) - dt / (Re * pow(hy, 2));
				B[i] = 1.0 + 2.0 * dt / (Re * pow(hy, 2));
				C[i] = dt * speedY[i][j] / (2.0 * hy) - dt / (Re * pow(hy, 2));
				F[i] = omega[i][j];*/
			}

			A[nY] = 0.0;
			B[nY] = 1.0;
			C[nY] = 0.0;
			F[nY] = -2.0 * (psi[nY][j] - psi[nY - 1][j]) / pow(hy, 2);

			Tom(A, B, C, F, nY, 0, e, nY);

			for (int k = 0; k <= nY; k++) {
				omega[k][j] = e[k];
			}
		}

		delete[] e;

		for (int i = 0; i <= nY; i++)
		{
			omega[i][nX] = omega[i][nX - 1];
		}
	}
	delete[] A;
	delete[] B;
	delete[] C;
	delete[] F;
	for (int i = nY; i >= 0; i--) {
		if (i % 1 == 0) {
			for (int j = 0; j <= nX; j++) {
				if (j % 1 == 0)
					cout << setprecision(3) << fixed << setw(3) << psi[i][j] << "   ";
			}
			cout << endl;
		}
	}cin.get();

	cout << "finish \n";

	for (int i = nY; i >= 0; i--) {
		for (int j = 0; j <= nX; j++) {
			out << j * hx << " " << i * hy << " " << psi[i][j] << endl;
		}
	}
	deleteArray(psi, nY + 1);
	out.close();

	ofstream f1("D:\\speedX.dat");
	for (int i = 0; i <= nY; i++) {
		f1 << i * hy << " " << speedX[i][j1/2] << " " << speedX[i][nX/2] << " " << speedX[i][nX - 1] << endl;
	}
	f1.close();
	return 0;
}


void Tom(double* A, double* B, double* C, double* F, int n, int nStart, double* Fun, int nAll)
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
		Fun[nStart + i] = ai[i] * Fun[nStart + i + 1] + bi[i];
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
	catch (std::exception) {
		return false;
	}
}

void borderConditionsPsi(double** arr, int nY, int nX, double hy, double hx, double xBarrierLenght, double yDistancesToBarrierDown, double yDistancesToBarrierUp, double xDistancesToBarrier, double** speedX, double** speedY)
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
	for (int i = 0; i < nY - 1; i++) {
		if (i != 0 && i * hy < yDistancesToBarrierDown || i * hy > yDistancesToBarrierUp) {
			speedX[i][0] = Ux;
			speedY[i][0] = 0.0;
		}
	}
	//граничные условия для скорости на верхней границе
	for (int j = 0; j < nX; j++) {
		speedX[nY - 1][j] = 0.0;
		speedY[nY - 1][j] = 0.0;
	}
	//граничные условия для скорости на правой
	for (int i = 0; i < nY; i++) {
		speedX[i][nX - 1] = 0.0;
		speedY[i][nX - 1] = 0.0;
	}
}