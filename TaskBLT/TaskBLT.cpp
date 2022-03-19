#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h> 
#include <direct.h>
#include <map>

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::fixed;
using std::setw;
using std::setprecision;
using std::ifstream;

void Tom(double* A, double* B, double* C, double* F, int n, int nXStart, double* Fun, int nAll);
double** getArray(int nY, int nX);
bool deleteArray(double**& arr, int nY);
void borderConditionsPsi(double** arr, int nY, int nX, double hy, double hx, double xBarrierLenght, double yDistancesToBarrierDown, double yDistancesToBarrierUp, double xDistancesToBarrier, double** speedX, double** speedY, double** teta);
double getMax(double** arr, int nY, int nX);
std::map<std::string, double> setParams();
std::string getPath();

int main()
{
	setlocale(LC_ALL, "Russian");

	std::map<std::string, double> mp = setParams();
	//число Рейнольдса
	double Re = (mp["Re"]) ? mp["Re"] : 0.0;
	//число Грасгофа
	double Gr = (mp["Gr"]) ? mp["Gr"] : 0.0;
	//число Прандтля
	double Pr = (mp["Pr"]) ? mp["Pr"] : 0.0;
	cout << "Re: " << Re << "\nPr: " << Pr << "\nGr: " << Gr << endl;

	/** Объект **/
	//размеры образца по x, y
	double xLenght = 3.0, yLenght = 3.0;
	//количество областей по x
	const int n = 3;
	//количество областей по y
	const int m = 3;
	//количество узлов по x
	int nX = (mp["nX"]) ? static_cast<int>(mp["nX"]) : 100;
	//количество узлов по y
	int nY = (mp["nY"]) ? static_cast<int>(mp["nY"]) : 100;

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
	double** teta = getArray(nY + 1, nX + 1);
	double** speedX = getArray(nY + 1, nX + 1);
	double** speedY = getArray(nY + 1, nX + 1);

	double** psi_n = getArray(nY + 1, nX + 1);
	double** omega_n = getArray(nY + 1, nX + 1);
	double** teta_n = getArray(nY + 1, nX + 1);


	//шаг по времени
	const double dt = 0.001;
	//шаг по x
	hx = xLenght / nX;
	//шаг по y
	hy = yLenght / nY;
	cin.get();
	borderConditionsPsi(psi, nY + 1, nX + 1, hy, hx, xBarrierLenght, yDistancesToBarrierDown, yDistancesToBarrierUp, xDistancesToBarrier, speedX, speedY, teta);
	if (mp["web"] == 1)
		for (int i = nY; i >= 0; i--) {
			for (int j = 0; j <= nX; j++) {
				cout << teta[i][j] << "   \t";
			}
			cout << endl;
		}

	int i1, i2, j1;

	i1 = static_cast<int>(yDistancesToBarrierDown / hy);

	i2 = static_cast<int>((yDistancesToBarrierUp) / hy) + 1;

	j1 = static_cast<int>((xDistancesToBarrier + xBarrierLenght) / hx);
	cout << "Test: " << j1 * hx << endl;
	cout << "i1:" << i1 << "; i2:" << i2 << "; j1:" << j1 << endl;
	double timer = 0.0;
	int k = 0;
	int max = std::max(nX, nY);
	double* A = new double[max + 1];
	double* B = new double[max + 1];
	double* C = new double[max + 1];
	double* F = new double[max + 1];
	while (timer < 10) {
		timer = dt * k++;
		cout << "Итерация: " << k << "; Время: " << timer << "; Точность: ";

		for (int i = 0; i <= nY; i++) {
			for (int j = 0; j <= nX; j++) {
				psi_n[i][j] = psi[i][j];
				omega_n[i][j] = omega[i][j];
				teta_n[i][j] = teta[i][j];
			}
		}

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

		if (mp["omega"] == 1) {
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
					F[j] = omega[i][j] - dt * Gr * (teta[i][j + 1] - teta[i][j - 1]) / (2.0 * hx * pow(Re, 2));
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
					F[j - j1] = omega[i][j] - dt * Gr * (teta[i][j + 1] - teta[i][j - 1]) / (2.0 * hx * pow(Re, 2));
				}

				A[nX - j1] = -1.0;
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
					F[j] = omega[i][j] - dt * Gr * (teta[i][j + 1] - teta[i][j - 1]) / (2.0 * hx * pow(Re, 2));
				}

				A[nX] = -1.0;
				B[nX] = 1.0;
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

			//прогонка по Y 3 область для вихрей
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

			for (int i = 0; i <= nY; i++)
			{
				omega[i][nX] = omega[i][nX - 1];
			}
		}

		if (mp["teta"] == 1) {
			//прогонка по X 1 область для температуры
			for (int i = 1; i < i1; i++) {
				A[0] = 0.0;
				B[0] = 1.0;
				C[0] = 0.0;
				F[0] = 0.0;

				for (int j = 1; j < nX; j++) {
					A[j] = -dt * (abs(speedX[i][j]) + speedX[i][j]) / (2.0 * hx) - dt / (Re * Pr * pow(hx, 2));
					B[j] = 1.0 + dt * abs(speedX[i][j]) / hx + 2.0 * dt / (Re * Pr * pow(hx, 2));
					C[j] = -dt * (abs(speedX[i][j]) - speedX[i][j]) / (2.0 * hx) - dt / (Re * Pr * pow(hx, 2));
					F[j] = teta[i][j];
				}

				A[nX] = -1.0;
				B[nX] = 1.0;
				C[nX] = 0.0;
				F[nX] = 0.0;
				Tom(A, B, C, F, nX, 0, teta[i], nX);
			}

			//прогонка по X 2 область для температуры
			for (int i = i1; i < i2; i++) {
				A[0] = 0.0;
				B[0] = 1.0;
				C[0] = 0.0;
				F[0] = teta[i][j1];

				for (int j = j1 + 1; j < nX; j++) {
					A[j - j1] = -dt * (abs(speedX[i][j]) + speedX[i][j]) / (2.0 * hx) - dt / (Re * Pr * pow(hx, 2));
					B[j - j1] = 1.0 + dt * abs(speedX[i][j]) / hx + 2.0 * dt / (Re * Pr * pow(hx, 2));
					C[j - j1] = -dt * (abs(speedX[i][j]) - speedX[i][j]) / (2.0 * hx) - dt / (Re * Pr * pow(hx, 2));
					F[j - j1] = teta[i][j];
				}

				A[nX - j1] = -1.0;
				B[nX - j1] = 1.0;
				C[nX - j1] = 0.0;
				F[nX - j1] = 0.0;

				Tom(A, B, C, F, nX - j1, j1, teta[i], nX);
			}

			//прогонка по X 3 область для температуры
			for (int i = i2; i < nY; i++) {
				A[0] = 0.0;
				B[0] = 1.0;
				C[0] = 0.0;
				F[0] = 0.0;

				for (int j = 1; j < nX; j++) {
					A[j] = -dt * (abs(speedX[i][j]) + speedX[i][j]) / (2.0 * hx) - dt / (Re * Pr * pow(hx, 2));
					B[j] = 1.0 + dt * abs(speedX[i][j]) / hx + 2.0 * dt / (Re * Pr * pow(hx, 2));
					C[j] = -dt * (abs(speedX[i][j]) - speedX[i][j]) / (2.0 * hx) - dt / (Re * Pr * pow(hx, 2));
					F[j] = teta[i][j];
				}

				A[nX] = -1.0;
				B[nX] = 1.0;
				C[nX] = 0.0;
				F[nX] = 0.0;

				Tom(A, B, C, F, nX, 0, teta[i], nX);
			}

			//прогонка по Y 1 область для температуры
			for (int j = 1; j <= j1; j++) {


				for (int k = 0; k <= nY; k++) {
					e[k] = teta[k][j];
				}

				A[0] = 0.0;
				B[0] = -1.0;
				C[0] = 1.0;
				F[0] = 0.0;

				for (int i = 1; i < i1; i++) {
					A[i] = -dt * (abs(speedY[i][j]) + speedY[i][j]) / (2.0 * hy) - dt / (Re * Pr * pow(hy, 2));
					B[i] = 1.0 + dt * abs(speedY[i][j]) / hy + 2.0 * dt / (Re * Pr * pow(hy, 2));
					C[i] = -dt * (abs(speedY[i][j]) - speedY[i][j]) / (2.0 * hy) - dt / (Re * Pr * pow(hy, 2));
					F[i] = teta[i][j];
				}

				A[i1] = 0.0;
				B[i1] = 1.0;
				C[i1] = 0.0;
				F[i1] = teta[i1][j];

				Tom(A, B, C, F, i1, 0, e, i1);

				//прогонка по Y 2 область для температуры
				A[0] = 0.0;
				B[0] = 1.0;
				C[0] = 0.0;
				F[0] = teta[i2 - 1][j];

				for (int i = i2; i < nY; i++) {
					A[i - i2 + 1] = -dt * (abs(speedY[i][j]) + speedY[i][j]) / (2.0 * hy) - dt / (Re * Pr * pow(hy, 2));
					B[i - i2 + 1] = 1.0 + dt * abs(speedY[i][j]) / hy + 2.0 * dt / (Re * Pr * pow(hy, 2));
					C[i - i2 + 1] = -dt * (abs(speedY[i][j]) - speedY[i][j]) / (2.0 * hy) - dt / (Re * Pr * pow(hy, 2));
					F[i - i2 + 1] = teta[i][j];
				}

				A[nY - i2 + 1] = -1.0;
				B[nY - i2 + 1] = 1.0;
				C[nY - i2 + 1] = 0.0;
				F[nY - i2 + 1] = 0.0;

				Tom(A, B, C, F, nY - i2 + 1, i2 - 1, e, nY);

				for (int k = 0; k <= nY; k++) {
					teta[k][j] = e[k];
				}
			}

			//прогонка по Y 3 область для температуры
			for (int j = j1 + 1; j < nX; j++) {

				A[0] = 0.0;
				B[0] = -1.0;
				C[0] = 1.0;
				F[0] = 0.0;

				for (int i = 1; i < nY; i++) {
					A[i] = -dt * (abs(speedY[i][j]) + speedY[i][j]) / (2.0 * hy) - dt / (Re * Pr * pow(hy, 2));
					B[i] = 1.0 + dt * abs(speedY[i][j]) / hy + 2.0 * dt / (Re * Pr * pow(hy, 2));
					C[i] = -dt * (abs(speedY[i][j]) - speedY[i][j]) / (2.0 * hy) - dt / (Re * Pr * pow(hy, 2));
					F[i] = teta[i][j];
				}

				A[nY] = -1.0;
				B[nY] = 1.0;
				C[nY] = 0.0;
				F[nY] = 0.0;

				Tom(A, B, C, F, nY, 0, e, nY);

				for (int k = 0; k <= nY; k++) {
					teta[k][j] = e[k];
				}
			}

			for (int i = 0; i <= nY; i++) {
				teta[i][0] = teta[i][1];
				teta[i][nY] = teta[i][nY - 1];
			}

			for (int j = 1; j < nX; j++) {
				teta[0][j] = teta[1][j];
				teta[nY][j] = teta[nY - 1][j];
			}
		}
		delete[] e;

		double** temp_psi = getArray(nY + 1, nX + 1);
		double** temp_omega = getArray(nY + 1, nX + 1);
		double** temp_teta = getArray(nY + 1, nX + 1);

		for (int i = 0; i <= nY; i++) {
			for (int j = 0; j <= nX; j++) {
				temp_psi[i][j] = abs(psi_n[i][j] - psi[i][j]);
				temp_omega[i][j] = abs(omega_n[i][j] - omega[i][j]);
				temp_teta[i][j] = abs(teta_n[i][j] - teta[i][j]);
			}
		}
		double eps = 0.00001;
		double maximum = std::max(getMax(temp_psi, nY + 1, nX + 1), std::max(getMax(temp_omega, nY + 1, nX + 1), getMax(temp_teta, nY + 1, nX + 1)));
		cout << maximum << endl;
		if (maximum < eps) {
			break;
		}
		deleteArray(temp_psi, nY + 1);
		deleteArray(temp_teta, nY + 1);
		deleteArray(temp_omega, nY + 1);
	}
	delete[] A;
	delete[] B;
	delete[] C;
	delete[] F;

	std::string path = "D:\\";
	mp["path"] = 1;
	if (mp["path"] == 1)
		path = getPath() + "\\\\";

	if (mp["web"] == 1)
		for (int i = nY; i >= 0; i--) {
			if (i % 1 == 0) {
				for (int j = 0; j <= nX; j++) {
					if (j % 1 == 0)
						cout << setprecision(3) << fixed << setw(3) << psi[i][j] << "   ";
				}
				cout << endl;
			}
		}
	{
		ofstream out(path + "data_psi.dat");
		for (int i = nY; i >= 0; i--) {
			for (int j = 0; j <= nX; j++) {
				out << j * hx << " " << i * hy << " " << psi[i][j] << endl;
			}
		}
		out.close();
	}
	{
		ofstream out(path + "data_teta.dat");
		for (int i = nY; i >= 0; i--) {
			for (int j = 0; j <= nX; j++) {
				out << j * hx << " " << i * hy << " " << teta[i][j] << endl;
			}
		}
		out.close();
	}
	{
		ofstream out(path + "data_omega.dat");
		for (int i = nY; i >= 0; i--) {
			for (int j = 0; j <= nX; j++) {
				out << j * hx << " " << i * hy << " " << omega[i][j] << endl;
			}
		}
		out.close();
	}
	deleteArray(psi, nY + 1);
	deleteArray(teta, nY + 1);
	deleteArray(omega, nY + 1);
	cin.get();
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

void borderConditionsPsi(double** arr, int nY, int nX, double hy, double hx, double xBarrierLenght, double yDistancesToBarrierDown, double yDistancesToBarrierUp, double xDistancesToBarrier, double** speedX, double** speedY, double** teta)
{
	/** Граничные условия **/
	double Ux = 1.0, Uy = 0.0, tetaB = 4.0;

	//граничные условия на нижней границе
	for (int j = 0; j < nX; j++) {
		arr[0][j] = 0.0;
		teta[0][j] = 0.0;
	}

	//граничные условия на левой границе
	for (int i = 0; i < nY; i++) {
		int nXStart;

		if (i * hy < yDistancesToBarrierDown) {
			nXStart = static_cast<int>(xDistancesToBarrier / hx);

			arr[i][nXStart] = i * hy * Ux;

			teta[i][nXStart] = 0.0;
		}
		if (i * hy >= yDistancesToBarrierDown && i * hy <= yDistancesToBarrierUp) {
			nXStart = static_cast<int>((xDistancesToBarrier + xBarrierLenght) / hx);
			arr[i][nXStart] = yDistancesToBarrierDown * Ux;
			teta[i][nXStart] = tetaB;
			if (abs(i * hy - yDistancesToBarrierDown) < 0.0000001 || abs(i * hy - yDistancesToBarrierUp) < 0.0000001) {
				for (int j = 0; j <= nXStart; j++) {
					arr[i][j] = yDistancesToBarrierDown * Ux;

					teta[i][j] = tetaB;
				}
			}
		}
		if (i * hy > yDistancesToBarrierUp) {
			nXStart = static_cast<int>(xDistancesToBarrier / hx);

			arr[i][nXStart] = Ux * (i * hy - yDistancesToBarrierUp) + Ux * yDistancesToBarrierDown;

			teta[i][nXStart] = 0.0;
		}
	}

	//граничные условия на верхней границе
	for (int j = 1; j < nX; j++) {
		arr[nY - 1][j] = arr[nY - 1][0];
		teta[nY - 1][j] = 0.0;
	}

	//граничные условия на правой границе
	for (int i = nY - 2; i >= 1; i--) {
		arr[i][nX - 1] = 0.0;
		teta[i][nX - 1] = 0.0;
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

double getMax(double** arr, int nY, int nX)
{
	double max = 0.0;
	for (int i = 0; i < nY; i++) {
		for (int j = 0; j < nX; j++) {
			if (max < arr[i][j]) {
				max = arr[i][j];
			}
		}
	}
	return max;
}

std::map<std::string, double> setParams()
{
	std::string name = _getcwd(NULL, 0);
	std::string newName;
	size_t len = name.length();
	for (int i = 0; i < len; i++) {
		if (name[i] == '\\') {
			newName += "\\";
		}
		newName += name[i];
	}
	ifstream in(std::string(newName + "\\\\Params.txt"));
	std::string buff;
	std::map<std::string, double> mp;
	while (in >> buff) {
		double num;
		in >> num;
		mp[buff] = num;
	}
	return mp;
}

std::string getPath()
{
	std::string name = _getcwd(NULL, 0);
	std::string newName;
	size_t len = name.length();
	for (int i = 0; i < len; i++) {
		if (name[i] == '\\') {
			newName += "\\";
		}
		newName += name[i];
	}
	return newName;
}