#include <iostream>
#include <fstream>

using namespace std;

void Tom(double* A, double* B, double* C, double* F, int n, int nXStart, double* Fun, int nAll);
double** getArray(int nY, int nX);
bool deleteArray(double**& arr, int nY);
void p();

int main()
{
	setlocale(LC_ALL, "Russian");

	/** Объект **/
	//размеры образца по x, y
	double xLenght = 1.0, yLenght = 0.5;
	//количество областей по x
	const int n = 3;
	//количество областей по y
	const int m = 3;
	//количество узлов по x
	const int nX = 10;
	//количество узлов по y
	const int nY = 5;

	double** psiN = getArray(nY + 1, nX + 1);
	double** psi = getArray(nY + 1, nX + 1);

	/** Препятствие **/
	//размеры препятствия по x, y
	double xBarrierLenght = 0.2, yBarrierLenght = 0.1;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ X
	double xDistancesToBarrier = 0.0;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ Y
	double yDistancesToBarrierDown = 0.2;
	double yDistancesToBarrierUp = yBarrierLenght + yDistancesToBarrierDown;

	/** Граничные условия **/
	//на левой границе
	double Ux = 1.0, Uy = 0.0;
	//на правой границе
	double Psi = 0.0;

	//шаг по времени
	const double dt = 0.01;
	//шаг по x
	double hx = xLenght / nX;
	//шаг по y
	double hy = yLenght / nY;

	ofstream out("D://datka.txt");
	//прогонка по X
	for (int i = 0; i <= nY; i++) {//Фиг пойми что тут ставить i=0 или i=1. Ну с 0 работает неплохо

		int nXStart;

		if (i * hy < yDistancesToBarrierDown) {
			nXStart = static_cast<int>(xDistancesToBarrier / hx);
			//cout << "Первая область";
		}
		if (i * hy >= yDistancesToBarrierDown && i * hy <= yDistancesToBarrierUp) {
			nXStart = static_cast<int>((xDistancesToBarrier + xBarrierLenght) / hx);
			//cout << "Вторая область";
		}
		if (i * hy > yDistancesToBarrierUp) {
			nXStart = static_cast<int>(xDistancesToBarrier / hx);
			//cout << "Третья область";
		}

		int size = nX - nXStart;

		double* A = new double[size + 1];
		double* B = new double[size + 1];
		double* C = new double[size + 1];
		double* F = new double[size + 1];

		A[0] = 0.0; //???
		B[0] = 1.0; //???
		C[0] = 0.0; //???
		F[0] = 0.0; //???

		for (int j = 1; j < size; j++) {
			A[j] = dt / pow(hx, 2);
			B[j] = -1.0 - 2.0 * dt / pow(hx, 2);
			C[j] = dt / pow(hx, 2);
			F[j] = -psi[i][j];
		}

		A[size] = 0.0; //???
		B[size] = 1.0; //???
		C[size] = 0.0; //???
		F[size] = 100.0; //???

		Tom(A, B, C, F, size, nXStart, psi[i], nX);
		delete[] A, B, C, F;
	}

	for (int i = 0; i <= nY; i++) {
		for (int j = 0; j <= nX; j++) {
			cout << psi[i][j] << "   ";
		}
		cout << endl;
	}

	for (int i = 0; i <= nY; i++)
		cout << i << "  " << i * hy << endl;

	//прогонка по Y
	for (int j = 0; j <= nX; j++) {//Фиг пойми что тут ставить i=0 или i=1. Ну с 0 работает неплохо

		int nYStart, nYEnd, size;
		double* A = new double[nY + 1]{ 0 };
		double* B = new double[nY + 1]{ 0 };
		double* C = new double[nY + 1]{ 0 };
		double* F = new double[nY + 1]{ 0 };

		double* psiN = new double[nY + 1];

		for (int k = 0; k <= nY; k++) {
			psiN[k] = psi[k][j];
		}

		if (j * hx <= xBarrierLenght + xDistancesToBarrier) {
			nYStart = 0;
			nYEnd = static_cast<int>(yDistancesToBarrierDown / hy);
			size = nYEnd - nYStart;

			A[0] = 0.0; //???
			B[0] = 1.0; //???
			C[0] = 0.0; //???
			F[0] = 0.0; //???

			for (int i = 0; i < size; i++) {
				A[i] = dt / pow(hy, 2);
				B[i] = -1.0 - dt / pow(hy, 2);
				C[i] = dt / pow(hy, 2);
				F[i] = -psi[i][j];
			}

			A[size] = 0.0; //???
			B[size] = 1.0; //???
			C[size] = 0.0; //???
			F[size] = 1.0; //???

			Tom(A, B, C, F, size, nYStart, psiN, nY);
		}
		if (j * hx <= xBarrierLenght + xDistancesToBarrier) {
			nYStart = static_cast<int>((yDistancesToBarrierDown + yBarrierLenght) / hy);
			nYEnd = nY;
			size = nYEnd - nYStart;

			A[0] = 0.0; //???
			B[0] = 1.0; //???
			C[0] = 0.0; //???
			F[0] = 0.0; //???

			for (int i = 0; i < size; i++) {
				A[i] = dt / pow(hy, 2);
				B[i] = -1.0 - dt / pow(hy, 2);
				C[i] = dt / pow(hy, 2);
				F[i] = -psi[i][j];
			}

			A[size] = 0.0; //???
			B[size] = 1.0; //???
			C[size] = 0.0; //???
			F[size] = 1.0; //???

			Tom(A, B, C, F, size, nYStart, psiN, nY);
		}
		if (j * hx > xBarrierLenght + xDistancesToBarrier) {
			nYStart = 0;
			nYEnd = nY;
			size = nYEnd - nYStart;

			A[0] = 0.0; //???
			B[0] = 1.0; //???
			C[0] = 0.0; //???
			F[0] = 0.0; //???

			for (int i = 0; i < size; i++) {
				A[i] = dt / pow(hy, 2);
				B[i] = -1.0 - dt / pow(hy, 2);
				C[i] = dt / pow(hy, 2);
				F[i] = -psi[i][j];
			}

			A[size] = 0.0; //???
			B[size] = 1.0; //???
			C[size] = 0.0; //???
			F[size] = 1.0; //???

			Tom(A, B, C, F, size, nYStart, psiN, nY);
		}

		for (int k = 0; k <= nY; k++) {
			psi[k][j] = psiN[k];
		}
	}



	for (int i = 0; i <= nY; i++) {
		for (int j = 0; j <= nX; j++) {
			cout << psi[i][j] << "   ";
		}
		cout << endl;
	}

	deleteArray(psiN, nY);
	deleteArray(psi, nY);
	out.close();

	return 0;
}

void p()
{

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
		//	cout << bi[i] << endl; cin.get();
	}

	Fun[nAll] = (F[n] - A[n] * bi[n - 1]) / (A[n] * ai[n - 1] + B[n]);

	//cout << nXStart + n << "   " << nAll << endl; cin.get();
	//cout << Fun[nAll] << "   " << Fun[nXStart + n] << "   " << A[n] << endl;
	//int size = nX - nXStart;

	for (int i = n - 1; i >= 0; i--) {
		//cout << Fun[nXStart + i] << "   " << A[i] << endl;
		Fun[nXStart + i] = ai[i] * Fun[nXStart + i + 1] + bi[i];
	}
	delete[] ai, bi;
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