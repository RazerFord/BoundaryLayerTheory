#include <iostream>

using namespace std;

void Tom(double* A, double* B, double* C, double* F);
double** getArray(int nY, int nX);
bool deleteArray(double**& arr, int nY);

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
	const int nX = 2000;
	//количество узлов по y
	const int nY = 20000;

	//double* ddd = new double[nX * nY]{ 0 };
	//cout << ddd[2000 * 2000 - 1];
	//cin.get();
	double** psiN = getArray(nY, nX);
	cout << psiN[nY - 1][nX - 1] << endl;
	cin.get();
	double** psi = getArray(nY, nX);
	cout << psi[nY - 1][nX - 1] << endl;

	/** Препятствие **/
	//размеры препятствия по x, y
	double xBarrierLenght = 0.2, yBarrierLenght = 0.1;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ X
	double xDistancesToBarrier = 0.0;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ Y
	double yDistancesToBarrierDown = 0.2;
	double yDistancesToBarrierUp = yBarrierLenght + yDistancesToBarrierDown;

	/** Граничные условия **/
	/*
	*нужно задать граничные условия
	*/

	//шаг по времени
	const double dt = 0.01;
	//шаг по x
	double hx = xLenght / nX;
	//шаг по y
	double hy = yLenght / nY;

	for (int i = 0; i <= nY; i++) {
		cout << i << "  :=  " << i * hy << "    ";
		int nXStart;
		if (i * hy < yDistancesToBarrierDown) {
			nXStart = static_cast<int>(xDistancesToBarrier / hx);
			cout << "Первая область";
		}
		if (i * hy >= yDistancesToBarrierDown && i * hy <= yDistancesToBarrierUp) {
			nXStart = static_cast<int>((xDistancesToBarrier + xBarrierLenght) / hx);
			cout << "Вторая область";
		}
		if (i * hy > yDistancesToBarrierUp) {
			nXStart = static_cast<int>(xDistancesToBarrier / hx);
			cout << "Третья область";
		}

		cout << "\n/*** " << nXStart << ", " << nXStart * hx << " ***/\n\n";

		//		double* A = new double[nX + 1 - nXStart];
		//		double* B = new double[nX + 1 - nXStart];
		//		double* C = new double[nX + 1 - nXStart];
		//		double* F = new double[nX + 1 - nXStart];

		for (int j = nXStart; j <= nX; j++) {
			//			A[j - nXStart] = dt / pow(hx, 2);
			//			B[j - nXStart] = -1.0 - dt / pow(hx, 2);
			//			C[j - nXStart] = dt / pow(hx, 2);
			//			F[j - nXStart] = -psi[i][j];
		}

		//		Tom(A, B, C, F);
	}

	deleteArray(psiN, nY);
	deleteArray(psi, nY);

	return 0;
}

void Tom(double* A, double* B, double* C, double* F)
{

}

double** getArray(int nY, int nX)
{
	double** arr = new double* [nY + 1];
	for (int i = 0; i <= nY; i++) {
		arr[i] = new double[nX + 1]{ 0 };
	}
	return arr;
}

bool deleteArray(double**& arr, int nY)
{
	try {
		for (int i = 0; i <= nY; i++) {
			delete[] arr[i];
		}
		delete[] arr;
		return true;
	}
	catch (exception) {
		return false;
	}
}

