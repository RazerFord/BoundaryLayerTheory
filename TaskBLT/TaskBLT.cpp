#include <iostream>

using namespace std;



int main()
{
	/**Объект**/
	//размеры образца по x, y
	double xLenght = 1.0, yLenght = 0.5;
	//количество областей по x
	const int n = 3;
	//количество областей по y
	const int m = 3;
	//количество узлов по x
	const int nX = 100;
	//количество узлов по y
	const int nY = 8;

	/**Препятствие**/
	//размеры препятствия по x, y
	double xBarrierLenght = 0.2, yBarrierLenght = 0.1;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ X
	double xDistancesToBarrier = 0.0;
	//расстояние от препятствия от точки (0,0) /Это левый нижний угол/ по ОСИ Y
	double yDistancesToBarrierDown = 0.2;
	double yDistancesToBarrierUp = yLenght - yBarrierLenght - yDistancesToBarrierDown;
	std::cout << yDistancesToBarrierUp << "    " << std::endl;

	/**Граничные условия**/
	/*
	*нужно задать граничные условия
	*/

	//шаг по x
	double hx = xLenght / nX;
	//шаг по y
	double hy = yLenght / nY;

	int correction = 0;

	/**Разбиение на области по оси X**/
	/*Первая область*/
	//Количество узлов в первой области по X занимаемой барьером
	int n_xDistancesToBarrier_1 = 0;
	//Количество узлов в первой области по X
	int nX1 = nX - static_cast<int>(xDistancesToBarrier / hx) - n_xDistancesToBarrier_1;
	//Количество узлов в первой области по Y
	int nY1 = static_cast<int>(yDistancesToBarrierDown / hy);
	cout << "   nY1 = " << nY1 << "   nY1*h1 = " << nY1 * hy << std::endl;
	if (nY1 * hy == yDistancesToBarrierDown) {
		nY1 -= 1;
		correction += 1;
	}

	/*Третья область*/
	//Количество узлов в первой области по X занимаемой барьером
	int n_xDistancesToBarrier_3 = 0;
	//Количество узлов в первой области по X
	int nX3 = nX - static_cast<int>(xDistancesToBarrier / hx) - n_xDistancesToBarrier_3;
	//Количество узлов в первой области по Y
	int nY3 = static_cast<int>(yDistancesToBarrierUp / hy);
	cout << "   nY3 = " << nY3 << "   nY3*h3 = " << nY3 * hy << std::endl;
	if (nY3 * hy == yDistancesToBarrierUp) {
		nY3 -= 1;
		correction += 1;
	}

	/*Вторая область*/
	//Количество узлов в первой области по X занимаемой барьером
	int n_xDistancesToBarrier_2 = static_cast<int>(xBarrierLenght / hx);
	//Количество узлов в второй области по X
	int nX2 = nX - static_cast<int>(xDistancesToBarrier / hx) - n_xDistancesToBarrier_2;
	//Количество узлов в второй области по Y
	int nY2 = static_cast<int>(yBarrierLenght / hy) + correction;

	for (int i = 0; i <= nY; i++) {
		std::cout << i * hy << std::endl;
	}

	cout << "nX1 = " << nX1 << "   nY1 = " << nY1 << "   nY1*h1 = " << nY1 * hy << std::endl;
	cout << "nX2 = " << nX2 << "   nY2 = " << nY2 << "   nY2*h2 = " << nY2 * hy << std::endl;
	cout << "nX3 = " << nX3 << "   nY3 = " << nY3 << "   nY3*h3 = " << nY3 * hy << std::endl;
	cout << nY1 * hy + nY2 * hy + nY3 * hy << endl;
	return 0;
}
