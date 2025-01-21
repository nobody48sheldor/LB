#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>

namespace obstSetup {

double distance(double x0, double y0, double x1, double y1) {
	double dist = std::sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
	return dist;
}

void obstacle_func_ball(int Nx, int Ny, double radius, double x_0,double y_0, std::vector<double> x, std::vector<double> y, std::vector<std::vector<bool>> &obs) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			if (distance(x[i], y[j], x_0, y_0) < radius) {
				obs[i][j] = true;
			}
		}
	}
}

}
