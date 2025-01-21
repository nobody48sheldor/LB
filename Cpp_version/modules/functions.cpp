#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>

namespace func {

void push_1_to_0(int Nx, int Ny, int velocities, std::vector<std::vector<std::vector<double>>> &V_0, std::vector<std::vector<std::vector<double>>> V_1) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			for (int d=0; d<velocities; d++) {
				V_0[i][j][d] = V_1[i][j][d];
			}
		}
	}
}

void write(int Nx, int Ny, int l, std::vector<std::vector<double>> array, int nx, int ny) {
	std::ofstream file("res/res_"+std::to_string(l)+".dat");
	if (!file.is_open()) {
		std::cerr << "Failed to open file for writing." << std::endl;
		return;
	}
	int coef_x = Nx/(nx+1);
	int coef_y = Ny/(ny+1);
	//file << nx << "|" << ny << "|" << "\n";
	file << Nx << "|" << Ny << "|" << "\n";
	//for (int i = 0; i < nx; i++) {
	for (int i = 0; i < Nx; i++) {
		//for (int j = 0; j < ny; j++) {
		for (int j = 0; j < Ny; j++) {
			//file << std::sqrt( (array[coef_x*i][coef_y*j][coef_z*k][0]*array[coef_x*i][coef_y*j][coef_z*k][0]) + (array[coef_x*i][coef_y*j][coef_z*k][1]*array[coef_x*i][coef_y*j][coef_z*k][1]) + (array[coef_x*i][coef_y*j][coef_z*k][2]*array[coef_x*i][coef_y*j][coef_z*k][2]) ) << "/";
			//file << array[nx/2 -1 + coef_x*i][ny/2 -1 + coef_y*j] << "/";
			file << array[i][j] << "/";
		}
	}
	std::cout << "data of time " << l << " written" << std::endl;

	file.close();
}

double sign(double val) {
	if (val >= 0) {
			return 1.0;
	}
	else {
			return -1.0;
	}
}

void curl(int Nx, int Ny, int Nz, std::vector<std::vector<std::vector<std::vector<double>>>> momentum, std::vector<std::vector<std::vector<double>>> &curl_array) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			for (int k=0; k<Nz; k++) {
				if ( (i==0) || (i==(Nx-1)) || (j==0) || (j==(Ny-1)) || (k==0) || (k==(Nz-1)) ) {
					curl_array[i][j][k] = 0.0;
				}
				else {
					curl_array[i][j][k] = sign(( ( momentum[i][j][k+1][0] - momentum[i][j][k-1][0]) - (momentum[i+1][j][k][2] - momentum[i-1][j][k][2]) )) * 0.5 * std::sqrt(
						( ( momentum[i][j+1][k][2] - momentum[i][j-1][k][2]) - (momentum[i][j][k+1][1] - momentum[i][j][k-1][1]) ) * ( ( momentum[i][j+1][k][2] - momentum[i][j-1][k][2]) - (momentum[i][j][k+1][1] - momentum[i][j][k-1][1]) )
						+ 
						( ( momentum[i][j][k+1][0] - momentum[i][j][k-1][0]) - (momentum[i+1][j][k][2] - momentum[i-1][j][k][2]) ) * ( ( momentum[i][j][k+1][0] - momentum[i][j][k-1][0]) - (momentum[i+1][j][k][2] - momentum[i-1][j][k][2]) )
						+
						( ( momentum[i+1][j][k][1] - momentum[i-1][j][k][1]) - (momentum[i][j+1][k][0] - momentum[i][j-1][k][0]) ) * ( ( momentum[i+1][j][k][1] - momentum[i-1][j][k][1]) - (momentum[i][j+1][k][0] - momentum[i][j-1][k][0]) )
					);
				}
			}
		}
	}
}

void velocity(int Nx, int Ny, std::vector<std::vector<std::vector<double>>> momentum, std::vector<std::vector<double>> &velocity_arr) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			velocity_arr[i][j] = std::sqrt( momentum[i][j][0]*momentum[i][j][0] + momentum[i][j][1]*momentum[i][j][1]);
		}
	}
}

void reset(int Nx, int Ny, int velocities, std::vector<std::vector<double>> &rho, std::vector<std::vector<std::vector<double>>> &momentum, std::vector<std::vector<std::vector<double>>> &boundary, std::vector<std::vector<std::vector<double>>> &boundary_inv) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			momentum[i][j][0] = 0.0;
			momentum[i][j][1] = 0.0;
			rho[i][j] = 0.0;
			for (int d=0; d<velocities; d++) {
				boundary[i][j][d] = 0.0;
				boundary_inv[i][j][d] = 0.0;
			}
		}
	}
}

void initialisation(int Nx, int Ny, int velocities, std::vector<std::vector<std::vector<double>>> &V_0, std::vector<std::vector<std::vector<double>>> &V_1, std::vector<std::vector<bool>> &obs, std::vector<std::vector<double>> &rho, std::vector<std::vector<std::vector<double>>> &momentum, std::vector<std::vector<std::vector<double>>> &Veq, std::vector<std::vector<double>> &curl_arr, std::vector<std::vector<double>> &velocity_arr, std::vector<std::vector<std::vector<double>>> &boundary,std::vector<std::vector<std::vector<double>>> &boundary_inv, double u) {

	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_real_distribution<double> distrib(-u, u);

	V_0.resize(Nx);
	V_1.resize(Nx);
	obs.resize(Nx);
	rho.resize(Nx);
	momentum.resize(Nx);
	Veq.resize(Nx);
	curl_arr.resize(Nx);
	velocity_arr.resize(Nx);
	boundary.resize(Nx);
	boundary_inv.resize(Nx);
	for (int i=0; i<Nx; i++) {
		V_0[i].resize(Ny);
		V_1[i].resize(Ny);
		obs[i].resize(Ny);
		rho[i].resize(Ny);
		momentum[i].resize(Ny);
		Veq[i].resize(Ny);
		curl_arr[i].resize(Ny);
		velocity_arr[i].resize(Ny);
		boundary[i].resize(Ny);
		boundary_inv[i].resize(Ny);
		for (int j=0; j<Ny; j++) {
			V_0[i][j].resize(velocities);
			V_1[i][j].resize(velocities);
			Veq[i][j].resize(velocities);
			boundary[i][j].resize(velocities);
			boundary_inv[i][j].resize(velocities);
			obs[i][j] = false;
			rho[i][j] = 0.0;
			momentum[i][j].resize(2);
			momentum[i][j][0] = 0.0; // momentum x
			momentum[i][j][1] = 0.0; // momentum y
			for (int d=0; d<velocities; d++) {
				V_0[i][j][d] = 1.0 + distrib(gen);
				V_1[i][j][d] = 1.0 + distrib(gen);
				Veq[i][j][d] = 1.0;
				boundary[i][j][d] = 1.0;
				boundary_inv[i][j][d] = 1.0;
			}
		}
	}

	std::cout << std::endl;
	std::cout << "~~ V, Veq, obs, rho, momentum, curl_arr, velocity_arr have been initialized ~~" << std::endl;
	std::cout << std::endl;
}

std::vector<double> linspace(double start, double end, int steps_number)
{
	std::vector<double> arr;
	arr.resize(steps_number);
	for (int i=0; i<steps_number; i++) {
		arr[i] = start + ( (double)i ) * (end-start)/( ( (double) steps_number) - 1.0 );
	}
	return arr;
}

}
