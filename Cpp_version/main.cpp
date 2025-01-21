#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
#include "modules/functions.cpp"
#include "modules/obstacle_setup.cpp"

void init_conditions(int Nx, int Ny, int velocities, std::vector<std::vector<std::vector<double>>> &V_0, std::vector<std::vector<std::vector<double>>> &V_1) {

	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			V_0[i][j][3] = 6.1; // vitesse vers la droite
			V_1[i][j][3] = 6.1; // vitesse vers la droite
		}
	}

	std::cout << std::endl;
	std::cout << "~~ initial conditions have been applied to V ~~" << std::endl;
	std::cout << std::endl;
}

void border(int Nx, int Ny, std::vector<std::vector<std::vector<double>>> V_0, std::vector<std::vector<std::vector<double>>> &V_1) {
	// i fixé
	for (int j=0; j<Ny; j++) {
		V_1[Nx-1][j][6] = V_0[Nx-2][j][6];
		V_1[Nx-1][j][7] = V_0[Nx-2][j][7];
		V_1[Nx-1][j][8] = V_0[Nx-2][j][8];
	}

	for (int j=0; j<Ny; j++) {
		V_1[0][j][2] = V_0[1][j][2];
		V_1[0][j][3] = V_0[1][j][3];
		V_1[0][j][4] = V_0[1][j][4];
	}

	// j fixé
	for (int i=0; i<Nx; i++) {
		V_1[i][0][8] = V_0[i][1][8];
		V_1[i][0][1] = V_0[i][1][1];
		V_1[i][0][2] = V_0[i][1][2];
	}

	for (int i=0; i<Nx; i++) {
		V_1[i][Ny-1][6] = V_0[i][Ny-2][6];
		V_1[i][Ny-1][5] = V_0[i][Ny-2][5];
		V_1[i][Ny-1][4] = V_0[i][Ny-2][4];

	}
}

void boundary_compute(int Nx, int Ny, int velocities, std::vector<std::vector<std::vector<double>>> V_0, std::vector<std::vector<bool>> obs, std::vector<std::vector<std::vector<double>>> &boundary) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			if (obs[i][j]) {
				for (int d=0; d<velocities; d++) {
					boundary[i][j][d] = V_0[i][j][d];
				}
			}
		}
	}
}

void boundary_inv_compute(int Nx, int Ny, std::vector<std::vector<std::vector<double>>> boundary, std::vector<std::vector<std::vector<double>>> &boundary_inv, std::vector<std::vector<bool>> obs) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			if (obs[i][j]) {
				// on inverse les vitesse de boundary (rebond)
				boundary_inv[i][j][0] = boundary[i][j][0];
				boundary_inv[i][j][1] = boundary[i][j][5];
				boundary_inv[i][j][2] = boundary[i][j][6];
				boundary_inv[i][j][3] = boundary[i][j][7];
				boundary_inv[i][j][4] = boundary[i][j][8];
				boundary_inv[i][j][5] = boundary[i][j][1];
				boundary_inv[i][j][6] = boundary[i][j][2];
				boundary_inv[i][j][7] = boundary[i][j][3];
				boundary_inv[i][j][8] = boundary[i][j][4];
			}
		}
	}
}

void bounce_obstacle_V(int Nx, int Ny, int velocities, std::vector<std::vector<std::vector<double>>> &V_1, std::vector<std::vector<bool>> obs, std::vector<std::vector<std::vector<double>>> boundary_inv) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			if (obs[i][j]) {
				for (int d=0;  d<velocities; d++) {
					V_1[i][j][d] = boundary_inv[i][j][d];
				}
			}
		}
	}
}

void reset_momentum_obs(int Nx, int Ny, std::vector<std::vector<bool>> obs, std::vector<std::vector<std::vector<double>>> &momentum) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			if (obs[i][j]) {
				momentum[i][j][0] = 0.0;
				momentum[i][j][1] = 0.0;
			}
		}
	}
}


void streaming(int Nx, int Ny, int velicities, std::vector<std::vector<double>> directions, std::vector<std::vector<std::vector<double>>> V_0, std::vector<std::vector<std::vector<double>>> &V_1, int axis) {

	if (axis == 0) {
	for (int d=0; d<velicities; d++) {
		int cx = directions[d][0];
		for (int i=0; i<Nx; i++) {
			for (int j=0; j<Ny; j++) {
				int i_shift = ((i-cx)%Nx + Nx) % Nx;
				V_1[i][j][d] = V_0[i_shift][j][d]; // axis 0
			}
		}
	}
	}
	
	if (axis == 1) {

	for (int d=0; d<velicities; d++) {
		int cy= directions[d][1];
		for (int i=0; i<Nx; i++) {
			for (int j=0; j<Ny; j++) {
				int j_shift = ((j-cy)%Ny + Ny) % Ny;
				V_1[i][j][d] = V_0[i][j_shift][d]; // axis 1
			}
		}
	}
	}


	if (axis == -1) {

		for (int d=0; d<velicities; d++) {
		
			int cx = directions[d][0];
			int cy = directions[d][1];

			std::cout << cx << " <- cx "<< std::endl;
			std::cout << cy << " <- cy "<< std::endl;

			for (int i=0; i<Nx; i++) {
				for (int j=0; j<Ny; j++) {
					int i_shift = ((i-cx)%Nx + Nx) % Nx;
					V_1[i][j][d] = V_0[i_shift][j][d]; // pour x
					//std::cout << i << " | " << i_shift << " <- i - ishift  | cx -> " << cx << std::endl;
				}
			}

			for (int i=0; i<Nx; i++) {
				for (int j=0; j<Ny; j++) {
					int j_shift = ((j-cy)%Ny + Ny) % Ny;
					V_1[i][j][d] = V_0[i][j_shift][d]; // pour y
					//std::cout << j << " | " << j_shift << " <- j - jshift  | cy -> " << cy << std::endl;
				}
			}

		}
	}
}


void sum_speed_and_momentum(int Nx, int Ny, int velocities, std::vector<std::vector<double>> directions, std::vector<std::vector<std::vector<double>>> V_0, std::vector<std::vector<double>> &rho, std::vector<std::vector<std::vector<double>>> &momentum) {

	// reset first to 0
	
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			rho[i][j] = 0.0;
		}
	}

	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			for (int d=0; d<velocities; d++) {
				momentum[i][j][0] = 0.0;
				momentum[i][j][1] = 0.0;
			}
		}
	}

	// sum values (to 0.0)

	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			for (int d=0; d<velocities; d++) {
				rho[i][j] += V_0[i][j][d];
			}
		}
	}

	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			for (int d=0; d<velocities; d++) {
				momentum[i][j][0] += (directions[d][0]*V_0[i][j][d]) / rho[i][j];
				momentum[i][j][1] += (directions[d][1]*V_0[i][j][d]) / rho[i][j];
			}
		}
	}
}

void compute_Veq(int Nx, int Ny, int velocities, std::vector<std::vector<double>> directions, std::vector<double> weights, std::vector<std::vector<double>> rho, std::vector<std::vector<std::vector<double>>> momentum, std::vector<std::vector<std::vector<double>>> &Veq) {

	// reset
	//
	for (int d=0; d<velocities; d++) {
		for (int i=0; i<Nx; i++) {
			for (int j=0; j<Ny; j++) {
				Veq[i][j][d] = 0.0;
			}
		} 
	}

	// compute Veq
	//
	for (int d=0; d<velocities; d++) {


		for (int i=0; i<Nx; i++) {
			for (int j=0; j<Ny; j++) {
				Veq[i][j][d] = rho[i][j] * weights[d] * ( 1.0 + ( 3.0 * ( (directions[d][0] * momentum[i][j][0]) + (directions[d][1] * momentum[i][j][1]) ) ) + ( 4.5 * ( (directions[d][0] * momentum[i][j][0]) + (directions[d][1] * momentum[i][j][1]) )*( (directions[d][0] * momentum[i][j][0]) + (directions[d][1] * momentum[i][j][1]) ) ) - ( 1.5 * ( (momentum[i][j][0]*momentum[i][j][0]) + (momentum[i][j][1]*momentum[i][j][1]) ) ));
			}
		}

	}

}


void compute_next_V(int Nx, int Ny, int velocities, std::vector<std::vector<std::vector<double>>> V_0, std::vector<std::vector<std::vector<double>>> &V_1, std::vector<std::vector<std::vector<double>>> Veq, double tau) {
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			for (int d=0; d<velocities; d++) {
				V_1[i][j][d] = V_0[i][j][d] - ( (V_0[i][j][d] - Veq[i][j][d]) / tau);
			}
		}
	}
}





int main() {
	int N = 799;

	int Nx = N;
	int Ny = N;

	int Nt = 10001;

	constexpr double L = 0.5; // [-L,L] 1m // just to give it a number, does not influence the result

	// double tau = 1.5; // viscosite dynamique
	constexpr double tau = 10.0; // viscosite dynamique
	constexpr double u = 0.01; // random speed deviation at the initial time to avoid having weird things sometime
	int velocities = 9; // nombre de vitesses discretises

	std::vector<double> x = func::linspace(-L, L, Nx);
	std::vector<double> y = func::linspace(-L, L, Ny);

	std::vector<std::vector<double>> directions = {
	{ 0.0,  0.0}, // 0
	{ 0.0,  1.0}, // 1
	{ 1.0,  1.0}, // 2
	{ 1.0,  0.0}, // 3
	{ 1.0, -1.0}, // 4
	{ 0.0, -1.0}, // 5
	{-1.0, -1.0}, // 6
	{-1.0,  0.0}, // 7
	{-1.0,  1.0}, // 8
	};

	std::vector<double> weights = {
	(4.0/9.0),	// 0
	(1.0/9.0),	// 1
	(1.0/36.0),	// 2
	(1.0/9.0),	// 3
	(1.0/36.0),	// 4
	(1.0/9.0),	// 5
	(1.0/36.0),	// 6
	(1.0/9.0),	// 7
	(1.0/36.0),	// 8
	};

	// init
	
	std::vector<std::vector<std::vector<double>>> V_0; // Nx, Ny, direction
	std::vector<std::vector<std::vector<double>>> V_1; // Nx, Ny, direction
	std::vector<std::vector<bool>> obs; // Nx, Ny
	std::vector<std::vector<double>> rho; // Nx, Ny
	std::vector<std::vector<std::vector<double>>> momentum; // Nx, Ny, (x,y,z)
	std::vector<std::vector<std::vector<double>>> Veq; // Nx, Ny, directions
	std::vector<std::vector<std::vector<double>>> boundary; // Nx, Ny, direction
	std::vector<std::vector<std::vector<double>>> boundary_inv; // Nx, Ny, direction
	std::vector<std::vector<double>> curl_arr; // Nx, Ny
	std::vector<std::vector<double>> velocity_arr; // Nx, Ny
	

	func::initialisation(Nx, Ny, velocities, V_0, V_1, obs, rho, momentum, Veq, curl_arr, velocity_arr, boundary, boundary_inv, u); // resize correctly the differents arrays
	init_conditions(Nx, Ny, velocities, V_0, V_1); // apply the inital conditions to V

	obstSetup::obstacle_func_ball(Nx, Ny, 0.1, -0.3, 0.0, x, y, obs); // obstacle est une boule de rayon 1

	for (int l=1; l<(Nt); l++) {
		std::cout << l << " / " << Nt << " ~~ ";
		// func::reset(Nx, Ny, velocities, rho, momentum, boundary, boundary_inv); // reset les array rho et momentum
		
		// boundaries conditions
		//
		border(Nx, Ny, V_0, V_1);
		func::push_1_to_0(Nx, Ny, velocities, V_0, V_1);


		// roll velocities
		//
		streaming(Nx, Ny, velocities, directions,  V_0, V_1, -1);
		func::push_1_to_0(Nx, Ny, velocities, V_0, V_1);


		// obstacle
		//
		boundary_compute(Nx, Ny, velocities, V_0, obs, boundary); // copy velocities in the obstacle to the variable "boundary"
		boundary_inv_compute(Nx, Ny, boundary, boundary_inv, obs); // on inverse les vitesses copiés de "boundary" vers "boundary_inv" toujours sans modifier V


		// compute rho and momentum
		//
		sum_speed_and_momentum(Nx, Ny, velocities, directions, V_0, rho, momentum); // compute rho and momentum


		// inverse velocities in obstacle
		//
		bounce_obstacle_V(Nx, Ny, velocities, V_1, obs, boundary_inv); // on remplie les vitesses inversés copiées dans boundary dans V
		func::push_1_to_0(Nx, Ny, velocities, V_0, V_1);
		//reset_momentum_obs(Nx, Ny, obs, momentum); // put 0 momentum in obstacle

		// collisions
		//
		compute_Veq(Nx, Ny, velocities, directions, weights, rho, momentum, Veq);
		compute_next_V(Nx, Ny, velocities, V_0, V_1, Veq, tau);
		func::push_1_to_0(Nx, Ny, velocities, V_0, V_1);

		//bounce_obstacle_V(Nx, Ny, velocities, V_1, obs, boundary_inv); // on remplie les vitesses inversés copiées dans boundary dans V
		// func::push_1_to_0(Nx, Ny, velocities, V_0, V_1);
		

		std::cout << "momentum_x = " << momentum[Nx/2][Ny/2][0] << " | momentum_y = " << momentum[Nx/2][Ny/2][1] << std::endl;

		// compute the speed
		//
		//func::curl(Nx, Ny, Nz, momentum, curl_arr);
		if ( ((l-1)%20) == 0 ) {
			func::velocity(Nx, Ny, momentum, velocity_arr);
			func::write(Nx, Ny, l, velocity_arr, 20, 20);
		}
	}
	return 0;
}
