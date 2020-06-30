#include <iostream>
#include "VTK.h"
#include "Massive3D.h"

constexpr int Tmax = 10000;
constexpr double g = 0.00001;
constexpr double tau = 1.0;
constexpr double x0 = 25;
constexpr double L = 10;
constexpr double rhomax = 1.0;
constexpr double rhomin = 1.0;
constexpr double PI = 3.1415926535;

using namespace std;


inline void initial_condition(vec & rho,vec& ux,vec& uy, vec& uz, vec& N0, vec& N1, vec& N2, vec& N3, vec& N4, vec& N5, vec& N6,
	vec& N7, vec& N8, vec& N9, vec& N10, vec& N11, vec& N12, vec& N13, vec& N14, vec& N15, vec& N16, vec& N17, vec& N18) {

	for (int k = 1; k <= Nz; k++) {
		for (int i = 1; i <= Nx; i++) {
			for (int j = 1; j <= Ny; j++) {
				ux[k][i][j] = 0.0;
				uy[k][i][j] = 0.0;
				uz[k][i][j] = 0.0;

				//if (i < Nx/2) {
				//	rho[k][i][j] = 1.;
				//}
				//else {
				//	rho[k][i][j] = 0.9;
				//}

					if (k <= x0 - L / 2) {
						rho[k][i][j] = rhomax;
					}
					if (k > x0 - L / 2 && k < L / 2 + x0) {
						rho[k][i][j] = (rhomin + (rhomax - rhomin) * (1 - sin((k - double(x0)) * PI / L)) / 2);
					}
					if (k >= x0 + L / 2) {
						rho[k][i][j] = rhomin;
					}


				N0[k][i][j] = (1. / 3) * rho[k][i][j] * (1 - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N1[k][i][j] = (1. / 18) * rho[k][i][j] * (1 + 3 * ux[k][i][j] + 4.5 * (ux[k][i][j]) * (ux[k][i][j])
					- 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N2[k][i][j] = (1. / 18) * rho[k][i][j] * (1 - 3 * ux[k][i][j] + 4.5 * (ux[k][i][j]) * (ux[k][i][j])
					- 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N3[k][i][j] = (1. / 18) * rho[k][i][j] * (1 + 3 * uy[k][i][j] + 4.5 * (uy[k][i][j]) * (uy[k][i][j])
					- 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N4[k][i][j] = (1. / 18) * rho[k][i][j] * (1 - 3 * uy[k][i][j] + 4.5 * (uy[k][i][j]) * (uy[k][i][j])
					- 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N5[k][i][j] = (1. / 18) * rho[k][i][j] * (1 + 3 * uz[k][i][j] + 4.5 * (uz[k][i][j]) * (uz[k][i][j])
					- 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N6[k][i][j] = (1. / 18) * rho[k][i][j] * (1 - 3 * uz[k][i][j] + 4.5 * (uz[k][i][j]) * (uz[k][i][j])
					- 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N7[k][i][j] = (1. / 36) * rho[k][i][j] * (1 + 3 * (ux[k][i][j] + uy[k][i][j]) + 4.5 * (ux[k][i][j] + uy[k][i][j]) * 
					(ux[k][i][j] + uy[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N8[k][i][j] = (1. / 36) * rho[k][i][j] * (1 + 3 * (-ux[k][i][j] + uy[k][i][j]) + 4.5 * (-ux[k][i][j] + uy[k][i][j]) *
					(-ux[k][i][j] + uy[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N9[k][i][j] = (1. / 36) * rho[k][i][j] * (1 + 3 * (ux[k][i][j] - uy[k][i][j]) + 4.5 * (ux[k][i][j] - uy[k][i][j]) *
					(ux[k][i][j] - uy[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N10[k][i][j] = (1. / 36) * rho[k][i][j] * (1 - 3 * (ux[k][i][j] + uy[k][i][j]) + 4.5 * (ux[k][i][j] + uy[k][i][j]) *
					(ux[k][i][j] + uy[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N11[k][i][j] = (1. / 36) * rho[k][i][j] * (1 + 3 * (ux[k][i][j] + uz[k][i][j]) + 4.5 * (ux[k][i][j] + uz[k][i][j]) *
					(ux[k][i][j] + uz[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N12[k][i][j] = (1. / 36) * rho[k][i][j] * (1 + 3 * (-ux[k][i][j] + uz[k][i][j]) + 4.5 * (-ux[k][i][j] + uz[k][i][j]) *
					(-ux[k][i][j] + uz[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N13[k][i][j] = (1. / 36) * rho[k][i][j] * (1 + 3 * (ux[k][i][j] - uz[k][i][j]) + 4.5 * (ux[k][i][j] - uz[k][i][j]) *
					(ux[k][i][j] - uz[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N14[k][i][j] = (1. / 36) * rho[k][i][j] * (1 - 3 * (ux[k][i][j] + uz[k][i][j]) + 4.5 * (ux[k][i][j] + uz[k][i][j]) *
					(ux[k][i][j] + uz[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N15[k][i][j] = (1. / 36) * rho[k][i][j] * (1 + 3 * (uy[k][i][j] + uz[k][i][j]) + 4.5 * (uy[k][i][j] + uz[k][i][j]) *
					(uy[k][i][j] + uz[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N16[k][i][j] = (1. / 36) * rho[k][i][j] * (1 + 3 * (-uy[k][i][j] + uz[k][i][j]) + 4.5 * (-uy[k][i][j] + uz[k][i][j]) *
					(-uy[k][i][j] + uz[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N17[k][i][j] = (1. / 36) * rho[k][i][j] * (1 + 3 * (uy[k][i][j] - uz[k][i][j]) + 4.5 * (uy[k][i][j] - uz[k][i][j]) *
					(uy[k][i][j] - uz[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

				N18[k][i][j] = (1. / 36) * rho[k][i][j] * (1 - 3 * (uy[k][i][j] + uz[k][i][j]) + 4.5 * (uy[k][i][j] + uz[k][i][j]) *
					(uy[k][i][j] + uz[k][i][j]) - 1.5 * (ux[k][i][j] * ux[k][i][j] + uy[k][i][j] * uy[k][i][j] + uz[k][i][j] * uz[k][i][j]));

			}
		}
	}
}

inline void equilibrium_func(double ux,double uy, double uz,double rho, vector<double>& N_eq) {
	N_eq[0] = (1. / 3) * rho * (1 - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[1] = (1. / 18) * rho * (1 + 3 * ux + 4.5 * ux * ux - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[2] = (1. / 18) * rho * (1 - 3 * ux + 4.5 * ux * ux - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[3] = (1. / 18) * rho * (1 + 3 * uy + 4.5 * uy * uy - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[4] = (1. / 18) * rho * (1 - 3 * uy + 4.5 * uy * uy - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[5] = (1. / 18) * rho * (1 + 3 * uz + 4.5 * uz * uz - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[6] = (1. / 18) * rho * (1 - 3 * uz + 4.5 * uz * uz - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[7] = (1. / 36) * rho * (1 + 3 * (ux + uy) + 4.5 * (ux + uy) * (ux + uy) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[8] = (1. / 36) * rho * (1 + 3 * (-ux + uy) + 4.5 * (-ux + uy) * (-ux + uy) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[9] = (1. / 36) * rho * (1 + 3 * (ux - uy) + 4.5 * (ux - uy) * (ux - uy) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[10] = (1. / 36) * rho * (1 - 3 * (ux + uy) + 4.5 * (ux + uy) * (ux + uy) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[11] = (1. / 36) * rho * (1 + 3 * (ux + uz) + 4.5 * (ux + uz) * (ux + uz) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[12] = (1. / 36) * rho * (1 + 3 * (-ux + uz) + 4.5 * (-ux + uz) * (-ux + uz) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[13] = (1. / 36) * rho * (1 + 3 * (ux - uz) + 4.5 * (ux - uz) * (ux - uz) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[14] = (1. / 36) * rho * (1 - 3 * (ux + uz) + 4.5 * (ux + uz) * (ux + uz) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[15] = (1. / 36) * rho * (1 + 3 * (uy + uz) + 4.5 * (uy + uz) * (uy + uz) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[16] = (1. / 36) * rho * (1 + 3 * (-uy + uz) + 4.5 * (-uy + uz) * (-uy + uz) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[17] = (1. / 36) * rho * (1 + 3 * (uy - uz) + 4.5 * (uy - uz) * (uy - uz) - 1.5 * (ux * ux + uy * uy + uz * uz));

	N_eq[18] = (1. / 36) * rho * (1 - 3 * (uy + uz) + 4.5 * (uy + uz) * (uy + uz) - 1.5 * (ux * ux + uy * uy + uz * uz));
}

inline void boundary_conditions(vec& N0, vec& N1, vec& N2, vec& N3, vec& N4, vec& N5, vec& N6,
	vec& N7, vec& N8, vec& N9, vec& N10, vec& N11, vec& N12, vec& N13, vec& N14, vec& N15, vec& N16, vec& N17, vec& N18) {

	/*periodic boundary condition*/
	for (int k = 1; k <= Nz; ++k) {
		for (int j = 1; j <= Ny; ++j) {
			N2[k][Nx + 1][j] = N2[k][1][j];
			N12[k][Nx + 1][j] = N12[k][1][j];
			N14[k][Nx + 1][j] = N14[k][1][j];
			N10[k][Nx + 1][j] = N10[k][1][j];
			N8[k][Nx + 1][j] = N8[k][1][j];

			N7[k][0][j] = N7[k][Nx][j];
			N1[k][0][j] = N1[k][Nx][j];
			N9[k][0][j] = N9[k][Nx][j];
			N11[k][0][j] = N11[k][Nx][j];
			N13[k][0][j] = N13[k][Nx][j];
		}
	}
	
	for (int k = 1; k <= Nz; k++) {
		for (int i = 1; i <= Nx; i++) {
			N8[k][i][0] = N8[k][i][Ny];
			N3[k][i][0] = N3[k][i][Ny];
			N7[k][i][0] = N7[k][i][Ny];
			N15[k][i][0] = N15[k][i][Ny];
			N17[k][i][0] = N17[k][i][Ny];

			N10[k][i][Ny + 1] = N10[k][i][1];
			N18[k][i][Ny + 1] = N18[k][i][1];
			N16[k][i][Ny + 1] = N16[k][i][1];
			N9[k][i][Ny + 1] = N9[k][i][1];
			N4[k][i][Ny + 1] = N4[k][i][1];
		}
	}

	/*for (int i = 1; i <= Nx; i++) {
		for (int j = 1; j <= Ny; j++) {
			N6[Nz + 1][i][j] = N6[1][i][j];
			N18[Nz + 1][i][j] = N18[1][i][j];
			N17[Nz + 1][i][j] = N17[1][i][j];
			N14[Nz + 1][i][j] = N14[1][i][j];
			N13[Nz + 1][i][j] = N13[1][i][j];

			N5[0][i][j] = N5[Nz][i][j];
			N15[0][i][j] = N15[Nz][i][j];
			N16[0][i][j] = N16[Nz][i][j];
			N11[0][i][j] = N11[Nz][i][j];
			N12[0][i][j] = N12[Nz][i][j];
		}
	}*/


	for (int k = 1; k <= Nz; k++) {

		N7[k][0][0] = N7[k][Nx][Ny];
		N9[k][0][Ny + 1] = N9[k][Nx][1];

		N8[k][Nx + 1][0] = N8[k][1][Ny];
		N10[k][Nx + 1][Ny + 1] = N10[k][1][1];

	}

	/*for (int i = 0; i <= Nx + 1; i++) {

		N18[Nz + 1][i][Ny + 1] = N18[1][i][1];
		N17[Nz + 1][i][0] = N17[1][i][Ny];

		N15[0][i][0] = N15[Nz][i][Ny];
		N16[0][i][Ny + 1] = N16[Nz][i][1];

	}*/

	/*for (int j = 0; j <= Ny + 1; j++) {

		N12[0][Nx + 1][j] = N12[Nz][1][j];
		N11[0][0][j] = N11[Nz][Nx][j];

		N14[Nz + 1][Nx + 1][j] = N14[1][1][j];
		N13[Nz + 1][0][j] = N13[1][Nx][j];

	}*/

	/*bounce-back boundary condition z = Nz + 1 and z = 0*/
	for (int i = 1; i <= Nx; i++) {
		for (int j = 1; j <= Ny; j++) {
			N5[0][i][j] = N6[1][i][j];
			N15[0][i][j - 1] = N18[1][i][j];
			N16[0][i][j + 1] = N17[1][i][j];
			N11[0][i - 1][j] = N14[1][i][j];
			N12[0][i + 1][j] = N13[1][i][j];

			N6[Nz + 1][i][j] = N5[Nz][i][j];
			N18[Nz + 1][i][j + 1] = N15[Nz][i][j];
			N17[Nz + 1][i][j - 1] = N16[Nz][i][j];
			N14[Nz + 1][i + 1][j] = N11[Nz][i][j];
			N13[Nz + 1][i - 1][j] = N12[Nz][i][j];
		}
	}


}

inline void solve() {

	initial_condition(rho, Ux, Uy, Uz, N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18);
	int t = 1;
	double impx = 0, impy = 0, impz = 0;
	double m = 0;
	while (t <= Tmax) {

		boundary_conditions(N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18);

		/*---------Streaming step---------*/
		for (int k = 1; k <= Nz; k++) {
			for (int i = 1; i <= Nx; i++) {
				for (int j = 1; j <= Ny; j++) {
					N1_tmp[k][i][j] = N1[k][i - 1][j];
					N2_tmp[k][i][j] = N2[k][i + 1][j];
					N3_tmp[k][i][j] = N3[k][i][j - 1];
					N4_tmp[k][i][j] = N4[k][i][j + 1];
					N5_tmp[k][i][j] = N5[k - 1][i][j];
					N6_tmp[k][i][j] = N6[k + 1][i][j];
					N7_tmp[k][i][j] = N7[k][i - 1][j - 1];
					N8_tmp[k][i][j] = N8[k][i + 1][j - 1];
					N9_tmp[k][i][j] = N9[k][i - 1][j + 1];
					N10_tmp[k][i][j] = N10[k][i + 1][j + 1];
					N11_tmp[k][i][j] = N11[k - 1][i - 1][j];
					N12_tmp[k][i][j] = N12[k - 1][i + 1][j];
					N13_tmp[k][i][j] = N13[k + 1][i - 1][j];
					N14_tmp[k][i][j] = N14[k + 1][i + 1][j];
					N15_tmp[k][i][j] = N15[k - 1][i][j - 1];
					N16_tmp[k][i][j] = N16[k - 1][i][j + 1];
					N17_tmp[k][i][j] = N17[k + 1][i][j - 1];
					N18_tmp[k][i][j] = N18[k + 1][i][j + 1];
				}
			}
		}
	
		N1_tmp.swap(N1);
		N2_tmp.swap(N2);
		N3_tmp.swap(N3);
		N4_tmp.swap(N4);
		N5_tmp.swap(N5);
		N6_tmp.swap(N6);
		N7_tmp.swap(N7);
		N8_tmp.swap(N8);
		N9_tmp.swap(N9);
		N10_tmp.swap(N10);
		N11_tmp.swap(N11);
		N12_tmp.swap(N12);
		N13_tmp.swap(N13);
		N14_tmp.swap(N14);
		N15_tmp.swap(N15);
		N16_tmp.swap(N16); 
		N17_tmp.swap(N17);
		N18_tmp.swap(N18);
		/*------------------*/
		for (int k = 1; k <= Nz; k++) {
			for (int i = 1; i <= Nx; i++) {
				for (int j = 1; j <= Ny; j++) {

					rho[k][i][j] = N0[k][i][j] + N1[k][i][j] + N2[k][i][j] + N3[k][i][j] + N4[k][i][j] + N5[k][i][j] + N6[k][i][j] +
						N7[k][i][j] + N8[k][i][j] + N9[k][i][j] + N10[k][i][j] + N11[k][i][j] + N12[k][i][j] + N13[k][i][j] + N14[k][i][j] +
						N15[k][i][j] + N16[k][i][j] + N17[k][i][j] + N18[k][i][j];

					Ux[k][i][j] = (N1[k][i][j] + N11[k][i][j] + N13[k][i][j] + N9[k][i][j] + N7[k][i][j] - N2[k][i][j] - N8[k][i][j] -
						N14[k][i][j] - N10[k][i][j] - N12[k][i][j]) / rho[k][i][j];

					Uy[k][i][j] = (N3[k][i][j] + N15[k][i][j] + N17[k][i][j] + N8[k][i][j] + N7[k][i][j] - N18[k][i][j] - N4[k][i][j] -
						N16[k][i][j] - N10[k][i][j] - N9[k][i][j]) / rho[k][i][j];

					Uz[k][i][j] = (N6[k][i][j] + N14[k][i][j] + N13[k][i][j] + N18[k][i][j] + N17[k][i][j] - N5[k][i][j] - N11[k][i][j] -
						N12[k][i][j] - N15[k][i][j] - N16[k][i][j]) / rho[k][i][j];

					m += rho[k][i][j];

					impx += rho[k][i][j] * Ux[k][i][j];
					impy += rho[k][i][j] * Uy[k][i][j];
					impz += rho[k][i][j] * Uz[k][i][j];


				}
			}
		}
		/*-----Collision step------*/
		for (int k = 1; k <= Nz; k++) {
			for (int i = 1; i <= Nx; i++) {
				for (int j = 1; j <= Ny; j++) {

					vector<double> Neq(19);

					equilibrium_func(Ux[k][i][j], Uy[k][i][j], Uz[k][i][j], rho[k][i][j], Neq);

					Ux[k][i][j] += g/2;

					vector<double> Neq_d(19);

					equilibrium_func(Ux[k][i][j] + g/2, Uy[k][i][j], Uz[k][i][j], rho[k][i][j], Neq_d);

					N0[k][i][j] += ((Neq[0] - N0[k][i][j]) / tau) + Neq_d[0] - Neq[0];
					N1[k][i][j] += ((Neq[1] - N1[k][i][j]) / tau) + Neq_d[1] - Neq[1];
					N2[k][i][j] += ((Neq[2] - N2[k][i][j]) / tau) + Neq_d[2] - Neq[2];
					N3[k][i][j] += ((Neq[3] - N3[k][i][j]) / tau) + Neq_d[3] - Neq[3];
					N4[k][i][j] += ((Neq[4] - N4[k][i][j]) / tau) + Neq_d[4] - Neq[4];
					N5[k][i][j] += ((Neq[5] - N5[k][i][j]) / tau) + Neq_d[5] - Neq[5];
					N6[k][i][j] += ((Neq[6] - N6[k][i][j]) / tau) + Neq_d[6] - Neq[6];
					N7[k][i][j] += ((Neq[7] - N7[k][i][j]) / tau) + Neq_d[7] - Neq[7];
					N8[k][i][j] += ((Neq[8] - N8[k][i][j]) / tau) + Neq_d[8] - Neq[8];
					N9[k][i][j] += ((Neq[9] - N9[k][i][j]) / tau) + Neq_d[9] - Neq[9];
					N10[k][i][j] += ((Neq[10] - N10[k][i][j]) / tau) + Neq_d[10] - Neq[10];
					N11[k][i][j] += ((Neq[11] - N11[k][i][j]) / tau) + Neq_d[11] - Neq[11];
					N12[k][i][j] += ((Neq[12] - N12[k][i][j]) / tau) + Neq_d[12] - Neq[12];
					N13[k][i][j] += ((Neq[13] - N13[k][i][j]) / tau) + Neq_d[13] - Neq[13];
					N14[k][i][j] += ((Neq[14] - N14[k][i][j]) / tau) + Neq_d[14] - Neq[14];
					N15[k][i][j] += ((Neq[15] - N15[k][i][j]) / tau) + Neq_d[15] - Neq[15];
					N16[k][i][j] += ((Neq[16] - N16[k][i][j]) / tau) + Neq_d[16] - Neq[16];
					N17[k][i][j] += ((Neq[17] - N17[k][i][j]) / tau) + Neq_d[17] - Neq[17];
					N18[k][i][j] += ((Neq[18] - N18[k][i][j]) / tau) + Neq_d[18] - Neq[18];
					
				}
			}
		}
		/*--------------*/
		for (int k = 1; k <= Nz; k++) {
			for (int i = 1; i <= Nx; i++) {
				for (int j = 1; j <= Ny; j++) {
					U[k][i][j] = sqrt(Ux[k][i][j] * Ux[k][i][j] + Uy[k][i][j] * Uy[k][i][j] + Uz[k][i][j] * Uz[k][i][j]);
				}
			}
		}
		if (t % 20 == 0) {
			SaveVTKFile(t, rho, U, Nx, Ny, Nz);
			cout << "t = " << t << "; m = " << m << " impx = " << impx << " impy = " << impy << " impz = " << impz << endl;
		}


		t++;
		m = 0;
		impx = 0;
		impy = 0;
		impz = 0;

	}
}

int main() {
	solve();
	return 0;
}