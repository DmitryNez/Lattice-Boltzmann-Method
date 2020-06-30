#pragma once
#include <sstream>
#include <vector>
#include <fstream>
using namespace std;
using vec = vector<vector<vector<double>>>;

void SaveVTKFile(int tStep, vec& p1, vec& u, int Nx, int Ny, int Nz)
{
	stringstream fname;
	fname << "VTK/adv_";
	if (tStep < 10) fname << "0";
	if (tStep < 100) fname << "0";
	if (tStep < 1000) fname << "0";
	if (tStep < 10000) fname << "0";
	if (tStep < 100000) fname << "0";
	if (tStep < 1000000) fname << "0";
	if (tStep < 10000000) fname << "0";
	fname << tStep << ".vtk";
	ofstream vtk_file(fname.str().c_str());
	vtk_file << "# vtk DataFile Version 3.0\n";
	vtk_file << "Phase field with advection\n";
	vtk_file << "ASCII\n";
	vtk_file << "DATASET RECTILINEAR_GRID\nDIMENSIONS " << Nx << " " <<
		Ny << " " << Nz << endl;

	vtk_file << "Z_COORDINATES " << Nz << " double\n";
	for (int i = 1; i <= Nz; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "X_COORDINATES " << Nx << " double\n";
	for (int i = 1; i <= Nx; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "Y_COORDINATES " << Ny << " double\n";
	for (int i = 1; i <= Ny; i++) vtk_file << i << " ";
	vtk_file << endl;


	vtk_file << "POINT_DATA " << Nx * Ny * Nz << endl;

	vtk_file << "SCALARS rho double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";

	for (int k = 1; k <= Nz; k++) {
		for (int j = 1; j <= Ny; j++) {
			for (int i = 1; i <= Nx; i++) {
				vtk_file << p1[k][i][j] << " ";
			}
		}
	}

	/*vtk_file << "SCALARS rho2 double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";

	for (int k = 1; k <= Nz; k++) {
		for (int j = 1; j <= Ny; j++) {
			for (int i = 1; i <= Nx; i++) {
				vtk_file << p2[k][i][j] << " ";
			}
		}
	}*/
	vtk_file << endl;
	vtk_file << "SCALARS u double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";

	for (int k = 1; k <= Nz; k++) {
		for (int j = 1; j <= Ny; j++) {
			for (int i = 1; i <= Nx; i++) {
				vtk_file << u[k][i][j] << " ";
			}
		}
	}

	vtk_file << endl;
	vtk_file.close();

}