#pragma once
#include <vector>
constexpr int Nx = 100;
constexpr int Ny = 50;
constexpr int Nz = 50;

using vec = vector<vector<vector<double>>>;

vec N0(Nx + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N1(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N2(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N3(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N4(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N5(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N6(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N7(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N8(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N9(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N10(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N11(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N12(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N13(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N14(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N15(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N16(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N17(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N18(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));

vec N0_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N1_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N2_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N3_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N4_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N5_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N6_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N7_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N8_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N9_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N10_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N11_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N12_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N13_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N14_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N15_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N16_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N17_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec N18_tmp(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));

vec rho(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec Ux(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec Uy(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec Uz(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));
vec U(Nz + 2, (vector<vector<double>>(Nx + 2, vector<double>(Ny + 2))));