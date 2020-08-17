#ifndef PARAM_H // header guard
#define PARAM_H

#include<algorithm>
#include<iostream>
#include<fstream>
#include<iterator>
#include<cmath>
#include<vector>
#include<cstddef>
struct Position
{
	double x;
	double y;
	double z;
};
const int max_particles = 1000;
std::vector<Position> Particle(max_particles);
void initial_config(int write_index);
void sc(int write_index); void fcc(int write_index); void write_config(int write_ind=0);
void bcc(int write_index);
int run_index , N , start_index ;
double Temp, rho ,box;
double potential(int i,double x,double y, double z);
double periodic(double pos);
double randnum();
double presssure();
#endif

