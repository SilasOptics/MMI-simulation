

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <lapacke.h>

#define _USE_MATH_DEFINES

#define NEUMANN   0x00
#define DIRICHLET 0x01
#define ANALYTIC  0x02


const double wavelength = 1.55; // microns
const double wavenumber = 2.0 * M_PI / wavelength;

const double n_core = 3.5;
const double n_clad = n_core - 2.0;

const double er_core = n_core * n_core;
const double er_clad = n_clad * n_clad;

const double W_in = 1.0; // microns
double W_sim;

const double kW_in = wavenumber * W_in;
double kW_sim;

unsigned int N_in;
unsigned int N_sim;

double h;
double od;
double gamma_analytic;
const double TOL = DBL_EPSILON * (1 << 15);

void waveguide_er(double *, unsigned int);
unsigned int helmholtz(double *, char, double **);

int main(int argc, char **argv)
{
	printf("%e\n", TOL);
	N_in = 1 << 8;
	h  = kW_in / (double) N_in;
	od = 1.0 / (h * h);
	gamma_analytic = exp(-sqrt(er_core - er_clad) * h);
	
	FILE *file = fopen("test.csv", "wt");
	if (file == NULL) exit(EXIT_FAILURE);
	
	for (N_sim = N_in; N_sim <= N_in << 2; N_sim += 2)
	{
		printf("%d\n", N_sim);
		kW_sim = (double) N_sim / (double) N_in * kW_in;
		W_sim  = kW_sim / wavenumber;
		
		double *er = malloc(sizeof(double) * N_sim);
		waveguide_er(er, N_in);
		double *er_eff;
		unsigned int num_modes = helmholtz(er, ANALYTIC, &er_eff);
		fprintf(file, "%f,", W_sim);
		for (unsigned int j = 0; j < num_modes; j++) fprintf(file, "%f,", sqrt(er_eff[j]));
		fprintf(file, "\n");
		
		free(er);
	}
	
	fclose(file);
	
	return 0;
}




void waveguide_er(double *er, unsigned int N_core)
{
	unsigned int N1 = (N_sim - N_core) / 2;
	unsigned int N2 = N1 + N_core;
	unsigned int i;
	for (i = 0;  i < N1;    i++) er[i] = er_clad;
	for (i = N1; i < N2;    i++) er[i] = er_core;
	for (i = N2; i < N_sim; i++) er[i] = er_clad;
}




unsigned int helmholtz(double *er, char boundary, double **er_eff)
{
	double *diag = malloc(sizeof(double) * N_sim);
	double *off_diag = malloc(sizeof(double) * (N_sim - 1));
	unsigned int i;
	for (i = 0; i < N_sim; i++) diag[i] = er[i] - 2.0 * od;
	for (i = 0; i < N_sim - 1; i++) off_diag[i] = od;
	
	double gamma;
	if (boundary == NEUMANN)   gamma = 0.0;
	if (boundary == DIRICHLET) gamma = 1.0;
	if (boundary == ANALYTIC)  gamma = gamma_analytic;
	diag[0] += gamma * od;
	diag[N_sim - 1] += gamma * od;
	
	int check = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', N_sim, diag, off_diag, NULL, N_sim);
	if (check) exit(EXIT_FAILURE);
	
	unsigned int num_modes = 0;
	for (i = 0; i < N_sim; i++) if (er_clad < diag[i] && diag[i] < er_core) num_modes++;
	*er_eff = (double *) malloc(sizeof(double) * num_modes);
	for (i = 0; i < num_modes; i++) (*er_eff)[i] = diag[N_sim - 1 - i];
	
	if (boundary == ANALYTIC)
	{
		for (int j = 0; j < num_modes; j++)
		{
			double inc;
			do
			{
				for (int i = 0; i < N_sim; i++) diag[i] = er[i] - 2.0 * od;
				gamma = exp(-sqrt((*er_eff)[j] - er_clad) * h);
				diag[0] += gamma * od;
				diag[N_sim - 1] += gamma * od;
				for (int i = 0; i < N_sim - 1; i++) off_diag[i] = od;
				
				check = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', N_sim, diag, off_diag, NULL, N_sim);
				if (check) exit(EXIT_FAILURE);
				
				inc = diag[N_sim - 1 - j] - (*er_eff)[j];
				(*er_eff)[j] = diag[N_sim - 1 - j];
			}
			while(fabs(inc) > TOL);
		}
	}
	
	free(diag);
	free(off_diag);
	
	return num_modes;
}


