

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

#define _USE_MATH_DEFINES


void symm_slab_dker(double, double, double *, int, int);
int slove_waveguide(double *, unsigned int, double **, double ***);
void simulate_MMI(double *, double, double **, double *, unsigned int, unsigned int);
double sudden_approx(double *, double **, unsigned int, unsigned int);




int main(int argc, char **argv)
{
	const double wavelength = 1.55; // microns
	const double k0 = 2.0 * M_PI / wavelength;
	
	const double ncore = 3.5;
//	const double nclad = 3.1693;
	
	const unsigned int N = 1 << 8;
	
	const double L = 8.0; // microns
	const double kL = k0 * L;
	const double inc = kL / (double) N; // step size
	
	const double dker_core = inc * inc * ncore * ncore; // relative permittivity normalised by the step size
//	const double dker_clad = inc * inc * nclad * nclad;
	
	
	FILE *file = fopen("quality vs inc_n.csv", "wt");
	
	double dker_inc = (dker_core - inc * inc) / 100.0;
	
	for (double dker_clad = inc * inc; dker_clad < dker_core; dker_clad += dker_inc)
	{
		
		double *dker = malloc(sizeof(double) * N); // relative permittivity of the waveguide
		symm_slab_dker(dker_core, dker_clad, dker, N, N >> 5);
		double *dker_eff;
		double **field;
		int num_modes = slove_waveguide(dker, N, &dker_eff, &field);
		
		printf("\n%d\n", num_modes);
		
		for (int j = 0; j < num_modes; j++) printf("%f\n", sqrt(dker_eff[j]) / inc);
	//	for (int i = 0; i < N; i++) printf("%f\n", field[0][i]);
		
		
		double *dker_MMI = malloc(sizeof(double) * N);
		symm_slab_dker(dker_core, dker_clad, dker_MMI, N, N >> 1);
		double *dker_eff_MMI;
		double **field_MMI;
		int num_modes_MMI = slove_waveguide(dker_MMI, N, &dker_eff_MMI, &field_MMI);
		
		printf("\n%d\n", num_modes_MMI);
		
		for (int j = 0; j < num_modes_MMI; j++) printf("%f\n", sqrt(dker_eff_MMI[j]) / inc);
	//	for (int i = 0; i < N; i++) printf("%f\n", field[0][i]);
		
		
		double quality = sudden_approx(field[0], field_MMI, N, num_modes_MMI);
		
		printf("\n%f	%f\n", sqrt(dker_clad) / inc, quality);
		fprintf(file, "%f,%f\n", sqrt(dker_clad) / inc, quality);
		
		
		free(dker);
		free(dker_eff);
		
		for (int i = 0; i < num_modes; i++) free(field[i]);
		free(field);
		
		free(dker_MMI);
		free(dker_eff_MMI);
		
		for (int i = 0; i < num_modes_MMI; i++) free(field_MMI[i]);
		free(field_MMI);
		
	}
	
	fclose(file);
	
//	simulate_MMI(field[0], dker_eff[0], field_MMI, dker_eff_MMI, N, num_modes_MMI);
	
	
	return 0;
}




void symm_slab_dker(double dker_core, double dker_clad, double *dker, int N, int N_core)
{
	int N1 = (N - N_core) / 2;
	int N2 = (N + N_core) / 2;
	
	for (int i = 0; i < N1; i++) dker[i] = dker_clad;
	for (int i = N1; i < N2; i++) dker[i] = dker_core;
	for (int i = N2; i < N; i++) dker[i] = dker_clad;
}




int slove_waveguide(double *dker, unsigned int N, double **dker_eff, double ***field)
{
	double dker_min = HUGE_VAL;
	double dker_max = -HUGE_VAL;
	
	for (int i = 0; i < N; i++)
	{
		dker_min = fmin(dker[i], dker_min);
		dker_max = fmax(dker[i], dker_max);
	}
	
	double V = dker_max - dker_min;
	
	double *diag = malloc(sizeof(double) * N);
	double *off_diag = malloc(sizeof(double) * (N - 1));
	
	for (int i = 0; i < N; i++) diag[i] = dker[i] - 2.0;
	
	double gL = exp(-sqrt(dker_max - dker[0]));
	double gR = exp(-sqrt(dker_max - dker[N - 1]));
	
	diag[0] += gL;
	diag[N - 1] += gR;
	
	for (int i = 0; i < N - 1; i++) off_diag[i] = 1.0;
	
	int check = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', N, diag, off_diag, NULL, N);
	if (check) exit(EXIT_FAILURE);
	
	
	int num_modes = 0;
	for (int i = 0; i < N; i++) if (dker_min < diag[i] && diag[i] < dker_max) num_modes++;
	
	*dker_eff = (double *) malloc(sizeof(double) * num_modes);
	*field = (double **) malloc(sizeof(double *) * num_modes);
	
	for (int i = 0; i < num_modes; i++) (*dker_eff)[i] = diag[N - 1 - i];
	
	
	for (int j = 0; j < num_modes; j++)
	{
		double tol;
		do
		{
			for (int i = 0; i < N; i++) diag[i] = dker[i] - 2.0;
			
			gL = exp(-sqrt((*dker_eff)[j] - dker[0]));
			gR = exp(-sqrt((*dker_eff)[j] - dker[N - 1]));
			
			diag[0] += gL;
			diag[N - 1] += gR;
			
			for (int i = 0; i < N - 1; i++) off_diag[i] = 1.0;
			
			check = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', N, diag, off_diag, NULL, N);
			if (check) exit(EXIT_FAILURE);
			
			
			tol = (diag[N - 1 - j] - (*dker_eff)[j]) / V;
			(*dker_eff)[j] = diag[N - 1 - j];
		}
		while(fabs(tol) > 1.0E-10);
	}
	
	
	double *work = malloc(sizeof(double) * N * N);
	
	for (int j = 0; j < num_modes; j++)
	{
		(*field)[j] = malloc(sizeof(double) * N);
		
		for (int i = 0; i < N; i++) diag[i] = dker[i] - 2.0;
		
		gL = exp(-sqrt((*dker_eff)[j] - dker[0]));
		gR = exp(-sqrt((*dker_eff)[j] - dker[N - 1]));
		
		diag[0] += gL;
		diag[N - 1] += gR;
		
		for (int i = 0; i < N - 1; i++) off_diag[i] = 1.0;
		
		check = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', N, diag, off_diag, work, N);
		if (check) exit(EXIT_FAILURE);
		
		
		(*dker_eff)[j] = diag[N - 1 - j];
		
		for (int i = 0; i < N; i++) (*field)[j][i] = work[N * (N - 1 - j) + i];
	}
	
//	for (int j = 0; j < num_modes; j++) printf("%f\n", (*dker_eff)[j]);
//	for (int i = 0; i < N; i++) printf("%f\n", field[1][i]);
	
	
	free(diag);
	free(off_diag);
	free(work);
	
	
	return num_modes;
}




double sudden_approx(double *field_in, double **field_MMI, unsigned int N, unsigned int num_modes_MMI)
{
	double *overlap = malloc(sizeof(double) * num_modes_MMI);
	double quality = 0.0;
	
	printf("\n");
	
	for (int j = 0; j < num_modes_MMI; j++)
	{
		overlap[j] = 0.0;
		
		for (int i = 0; i < N; i++) overlap[j] += field_in[i] * field_MMI[j][i];
		
		printf("%f\n", overlap[j]);
		quality += overlap[j] * overlap[j];
		
		for (int i = 0; i < N; i++) field_MMI[j][i] *= overlap[j];
	}
	
	printf("\n%f\n", quality);
	
	free(overlap);
	
	
	return quality;
}




void simulate_MMI(double *field_in, double dker_eff, double **field_MMI, double *dker_eff_MMI, unsigned int N, unsigned int num_modes_MMI)
{
	sudden_approx(field_in, field_MMI, N, num_modes_MMI);
	
	double **heightmap = malloc(sizeof(double *) * N);
	for (int i = 0; i < N; i++) heightmap[i] = malloc(sizeof(double) * N);
	
	double complex field;
	double intensity;
	double intensity_max = -HUGE_VAL;
	
	for (int k = 0; k < N; k++)
	for (int i = 0; i < N; i++)
	{
		field = 0.0;
		for (int j = 0; j < num_modes_MMI; j++)
		{
			field += field_MMI[j][i] * cexp(I * sqrt(dker_eff_MMI[j]) * k * 5.0);
		}
		
		intensity = field * conj(field);
		intensity_max = fmax(intensity_max, intensity);
		
		heightmap[k][i] = intensity;
	}
	
	for (int k = 0; k < N; k++)
	for (int i = 0; i < N; i++)
	{
		heightmap[k][i] /= intensity_max;
	}
	
	
	unsigned char pix;
	
	FILE *file = fopen("debug 2.pgm", "wb");
	fprintf(file, "P5\n");
	fprintf(file, "%d %d\n", N, N);
	fprintf(file, "%d\n", UCHAR_MAX);
	
	for (int k = 0; k < N; k++)
	for (int i = 0; i < N; i++)
	{
		pix = UCHAR_MAX;
		pix *= heightmap[k][i];
		
		fwrite(&pix, sizeof(pix), 1, file);
	}
	
	fclose(file);
	
	
	for (int i = 0; i < N; i++) free(heightmap[i]);
	free(heightmap);
}

