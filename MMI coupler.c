/*
 * MMI coupler.c
 * 
 * Copyright 2025 Foundation <Foundation@DESKTOP-MLC085A>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

// gcc -Wall -o "%e" "%f" -IC:\mingw64\include\openblas -LC:\mingw64\lib -lopenblas -llapack -lblas

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <lapacke.h>

#define _USE_MATH_DEFINES


const double wavelength = 1.55; // microns
const double wavenumber = 2.0 * M_PI / wavelength;

const double n_core = 3.5;
const double n_clad = 3.1693;
//const double n_clad = n_core - 0.0625;

const unsigned int N_sim = 1<<9;
const unsigned int N_in  = N_sim>>4;
const unsigned int N_MMI = N_sim>>1;
const unsigned int N_sep = N_sim>>5;
//unsigned int N_in;

const double uwidth_sim = (double) N_sim / (double) N_MMI;
const double uwidth_in  = (double) N_in  / (double) N_MMI;
const double uwidth_sep = (double) N_sep / (double) N_MMI;
const double uwidth_MMI = 1.0;
//double uwidth_in;

const double width_in = 0.5; // microns
//const double width_MMI = 3.0; // microns
const double width_MMI = width_in  * uwidth_MMI / uwidth_in;
//const double width_in = width_MMI * uwidth_in / uwidth_MMI;
const double width_sim = width_MMI * uwidth_sim / uwidth_MMI;
const double width_sep = width_MMI * uwidth_sep / uwidth_MMI;
//double width_in;

const double kwidth_in = wavenumber * width_in;
const double kwidth_MMI = wavenumber * width_MMI;
const double kwidth_sim = wavenumber * width_sim;
const double kwidth_sep = wavenumber * width_sep;
//double kwidth_in;

const double kx_inc = kwidth_MMI / (double) N_MMI;
const double dker_core = kx_inc * kx_inc * n_core * n_core;
const double dker_clad = kx_inc * kx_inc * n_clad * n_clad;

//const double kz_inc = kx_inc * n_core * kwidth_MMI / 2.0 / M_PI * 0.9;
//const double kz_inc = kx_inc * n_core * kwidth_MMI / 2.0 / M_PI * 3.0;
//const double kz_inc = kx_inc * n_core * kwidth_MMI / 2.0 / M_PI * 32.0;
//const double kz_inc = kx_inc * n_core * kwidth_MMI / 2.0 / M_PI * 4000.0;

const double kz_inc = kx_inc * n_core * kwidth_MMI / 2.0 / M_PI * 0.7;
//const double kz_inc = kx_inc * n_core * kwidth_MMI / 2.0 / M_PI * 2.3;
//const double kz_inc = kx_inc * n_core * kwidth_MMI / 2.0 / M_PI * 25.0;
//const double kz_inc = kx_inc * n_core * kwidth_MMI / 2.0 / M_PI * 3000.0;

void symm_slab_dker(double *, int);
int slove_waveguide(double *, double **, double ***);
void field_to_csv(unsigned int, double **);
double simulate_MMI(unsigned int, double *, double **, double *);
void offset_slab_dker(double *, unsigned int, unsigned int);
void coupler_dker(double *, double *, unsigned int);
void coupler3_dker(double *, double *, unsigned int);

void colour_gradient();


int main(int argc, char **argv)
{
	colour_gradient();
	
	printf("n_core = %f\n", n_core);
	printf("n_clad = %f\n", n_clad);
	printf("width_in  = %f microns\n", width_in);
	printf("width_MMI = %f microns\n", width_MMI);
	printf("width_sim = %f microns\n", width_sim);
	printf("width_sep = %f microns\n", width_sep);
	printf("approx reimage length = %f microns\n", kz_inc / wavenumber * N_MMI);
	
	double *dker_in  = malloc(sizeof(double) * N_sim);
	double *dker_MMI = malloc(sizeof(double) * N_sim);
	
//	symm_slab_dker(dker_in, N_in);
//	offset_slab_dker(dker_in, N_in, (N_sim-N_MMI)/2 + N_MMI/3);
	coupler_dker(dker_in, dker_MMI, N_sep);
	double *dker_eff_in;
	double **field_in;
	int num_modes_in = slove_waveguide(dker_in, &dker_eff_in, &field_in);
	printf("num_modes_in = %d\n", num_modes_in);
	printf("n_eff_in\n");
	for (int j = 0; j < num_modes_in; j++) printf("%f\n", sqrt(dker_eff_in[j]) / kx_inc);
//	field_to_csv(num_modes_in, field_in);
	
	
//	double *dker_MMI = malloc(sizeof(double) * N_sim);
//	symm_slab_dker(dker_MMI, N_MMI);
	double *dker_eff_MMI;
	double **field_MMI;
	int num_modes_MMI = slove_waveguide(dker_MMI, &dker_eff_MMI, &field_MMI);
	printf("num_modes_MMI = %d\n", num_modes_MMI);
	printf("n_eff_MMI\n");
	for (int j = 0; j < num_modes_MMI; j++) printf("%f\n", sqrt(dker_eff_MMI[j]) / kx_inc);
//	field_to_csv(num_modes_MMI, field_MMI);
	
	simulate_MMI(num_modes_MMI, dker_eff_MMI, field_MMI, field_in[0]);
	
	free(dker_in);
	free(dker_eff_in);
	free(field_in);
	free(dker_MMI);
	free(dker_eff_MMI);
	free(field_MMI);
	
	return 0;
}




void symm_slab_dker(double *dker, int N_core)
{
	int N1 = (N_sim - N_core) / 2;
	int N2 = (N_sim + N_core) / 2;
	
	for (int i = 0;  i < N1;    i++) dker[i] = dker_clad;
	for (int i = N1; i < N2;    i++) dker[i] = dker_core;
	for (int i = N2; i < N_sim; i++) dker[i] = dker_clad;
}




void offset_slab_dker(double *dker, unsigned int N_core, unsigned int offset)
{
	for (int i = 0; i < offset; i++) dker[i] = dker_clad;
	for (int i = offset; i < offset + N_core; i++) dker[i] = dker_core;
	for (int i = offset + N_core; i < N_sim; i++) dker[i] = dker_clad;
}




void coupler_dker(double *dker_in, double *dker_mmi, unsigned int separation)
{
	unsigned int N1 = (N_sim - separation)/2 - N_in;
	unsigned int N2 = (N_sim - separation)/2;
	unsigned int N3 = (N_sim + separation)/2;
	unsigned int N4 = (N_sim + separation)/2 + N_in;
	
	for (int i = 0;  i < N1;    i++) dker_in[i] = dker_clad;
	for (int i = N1; i < N2;    i++) dker_in[i] = dker_core;
	for (int i = N2; i < N_sim; i++) dker_in[i] = dker_clad;
	
	for (int i = 0;  i < N1;    i++) dker_mmi[i] = dker_clad;
	for (int i = N1; i < N2;    i++) dker_mmi[i] = dker_core;
	for (int i = N2; i < N3;    i++) dker_mmi[i] = dker_clad;
	for (int i = N3; i < N4;    i++) dker_mmi[i] = dker_core;
	for (int i = N4; i < N_sim; i++) dker_mmi[i] = dker_clad;
}




void coupler3_dker(double *dker_in, double *dker_mmi, unsigned int separation)
{
	unsigned int N1 = N_sim/2 - N_in/2 - separation - N_in;
	unsigned int N2 = N_sim/2 - N_in/2 - separation;
	unsigned int N3 = N_sim/2 - N_in/2;
	unsigned int N4 = N_sim/2 + N_in/2;
	unsigned int N5 = N_sim/2 + N_in/2 + separation;
	unsigned int N6 = N_sim/2 + N_in/2 + separation + N_in;
	
	for (int i = 0;  i < N1;    i++) dker_in[i] = dker_clad;
	for (int i = N1; i < N2;    i++) dker_in[i] = dker_core;
	for (int i = N2; i < N_sim; i++) dker_in[i] = dker_clad;
	
	for (int i = 0;  i < N1;    i++) dker_mmi[i] = dker_clad;
	for (int i = N1; i < N2;    i++) dker_mmi[i] = dker_core;
	for (int i = N2; i < N3;    i++) dker_mmi[i] = dker_clad;
	for (int i = N3; i < N4;    i++) dker_mmi[i] = dker_core;
	for (int i = N4; i < N5;    i++) dker_mmi[i] = dker_clad;
	for (int i = N5; i < N6;    i++) dker_mmi[i] = dker_core;
	for (int i = N6; i < N_sim; i++) dker_mmi[i] = dker_clad;
}




int slove_waveguide(double *dker, double **dker_eff, double ***field)
{
	double Vsq = dker_core - dker_clad;
	
	double *diag = malloc(sizeof(double) * N_sim);
	double *off_diag = malloc(sizeof(double) * (N_sim - 1));
	
	for (int i = 0; i < N_sim; i++) diag[i] = dker[i] - 2.0;
	double gamma = exp(-sqrt(dker_core - dker_clad));
	diag[0] += gamma;
	diag[N_sim - 1] += gamma;
	
	for (int i = 0; i < N_sim - 1; i++) off_diag[i] = 1.0;
	
	int check = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', N_sim, diag, off_diag, NULL, N_sim);
	if (check) exit(EXIT_FAILURE);
	
	
	int num_modes = 0;
	for (int i = 0; i < N_sim; i++) if (dker_clad < diag[i] && diag[i] < dker_core) num_modes++;
	
	*dker_eff = (double *) malloc(sizeof(double) * num_modes);
	*field = (double **) malloc(sizeof(double *) * num_modes);
	
	for (int i = 0; i < num_modes; i++) (*dker_eff)[i] = diag[N_sim - 1 - i];
	
	
	for (int j = 0; j < num_modes; j++)
	{
		double tol;
		do
		{
			for (int i = 0; i < N_sim; i++) diag[i] = dker[i] - 2.0;
			gamma = exp(-sqrt((*dker_eff)[j] - dker_clad));
			diag[0] += gamma;
			diag[N_sim - 1] += gamma;
			for (int i = 0; i < N_sim - 1; i++) off_diag[i] = 1.0;
			
			check = LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', N_sim, diag, off_diag, NULL, N_sim);
			if (check) exit(EXIT_FAILURE);
			
			tol = (diag[N_sim - 1 - j] - (*dker_eff)[j]) / Vsq;
			(*dker_eff)[j] = diag[N_sim - 1 - j];
		}
		while(fabs(tol) > 1.0E-10);
	}
	
	double *work = malloc(sizeof(double) * N_sim * N_sim);
	
	for (int j = 0; j < num_modes; j++)
	{
		(*field)[j] = malloc(sizeof(double) * N_sim);
		
		for (int i = 0; i < N_sim; i++) diag[i] = dker[i] - 2.0;
		gamma = exp(-sqrt((*dker_eff)[j] - dker_clad));
		diag[0] += gamma;
		diag[N_sim - 1] += gamma;
		for (int i = 0; i < N_sim - 1; i++) off_diag[i] = 1.0;
		
		check = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', N_sim, diag, off_diag, work, N_sim);
		if (check) exit(EXIT_FAILURE);
		
		(*dker_eff)[j] = diag[N_sim - 1 - j];
		for (int i = 0; i < N_sim; i++) (*field)[j][i] = work[N_sim * (N_sim - 1 - j) + i];
	}
	
	
	free(diag);
	free(off_diag);
	free(work);
	
	
	return num_modes;
}




void field_to_csv(unsigned int num_modes, double **field)
{
	FILE *file = fopen("field coupler.csv", "wt");
	
	for (int i = 0; i < N_sim; i++)
	{
		fprintf(file, "%f,", width_sim * (double)i / (double)N_sim);
		
		for (int j = 0; j < num_modes; j++)
		{
			fprintf(file, "%f,", field[j][i]);
		}
		fprintf(file, "\n");
	}
	
	fclose(file);
}




void MMI_picture(unsigned int, double *, double **);
void MMI_picture_colour(unsigned int, double *, double **);
void field_MMI_depth(double, int, double *, double **, double complex *);


double simulate_MMI(unsigned int num_modes, double *dker_eff_MMI, double **field_MMI, double *field_in)
{
	double *overlap = malloc(sizeof(double) * num_modes);
	double **field_scaled = malloc(sizeof(double *) * num_modes);
	for (int j = 0; j < num_modes; j++) field_scaled[j] = malloc(sizeof(double) * N_sim);
	double *field_start = malloc(sizeof(double) * N_sim);
	double efficiency = 0.0;
	
	for (int i = 0; i < N_sim; i++) field_start[i] = 0.0;
	
	for (int j = 0; j < num_modes; j++)
	{
		overlap[j] = 0.0;
		
		for (int i = 0; i < N_sim; i++)
		{
			field_scaled[j][i] = field_MMI[j][i];
			overlap[j] += field_in[i] * field_MMI[j][i];
		}
		
		efficiency += overlap[j] * overlap[j];
		
		for (int i = 0; i < N_sim; i++)
		{
			field_scaled[j][i] *= overlap[j];
			field_start[i] += field_scaled[j][i];
		}
	}
	
	printf("efficiency = %f\n", efficiency);
	double objective = 0.0;
	for (int i = 0; i < N_sim; i++) objective += field_start[i] * field_start[i];
	printf("match = %f\n", objective);
	
	double *n_eff_MMI = malloc(sizeof(double) * num_modes);
	for (int j = 0; j < num_modes; j++) n_eff_MMI[j] = sqrt(dker_eff_MMI[j]) / kx_inc;
	MMI_picture(num_modes, n_eff_MMI, field_scaled);
	MMI_picture_colour(num_modes, n_eff_MMI, field_scaled);
	
	
	double complex *field = malloc(sizeof(double complex) * N_sim);
	double objective_max = 0.0;
	double klength_reimage;
	
	for (int k = N_MMI; k < 2 * N_MMI; k++)
	{
		objective = 0.0;
		field_MMI_depth(kz_inc * k, num_modes, n_eff_MMI, field_scaled, field);
		for (int i = 0; i < N_sim; i++) objective += fabs(field_start[i]) * cabs(field[i]);
		
		if (objective > objective_max)
		{
			klength_reimage = kz_inc * k;
			objective_max = objective;
		}
	}
	
	for (double dkz = kz_inc / 2.0; dkz > 1.0E-10; dkz /= 2.0)
	{
		field_MMI_depth(klength_reimage + dkz, num_modes, n_eff_MMI, field_scaled, field);
		objective = 0.0;
		for (int i = 0; i < N_sim; i++) objective += fabs(field_start[i]) * cabs(field[i]);
		
		while (objective > objective_max)
		{
			klength_reimage += dkz;
			objective_max = objective;
			
			field_MMI_depth(klength_reimage + dkz, num_modes, n_eff_MMI, field_scaled, field);
			objective = 0.0;
			for (int i = 0; i < N_sim; i++) objective += fabs(field_start[i]) * cabs(field[i]);
		}
		
		field_MMI_depth(klength_reimage - dkz, num_modes, n_eff_MMI, field_scaled, field);
		objective = 0.0;
		for (int i = 0; i < N_sim; i++) objective += fabs(field_start[i]) * cabs(field[i]);
		
		while (objective > objective_max)
		{
			klength_reimage -= dkz;
			objective_max = objective;
			
			field_MMI_depth(klength_reimage - dkz, num_modes, n_eff_MMI, field_scaled, field);
			objective = 0.0;
			for (int i = 0; i < N_sim; i++) objective += fabs(field_start[i]) * cabs(field[i]);
		}
	}
	
	printf("reimage length = %f\n", klength_reimage / wavenumber);
//	printf("match = %f\n", objective_max);
	
	
	double klength_2split = klength_reimage / 2.0;
	
	field_MMI_depth(klength_2split, num_modes, n_eff_MMI, field_scaled, field);
	
	double intensity_max;
	int peak_index;
	intensity_max = -HUGE_VAL;
	for (int i = 0; i < N_sim/2; i++)
	{
		if (cabs(field[i]) > intensity_max)
		{
			intensity_max = cabs(field[i]);
			peak_index = i;
		}
	}
	
	double split2_efficiency = 0.0;
	for (int i = 0; i < N_sim/2; i++) split2_efficiency += cabs(field[i]) * fabs(field_in[i + N_sim/2 - peak_index]);
	split2_efficiency *= split2_efficiency;
	printf("\n2 split efficiency %f\n", split2_efficiency);
	
	free(overlap);
	free(n_eff_MMI);
	for (int j = 0; j < num_modes; j++) free(field_scaled[j]);
	free(field_scaled);
	free(field_start);
	
	
	return split2_efficiency;
}




void field_MMI_depth(double kz, int num_modes, double *n_eff_MMI, double **field_scaled, double complex *field)
{
	double complex *phaser = malloc(sizeof(double complex) * num_modes);
	for (int j = 0; j < num_modes; j++) phaser[j] = cexp(I * n_eff_MMI[j] * kz);
	
	for (int i = 0; i < N_sim; i++)
	{
		field[i] = 0.0;
		for (int j = 0; j < num_modes; j++) field[i] += field_scaled[j][i] * phaser[j];
	}
	
	free(phaser);
}




void MMI_picture(unsigned int num_modes, double *n_eff_MMI, double **field_scaled)
{
	double complex *field = malloc(sizeof(double complex) * N_sim);
	double **heightmap = malloc(sizeof(double *) * 2 * N_MMI);
	for (int k = 0; k < 2 * N_MMI; k++) heightmap[k] = malloc(sizeof(double) * N_sim);
	
	double height_max = -HUGE_VAL;
	
	for (int k = 0; k < 2 * N_MMI; k++)
	{
		field_MMI_depth(kz_inc * k, num_modes, n_eff_MMI, field_scaled, field);
		for (int i = 0; i < N_sim; i++)
		{
		//	heightmap[k][i] = cabs(field[i]);
			heightmap[k][i] = field[i] * conj(field[i]);
			height_max = fmax(height_max, heightmap[k][i]);
		}
	}
	
	unsigned char *picture = malloc(sizeof(unsigned char) * N_sim * N_sim);
	
	for (int k = 0; k < 2 * N_MMI; k++)
	for (int i = 0; i < N_sim; i++)
	{
		heightmap[k][i] /= height_max;
		picture[N_sim * k + i] = UCHAR_MAX * heightmap[k][i];
	}
	
	
	free(field);
	for (int k = 0; k < 2 * N_MMI; k++) free(heightmap[k]);
	free(heightmap);
	
	FILE *file = fopen("MMI coupler.pgm", "wb");
	fprintf(file, "P5\n%d %d\n%d\n", N_sim, 2 * N_MMI, UCHAR_MAX);
	fwrite(picture, sizeof(unsigned char), 2 * N_MMI * N_sim, file);
	fclose(file);
	free(picture);
}




/*
const double
 r = 6.0/7.0,
  y = 5.0/7.0,
   g = 4.0/7.0,
    c = 3.0/7.0,
     b = 2.0/7.0,
      m = 1.0/7.0;
*/

const double
 r = 1.0/2.0,
  y = 1.0/3.0,
   g = 1.0/4.0,
    c = 1.0/5.0,
     b = 1.0/6.0,
      m = 1.0/7.0;


unsigned char red(double x)
{
	double ans;
	if      (x < m) ans = (x - 0)/(m - 0);
	else if (x < b) ans = (x - b)/(m - b);
	else if (x < c) ans = 0.0;
	else if (x < g) ans = 0.0;
	else if (x < y) ans = (x - g)/(y - g);
	else if (x < r) ans = 1.0;
	else            ans = 1.0;
	return UCHAR_MAX * fabs(ans);
}

unsigned char green(double x)
{
	double ans;
	if      (x < m) ans = 0.0;
	else if (x < b) ans = 0.0;
	else if (x < c) ans = (x - b)/(c - b);
	else if (x < g) ans = 1.0;
	else if (x < y) ans = 1.0;
	else if (x < r) ans = (x - r)/(y - r);
	else            ans = (x - r)/(1 - r);
	return UCHAR_MAX * fabs(ans);
}

unsigned char blue(double x)
{
	double ans;
	if      (x < m) ans = (x - 0)/(m - 0);
	else if (x < b) ans = 1.0;
	else if (x < c) ans = 1.0;
	else if (x < g) ans = (x - g)/(c - g);
	else if (x < y) ans = 0.0;
	else if (x < r) ans = 0.0;
	else            ans = (x - r)/(1 - r);
	return UCHAR_MAX * fabs(ans);
}


void MMI_picture_colour(unsigned int num_modes, double *n_eff_MMI, double **field_scaled)
{
	double complex *field = malloc(sizeof(double complex) * N_sim);
	double **heightmap = malloc(sizeof(double *) * 2 * N_MMI);
	for (int k = 0; k < 2 * N_MMI; k++) heightmap[k] = malloc(sizeof(double) * N_sim);
	
	double height_max = -HUGE_VAL;
	
	for (int k = 0; k < 2 * N_MMI; k++)
	{
		field_MMI_depth(kz_inc * k, num_modes, n_eff_MMI, field_scaled, field);
		for (int i = 0; i < N_sim; i++)
		{
		//	heightmap[k][i] = cabs(field[i]);
			heightmap[k][i] = field[i] * conj(field[i]);
			height_max = fmax(height_max, heightmap[k][i]);
		}
	}
	
	struct pixel { unsigned char r, g, b; };
	struct pixel *picture = malloc(sizeof(struct pixel) * N_sim * N_sim);
	
	for (int k = 0; k < 2 * N_MMI; k++)
	for (int i = 0; i < N_sim; i++)
	{
		heightmap[k][i] /= height_max;
		picture[N_sim * k + i].r =   red(heightmap[k][i]);
		picture[N_sim * k + i].g = green(heightmap[k][i]);
		picture[N_sim * k + i].b =  blue(heightmap[k][i]);
	}
	
	
	free(field);
	for (int k = 0; k < 2 * N_MMI; k++) free(heightmap[k]);
	free(heightmap);
	
	FILE *file = fopen("MMI coupler colour.ppm", "wb");
	fprintf(file, "P6\n%d %d\n%d\n", N_sim, 2 * N_MMI, UCHAR_MAX);
	fwrite(picture, sizeof(struct pixel), 2 * N_MMI * N_sim, file);
	fclose(file);
	free(picture);
}




void MMI_sect_picture_colour(double kz, double scale, unsigned int num_modes, double *n_eff_MMI, double **field_scaled)
{
	double complex *field = malloc(sizeof(double complex) * N_sim);
	double *heightmap = malloc(sizeof(double) * N_sim);
	double height_max = -HUGE_VAL;
	field_MMI_depth(kz/scale, num_modes, n_eff_MMI, field_scaled, field);
	
	for (int i = 0; i < N_sim; i++)
	{
		heightmap[i] = field[i] * conj(field[i]);
		height_max = fmax(height_max, heightmap[i]);
	}
	
	struct pixel { unsigned char r, g, b; };
	struct pixel *picture = malloc(sizeof(struct pixel) * N_sim * N_sim);
	
	for (int i = 0; i < N_sim; i++)
	{
		heightmap[i] /= height_max * scale;
		picture[i].r =   red(heightmap[i]);
		picture[i].g = green(heightmap[i]);
		picture[i].b =  blue(heightmap[i]);
	}
	
	free(field);
	free(heightmap);
	
	FILE *file = fopen("MMI coupler sect colour.ppm", "wb");
	fprintf(file, "P6\n%d %d\n%d\n", N_sim, 1, UCHAR_MAX);
	fwrite(picture, sizeof(struct pixel), N_sim, file);
	fclose(file);
	free(picture);
}




void colour_gradient()
{
	unsigned int L = 256;
	double *lerp = malloc(sizeof(double) * L);
	
	for (int i = 0; i < L; i++)
	{
		lerp[i] = (double)i / (double)L;
	}
	
	struct pixel { unsigned char r, g, b; };
	struct pixel *grad = malloc(sizeof(struct pixel) * L);
	
	for (int i = 0; i < L; i++)
	{
		grad[i].r =   red(lerp[i]);
		grad[i].g = green(lerp[i]);
		grad[i].b =  blue(lerp[i]);
	}
	
	FILE *file = fopen("colour gradient.ppm", "wb");
	fprintf(file, "P6\n%d %d\n%d\n", L, 1, UCHAR_MAX);
	fwrite(grad, sizeof(struct pixel), L, file);
	fclose(file);
	free(lerp);
	free(grad);
}


