/*
 * MMI.c
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <lapacke.h>

#define _USE_MATH_DEFINES


const double wavelength = 1.55; // microns
const double wavenumber = 2.0 * M_PI / wavelength;

const double n_core = 3.5;
//const double n_clad = 3.1693;
const double n_clad = n_core - 2.0;

const unsigned int N_sim = 1 << 9;
const unsigned int N_in  = N_sim / 10;
const unsigned int N_MMI = N_sim >> 1;
//unsigned int N_in;

const double uwidth_sim = (double) N_sim / (double) N_MMI;
const double uwidth_in  = (double) N_in  / (double) N_MMI;
const double uwidth_MMI = 1.0;
//double uwidth_in;

//const double width_in = 1.0; // microns
const double width_MMI = 10.0; // microns
//const double width_MMI  = width_in  * uwidth_MMI / uwidth_in;
const double width_in = width_MMI * uwidth_in / uwidth_MMI;
const double width_sim  = width_MMI * uwidth_sim / uwidth_MMI;
//double width_in;

const double kwidth_in  = wavenumber * width_in;
const double kwidth_MMI = wavenumber * width_MMI;
const double kwidth_sim = wavenumber * width_sim;
//double kwidth_in;

const double kx_inc    = kwidth_MMI / (double) N_MMI;
const double dker_core = kx_inc * kx_inc * n_core * n_core;
const double dker_clad = kx_inc * kx_inc * n_clad * n_clad;

const double kz_inc = kx_inc * n_core * kwidth_MMI / 2.0 / M_PI;

void symm_slab_dker(double *, int);
int slove_waveguide(double *, double **, double ***);
void field_to_csv(unsigned int, double **);
double simulate_MMI(unsigned int, double *, double **, double *);
void offset_slab_dker(double *, unsigned int, unsigned int);


int main(int argc, char **argv)
{
	printf("n_core = %f\n", n_core);
	printf("n_clad = %f\n", n_clad);
	printf("width_in = %f microns\n", width_in);
	printf("width_MMI = %f microns\n", width_MMI);
	printf("width_sim = %f microns\n", width_sim);
	printf("approx reimage length = %f microns\n", kz_inc / wavenumber * N_MMI);
	
	double *dker_in = malloc(sizeof(double) * N_sim);
	
	symm_slab_dker(dker_in, N_in);
//	offset_slab_dker(dker_in, N_in, 100);
	double *dker_eff_in;
	double **field_in;
	int num_modes_in = slove_waveguide(dker_in, &dker_eff_in, &field_in);
	printf("num_modes_in = %d\n", num_modes_in);
	printf("n_eff_in\n");
	for (int j = 0; j < num_modes_in; j++) printf("%f\n", sqrt(dker_eff_in[j]) / kx_inc);
//	field_to_csv(num_modes_in, field_in);
	
	
	double *dker_MMI = malloc(sizeof(double) * N_sim);
	symm_slab_dker(dker_MMI, N_MMI);
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
	FILE *file = fopen("field.csv", "wt");
	
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
//	MMI_picture(num_modes, n_eff_MMI, field_scaled);
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
		picture[N_sim * k + i] = UCHAR_MAX * (1.0 - heightmap[k][i]);
	}
	
	
	free(field);
	for (int k = 0; k < 2 * N_MMI; k++) free(heightmap[k]);
	free(heightmap);
	
	FILE *file = fopen("MMI.pgm", "wb");
	fprintf(file, "P5\n%d %d\n%d\n", N_sim, 2 * N_MMI, UCHAR_MAX);
	fwrite(picture, sizeof(unsigned char), 2 * N_MMI * N_sim, file);
	fclose(file);
	free(picture);
}




unsigned char red(double x)
{
	double y;
	if      (x < 0.20) y = 0.0;
	else if (x < 0.25) y = 0.0;
	else if (x < 0.33) y = (x - 0.25)/(0.33 - 0.25);
	else if (x < 0.50) y = 1.0;
	else               y = 1.0;
	return UCHAR_MAX * fabs(y);
}

unsigned char green(double x)
{
	double y;
	if      (x < 0.20) y = 0.0;
	else if (x < 0.25) y = (x - 0.20)/(0.25 - 0.20);
	else if (x < 0.33) y = 1.0;
	else if (x < 0.50) y = (x - 0.50)/(0.33 - 0.50);
	else               y = (x - 0.50)/(1.00 - 0.50);
	return UCHAR_MAX * fabs(y);
}

unsigned char blue(double x)
{
	double y;
	if      (x < 0.20) y = (x - 0.00)/(0.20 - 0.00);
	else if (x < 0.25) y = 1.0;
	else if (x < 0.33) y = (x - 0.33)/(0.25 - 0.33);
	else if (x < 0.50) y = 0.0;
	else               y = (x - 0.50)/(1.00 - 0.50);
	return UCHAR_MAX * fabs(y);
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
	
	FILE *file = fopen("MMI.ppm", "wb");
	fprintf(file, "P6\n%d %d\n%d\n", N_sim, 2 * N_MMI, UCHAR_MAX);
	fwrite(picture, sizeof(unsigned char) * 3, 2 * N_MMI * N_sim, file);
	fclose(file);
	free(picture);
}


