

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

#define _USE_MATH_DEFINES


void symm_slab_dker(double, double, double *, int, int);
int slove_waveguide(double *, unsigned int, double **, double ***);
double simulate_MMI(double *, double, double, double, int, double **, double *, unsigned int, unsigned int);




int main(int argc, char **argv)
{
	const double wavelength = 1.55; // microns
	const double k0 = 2.0 * M_PI / wavelength;
	
	const double ncore = 3.5;
	const double nclad = 3.1693;
//	const double nclad = ncore - 2.5;
	
	const unsigned int N = 1 << 9;
//	const unsigned int N_core_in = N >> 4;
	const unsigned int N_core_MMI = N >> 1;
	
	const double W_MMI = 4.0; // microns
	const double L = W_MMI * (double) N / (double) N_core_MMI; // microns
	
	const double kW_MMI = k0 * W_MMI;
	const double kL = k0 * L;
	
	const double inc = kL / (double) N; // step size
	
	const double kdepth_repeat = ncore * kW_MMI * kW_MMI / 2.0 / M_PI;
//	printf("analytic repeat length %f\n", kdepth_repeat / k0);
	
	const double dker_core = inc * inc * ncore * ncore; // relative permittivity normalised by the step size
	const double dker_clad = inc * inc * nclad * nclad;
	
	
	double *dker_in = malloc(sizeof(double) * N); // relative permittivity of the waveguide
	double *dker_MMI = malloc(sizeof(double) * N);
	
	symm_slab_dker(dker_core, dker_clad, dker_MMI, N, N_core_MMI);
	double *dker_eff_MMI;
	double **field_MMI;
	int num_modes_MMI = slove_waveguide(dker_MMI, N, &dker_eff_MMI, &field_MMI);
	printf("%d\n", num_modes_MMI);
	for (int j = 0; j < num_modes_MMI; j++) printf("%f\n", sqrt(dker_eff_MMI[j]) / inc);
	
	
	FILE *file = fopen("power split vs in width.csv", "wt");
	
	for (int N_core_in = 1<<1; N_core_in < N_core_MMI/2; N_core_in+=2)
	{
		double W = W_MMI * (double) N_core_in / (double) N_core_MMI; // microns
		double kW = k0 * W;
		symm_slab_dker(dker_core, dker_clad, dker_in, N, N_core_in);
		double *dker_eff_in;
		double **field_in;
		
		int num_modes_in = slove_waveguide(dker_in, N, &dker_eff_in, &field_in);
		printf("%d\n", num_modes_in);
		for (int j = 0; j < num_modes_in; j++) printf("%f\n", sqrt(dker_eff_in[j]) / inc);
		
		
		fprintf(file, "%f,%f\n", kW / k0, simulate_MMI(field_in[0], dker_eff_in[0], inc, kdepth_repeat, N_core_MMI, field_MMI, dker_eff_MMI, N, num_modes_MMI));
	//	printf("\nckeck %d\n", N_core_in);
		
		free(dker_eff_in);
		for (int i = 0; i < num_modes_in; i++) free(field_in[i]);
		free(field_in);
	}
	
	
	free(dker_in);
	free(dker_MMI);
	free(dker_eff_MMI);
	for (int i = 0; i < num_modes_MMI; i++) free(field_MMI[i]);
	free(field_MMI);
	
	
	fclose(file);
	
	return 0;
}




void symm_slab_dker(double dker_core, double dker_clad, double *dker, int N, int N_core)
{
	int N1 = (N - N_core) / 2;
	int N2 = (N + N_core) / 2;
	
	for (int i = 0;  i < N1; i++) dker[i] = dker_clad;
	for (int i = N1; i < N2; i++) dker[i] = dker_core;
	for (int i = N2; i < N;  i++) dker[i] = dker_clad;
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
	
	
	free(diag);
	free(off_diag);
	free(work);
	
	
	return num_modes;
}




double field_MMI_depth(double, int, double *, int, double complex *, double **, double *);


double simulate_MMI(double *field_in, double dker_eff_in, double inc, double kdepth_repeat, int N_core_MMI, double **field_MMI, double *dker_eff_MMI, unsigned int N, unsigned int num_modes_MMI)
{
	double *overlap = malloc(sizeof(double) * num_modes_MMI);
	double **field_MMI_scaled = malloc(sizeof(double *) * num_modes_MMI);
	for (int j = 0; j < num_modes_MMI; j++) field_MMI_scaled[j] = malloc(sizeof(double) * N);
	double efficiency = 0.0;
	
	printf("\n");
	
	for (int j = 0; j < num_modes_MMI; j++)
	{
		overlap[j] = 0.0;
		
		for (int i = 0; i < N; i++) overlap[j] += field_in[i] * field_MMI[j][i];
		
		printf("%f\n", overlap[j]);
		efficiency += overlap[j] * overlap[j];
		
		for (int i = 0; i < N; i++) field_MMI_scaled[j][i] = field_MMI[j][i] * overlap[j];
	}
	
	printf("%f\n", efficiency);
	
	
	double *field_start = malloc(sizeof(double) * N);
	
	for (int i = 0; i < N; i++)
	{
		field_start[i] = 0.0;
		
		for (int j = 0; j < num_modes_MMI; j++)
		{
			field_start[i] += field_MMI_scaled[j][i];
		}
	}
	
	free(overlap);
	
	
	double complex *fieldmap = malloc(sizeof(double complex) * N);
	double length_adj = kdepth_repeat * 2.0 / inc;
	double objective;
	struct { double obj; double z; } find_repeat;
	find_repeat.obj = HUGE_VAL;
	
	for (int k = N/2; k < N; k++)
	{
		double z = (double) k / (double) N;
		z *= length_adj;
		objective = field_MMI_depth(z, num_modes_MMI, dker_eff_MMI, N, fieldmap, field_MMI_scaled, field_in);
		
		if (objective < find_repeat.obj)
		{
			find_repeat.obj = objective;
			find_repeat.z = z;
		}
	}
	
	
	for (double dz = length_adj / (double)N; dz > 1E-20; dz /= 2.0)
	{
		objective = field_MMI_depth(find_repeat.z + dz, num_modes_MMI, dker_eff_MMI, N, fieldmap, field_MMI_scaled, field_in);
		
		while (objective < find_repeat.obj)
		{
			find_repeat.obj = objective;
			find_repeat.z += dz;
			objective = field_MMI_depth(find_repeat.z + dz, num_modes_MMI, dker_eff_MMI, N, fieldmap, field_MMI_scaled, field_in);
		}
		
		objective = field_MMI_depth(find_repeat.z - dz, num_modes_MMI, dker_eff_MMI, N, fieldmap, field_MMI_scaled, field_in);
		
		while (objective < find_repeat.obj)
		{
			find_repeat.obj = objective;
			find_repeat.z -= dz;
			objective = field_MMI_depth(find_repeat.z - dz, num_modes_MMI, dker_eff_MMI, N, fieldmap, field_MMI_scaled, field_in);
		}
	}
	
	FILE *file = fopen("repeat length.csv", "at");
	fprintf(file, "%f\n", find_repeat.z / length_adj);
	fclose(file);
	
	field_MMI_depth(find_repeat.z / 2.0, num_modes_MMI, dker_eff_MMI, N, fieldmap, field_MMI_scaled, field_in);
	
//	FILE *file = fopen("dump.csv", "wt");
//	for (int i = 0; i < N; i++) fprintf(file, "%f\n", cabs(fieldmap[i]));
//	fclose(file);
	
	
	struct { double intensity; int indx; } peak;
	peak.intensity = -HUGE_VAL;
	
	for (int i = 0; i < N/2; i++)
	{
		if (cabs(fieldmap[i]) > peak.intensity)
		{
			peak.intensity = cabs(fieldmap[i]);
			peak.indx = i;
		}
	}
	
	
	if (cabs(fieldmap[peak.indx-1]) > cabs(fieldmap[peak.indx+1]))
	{
		peak.indx -= 1;
		peak.intensity = cabs(fieldmap[peak.indx]);
	}
	
	struct { double f1, df1, f2, df2, h; } peakinfo;
	peakinfo.f1 = cabs(fieldmap[peak.indx]);
	peakinfo.f2 = cabs(fieldmap[peak.indx+1]);
	peakinfo.df1 = (cabs(fieldmap[peak.indx+1]) - cabs(fieldmap[peak.indx-1])) / 2.0;
	peakinfo.df2 = (cabs(fieldmap[peak.indx+2]) - cabs(fieldmap[peak.indx])) / 2.0;
	
	peakinfo.h = (peakinfo.f1 - peakinfo.f2 - peakinfo.df2 / (double)N) / (peakinfo.df1 - peakinfo.df2);
	
	file = fopen("split width.csv", "at");
	fprintf(file, "%f\n", (double) (N - 2*peak.indx + peakinfo.h) / (double) N);
	fclose(file);
	
	
	double complex split2_quality = 0.0;
	
	for (int i = 0; i < N/2; i++)
	{
		split2_quality += fieldmap[i] * field_in[i + N/2 - peak.indx];
	}
	
	printf("\n2 split quality %f\n", cabs(split2_quality * conj(split2_quality)));
	
	
	free(fieldmap);
	for (int i = 0; i < num_modes_MMI; i++) free(field_MMI_scaled[i]);
	free(field_MMI_scaled);
	free(field_start);
	
	
	return cabs(split2_quality * conj(split2_quality));
}




double field_MMI_depth(double z, int num_modes_MMI, double *dker_eff_MMI, int N, double complex *field, double **field_MMI, double *field_in)
{
	double objective = 0.0;
	
	double complex *phaser = malloc(sizeof(double complex) * num_modes_MMI);
	
	for (int j = 0; j < num_modes_MMI; j++)
	{
		phaser[j] = cexp(I * sqrt(dker_eff_MMI[j]) * z);
	}
	
	
	for (int i = 0; i < N; i++)
	{
		field[i] = 0.0;
		for (int j = 0; j < num_modes_MMI; j++) field[i] += field_MMI[j][i] * phaser[j];
		objective += (fabs(field_in[i]) - cabs(field[i])) * (fabs(field_in[i]) - cabs(field[i]));
	}
	
	free(phaser);
	
	
	return objective;
}

