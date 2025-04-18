//	g++ -o "%e" "%f" -IC:\mingw64\include\openblas -LC:\mingw64\lib -lopenblas -llapack -lblas

#include <iostream>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <lapacke.h>

#define _USE_MATH_DEFINES

const double k0 = 2.0 * M_PI / 1.55;
const int N = 96;
std::array<double, N> rel_per, diag;
std::array<double, N-1> off_diag;
std::array<std::array<double, N>, N> field;


double w = 2.0/17.0;
double indx(double x)
{
	if (x > 0.5 - w/2.0 and x < 0.5 + w/2.0) return 3.5;
	else return 3.1693;
}


int main(int argc, char **argv)
{
	double inc = 0.4 / w / (double) N;
	inc *= k0;
	inc *= inc;
	
	std::ofstream file("rel_per.csv");
	if (!file.is_open()) return 1;
	
	for (int i = 0; i < N; i++) rel_per[i] = i / (double) N;
	for (auto &it : rel_per)
	{
		it = indx(it);
		file << it << std::endl;
		it *= it;
		it *= inc;
	}
	
	file.close();
	
	
	double neff;
	neff = *std::max_element(rel_per.begin(), rel_per.end());
	std::cout << sqrt(neff / inc) << std::endl;
	
	double gL, gR, diff;
	
	do
	{
		
	//	gL = gR = 1.0;
		gL = exp(-sqrt(neff - rel_per.front()));
		gR = exp(-sqrt(neff - rel_per.back()));
		
		diag = rel_per;
		for (auto &it : diag) it -= 2.0;
		diag.front() += gL;
		diag.back() += gR;
		
		std::fill(off_diag.begin(), off_diag.end(), 1.0);
		
		for (auto &it : field) std::fill(it.begin(), it.end(), 0.0);
		
		
		LAPACKE_dstev(LAPACK_ROW_MAJOR, 'V', N, diag.data(), off_diag.data(), field[0].data(), N);
		
		
		diff = abs(1.0 - neff / diag[N-1]);
		
		 neff = diag[N-1];
		 std::cout << sqrt(neff / inc) << std::endl;
		
	}
	while (diff > 1e-5);
	
	file.open("output TE.csv");
	if (!file.is_open()) return 1;
	
	for (const auto &it : diag) file << sqrt(it / inc) << ",";
	file << std::endl;
	
	for (int i = 0; i < field.size(); i++)
	{
		for (const auto &jt : field[i])
		{
			file << jt << ",";
		}
		file << std::endl;
	}
	
	file.close();
	
	
	return 0;
}

