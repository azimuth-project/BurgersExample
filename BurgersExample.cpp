//============================================================================
// Name        : BurgersExample.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <fstream>
#include <iostream>
#include <math.h>

#include "azimuthmath/specfunc/bessel/bessel.h"

using namespace std;
using namespace azimuthmath;

/**
 * This class computes the solution to the 1D Burgers equation
 * with periodic boundary conditions on the interval [0, 2 PI].
 * This is formula (8.7.7) in
 * <ul> Nokenath Debnath: "Nonlinear Partial Differential Equations For Scientists and Engineers",
 * 2nd edition, Birkhäuser Boston, 2005
 */
class BurgersSolution {
private:
	double u_zero;

	double viscosity;

	double bessel_argument;

	static int const bessel_evaluations = 30;

	double bessel_values[bessel_evaluations];

public:
	BurgersSolution(double const u_zero, double const viscosity);
	double evaluate(double xValue, double tValue);
};

BurgersSolution::BurgersSolution(double const u_zero, double const viscosity) {

	this->u_zero = u_zero;
	this->viscosity = viscosity;
	bessel_argument = u_zero / viscosity;

	int index;
	for (index = 0; index < bessel_evaluations; index++) {
		bessel_values[index] = azimuthmath::specfunc::bessel::bessel_i(index,
				bessel_argument);
	}
}

double BurgersSolution::evaluate(double xValue, double tValue) {

	double exp_values[bessel_evaluations];

	int index;

	for (index = 0; index < bessel_evaluations; index++) {
		double indexAsDouble = (double) index;
		exp_values[index] = exp(
				-indexAsDouble * indexAsDouble * viscosity * tValue * 0.25);
	}

	double nominator = 0;

	// calculating the nominator, note that we have to start at index 1 here
	// because the sum starts with the I_1 term.
	for (index = 1; index < bessel_evaluations; index++) {

		double indexAsDouble = (double) index;

		double sinValue = sin(0.5 * indexAsDouble * xValue);

		double summand = indexAsDouble * bessel_values[index]
				* exp_values[index] * sinValue;

		nominator += summand;
	}

	nominator *= (2.0 * viscosity);

	double denomiator = 0.0;

	for (index = 1; index < bessel_evaluations; index++) {

		double indexAsDouble = (double) index;

		double cosValue = cos(0.5 * indexAsDouble * xValue);

		double summand = bessel_values[index] * exp_values[index] * cosValue;

		denomiator += summand;
	}

	denomiator *= 2.0;

	denomiator += bessel_values[0];

	double result = nominator / denomiator;

	return result;
}

int main() {
	cout << "!!!Hello World!!!" << endl;

	BurgersSolution burgersSolution(1.0, 0.5);

	ofstream fout("burgersEquationExactSolution.txt");

	fout << "# TODO: write header \n";

	double t = 0.0;
	double x;

	/*
	 for (x = 0.0; x < M_TWOPI; x += 0.1) {
	 fout << x << " " << burgersSolution.evaluate(x, t) << "\n";
	 }
	 */

	for (t = 0.0; t < 20.0; t += 0.1) {
		for (x = 0.0; x < M_TWOPI; x += 0.1) {
			fout << t << " " << x << " " << burgersSolution.evaluate(x, t)
					<< "\n";
		}
	}

	/*
	 double x;

	 for (x = -3.0; x < 3.0; x += 0.1) {
	 fout << x << " " << azimuthmath::specfunc::bessel::bessel_i(0, x)
	 << "\n";
	 }

	 for (x = -3.0; x < 3.0; x += 0.1) {
	 fout << x << " " << azimuthmath::specfunc::bessel::bessel_i(2, x)
	 << "\n";
	 }
	 */

	fout.close();

	cout << "End of main" << endl;

	return 0;
}
