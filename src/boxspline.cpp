#include "boxspline.h"

/**
 * @brief Function for evaluating base 1D monotonic box-spline
 *
 * This function evaluates a single 1D monotonic box-spline in a single point
 * It is a utility function, and is rarely used on it's own.
 * 
 * For k=0 (s={}) and i=0 it's the Heaviside step function
 * For k=0 (s={}) it is the (-i)-th integral of the Heaviside step function
 * For k>0 and i=0 it's the Heaviside step function, k-times convoluted with an interval kernel 
 * 
 * @code
 * double s[2] = {0.1, 0.1};
 * printf("Value: %lg\n", bxs_base(s,k,3.14,0);
 * @endcode
 * The "derivative leval" argument is the order of derivative (or integral in the negatives)
 * to calculate. Zero means evaluation of function, one evaluation of it's derivative, and
 * minus one means the evaluation of integral.
 *
 * @param s Intervals lengths of convolution kernels
 * @param k Number of convolution kernels
 * @param x Point of evaluation
 * @param i Derivative level (defalut: 0)
 * @return The value of base 1D box-spline
 * @note See the general documentation for definitions
 * @warning Does not check input
 */

double bxs_base(double *s, int order, double x, int i) {
	if (order > 0) {
		if (s[0] == 0) {
			return bxs_base(s + 1, order - 1, x, i);
		}
		else {
			return (bxs_base(s + 1, order - 1, x + s[0] / 2, i - 1) - bxs_base(s + 1, order - 1, x - s[0] / 2, i - 1)) / s[0];
		}
	}
	else {
		if (x < 0) return 0;
		if (i > 0) return 0;
		double ret = 1;
		for (; i < 0; i++) ret = -ret * x / i;
		return ret;
	}
}

double bxspline(double x, double* s, int order, double* tab, int Pars, bool per, int d) {
	double dx = 1.0 / Pars; // Divide the overall [a,b] interval into Pars parts
	int w = floor(x / dx); // Calculate in which small interval we are.
	double ret = 0;
	if (d <= 0) ret += tab[mod(w, Pars, per)];
	int i = w;
	for (i = w + 1; i < 100 * Pars; i++) {
		double c = bxs_base(s, order, x - i * dx, d);
		if (c == 0) break;
		ret = ret + c * (tab[mod(i, Pars, per)] - tab[mod(i - 1, Pars, per)]);
	}
	for (i = w; i > -100 * Pars; i--) {
		double c = bxs_base(s, order, i * dx - x, d);
		if (d % 2 == 1) c = -c;
		if (c == 0) break;
		ret = ret + c * (tab[mod(i - 1, Pars, per)] - tab[mod(i, Pars, per)]);
	}
	return ret;
}

void vbxspline(double x, double* s, int order, int Pars, bool per, int d, double* bxsrow){
	double dx = 1.0 / Pars; // Divide the overall [a,b] interval into Pars parts
	int w = floor(x / dx); // Calculate in which small interval we are.
	int i = w;
	for (i = 0; i < Pars; i++) bxsrow[i] = 0;
	if (d <= 0) bxsrow[mod(w, Pars, per)] += 1;
	for (i = w + 1; i < 100 * Pars; i++) {
		double c = bxs_base(s, order, x - i * dx, d);
		if (c == 0) break;
		bxsrow[mod(i, Pars, per)] += c;
		bxsrow[mod(i - 1, Pars, per)] += -c;
	}
	for (i = w; i > -100 * Pars; i--) {
		double c = bxs_base(s, order, i * dx - x, d);
		if (d % 2 == 1) c = -c;
		if (c == 0) break;
		bxsrow[mod(i - 1, Pars, per)] += c;
		bxsrow[mod(i, Pars, per)] += -c;
	}
}


inline int mod(int i, int Pars, bool per) {
	//option for nonperiodic movement, but ends are flat
	//if (i <= 0 && !per) return 0;
	//if (i >= Pars && !per) return Pars-1;
	return ((i % Pars) + Pars) % Pars;
}



