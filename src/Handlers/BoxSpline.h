#ifndef BOXSPLINE_H
#define BOXSPLINE_H

#include "../CommonHandler.h"

#include "vHandler.h"
#include "Callback.h"
#include "Design.h"

class  BoxSpline  : public  Design  {
	int Pars;
	int Pars2;
	double * tab2;
	double lower, upper;
	FILE * f;
	Handler * hand;
	bool per;
	int order;
public:
	static std::string xmlname;
int Init ();
int Finish ();
int NumberOfParameters ();
double Pos (int j);
int Parameters (int type, double * tab);
};

#endif // BOXSPLINE_H
