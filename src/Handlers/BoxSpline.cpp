#include "BoxSpline.h"
std::string BoxSpline::xmlname = "BoxSpline";
#include "../HandlerFactory.h"

int BoxSpline::Init () {
		Pars = -1;
		pugi::xml_attribute attr;
                pugi::xml_node par = node.first_child();
		if (! par) {
			ERROR("%s needs exacly one child! - none found\n", node.name());
			return -1;
		}
		hand = new Handler(par, solver);
		par = par.next_sibling();
		if (par) {
			ERROR("%s needs exacly one child! - more found\n", node.name());
			return -1;
		}
		if (! (*hand)) return -1;
		if ((*hand).Type() != HANDLER_DESIGN) {
			ERROR("%s needs child of design type!\n", node.name());
			return -1;
		}

		Pars2 = (*hand)->NumberOfParameters();
		attr = node.attribute("nodes");
		if (attr) {
			Pars = attr.as_int();
		} else {
			Pars = 10;
			notice("number of modes not set in %s - setting to %d\n",node.name(), Pars);
		}
		output("Lenght of time-resolved control: %d\n", Pars2);
		tab2 = new double[Pars2];
		output("Lenght of boxspline control: %d\n", Pars);

		attr = node.attribute("periodic");
		if (attr) {
			per=true;
		} else {
			per=false;
		}

		attr = node.attribute("order");
		if (attr) {
			order=attr.as_int();
		} else {
			order=3;
		}
		
		s = new double[4];
		attr = node.attribute("s0");
		if (attr) {
			s[0]=attr.as_double();
		} else {
			s[0]=0;
		}
		
		attr = node.attribute("s1");
		if (attr) {
			s[1]=attr.as_double();
		} else {
			s[1]=0;
		}
		
		attr = node.attribute("s2");
		if (attr) {
			s[2]=attr.as_double();
		} else {
			s[2]=0;
		}
		
		attr = node.attribute("s3");
		if (attr) {
			s[3]=attr.as_double();
		} else {
			s[3]=0;
		}

		attr = node.attribute("lower");
		if (attr) {
			lower = solver->units.alt(attr.value());
		} else {
			notice("lower bound not set in %s - setting to -1\n",node.name());
			lower = -1;
		}
		attr = node.attribute("upper");
		if (attr) {
			upper = solver->units.alt(attr.value());
		} else {
			notice("upper bound not set in %s - setting to 1\n",node.name());
			upper = 1;
		}
		return Design::Init();
	};


int BoxSpline::Finish () {
		delete hand;
		return 0;
	}


int BoxSpline::NumberOfParameters () {
		return Pars;
	};


double BoxSpline::Pos (int j) {
		return (1.0*j) / Pars2;
	};


int BoxSpline::Parameters (int type, double * tab) {
	if (solver->mpi_rank == 0) {
		switch(type) {
		case PAR_GET:
			output("Getting the params and making fourier decomposition\n");
			(*hand)->Parameters(type, tab2);
			{
							double * mat = (double*)malloc(sizeof(double)*Pars*Pars);
							double * Y = (double*)malloc(sizeof(double)*Pars);
							double * bxsrow = (double*)malloc(sizeof(double)*Pars);

							for (int i = 0; i < Pars; i++) {
								Y[i] = 0.0;
								for (int k = 0; k < Pars; k++) mat[i + Pars * k] = 0.0;
							}
							for (int j = 0; j < Pars2; j++) {
								double x = Pos(j);
								vbxspline(x, s, order, Pars, per, 0, bxsrow);
								for (int i = 0; i < Pars; i++) {
									Y[i] += bxsrow[i] * tab2[j];
									for (int k = 0; k < Pars; k++) mat[i + Pars * k] += bxsrow[i] * bxsrow[k];
								}
							}
							for (int i=0; i<Pars;i++) {
                                    for (int k=0; k<Pars;k++) printf("%5lf ", mat[i+Pars*k]);
									printf("| %5lf\n", Y[i]);
                            }
			    
			    GaussSolve (mat, Y, tab, Pars);
                            free(mat);
							free(bxsrow);
                            free(Y);
			}
			return 0;
		case PAR_SET:
							double * bxsrow = (double*)malloc(sizeof(double)*Pars);
							printf("Setting the params with a fourier series\n");
							for (int j = 0; j < Pars2; j++) {
								tab2[j] = 0;
								double x = Pos(j);
								vbxspline(x, s, order, Pars, per, 0, bxsrow);
								for (int i = 0; i < Pars; i++) tab2[j] += bxsrow[i] * tab[i];
							}

							free(bxsrow);
			(*hand)->Parameters(type, tab2);
			return 0;
		case PAR_GRAD:
			output("Getting gradient and making fourier decomposition\n");
			(*hand)->Parameters(type, tab2);
			for (int i=0; i<Pars;i++) tab[i] = 0;
			for (int j=0; j<Pars2;j++) {
				double x = Pos(j);
				vbxspline(x, s, order, Pars, per, 0, bxsrow);
				for (int i=0; i<Pars; j++) tab[i] += bxsrow[i] * tab2[j];
			}
			return 0;
		case PAR_UPPER:
			for (int i=0;i<Pars;i++) {
				tab[i] = upper;
			}
			return 0;
		case PAR_LOWER:
			for (int i=0;i<Pars;i++) {
				tab[i] = lower;
			}
			return 0;
		default:
			ERROR("Unknown type %d in call to Parameters in %s\n", type, node.name());
			exit(-1);
		}
	} else {
		return (*hand)->Parameters(type, NULL);
	}
	};


// Register the handler (basing on xmlname) in the Handler Factory
template class HandlerFactory::Register< GenericAsk< BoxSpline > >;
