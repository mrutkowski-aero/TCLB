#ifndef ACTHRESHOLDNOW_H
#define ACTHRESHOLDNOW_H

#include "../CommonHandler.h"

#include "vHandler.h"
#include "Action.h"
#include "GenericAction.h"

class  acThresholdNow  : public  GenericAction  {
	int par;
	double level;
	public:
int Init ();
};

#endif // ACTHRESHOLDNOW_H
