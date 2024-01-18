#pragma once

#include <math.h>
#include "const_param.h"

class atom
{
public:
	position r;
	velocity v;			// center-of-mass motion velocity
	velocity v_pre;		// center-of-mass motion velocity
	state s;
	double radius;
	double phi;
	double acc_x, acc_y;
	double l_rot = 0.0;		// orbital angular momentum around beam axis [kg*m^2/s]
	double E_kin = 0.0;

	atom();
	atom(position r0, velocity v0, state s0);	
};

