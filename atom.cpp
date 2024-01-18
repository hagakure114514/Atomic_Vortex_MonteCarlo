#include "atom.h"

atom::atom() {
	r.x = 100.0e-6;
	r.y = 100.0e-6;
	r.z = 0.0;
	v.vx = 0.5e-3;
	v.vy = 0.0;
	v.vz = 0.0;
	v_pre=v;
	radius = sqrt(r.x * r.x + r.y * r.y);
	phi = (r.x == 0.0 && r.y == 0.0) ? 0.0 : atan2(r.y, r.x);

	acc_x = 0.0, acc_y = 0.0;
	s = state::d1;
}


atom::atom(position r0, velocity v0, state s0) {
	r = r0;
	v = v0;
	s = s0;
	v_pre= v;
	radius = sqrt(r.x * r.x + r.y * r.y);
	phi = (r.x == 0.0 && r.y == 0.0) ? 0.0 : atan2(r.y, r.x);

	acc_x = 0.0, acc_y = 0.0;
}


