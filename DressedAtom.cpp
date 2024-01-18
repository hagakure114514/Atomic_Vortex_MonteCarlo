#include "DressedAtom.h"

DressedAtom::DressedAtom()
{
	count_sp = 0;
	detuning = detuning0;			//detuning between |e> and |g1>
	srand((unsigned)time(NULL));	//a magic phrase to make random numbers be based on time

}

///////////////////////////////////// Repump beam //////////////////////////////////////////////

// Repump
void DressedAtom::process_repump(atom* obj)
{
		if (repump_emission(obj) == 1) {		// spontaneous emission occur in Repump beam
		recoil_diss_pm(obj);
		count_sp++;
	}
}

// spontaneou semission hanbetu
bool DressedAtom::repump_emission(atom* obj)
{
	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	bool p_pm;
	double psp_pm = dist(rand_src);

	if (obj->s == state::d1) {
		p_pm = 0;

	}else if (obj->s == state::d2) {
		double p2 = branch * (1.0 - exp( - gamma * s2_pm(obj->radius) * dt / 2.0));
		p_pm = psp_pm<=p2 ? 1 : 0;
	}
	else {
		double p_diss = dist(rand_src);
		if (p_diss <= branch) {
			obj->s = state::d1;
			p_pm=0;
		}
		else {
			obj->s = state::d2;
			double p2 = branch * (1.0 - exp( - gamma * s2_pm(obj->radius) * dt / 2.0));
			p_pm = psp_pm <=p2 ? 1 : 0;
		}
	}
	return p_pm;
}

// recoil
void DressedAtom::recoil_diss_pm(atom* obj)
{
	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	
	double sp_psi_pm = 2.0 * M_PI * dist(rand_src);
	double sp_theata_pm = M_PI * dist(rand_src);

	if(flag_sp == 2) {
		// dipole radiation direction
		dip_sin2(&sp_theata_pm);
		obj->l_rot += hbar * (-(double)l_pm);
		obj->v.vx += hbar * k_wave * sin(sp_theata_pm) * cos(sp_psi_pm) / mass;
		obj->v.vy += hbar * k_wave * sin(sp_theata_pm) * sin(sp_psi_pm) / mass;
		obj->v.vz += - hbar * k_wave / mass + hbar * k_wave * cos(sp_theata_pm) / mass;
	}else{
		// default
		obj->l_rot += hbar * (-(double)l_pm);
		obj->v.vx += hbar * k_wave * sin(sp_theata_pm) * cos(sp_psi_pm) / mass;
		obj->v.vy += hbar * k_wave * sin(sp_theata_pm) * sin(sp_psi_pm) / mass;
		obj->v.vz += - hbar * k_wave / mass + hbar * k_wave * cos(sp_theata_pm) / mass;
	}

	obj->s = state::d1;
}




///////////////////////////////////// Guide beam //////////////////////////////////////////////


// Optcal potential
void DressedAtom::process_dipole(atom* obj)
{
	detuning_doppler(obj);
	force_dip(obj);
}

// dissipative force
void DressedAtom::process_diss(atom* obj)
{
	if (spontaneous_emission(obj) == 1) {		// spontaneous emission occur
		recoil_diss(obj);
		count_sp++;
	}

}





///////////////////////////////////// Step motion //////////////////////////////////////////////

// motion within a time step dt
void DressedAtom::step_motion(atom* obj)
{
	obj->r.x += ( 1.0/2.0 * (obj->v.vx + obj->v_pre.vx) - obj->l_rot * sin(obj->phi) / (obj->radius * mass) ) * dt + 1.0/2.0 * obj->acc_x * dt * dt;
	obj->r.y += ( 1.0/2.0 * (obj->v.vy + obj->v_pre.vy) + obj->l_rot * cos(obj->phi) / (obj->radius * mass) ) * dt + 1.0/2.0 * obj->acc_y * dt * dt;
	obj->r.z += 1.0/2.0 * (obj->v.vz + obj->v_pre.vz) * dt - 1.0/2.0 * G * dt * dt;

	obj->v.vx += obj->acc_x * dt;
	obj->v.vy += obj->acc_y * dt;
	obj->v.vz += - G * dt;
	obj->v_pre = obj->v;

	obj->radius = sqrt(obj->r.x * obj->r.x + obj->r.y * obj->r.y);
	obj->phi = (obj->r.x == 0.0 && obj->r.y == 0.0) ? obj->phi : atan2(obj->r.y, obj->r.x);		//An azimuth need not be defined when r=0 (i.e. x=0 and y=0) but I defined it as it becomes continuous.

	obj->acc_x = 0.0;
	obj->acc_y = 0.0;	
}



///////////////////////////////////// Functions for a Guide beam //////////////////////////////////////////////

// intensity function of Optical Vortex
double DressedAtom::intensity(double x)
{
	// Optical Vortex Electric Field: E0*sqrt(2)*r/w0 *exp(-r^2/w0^2)
	double ints = 2.0*I0*x*x/(w0*w0)*exp(-2.0*x * x / (w0 * w0));			// intensity [V/m]
	return ints;
}

// detuning doppler shift
void DressedAtom::detuning_doppler(atom* obj)
{
	double vphi = -(obj->v_pre.vx) * sin(obj->phi) + (obj->v_pre.vy) * cos(obj->phi) + obj->l_rot / (obj->radius * mass);
	detuning = detuning0 - k_wave * (obj->v_pre.vz) + l * vphi / (obj->radius);
}

// saturation parameter between |g1> and |e>
double DressedAtom::s1(double x)
{
	double s = intensity(x) / (I_s1 * (1.0 + 4.0 * detuning * detuning / (gamma1 * gamma1)));
	
	if (s > 0.1) {
		printf("s parameter1 is error %e\n", s);
	}
	
	return s;
}

double DressedAtom::s2(double x)
{
	double s = intensity(x) / (I_s2 * (1.0 + 4.0 * (detuning + delta_hfs) * (detuning + delta_hfs) / (gamma2 * gamma2)));

	if (s > 0.1) {
		printf("s parameter2 is error %e\n", s);
	}

	return s;
}

double DressedAtom::grad_s1(double x)
{
	double grad_int = 4.0 * I0/(w0*w0) *(x - 2.0*x * x*x / (w0 * w0)) * exp(-2 * x * x / (w0 * w0));
	double gs = grad_int /(I_s1 *(1 + 4 * detuning * detuning / (gamma1*gamma1)));
	return gs;
}

double DressedAtom::grad_s2(double x)
{
	double grad_int = 4.0 * I0 / (w0 * w0) * (x - 2.0 * x * x * x / (w0 * w0)) * exp(-2 * x * x / (w0 * w0));
	double gs = grad_int /(I_s2* (1 + 4 * (detuning + delta_hfs) * (detuning + delta_hfs) / (gamma2*gamma2)));
	return gs;
}

// accel change in Optical Potential
void DressedAtom::force_dip(atom* obj)
{
	double f_dip = 0.0;
	if (obj->s == state::d1) {
		f_dip = - 2.0/3.0 * hbar * detuning / 2.0 * grad_s1(obj->radius)/(1.0+s1(obj->radius));
	}
	else if (obj->s == state::d2) {
		f_dip = - 2.0/3.0 * hbar * (detuning + delta_hfs) / 2.0 * grad_s2(obj->radius) / (1.0+s2(obj->radius));
	}
	else {
		std::mt19937 rand_src(std::random_device{}());
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		double p_diss = dist(rand_src);

		if (p_diss <= branch) {
			f_dip = - 2.0/3.0 * hbar * detuning / 2.0 * grad_s1(obj->radius)/(1.0+s1(obj->radius));
			obj->s = state::d1;
		}
		else {
			f_dip = - 2.0/3.0 * hbar * (detuning + delta_hfs) / 2.0 * grad_s2(obj->radius) / (1.0+s2(obj->radius));
			obj->s = state::d2;
		}
	}

	if ( f_dip > 0 ) {
		f_dip = 0;
		printf("dipole force error\n");
	}

	obj->acc_x = f_dip * cos(obj->phi) / mass;
	obj->acc_y = f_dip * sin(obj->phi) / mass;
}



// spontaneou semission hanbetu
bool DressedAtom::spontaneous_emission(atom* obj)
{
	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	bool p;
	double psp = dist(rand_src);

	if (obj->s == state::d1) {
		double p1 = 1.0 - exp(-gamma * s1(obj->radius) * dt / 2.0);
		p = psp<=p1 ? 1 : 0;
	}
	else if (obj->s == state::d2) {
		double p2 = 1.0 - exp(-gamma * s2(obj->radius) * dt / 2.0);
		p = psp<=p2 ? 1 : 0;
	}
	else {
		double p_diss = dist(rand_src);
		if (p_diss <= branch) {
			obj->s = state::d1;
			double p1 = 1.0 - exp(-gamma * s1(obj->radius) * dt / 2.0);
			p = psp<=p1 ? 1 : 0;
		}
		else {
			obj->s = state::d2;
			double p2 = 1.0 - exp(-gamma * s2(obj->radius) * dt / 2.0);
			p = psp<=p2 ? 1 : 0;
		}
	}
	return p;
}

// recoil
void DressedAtom::recoil_diss(atom* obj)
{
	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	double p_recoil = dist(rand_src);
	double sp_psi = 2.0 * M_PI * dist(rand_src);
	double sp_theata = M_PI * dist(rand_src);

	if(flag_sp == 1){
		//no OAM mode
		obj->l_rot = 0.0;
		obj->v.vx += - hbar * k_wave * sin(sp_theata) * cos(sp_psi) / mass;
		obj->v.vy += - hbar * k_wave * sin(sp_theata) * sin(sp_psi) / mass;
		obj->v.vz += hbar * k_wave / mass - hbar * k_wave * cos(sp_theata) / mass;
	}else if(flag_sp == 2){
		// dipole radiation direction (azimuthal polarization)
		// sp_theata	: dipole radiation angle 0-pi
		// sp_psi		: angle around azimuthal axis
		dip_sin2(&sp_theata);
		obj->l_rot += hbar * (double)l;
		obj->v.vx += -hbar * k_wave * cos(sp_theata + obj->phi - M_PI / 2.0) * cos(sp_psi) / mass;
		obj->v.vy += -hbar * k_wave * sin(sp_theata + obj->phi - M_PI / 2.0) * cos(sp_psi) / mass;
		obj->v.vz += hbar * k_wave / mass - hbar * k_wave * sin(sp_psi) / mass;
	}else if(flag_sp == 4){
		// all directio (coherent OAM radiation)
		obj->l_rot += hbar * (double)l * 2.0;
		obj->v.vx += - hbar * k_wave * sin(sp_theata) * cos(sp_psi) / mass;
		obj->v.vy += - hbar * k_wave * sin(sp_theata) * sin(sp_psi) / mass;
		obj->v.vz += hbar * k_wave / mass - hbar * k_wave * cos(sp_theata) / mass;
	}else{
		// default
		obj->l_rot +=  hbar * (double)l;
		obj->v.vx += - hbar * k_wave * sin(sp_theata) * cos(sp_psi) / mass;
		obj->v.vy += - hbar * k_wave * sin(sp_theata) * sin(sp_psi) / mass;
		obj->v.vz += hbar * k_wave / mass - hbar * k_wave * cos(sp_theata) / mass;
	}

	if (obj->s == state::d1) {
		if (p_recoil <= branch) {
			// To |1,n-1>
			obj->s = state::d1;
		}else{
			// To |2,n-1>
			// Sisyphus cooling
			obj->s = state::d2;

			double uopt1 = 2.0/3.0 * hbar * detuning / 2.0 * log(1.0+s1(obj->radius));
			double uopt2 = 2.0/3.0 * hbar * (detuning + delta_hfs) / 2.0 * log(1.0+s2(obj->radius));
			double dK = uopt1 - uopt2;
			double K_r = 1.0/2.0 * mass * ((obj->v.vx*cos(obj->phi))*(obj->v.vx*cos(obj->phi))+(obj->v.vy*sin(obj->phi))*(obj->v.vy*sin(obj->phi)));
			double tau = K_r-dK > 0? (K_r-dK)/K_r : 0.0;

			// Kinetic energy loss
			obj->v.vx *= (cos(obj->phi)*cos(obj->phi)*tau+sin(obj->phi)*sin(obj->phi));
			obj->v.vy *= (cos(obj->phi)*cos(obj->phi)+sin(obj->phi)*sin(obj->phi)*tau);
		}
	}else{
		if (p_recoil <= branch) {
			// To |1,n-1>
			obj->s = state::d1;
		}else{
			// To |2,n-1>
			obj->s = state::d2;
		}
	}
}


// Calculate atom energy
void DressedAtom::calc_energy(atom* obj)
{
	double Uopt = ( obj->s == state::d1 )? 2.0/3.0 * hbar * detuning / 2.0 *  log(1.0+s1(obj->radius)): 2.0/3.0 * hbar * (detuning + delta_hfs) / 2.0 *  log(1.0+s2(obj->radius));
	obj->E_kin = ( 1.0 / 2.0 * mass *( obj->v.vx * obj->v.vx + obj->v.vy * obj->v.vy + obj->v.vz * obj->v.vz) + Uopt + mass * G * obj->r.z ) * 2.0 / (3.0 * k_b);

}



// random(dipole-radiation)
void DressedAtom::dip_sin2(double* x) {

	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	double p_uni, giji_ransu;
	do {
		p_uni = dist(rand_src);
		giji_ransu = M_PI * dist(rand_src);

	} while (p_uni > sin(giji_ransu) * sin(giji_ransu));

	*x = giji_ransu;
}



///////////////////////////////////// Functions for a repump beam //////////////////////////////////////////////

// saturation parameter between |g2> and |e>
double DressedAtom::s2_pm(double x)
{
	// Gaussian Electric Field: E0 *exp(-r^2/w0_pm^2)
	double ints_pm = 2.0*I0_pm*x*x/(w0_pm * w0_pm)*exp(-2.0*x * x / (w0_pm * w0_pm));			// intensity [V/m]
	double s_pm = ints_pm / (I_s2 * (1.0 + 4.0 * detuning_pm * detuning_pm / (gamma2 * gamma2)));
	
	return s_pm;
}




