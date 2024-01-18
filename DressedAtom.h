#pragma once

#include <random>	//necessary when calling the current time
#include "const_param.h"
#include "atom.h"

class DressedAtom
{
public:
	void process_repump(atom *obj);
	void process_dipole(atom *obj);
	void process_diss(atom* obj);
	void step_motion(atom* obj);
	void calc_energy(atom* obj);
	
	DressedAtom();
	int flag_sp;
	int count_sp;

private:
	double detuning;			//detuning between |e> and |g1> [rad/s]
	double k_wave = 2.0 * M_PI / lambda;								//wavenumber[/m]
	const double I0 = 2.0 * beam_power / (M_PI * w0 * w0);	//beam intensity coefficient [V/m]
	const double I0_pm = 2.0 * beam_power_pm / (M_PI * w0_pm * w0_pm);	//beam intensity coefficient [V/m]
	const double gamma1 = gamma * branch;									//spontaneous decay rate between |e> and |1>[rad/s]
	const double gamma2 = gamma * (1.0 - branch);							//spontaneous decay rate between |e> and |2>[rad/s]
	const double lambda1 = 780.2e-9;								//wavelength between |e> and |g1>[m]
	const double lambda2 = 780.2e-9;								//wavelength between |e> and |g2>[m]
	const double I_s1 = M_PI * h * c * gamma1 / (3.0 * lambda1 * lambda1 * lambda1);	//saturation intensity (|1>)[W/m^2]
	const double I_s2 = M_PI * h * c * gamma2 / (3.0 * lambda2 * lambda2 * lambda2);	//saturation intensity (|2>)[W/m^2]
	
	double intensity(double x);										// intensity of Optical Vortex
	void detuning_doppler(atom* obj);								// detuning doppler shift [rad/s]
	double s1(double x);											// saturation parameter between |g1> and |e>
	double s2(double x);											// saturation parameter between |g2> and |e>
	double grad_s1(double x);										// gradient of saturation parameter between |g1> and |e>
	double grad_s2(double x);										// gradient of saturation parameter between |g2> and |e>
	void force_dip(atom* obj);										// Optical potential
	bool spontaneous_emission(atom* obj);							// spontaneous emission
	void recoil_diss(atom* obj);									// recoil of spontaneous emission
	bool repump_emission(atom* obj);							// spontaneous emission
	void recoil_diss_pm(atom* obj);									// recoil of spontaneous emission

	double s2_pm(double x);											// saturation parameter between |g2> and |e> of repump

	void dip_sin2(double* x);										// random sin^2 x
};

