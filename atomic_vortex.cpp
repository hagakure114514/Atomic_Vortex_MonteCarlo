// atomic_vortex.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>
#include "const_param.h"
#include "atom.h"
#include "DressedAtom.h"
#include <fstream>

using std::endl;
using std::ofstream;

int main()
{
	int flag_mode = 0;
	int flag_sp_tmp = 3;

	printf("select simulation mode?\n 0:trajectory of an atom, 1:molasess, 2:multi execution \n");
	std::cin >> flag_mode;

	// spontaneous emission mode  
	printf("select spontaneous emission mode?\n 1:no OAM, 2:dipole-radiation, 3:all direction\n");
	std::cin >> flag_sp_tmp;

	if (flag_mode == 0) {

		redo:

		char fname[30];
		sprintf_s(fname, "traj_p50um_v1cm.csv");
		ofstream ofs(fname);        // ファイルパスを指定する

		// オブジェクトのコンストラクタ
		position r0 = { 100.0e-6 / sqrt(2), 100.0e-6 / sqrt(2), 0.0 };
		velocity v0 = { 0.0, 0.0, 0.0 };
		atom Rb87(r0, v0, state::d1);		// 原子オブジェクト
		atom* rb87 = &Rb87;
		DressedAtom OV1;			// dressed-atom状態オブジェクト

		// 配列の定義
		double x[jloop + 1] = {}, y[jloop + 1] = {}, z[jloop + 1] = {};
		double  vz[jloop + 1] = {}, E[jloop + 1] = {};											// 運動エネルギー　+　光ポテンシャル + 位置エネルギー

		// spontaneous emission mode  
		OV1.flag_sp = flag_sp_tmp;
		printf("selected sp mode:%d \n", OV1.flag_sp);

		// 時間ステップごとの運動
		for (int i = 0; i <= jloop; i++) {
			OV1.process_repump(rb87);
			OV1.process_dipole(rb87);
			OV1.process_diss(rb87);
			OV1.step_motion(rb87);
			OV1.calc_energy(rb87);

			x[i] = rb87->r.x;  y[i] = rb87->r.y;  z[i] = rb87->r.z;
			vz[i] = rb87->v.vz; E[i] = rb87->E_kin;

			ofs << x[i] << ", " << y[i] << ", " << z[i] << ", " << E[i] << "," << vz[i] << endl;

			if (z[i] < -0.26) {
				printf("z potision under -26 cm with %d processes (%e s)\n", i, 5.0e-5 * i);
				double angVf = (-(rb87->v.vx) * sin(rb87->phi) + (rb87->v.vy) * cos(rb87->phi)) + rb87->l_rot / (mass * rb87->radius);
				if (angVf > 0) printf("The rotate direction is untied to OV\n");
				break;
			}

			if (rb87->radius > w0 / sqrt(2.0)) {
				printf("conf radius out at z=%e cm\n", z[i] * 1e2);
				break;
			}
		}

		if (rb87->r.z > -0.26) {
			goto redo;
		}

		printf("spontaneous emission %d times\n", OV1.count_sp);
		return 0;

	}
	else if (flag_mode == 1) {
		// 変数の定義
		int sum_sp = 0;
		int count_guide = 0;
		int count_vphi = 0;

		// 配列の定義
		double x0[SAMPLE + 1] = {}, y0[SAMPLE + 1] = {}, z0[SAMPLE + 1] = {};
		double vx0[SAMPLE + 1] = {}, vy0[SAMPLE + 1] = {}, vz0[SAMPLE + 1] = {};
		double T_final = 0.0;

		rm_position(x0, y0, z0);
		rm_velocity(vx0, vy0, vz0);

		printf("simulation execution\n");

		char fname[30];
		sprintf_s(fname, "stat_azimuthal_vel.csv");
		ofstream ofs(fname);        // ファイルパスを指定する

		for (int ii = 0; ii < SAMPLE; ii++) {

			position r0 = { x0[ii], y0[ii], z0[ii] };
			velocity v0 = { vx0[ii], vy0[ii], vz0[ii] };

			// オブジェクトのコンストラクタ
			atom Rb87(r0, v0, state::d1);		// 原子オブジェクト
			atom* rb87 = &Rb87;
			DressedAtom OV1;			// dressed-atom状態オブジェクト
			OV1.flag_sp = flag_sp_tmp;

			// 時間ステップごとの運動
			for (int i = 0; i <= jloop; i++) {
				OV1.process_repump(rb87);
				OV1.process_dipole(rb87);
				OV1.process_diss(rb87);
				OV1.step_motion(rb87);

				if (rb87->r.z < -0.26) {
					count_guide++;
					sum_sp += OV1.count_sp;
					T_final += 1.0 / 2.0 * mass * (rb87->v.vx * rb87->v.vx + rb87->v.vy * rb87->v.vy + rb87->v.vz * rb87->v.vz);
					double angVf = (-(rb87->v.vx) * sin(rb87->phi) + (rb87->v.vy) * cos(rb87->phi)) * rb87->radius + rb87->l_rot / (mass);
					ofs << angVf << endl;
					if (angVf > 0) count_vphi++;
					break;
				}

				if (rb87->radius > w0 / sqrt(2.0)) {
					break;
				}
			}
		}

		printf("avarage spontaneous emission %d/%d times, ", sum_sp, count_guide);
		printf("guide efficiency %lf, ", (double)count_guide / (double)SAMPLE * 100);
		printf("unity of azimuthal direction %lf \n", (2.0 * (double)count_vphi - (double)count_guide) / (double)count_guide * 100);			// R-L/R+L [%]
		printf("avarage Temperature %lf \n mK ", 2.0 / (3.0 * k_b) * T_final / (double)count_guide * 1e3);

		ofs << (double)count_guide / (double)SAMPLE << ", " << (double)sum_sp / (double)count_guide << ", " << (2.0 * (double)count_vphi - (double)count_guide) / (double)count_guide << ", " << 2.0 / (3.0 * k_b) * T_final / (double)count_guide << endl;
		return 0;

	}
	else if (flag_mode == 2) {

		char fname[30];
		sprintf_s(fname, "multi_exec_result.csv");
		ofstream ofs(fname);        // ファイルパスを指定する

		// 変数の定義
		double ratio_guide_max = -100.0, ratio_rotation_max = -100.0, temp_max = 0.0;
		double ratio_guide_min = 100.0, ratio_rotation_min = 100.0, temp_min = 1e5;

		for (int jj = 0; jj < 5; jj++) {

			printf("simulation execution\n");

			// 変数の定義
			int sum_sp = 0;
			int count_guide = 0;
			int count_vphi = 0;
			double T_final = 0.0;

			// 配列の定義
			double x0[SAMPLE + 1] = {}, y0[SAMPLE + 1] = {}, z0[SAMPLE + 1] = {};
			double vx0[SAMPLE + 1] = {}, vy0[SAMPLE + 1] = {}, vz0[SAMPLE + 1] = {};

			rm_position(x0, y0, z0);
			rm_velocity(vx0, vy0, vz0);


			for (int ii = 0; ii < SAMPLE; ii++) {

				position r0 = { x0[ii], y0[ii], z0[ii] };
				velocity v0 = { vx0[ii], vy0[ii], vz0[ii] };

				// オブジェクトのコンストラクタ
				atom Rb87(r0, v0, state::d1);		// 原子オブジェクト
				atom* rb87 = &Rb87;
				DressedAtom OV1;			// dressed-atom状態オブジェクト
				OV1.flag_sp = flag_sp_tmp;

				// 時間ステップごとの運動
				for (int i = 0; i <= jloop; i++) {
					OV1.process_repump(rb87);
					OV1.process_dipole(rb87);
					OV1.process_diss(rb87);
					OV1.step_motion(rb87);

					if (rb87->r.z < -0.26) {
						count_guide++;
						sum_sp += OV1.count_sp;
						T_final += 1.0 / 2.0 * mass * (rb87->v.vx * rb87->v.vx + rb87->v.vy * rb87->v.vy + rb87->v.vz * rb87->v.vz);
						double angVf = (-(rb87->v.vx) * sin(rb87->phi) + (rb87->v.vy) * cos(rb87->phi)) * rb87->radius + rb87->l_rot / (mass);
						if (angVf > 0) count_vphi++;
						break;
					}

					if (rb87->radius > w0 / sqrt(2.0)) {
						break;
					}
				}
			}

			ratio_guide_max = std::max((double)count_guide / (double)SAMPLE, ratio_guide_max);
			ratio_rotation_max = std::max((2.0 * (double)count_vphi - (double)count_guide) / (double)count_guide, ratio_rotation_max);
			temp_max = std::max(2.0 / (3.0 * k_b) * T_final / (double)count_guide, temp_max);

			ratio_guide_min = std::min((double)count_guide / (double)SAMPLE, ratio_guide_min);
			ratio_rotation_min = std::min((2.0 * (double)count_vphi - (double)count_guide) / (double)count_guide, ratio_rotation_min);
			temp_min = std::min(2.0 / (3.0 * k_b) * T_final / (double)count_guide, temp_min);

		}
		
		printf("Max guide efficiency: %lf\n, ", ratio_guide_max  * 100);
		printf("Max unity of azimuthal direction: %lf\n", ratio_rotation_max * 100);			// R-L/R+L [%]
		printf("Max avarage Temperature: %lf mK \n ", temp_max * 1e3);

		printf("Minimum guide efficiency: %lf\n", ratio_guide_min * 100);
		printf("Minimum unity of azimuthal direction: %lf \n", ratio_rotation_min * 100);			// R-L/R+L [%]
		printf("Minimum avarage Temperature: %lf mK \n", temp_min * 1e3);

		ofs << ratio_guide_max << ", " << ratio_guide_min << ", " << ratio_rotation_max << ", " << ratio_rotation_min << ", " << temp_max  << ", " << temp_min << endl;

		return 0;

	}
	else{

		return 0;

	}
}


void gauss(double* x) {
	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	double p_uni, giji_ransu;
	do {
		p_uni = dist(rand_src);
		giji_ransu = 2.0 * 3.0 * r0 / 3.0 * dist(rand_src) - 3.0 * r0 / 3.0;
	} while (p_uni > exp(-giji_ransu * giji_ransu / (2.0 * r0 / 3.0 * r0 / 3.0)));
	*x = giji_ransu;
}

void max_boltz(double* vx) {
	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	double p_uni, giji_ransu;
	do {
		p_uni = dist(rand_src);
		giji_ransu = v_max * dist(rand_src);
	} while (p_uni > mass * giji_ransu * giji_ransu / (2.0 * k_b * temp) *exp( 1.0 -mass*giji_ransu * giji_ransu / (2.0 *k_b*temp)));
	*vx = giji_ransu;
}

void rm_position(double* x, double* y, double* z){
	std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

	for (int j = 0; j < SAMPLE; j++) {
		double r0;
		gauss(&r0);

		double phi = 2.0 * M_PI * dist(rand_src);
    	double psi = M_PI * dist(rand_src);

		x[j] = r0* cos(phi) * sin(psi);
		y[j] = r0* sin(phi) * sin(psi);
		z[j] = r0* cos(psi);
	}
}

void rm_velocity(double* vx, double* vy, double* vz){
	std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

	for (int j = 0; j < SAMPLE; j++) {
		double v0;
		max_boltz(&v0);

		double phi = 2.0 * M_PI * dist(rand_src);
    	double psi = M_PI * dist(rand_src);

		vx[j] = v0* cos(phi) * sin(psi);
		vy[j] = v0* sin(phi) * sin(psi);
		vz[j] = - v0 * abs(cos(psi));
	}
}