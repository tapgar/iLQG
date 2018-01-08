/*
 * Environment.h
 *
 *  Created on: Jul 5, 2017
 *      Author: tapgar
 */

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_


#include <stdbool.h>
#include "mujoco.h"
#include "Visualizer.h"
#include <Eigen/Dense>
#include "Common_Structs.h"

#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <vector>
#include <memory>

using namespace Eigen;
using namespace std;

class Environment {
public:
	Environment() { nMPCDir = 0;};
	virtual ~Environment() { };

	void Step(Matrix<double, nX, 1> targ_x, Matrix<double, nU, 1> u, Matrix<double, nU, nX> K);

	mjModel* getModel() { return mj_Model; };
	mjData* getData() { return mj_Data; };

	//calculate one step cost according to cost/reward function
	virtual void GetStateActionCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, Matrix<double, nU, 1> mU, double* l, bool bTerminalState = false) { };
	virtual void GetStateCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, double* l, bool bTerminalState = false) { };
	//update cost derivatives
	virtual void UpdateCostDerivatives(const mjModel* m, const mjData* dmain, Cost* C, bool bTerminal) { };
	//extracts the relevant states from the full mujoco state
	virtual void GetStateFromMujoco(const mjData* d, Matrix<double, nX, 1>* mX) { };
	//adds relevant states to A for x = f(x_dot)
	virtual void IntegrateJacobians(const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B) { };

	virtual void DrawTargetTrajectory(vector<StateMatrix> xnew) { };

	//performs finite differences on state and actions to compute the first fx and fu derivates and the full cost derivatives
	virtual void Jacobian(const mjModel* m, const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B, Cost* C, bool bTerminalState) {};
	virtual void JacobianControl(const mjModel* m, const mjData* dmain, Matrix<double, nL, nU>* B, int dist_idx = -1) {};

	virtual void GetBounds(Matrix<double, nU, 2>* bounds) { };

	virtual void GetInitU(vector<ActionMatrix>* u_init) { };

	double smooth_abs(double x, double alpha)
	{
		return sqrt(pow(x,2) + pow(alpha,2)) - alpha;
	}

	double smooth_abs_deriv(double x, double alpha)
	{
		return x/sqrt(pow(x,2) + pow(alpha,2));
	}

	double smooth_abs_second_deriv(double x, double alpha)
	{
		return 1/sqrt(pow(x,2) + pow(alpha,2)) - pow(x,2)*pow(pow(x,2) + pow(alpha,2), -1.5);
	}


	virtual int getMotorPosIdx(int u) { return 0; }
	virtual int getMotorVelIdx(int u) { return 0; }

	mjModel *mj_Model;
	mjData *mj_Data;
	bool mujoco_initialized;

	int nMPCDir;

	virtual void Init() { };

	virtual void Draw(mjData* d) { m_Vis->Draw(d); };

	virtual void UpdateControl(mjData* d, Matrix<double, nU, 1> u) {
		for (int i = 0; i < nU; i++)
			d->ctrl[i] = u(i,0);
	}

#if LOGGING

	void LogEnvVars(const mjModel* m, const mjData* dmain) { LogClassSpecificVars(m, dmain, &fileEnvVars); };

	virtual void LogClassSpecificVars(const mjModel* m, const mjData* dmain, ofstream* file) {};

	void MPC_End();
	void MPC_Start(int mpc_iter);

	void MPC_DirExist();

	void LogBackwardPass(Matrix<double, nX, 1> Vx, Matrix<double, nX, nX> Vxx, Matrix<double, nX, 1> Qx, Matrix<double, nX, nX> Qxx,
			Matrix<double, nU, nX> Qux, Matrix<double, nU, nU> Quu, Matrix<double, nU, 1> Qu, Matrix<double, nU, nX> fb, Matrix<double, nU, 1> ff, int back_iter);

	void LogForwardPass(Matrix<double, nX, nX> A, Matrix<double, nX, nU> B, Matrix<double, nX, 1> lx, Matrix<double, nX, nX> lxx,
			Matrix<double, nU, nX> lux, Matrix<double, nU, nU> luu, Matrix<double, nU, 1> lu);

	void LogTrajectory(Matrix<double, nX, 1> x, Matrix<double, nU, 1> u, double c);
	void LogEvalTrajectory(Matrix<double, nX, 1> x, Matrix<double, nU, 1> u, double c, int numUpdates);
	void LogHyperParams(double lambda, double alpha);

	void LogPI2Trajectory(Matrix<double, nX, 1> x, Matrix<double, nU, 1> u, double c, double S, double P, int iter);

	void UpdateFN(char* buff, std::string orig_fn);

	void UpdateFN(char* buff, std::string orig_fn, int iter);

	ofstream fileTrajx;
	ofstream fileTraju;
	ofstream fileTrajc;
	ofstream fileA;
	ofstream fileB;
	ofstream filelx;
	ofstream filelxx;
	ofstream filelu;
	ofstream filelux;
	ofstream fileluu;
	ofstream fileQx;
	ofstream fileQxx;
	ofstream fileQu;
	ofstream fileQux;
	ofstream fileQuu;
	ofstream fileVx;
	ofstream fileVxx;
	ofstream filefb;
	ofstream fileff;
	ofstream fileEnvVars;

	ofstream fileFollow;

	ofstream fileHP;
	ofstream fileEvalTrajx;
	ofstream fileEvalTraju;
	ofstream fileEvalTrajc;



	const std::string fileTrajx_name = "fileTrajx.csv";
	const std::string  fileTraju_name= "fileTraju.csv";
	const std::string  fileTrajc_name= "fileTrajc.csv";
	const std::string fileEvalTrajx_name = "fileEvalTrajx.csv";
	const std::string  fileEvalTraju_name= "fileEvalTraju.csv";
	const std::string  fileEvalTrajc_name= "fileEvalTrajc.csv";
	const std::string fileHP_name = "fileHP.csv";
	const std::string  fileA_name= "fileA.csv";
	const std::string  fileB_name= "fileB.csv";
	const std::string  filelx_name= "filelx.csv";
	const std::string  filelxx_name= "filelxx.csv";
	const std::string  filelu_name= "filelu.csv";
	const std::string  filelux_name= "filelux.csv";
	const std::string  fileluu_name= "fileluu.csv";
	const std::string  fileQx_name= "fileQx.csv";
	const std::string  fileQxx_name= "fileQxx.csv";
	const std::string  fileQu_name= "fileQu.csv";
	const std::string  fileQux_name= "fileQux.csv";
	const std::string  fileQuu_name= "fileQuu.csv";
	const std::string  fileVx_name= "fileVx.csv";
	const std::string  fileVxx_name= "fileVxx.csv";
	const std::string  fileff_name= "fileff.csv";
	const std::string  filefb_name= "filefb.csv";
	const std::string  fileEnvVars_name = "fileEnvVars.csv";
	const std::string  fileFollow_name = "fileFollow.csv";

	vector<shared_ptr<ofstream> > filepi2xs;
	const std::string  filepi2x_name= "pi2x";
	vector<shared_ptr<ofstream> > filepi2us;
	const std::string  filepi2u_name= "pi2u";
	vector<shared_ptr<ofstream> > filepi2cs;
	const std::string  filepi2c_name= "pi2c";
	vector<shared_ptr<ofstream> > filepi2Ss;
	const std::string  filepi2S_name= "pi2S";
	vector<shared_ptr<ofstream> > filepi2Ps;
	const std::string  filepi2P_name= "pi2P";

#endif

	Visualizer* m_Vis;
};


#endif /* ENVIRONMENT_H_ */
