/*
 * iLQG.h
 *
 *  Created on: Jun 21, 2017
 *      Author: tapgar
 */

#ifndef ILQG_H_
#define ILQG_H_

#include <Eigen/Dense>
#include <vector>
#include "Common_Structs.h"
#include "mujoco.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "math.h"
#include <stdio.h>

#include "eigenmvn.h"

#include <unistd.h>

#include "Visualizer.h"
#include "Environment.h"
#include <numeric>      // std::accumulate

#define debugBackward true
#define useEnergyState false


using namespace std;
using namespace Eigen;


class iLQG {
public:
	iLQG(Environment* env);
	virtual ~iLQG();

	void Run(const mjModel* m, const mjData* dmain);
	void RunMPC(const mjModel* m, const mjData* dmain, vector<ActionMatrix> u);

	void GetNextAction(Matrix<double, nU, 1>* u) { (*u) = unew[0]; };
	void GetTrajectory(vector<StateMatrix>* x, vector<ActionMatrix>* u) { (*u) = unew; (*x) = xnew; };
	void GetController(Matrix<double, nU, nX>* fb) { (*fb) = m_LastController.fb; }

private:

	//copy relevant data between mjData structs... apparently much faster than copying the full thing
	void CopyMJData(const mjModel* m, const mjData* dmain, mjData* dest);

	//simulate PD controller
	void TargetJointAngleToTorque(const mjData* d, Matrix<double, nU, 1> targ_q, Matrix<double, nU, 1>* u);

	//performs rollout for N iterations and records mujoco states in state_hist (pre allocated so it should be fast?)
	void Rollout(const mjModel* m, const mjData* dmain, vector<StateMatrix> x, vector<ActionMatrix> u, vector<StateMatrix>* x_res, vector<ActionMatrix>* u_res, vector<double>* c_res, vector<Gains> cntrl, double alpha, bool bCalcDerivs);

	//performs finite differences to calc first and second derivatives of cost w.r.t state and action
	void UpdateCost(Matrix<double, nX, 1> mX, Matrix<double, nU, 1> mU, Cost* C, bool bTerminalState);

	void SmoothAbsDeriv(double r, double alpha, double *dfdr, double *d2fdr2);

	bool BackwardPass(const mjModel* m, vector<ActionMatrix> u, double lambda, double *dV1, double *dV2m, int nBackPassNum );

	void GetBodyPosition(const mjData* dmain, int body_id, Vector3d* res);

	void GetCOMPosition(const mjModel* m, const mjData* dmain, Vector3d* com);

	double smooth_abs(double x, double alpha);

	void CalcExtVars(const mjModel* m, const mjData* d, double* cp_offset, double* cg_offsetXY, double* cg_offsetZ);

	void UpdateDerivatives(const mjModel* m, const mjData* dmain, int iter);

	bool BoxQP(Matrix<double, nU, nU> H, Matrix<double, nU, 1> g, Matrix<double, nU, 1> next_kff, Matrix<double, nU, 2> bounds, Matrix<double, nU, 1>* kff, MatrixXd* Hfree, vector<bool>* clamped);

	void ClampAction(Matrix<double, nU, 2> bounds, Matrix<double, nU, 1>* out);
	void IsClamped(Matrix<double, nU, 1> kff, Matrix<double, nU, 1> grad, Matrix<double, nU, 2> bounds, vector<bool>* isClamped);

	bool IncreaseLambda(double* lambda, double *dlambda);
	bool DecreaseLambda(double* lambda, double *dlambda);

	static constexpr double lambda_max = 1e10;
	static constexpr double lambda_min = 1e-6;
	static constexpr double lambda_factor = 1.6;

	vector<StateMatrix> xnew;
	vector<ActionMatrix> unew;
	vector<ActionMatrix> tau_log;
	vector<ActionCovMatrix> unew_std;
	vector<double> cnew;

	vector<StateMatrix> x_init;
	vector<ActionMatrix> u_init;

	vector<Gains> m_Controller;

	Gains m_LastController;

	vector<Derivs> m_Derivs;

	Matrix<double, nX, 1> Qx;
	Matrix<double, nU, 1> Qu;
	Matrix<double, nX, nX> Qxx;
	Matrix<double, nU, nU> Quu;
	Matrix<double, nU, nU> QuuF;
	Matrix<double, nU, nX> Qux;

	Matrix<double, nX, 1> Vx;
	Matrix<double, nX, nX> Vxx;

	Matrix<double, nX, 1> m_EndTargetPos;

	Environment* m_pEnv;

};

#endif /* ILQG_H_ */
