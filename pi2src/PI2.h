/*
 * PI2.h
 *
 *  Created on: Jul 13, 2017
 *      Author: tapgar
 */

#ifndef PI2_H_
#define PI2_H_

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
#include <sys/time.h>

#include "eigenmvn.h"

#include <unistd.h>

#include "Visualizer.h"
#include "Environment.h"
#include <numeric>      // std::accumulate

using namespace Eigen;
using namespace std;

typedef Matrix<double, nU, nP> BasisMatrix;

struct PI_Derivs {
	Matrix<double, nL, nU> B;
	Matrix<double, nL, 1> Phi;
};

struct Path {

	Matrix<double, nP, 1> cntrl; //noisy controls
	vector<StateMatrix> state; //dont really need this

	vector<double> cost; //cost at a given time step

	vector<double> S; //path integral cost
	vector<double> P; //weighted path integral cost

	vector<PI_Derivs> Derivs;

	//avoid calculating multiple times
	vector<ControllableMatrix> H;
	vector<ControllableMatrix> Hinv;

};

class PI2 {

public:

	PI2(Environment *pEnv);

	void RunMPC(const mjModel* m, const mjData* dmain, vector<ActionMatrix> u);

	void GetNextAction(Matrix<double, nU, 1>* u) { (*u) = unew[0]; };
	void GetTrajectory(vector<StateMatrix>* x, vector<ActionMatrix>* u) { (*u) = unew; (*x) = xnew; };

private:

	Environment *m_pEnv;

	vector<StateMatrix> xnew;
	vector<ActionMatrix> unew;
	vector<double> cnew;

	vector<Path> m_Paths;

//	Matrix<double, nP, nP> CovMat;
//
//	Matrix<double, nP, nP> R;
//	Matrix<double, nP, nP> Rinv;

	MatrixXd CovMat;

	MatrixXd R;
	MatrixXd Rinv;

	vector<BasisMatrix> G;

	double lambda;
	double noise_factor;

	static constexpr double h = 10.0;

	void NominalRollout(const mjModel* m, const mjData* dmain);

	void Rollouts(const mjModel* m, const mjData* dmain, vector<ActionMatrix> init_cntrl);
	void Rollout(const mjModel* m, const mjData* dmain, int path_idx);
	void UpdateDerivatives(const mjModel* m, const mjData* dmain, int path_idx, int time_idx);

	void CalculateP(int path_idx, int time_idx, double minS, double maxS);
	void CalculateS(int path_idx, int time_idx);

	void UpdateNominal();

	void CovarianceMatrixAdaptation();


};



#endif /* PI2_H_ */
