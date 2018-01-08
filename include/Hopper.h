/*
 * Hopper.h
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#ifndef Hopper_H_
#define Hopper_H_

#include <stdbool.h>
#include "mujoco.h"
#include "Visualizer.h"
#include <Eigen/Dense>
#include "Common_Structs.h"
#include "Environment.h"

using namespace Eigen;

class Hopper : public Environment {
public:
	Hopper(bool bSaveVid);
	virtual ~Hopper();

	//calculate one step cost according to cost/reward function
	void GetStateActionCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, Matrix<double, nU, 1> mU, double* l, bool bTerminalState = false);
	void GetStateCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, double* l, bool bTerminalState);

	//update cost derivatives
	void UpdateCostDerivatives(const mjModel* m, const mjData* dmain, Cost* C, bool bTerminal);

	//extracts the relevant states from the full mujoco state
	void GetStateFromMujoco(const mjData* d, Matrix<double, nX, 1>* mX);

	//adds relevant states to A for x = f(x_dot)
	void IntegrateJacobians(const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B);

	//adds trajectory points to visualization
	void DrawTargetTrajectory(vector<StateMatrix> xnew);

	//performs finite differences on state and actions to compute the first fx and fu derivates and the full cost derivatives
	void Jacobian(const mjModel* m, const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B, Cost* C, bool bTerminalState);

	//performs finite differences on actions to compute the first fu derivate
	void JacobianControl(const mjModel* m, const mjData* dmain, Matrix<double, nL, nU>* B, int dist_idx = -1);

	void GetInitU(vector<ActionMatrix>* u_init);

private:

	void Init();

	static constexpr double m_dCostCoeff_cmx = 5e1;
	static constexpr double m_dCostCoeff_cmz = 1e2;
	static constexpr double m_dCostCoeff_u = 1e-4;
	static const int m_nFootBodyID = 4;
	static const int m_nTorsoBodyID = 1;

	Matrix<double, nX, 1> dCMxdx;
	Matrix<double, nX, 1> dCMzdx;

	vector<int> toMJvel;

	void GetExtVars(const mjModel* m, const mjData* d, double* cp_offsetx, double* cp_offsetz);
	void GetBodyPosition(const mjData* dmain, int body_id, Vector3d* res);
	void GetCOMPosition(const mjModel* m, const mjData* dmain, Vector3d* com);

};

#endif /* Hopper_H_ */
