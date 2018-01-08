/*
 * DIP.h
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#ifndef DIP_H_
#define DIP_H_

#include <stdbool.h>
#include "mujoco.h"
#include "Visualizer.h"
#include <Eigen/Dense>
#include "Common_Structs.h"
#include "Environment.h"

using namespace Eigen;

class DIP : public Environment {
public:
	DIP(bool bSaveVid);
	virtual ~DIP();

	//calculate one step cost according to cost/reward function
	void GetStateActionCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, Matrix<double, nU, 1> mU, double* l, bool bTerminalState = false);

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

private:

	void Init();

	static constexpr double m_dCostCoeff_x = 0.5;
	static constexpr double m_dCostCoeff_xc = 0.0;
	static constexpr double m_dCostCoeff_y = 1.0;
	static constexpr double m_dCostCoeff_vx = 1.0;
	static constexpr double m_dCostCoeff_vy = 1.0;

	static constexpr double m_dCostCoeff_u = 1e-1;

	static constexpr double m_dCostCoeff_e = 1.0;

};

#endif /* DIP_H_ */
