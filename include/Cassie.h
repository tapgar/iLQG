/*
 * Cassie.h
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#ifndef Cassie_H_
#define Cassie_H_

#include <stdbool.h>
#include "mujoco.h"
#include "Visualizer.h"
#include <Eigen/Dense>
#include "Common_Structs.h"
#include "Environment.h"

using namespace Eigen;

class Cassie : public Environment {
public:
	Cassie(bool bSaveVid);
	virtual ~Cassie();

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

	void Draw(mjData* d) {

        Vector3d left_foot;
        Vector3d right_foot;
        Vector3d com;
        Vector3d torso;

        GetBodyPosition(d, m_nLeftFootBodyID, &left_foot);
        GetBodyPosition(d, m_nRightFootBodyID, &right_foot);
        GetBodyPosition(d, m_nTorsoBodyID, &torso);
        GetCOMPosition(mj_Model, d, &com);

        Vector3d ave_foot = (left_foot + right_foot)*0.5;
        Vector3d com_ver = torso;
        ave_foot(2) = 0.0;
        com(2) = 0.0;

        vector<Vector3d> vec;
        vec.push_back(ave_foot);
        vec.push_back(com);
        vec.push_back(com_ver);
        com_ver(2) = 1.0;
        vec.push_back(com_ver);
        m_Vis->DrawWithPoints(d, vec);
	};

	void GetInitU(vector<ActionMatrix>* u_init);

	void UpdateControl(mjData* d, Matrix<double, nU, 1> targ_q);

	void GetBounds(Matrix<double, nU, 2>* bounds) {
		//for direct torque control
		(*bounds) = TorqueBounds;
	}

	int getMotorPosIdx(int u) { return toMotorposidx[u]; }
	int getMotorVelIdx(int u) { return toMotorvelidx[u]; }

	void LogClassSpecificVars(const mjModel* m, const mjData* dmain, ofstream* outfile);

private:

	void Init();

	static constexpr double m_dCostCoeff_xd = 1e0;
	static constexpr double m_dCostCoeff_q = 1e0;
	static constexpr double m_dCostCoeff_foot = 1e1;
	static constexpr double m_dCostCoeff_torsoH = 5e-10;
	static constexpr double m_dCostCoeff_torsoV = 1e2;
	static constexpr double m_dCostCoeff_footYaw = 1e-10;

	static constexpr double m_dAlpha_xd = 0.1;
	static constexpr double m_dAlpha_q = 2.0;
	static constexpr double m_dAlpha_foot = 2.0;
	static constexpr double m_dAlpha_torsoH = 0.1;
	static constexpr double m_dAlpha_torsoV = 2.0;

	static constexpr double m_dTargetVel_mps = 0.5;

	static constexpr double m_dCostCoeff_u[] = {5e-5, 5e-5, 1e-6, 1e-6, 5e-4, 5e-5, 5e-5, 1e-6, 1e-6, 5e-4};

	static const int m_nLeftFootBodyID = 8;
	static const int m_nRightFootBodyID = 20;
	static const int m_nTorsoBodyID = 1;


	Matrix<double, nX, 1> dCMxydx;
	Matrix<double, nX, 1> dCMzdx;
	Matrix<double, nX, 1> dCPdx;

	Matrix<double, nU, 2> TorqueBounds;

    int toMJposidx[nQ];
    int toMJvelidx[nQd];

    int toMotorposidx[nU];
    int toMotorvelidx[nU];

    Matrix<double, nU, 1, DontAlign> Kp;
    Matrix<double, nU, 1, DontAlign> Kd;
    Matrix<double, nU, 1, DontAlign> GearN;

	void GetExtVars(const mjModel* m, const mjData* d, double* cp_offset, double* cg_offsetXY, double* cg_offsetZ);
	void GetBodyPosition(const mjData* dmain, int body_id, Vector3d* res);
	void GetCOMPosition(const mjModel* m, const mjData* dmain, Vector3d* com);
	void IntegrateQuaternion(const mjData* dmain, Matrix<double, nX, nX>* A, int start_pos_idx, int start_w_idx);

};

#endif /* Cassie_H_ */
