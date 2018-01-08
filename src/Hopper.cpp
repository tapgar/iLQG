/*
 * Hopper.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#include "Hopper.h"
#include <iostream>

using namespace Eigen;
using namespace std;

Hopper::Hopper(bool save_vid) {
	mujoco_initialized = false;
	Init();
	m_Vis = new Visualizer(mj_Model, save_vid, "Hopper");
}

Hopper::~Hopper() {
	// TODO Auto-generated destructor stub
}

void Hopper::Init() {

	// Activate mujoco and load the model if this is the first instance
	if (!mujoco_initialized) {
		mj_activate("/home/tapgar/.mujoco/mjkey.txt");
		char error[1000] = "Could not load binary model";
		mj_Model = mj_loadXML("/home/tapgar/DeepRLcourse/rllab/vendor/mujoco_models/hopper.xml", 0, error, 1000);
		if (!mj_Model) {
			mju_error_s("Load model error: %s", error);
			return;
		}
		mujoco_initialized = true;
	}

	mj_Model->opt.timestep = m_dDeltaSimTime_s;

	// Initialize mjData
	mj_Data = mj_makeData(mj_Model);

	mj_Data->qpos[2] = 0.5;
	mj_Data->qpos[3] = 0.707;

	mj_forward(mj_Model, mj_Data);

	for (int i = 0; i < nL; i++)
		toMJvel.push_back(i+3);

//	for (int i = 0; i < mj_Model->nbody; i++)
//	{
//		printf("Body: %d|\t%f\t%f\t%f\n", i, mj_Data->xpos[i*3], mj_Data->xpos[i*3 + 1], mj_Data->xpos[i*3 + 2]);
//	}

}

void Hopper::GetStateActionCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, Matrix<double, nU, 1> mU, double* l, bool bTerminalState)
{
	*l = 0;
	double cp_offsetx = 0.0;
	double cp_offsetz = 0.0;
	GetExtVars(m, dmain, &cp_offsetx, &cp_offsetz);
	*l += m_dCostCoeff_cmx*pow(cp_offsetx,2.0);
	*l += m_dCostCoeff_cmz*pow(0.8-cp_offsetz,2.0);
	if (bTerminalState)
		return;
	*l += m_dCostCoeff_u*1.0*mU.transpose()*mU;
}

void Hopper::GetStateCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, double* l, bool bTerminalState)
{
	*l = 0;
	double cp_offsetx = 0.0;
	double cp_offsetz = 0.0;
	GetExtVars(m, dmain, &cp_offsetx, &cp_offsetz);
	*l += m_dCostCoeff_cmx*pow(cp_offsetx,2.0);
	*l += m_dCostCoeff_cmz*pow(0.8-cp_offsetz,2.0);
	for (int i = 0; i < dmain->ncon; i++)
	{
		int body1 = m->geom_bodyid[dmain->contact[i].geom1];
		int body2 = m->geom_bodyid[dmain->contact[i].geom2];

//		printf("Contact %d \t %d|%d\n", i, body1, body2);

		if ((body1 == 0 && body2 == m_nFootBodyID) || (body2 == 0 && body1 == m_nFootBodyID))
			continue;

//		m_Vis->Draw((mjData*)dmain);
//		usleep(10E6);
		*l += 100.0;
	}
	if (bTerminalState)
		return;
	//*l *= m_dControlTime_s;
}

void Hopper::UpdateCostDerivatives(const mjModel* m, const mjData* dmain, Cost* C, bool bTerminal)
{

	//dCMdx needs to be calculated before calling this

	C->lu = Matrix<double, nU, 1>::Zero();
	C->luu = Matrix<double, nU, nU>::Zero();
	C->lux = Matrix<double, nU, nX>::Zero();

	double cp_offsetx = 0.0;
	double cp_offsetz = 0.0;
	GetExtVars(m, dmain, &cp_offsetx, &cp_offsetz);

	C->lx = -m_dCostCoeff_cmx*2.0*(0.0 - cp_offsetx)*dCMxdx;
	C->lx += -m_dCostCoeff_cmz*2.0*(0.8 - cp_offsetz)*dCMzdx;

	C->lxx = m_dCostCoeff_cmx*2.0*dCMxdx*dCMxdx.transpose();
	C->lxx += m_dCostCoeff_cmz*2.0*dCMzdx*dCMzdx.transpose();

	if (bTerminal)
		return;

	for (int i = 0; i < nU; i++)
	{
		C->lu(i,0) = 2.0*m_dCostCoeff_u*dmain->ctrl[i];
		C->luu(i,i) = 2.0*m_dCostCoeff_u;
	}
}

void Hopper::IntegrateJacobians(const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B)
{
	//rotational velocities
	for (int i = 0; i < nQd; i++)
		(*A)(i,nQ+i) = 1.0;

	(*A) = (*A) * m_dControlTime_s;
	(*B) = (*B) * m_dControlTime_s;

	for (int i = 0; i < nX; i++)
		(*A)(i,i) += 1.0;

#if PAUSE_VIS

	cout << "A matrix" << "\n" <<  (*A) << "\n";
	cout << "B matrix" << "\n" <<  (*B) << "\n";
	printf("Position: \n");
	for (int i = 0; i < nQ; i++)
		printf("%f\t", dmain->qpos[i]);
	printf("\nVelocity: \n");
	for (int i = 0; i < nQd; i++)
		printf("%f\t", dmain->qvel[i]);
	printf("\nTorques: \n");
	for (int i = 0; i < nU; i++)
		printf("%f\t", dmain->ctrl[i]);
	printf("\n");

	double cp_offsetx = 0.0;
	double cp_offsetz = 0.0;
	GetExtVars(mj_Model, dmain, &cp_offsetx, &cp_offsetz);

	printf("CP offset: %f\t%f\n", cp_offsetx, cp_offsetz);

	cout << "dCMxdx matrix" << "\n" <<  dCMxdx << "\n";
	cout << "dCMzdx matrix" << "\n" <<  dCMzdx << "\n";

	m_Vis->Draw((mjData*)dmain);
#endif

}


void Hopper::GetStateFromMujoco(const mjData* d, Matrix<double, nX, 1>* mX)
{
	for (int i = 0; i < nQ; i++)
		(*mX)(i,0) = d->qpos[i];
	for (int i = 0; i < nQd; i++)
		(*mX)(nQ+i,0) = d->qvel[i];
}

void Hopper::DrawTargetTrajectory(vector<StateMatrix> xnew)
{
	vector<Vector3d> pts;
	for (int i = 0; i < N; i++)
	{
		Vector3d pt = Vector3d::Zero();
		pt(2) = xnew[i](0,0);
		pt(0) = xnew[i](1,0);
		pts.push_back(pt);
	}
	m_Vis->SetTrajPoints(pts);
}

void Hopper::Jacobian(const mjModel* m, const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B, Cost* C, bool bTerminalState)
{

	/*
	 * This function performs finite differences about the current state x and action u.
	 *
	 * A (12 x 12)
	 * B (12 x 1)
	 *
	 * The output space u is actually reference joint positions (10) so the resulting actuator forces
	 * will be obtained by applying a PD controller about each of the joints.  To differentiate x_dot
	 * w.r.t control inputs u we will simulate the PD controller to find the equivalent motor torque
	 * then perform finite differences.  Quaternion differentiation is performed following mujoco
	 * derivative example.
	 *
	 * In the example derivative.cpp code they also allowed a couple iterations to improve qacc estimate
	 * which they termed warmup iterations.  The mujoco states passed to this function were collected
	 * during the rollout so the qacc value should already be "warmed up".
	 *
	 */

	mjData* d = mj_makeData(m);
	mjMARKSTACK
	mjtNum* center = mj_stackAlloc(d, m->nv);
    mjtNum* warmstart = mj_stackAlloc(d, m->nv);

	d->time = dmain->time;
	mju_copy(d->qpos, dmain->qpos, m->nq);
	mju_copy(d->qvel, dmain->qvel, m->nv);
	mju_copy(d->qacc, dmain->qacc, m->nv);
	mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, m->nv);
	mju_copy(d->userdata, dmain->userdata, m->nuserdata);
	mju_copy(d->qfrc_applied, dmain->qfrc_applied, m->nv);
	mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
	mju_copy(d->ctrl, dmain->ctrl, m->nu);

	mj_forward(m, d);

	for( int rep=1; rep<3; rep++ )
	   mj_forwardSkip(m, d, mjSTAGE_VEL, 1);

	// select output from forward dynamics
    mjtNum* output = d->qacc;
    double orig_cp_offsetx = 0.0;
    double orig_cp_offsetz = 0.0;
    GetExtVars(m, d, &orig_cp_offsetx, &orig_cp_offsetz);

    // save output for warmstart
    mju_copy(center, output, m->nv);
    mju_copy(warmstart, d->qacc_warmstart, m->nv);

    dCMxdx = Matrix<double, nX, 1>::Zero();
    dCMzdx = Matrix<double, nX, 1>::Zero();

	// finite-difference over reference angle: skip = mjSTAGE_VEL
	for( int i=0; i<nU; i++ )
	{

		d->ctrl[i] += m_deps;

		mju_copy(d->qacc_warmstart, warmstart, m->nv);
		mj_forwardSkip(m, d, mjSTAGE_VEL, 1);

		d->ctrl[i] -= m_deps;

		for (int j = 0; j < nQd; j++)
			(*B)(nQ+j,i) = (output[j] - center[j])/(m_deps);
	}

	// finite-difference over velocity: skip = mjSTAGE_POS
	for( int i=0; i < nQd; i++ )
	{
		//first subtract eps
		d->qvel[i] += m_deps;
		mju_copy(d->qacc_warmstart, warmstart, m->nv);
		mj_forwardSkip(m, d, mjSTAGE_POS, 1);

		// undo perturbation
		d->qvel[i] = dmain->qvel[i];

		for (int j = 0; j < nQd; j++)
		{
			(*A)(nQ+j,nQ+i) = (output[j] - center[j])/(m_deps);
		}
	}

	// finite-difference over position: skip = mjSTAGE_NONE
    for( int i=0; i<nQ; i++ )
    {
    	// get joint id for this dof
	    int jid = m->dof_jntid[i];

		// get quaternion address and dof position within quaternion (-1: not in quaternion)
		int quatadr = -1, dofpos = 0;
		if( m->jnt_type[jid]==mjJNT_BALL )
		{
			quatadr = m->jnt_qposadr[jid];
			dofpos = i - m->jnt_dofadr[jid];
		}
		else if( m->jnt_type[jid]==mjJNT_FREE && i>=m->jnt_dofadr[jid]+3 )
		{
			quatadr = m->jnt_qposadr[jid] + 3;
			dofpos = i - m->jnt_dofadr[jid] - 3;
		}

		// apply quaternion or simple perturbation
		if( quatadr>=0 )
		{
			mjtNum angvel[3] = {0,0,0};
			angvel[dofpos] = m_deps;
			mju_quatIntegrate(d->qpos+quatadr, angvel, 1);
		}
		else
			d->qpos[i] += m_deps;

		// evaluate dynamics, with center warmstart
		mju_copy(d->qacc_warmstart, warmstart, m->nv);
		mj_forwardSkip(m, d, mjSTAGE_NONE, 1);

		double cp_offsetx = 0.0;
		double cp_offsetz = 0.0;
		GetExtVars(m, d, &cp_offsetx, &cp_offsetz);

		//printf("%d\t%f\t%f\n", i, orig_cp_offset, cp_offset);

		dCMxdx(i,0) = (cp_offsetx - orig_cp_offsetx)/m_deps;
		dCMzdx(i,0) = (cp_offsetz - orig_cp_offsetz)/m_deps;

		// undo perturbation
		mju_copy(d->qpos, dmain->qpos, m->nq);

        for (int j = 0; j < nQd; j++)
		{
			(*A)(nQ+j,i) = (output[j] - center[j])/(m_deps);
		}
    }

    UpdateCostDerivatives(m, dmain, C, bTerminalState);

    mjFREESTACK

    mj_deleteData(d);

}


void Hopper::JacobianControl(const mjModel* m, const mjData* dmain, Matrix<double, nL, nU>* B, int disturb_idx)
{
//	for (int i = 0; i < nU; i++)
//		(*B)(i,i) = 1.0;
//	return;
	mjData* d = mj_makeData(m);
	mjMARKSTACK
	mjtNum* center = mj_stackAlloc(d, m->nv);
    mjtNum* warmstart = mj_stackAlloc(d, m->nv);

	d->time = dmain->time;
	mju_copy(d->qpos, dmain->qpos, m->nq);
	mju_copy(d->qvel, dmain->qvel, m->nv);
	mju_copy(d->qacc, dmain->qacc, m->nv);
	mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, m->nv);
	mju_copy(d->userdata, dmain->userdata, m->nuserdata);
	mju_copy(d->qfrc_applied, dmain->qfrc_applied, m->nv);
	mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
	mju_copy(d->ctrl, dmain->ctrl, m->nu);

	if (disturb_idx >= 0)
		d->qvel[toMJvel[disturb_idx]] += m_deps;

	mj_forward(m, d);

	for( int rep=1; rep<3; rep++ )
	   mj_forwardSkip(m, d, mjSTAGE_VEL, 1);

	// select output from forward dynamics
    mjtNum* output = d->qacc;

    // save output for warmstart
    mju_copy(center, output, m->nv);
    mju_copy(warmstart, d->qacc_warmstart, m->nv);

	// finite-difference over control input: skip = mjSTAGE_VEL
	for( int i=0; i<nU; i++ )
	{
		d->ctrl[i] += m_deps;

		mju_copy(d->qacc_warmstart, warmstart, m->nv);
		mj_forwardSkip(m, d, mjSTAGE_VEL, 1);

		d->ctrl[i] -= m_deps;

		for (int j = 0; j < nL; j++)
			(*B)(j,i) = (output[toMJvel[j]] - center[toMJvel[j]])/(m_deps);
	}

    mjFREESTACK

    mj_deleteData(d);

}

void Hopper::GetExtVars(const mjModel* m, const mjData* d, double* cp_offsetx, double* cp_offsetz)
{
	Vector3d foot_pos;
	Vector3d torso_pos;
	Vector3d com;
	GetBodyPosition(d, m_nFootBodyID, &foot_pos);
	GetBodyPosition(d, m_nTorsoBodyID, &torso_pos);
	GetCOMPosition(m, d, &com);

	//printf("COM/FOOT: %f\t%f\t%f\t|%f\t%f\t%f\n", com(0),com(1),com(2),foot_pos(0),foot_pos(1),foot_pos(2));

	(*cp_offsetx) = com(0) - foot_pos(0);
	(*cp_offsetz) = com(2) - foot_pos(2);
}

void Hopper::GetBodyPosition(const mjData* dmain, int body_id, Vector3d* res)
{
	for (int i = 0; i < 3; i++)
		(*res)(i) = dmain->xipos[body_id*3 + i];
}

void Hopper::GetCOMPosition(const mjModel* m, const mjData* dmain, Vector3d* com)
{
	double dTotalMass = 0;

	for (int i = 0; i < 3; i++)
		(*com)(i) = 0.0;

	for (int i = 1; i <= m_nFootBodyID; i++)
	{
		dTotalMass += m->body_mass[i];
		(*com)(0) += dmain->xipos[i*3]*m->body_mass[i];
		(*com)(1) += dmain->xipos[i*3 + 1]*m->body_mass[i];
		(*com)(2) += dmain->xipos[i*3 + 2]*m->body_mass[i];
	}

	for (int i = 0; i < 3; i++)
		(*com)(i) /= dTotalMass;
}

void Hopper::GetInitU(vector<ActionMatrix>* init_u)
{
	ActionMatrix null_u = ActionMatrix::Zero();

	for (int i = 0; i < nU; i++)
		null_u(i,0) = 0.0;

	for (int i =0 ; i < N; i++)
	{
		(*init_u).push_back(null_u);
	}
}

