/*
 * DIP.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#include "DIP.h"

using namespace Eigen;

DIP::DIP(bool save_vid) {
	mujoco_initialized = false;
	Init();
	m_Vis = new Visualizer(mj_Model, save_vid, "Double Inverted Pendulum");
}

DIP::~DIP() {
	// TODO Auto-generated destructor stub
}

void DIP::Init() {

	// Activate mujoco and load the model if this is the first instance
	if (!mujoco_initialized) {
		mj_activate("/home/tapgar/.mujoco/mjkey.txt");
		char error[1000] = "Could not load binary model";
		mj_Model = mj_loadXML("/home/tapgar/DeepRLcourse/rllab/vendor/mujoco_models/inverted_double_pendulum.xml", 0, error, 1000);
		if (!mj_Model) {
			mju_error_s("Load model error: %s", error);
			return;
		}
		mujoco_initialized = true;
	}

	mj_Model->opt.timestep = m_dDeltaSimTime_s;

	// Initialize mjData
	mj_Data = mj_makeData(mj_Model);

	//double qpos[3] = {0.0, 3.1415, 0.0};
	double qpos[3] = {0.0, 0.4, 0.0};
	mju_copy(mj_Data->qpos, qpos, 3);

}

void DIP::GetStateActionCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, Matrix<double, nU, 1> mU, double* l, bool bTerminalState)
{
	double l1 = 0.6;
	double l2 = l1;

	*l = 0;
//
//	if (bTerminalState)
//	{
		double xpos = mX(0,0) + l1*sin(mX(1,0)) + l2*sin(mX(1,0) + mX(2,0));
		double xposc = mX(0,0);
		double ypos = l1*cos(mX(1,0)) + l2*cos(mX(1,0) + mX(2,0));
		double xvel = mX(3,0) + mX(4,0)*l1*cos(mX(1,0)) + (mX(4,0) + mX(5,0))*cos(mX(1,0) + mX(2,0));
		double yvel = -mX(4,0)*l1*sin(mX(1,0)) - (mX(4,0) + mX(5,0))*sin(mX(1,0) + mX(2,0));

		*l += m_dCostCoeff_x*pow(xpos,2.0);
		*l += m_dCostCoeff_y*pow(l1+l2-ypos,2.0);
		*l += m_dCostCoeff_vx*pow(xvel,2.0);
		*l += m_dCostCoeff_vy*pow(yvel,2.0);
		*l += m_dCostCoeff_xc*pow(xpos-xposc,2.0);
		//printf("End State X: %f\t Y: %f\tVx: %f\tVy: %f\tEnd Cost: %f\n", xpos, ypos, xvel, yvel, *l);
//		return;
//	}

	*l += m_dCostCoeff_u*pow(mU(0,0),2.0);
}

void DIP::UpdateCostDerivatives(const mjModel* m, const mjData* dmain, Cost* C, bool bTerminal)
{

//	if (bTerminal)
//	{
		Matrix<double, nX, 1> dr1dx = Matrix<double, nX, 1>::Zero();
		Matrix<double, nX, 1> dr2dx = Matrix<double, nX, 1>::Zero();
		Matrix<double, nX, 1> dr3dx = Matrix<double, nX, 1>::Zero();
		Matrix<double, nX, 1> dr4dx = Matrix<double, nX, 1>::Zero();
		Matrix<double, nX, 1> dr5dx = Matrix<double, nX, 1>::Zero();

		double l1 = 0.6;
		double l2 = l1;

		dr1dx(0,0) = 1.0; //dx/dq0
		dr1dx(1,0) = l1*cos(dmain->qpos[1]) + l2*cos(dmain->qpos[1] + dmain->qpos[2]); //dx/dq1
		dr1dx(2,0) = l2*cos(dmain->qpos[1] + dmain->qpos[2]); //dx/dq2

		dr2dx(1,0) = -l1*sin(dmain->qpos[1]) - l2*sin(dmain->qpos[1] + dmain->qpos[2]); //dy/dq1
		dr2dx(2,0) = -l2*sin(dmain->qpos[1] + dmain->qpos[2]); //dy/dq2

		dr3dx(3,0) = 1.0; //dxd/dqd0
		dr3dx(4,0) = l1*cos(dmain->qpos[1]) + l2*cos(dmain->qpos[1] + dmain->qpos[2]); //dx/dq1
		dr3dx(5,0) = l2*cos(dmain->qpos[1] + dmain->qpos[2]); //dx/dq2

		dr4dx(4,0) = -l1*sin(dmain->qpos[1]) - l2*sin(dmain->qpos[1] + dmain->qpos[2]); //dy/dq1
		dr4dx(5,0) = -l2*sin(dmain->qpos[1] + dmain->qpos[2]); //dy/dq2

		dr5dx(1,0) = l1*cos(dmain->qpos[1]) + l2*cos(dmain->qpos[1] + dmain->qpos[2]); //dx/dq1
		dr5dx(2,0) = l2*cos(dmain->qpos[1] + dmain->qpos[2]); //dx/dq2

		C->lux = Matrix<double, nU, nX>::Zero();
		C->lx = Matrix<double, nX, 1>::Zero();

		Matrix<double, nX, 1> mX = Matrix<double, nX, 1>::Zero();
		GetStateFromMujoco(dmain, &mX);

		double xpos = mX(0,0) + l1*sin(mX(1,0)) + l2*sin(mX(1,0) + mX(2,0));
		double xposc = mX(0,0);
		double ypos = l1*cos(mX(1,0)) + l2*cos(mX(1,0) + mX(2,0));
		double xvel = mX(3,0) + mX(4,0)*l1*cos(mX(1,0)) + (mX(4,0) + mX(5,0))*cos(mX(1,0) + mX(2,0));
		double yvel = -mX(4,0)*l1*sin(mX(1,0)) - (mX(4,0) + mX(5,0))*sin(mX(1,0) + mX(2,0));

		C->lx = -m_dCostCoeff_x*2.0*(0.0 - xpos)*dr1dx;
		C->lx += -m_dCostCoeff_y*2.0*(l1+l2 - ypos)*dr2dx;
		C->lx += -m_dCostCoeff_vx*2.0*(0.0 - xvel)*dr3dx;
		C->lx += -m_dCostCoeff_vy*2.0*(0.0 - yvel)*dr4dx;
		C->lx += -m_dCostCoeff_xc*2.0*(xpos - xposc)*dr5dx;

		C->lxx = m_dCostCoeff_x*2.0*dr1dx*dr1dx.transpose();
		C->lxx += m_dCostCoeff_y*2.0*dr2dx*dr2dx.transpose();
		C->lxx += m_dCostCoeff_vx*2.0*dr3dx*dr3dx.transpose();
		C->lxx += m_dCostCoeff_vy*2.0*dr4dx*dr4dx.transpose();
		C->lxx += m_dCostCoeff_xc*2.0*dr5dx*dr5dx.transpose();
//		return;
//	}
//	C->lxx = Matrix<double, nX, nX>::Zero();
	C->lux = Matrix<double, nU, nX>::Zero();
//	C->lx = Matrix<double, nX, 1>::Zero();
	C->lu(0,0) = 2.0*m_dCostCoeff_u*dmain->ctrl[0];
	C->luu(0,0) = 2.0*m_dCostCoeff_u;
}

void DIP::IntegrateJacobians(const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B)
{
	//rotational velocities
	for (int i = 0; i < nQd; i++)
		(*A)(i,nQ+i) = 1.0;

	(*A) = (*A) * m_dControlTime_s;
	(*B) = (*B) * m_dControlTime_s;

	for (int i = 0; i < nX; i++)
		(*A)(i,i) += 1.0;
}


void DIP::GetStateFromMujoco(const mjData* d, Matrix<double, nX, 1>* mX)
{
	for (int i = 0; i < nQ; i++)
		(*mX)(i,0) = d->qpos[i];
	for (int i = 0; i < nQd; i++)
		(*mX)(nQ+i,0) = d->qvel[i];
}

void DIP::DrawTargetTrajectory(vector<StateMatrix> xnew)
{
	vector<Vector3d> pts;
	for (int i = 0; i < N; i++)
	{
		Vector3d pt = Vector3d::Zero();
		pt(0) = xnew[i](0,0) + 0.6*sin(xnew[i](1,0)) + 0.6*sin(xnew[i](1,0) + xnew[i](2,0));
		pt(2) = 0.6*cos(xnew[i](1,0)) + 0.6*cos(xnew[i](1,0) + xnew[i](2,0));
		pts.push_back(pt);
	}
	m_Vis->SetTrajPoints(pts);
}

void DIP::Jacobian(const mjModel* m, const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B, Cost* C, bool bTerminalState)
{

	/*
	 * This function performs finite differences about the current state x and action u.
	 *
	 * A (6 x 6)
	 * B (6 x 1)
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

    // save output for warmstart
    mju_copy(center, output, m->nv);
    mju_copy(warmstart, d->qacc_warmstart, m->nv);

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
