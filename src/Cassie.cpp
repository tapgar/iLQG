/*
 * Cassie.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#include "Cassie.h"
#include <iostream>

using namespace Eigen;
using namespace std;

constexpr double Cassie::m_dCostCoeff_u[];

Cassie::Cassie(bool save_vid) {
	mujoco_initialized = false;
	Init();
	m_Vis = new Visualizer(mj_Model, save_vid, "Cassie");

        for (int j = 13; j < 16; j++)
        {
    for (int i = 0; i < 5; i++)
    {

           mj_Data->qpos[j] = -M_PI_2 + M_PI*double(i)/4.0;
           mj_forward(mj_Model, mj_Data);
           m_Vis->Draw(mj_Data);
        }
    }

}

Cassie::~Cassie() {
	// TODO Auto-generated destructor stub
}

void Cassie::Init() {

	// Activate mujoco and load the model if this is the first instance
	if (!mujoco_initialized) {
		mj_activate("/home/tapgar/.mujoco/mjkey.txt");
		char error[1000] = "Could not load binary model";
		mj_Model = mj_loadXML("/home/tapgar/libcassie-master/cassie.xml", 0, error, 1000);
		if (!mj_Model) {
			mju_error_s("Load model error: %s", error);
			return;
		}
		mujoco_initialized = true;
	}

	mj_Model->opt.timestep = m_dDeltaSimTime_s;
	mj_Model->opt.tolerance = 0.0;
	mj_Model->opt.iterations = 100;

	// Initialize mjData
	mj_Data = mj_makeData(mj_Model);

    double qpos_init[] =
		{-0.0305, 0, 0.4973, -1.1997, 0, 1.4267, -1.5968,
		 -1.5244, .6472, 0, 0.9785, -0.0164, 0.01787, -0.2049,
		 -0.0305, 0, 0.4973, -1.1997, 0, 1.4267, -1.5968,
		 -1.5244, 0.6472, 0, 0.9786, 0.00386, -0.01524, -0.2051};

    mju_copy(&mj_Data->qpos[7], qpos_init, mj_Model->nq-7);

	mj_forward(mj_Model, mj_Data);

	//35 qpos states... get rid of qws
	for (int i = 0; i < nQ; i++)
	{
		toMJposidx[i] = i;
		if (i >= 3)
			toMJposidx[i]++;
		if (i >= 16)
			toMJposidx[i]++;
		if (i >= 29)
			toMJposidx[i]++;
	}
	for (int i = 0; i < nQd; i++)
		toMJvelidx[i] = i;

	for (int i = 0; i < nQ; i++)
		printf("%d\t%d\n",toMJposidx[i], toMJvelidx[i]);

	for (int i = 0; i < 4; i++)
	{
		toMotorposidx[i] = i + 6;
		toMotorvelidx[i] = nQ + i + 6;
	}
	toMotorposidx[4] = 12;
	toMotorvelidx[4] = nQ + 12;
	for (int i = 0; i < 4; i++)
	{
		toMotorposidx[i+5] = i + 19;
		toMotorvelidx[i+5] = nQ + i + 19;
	}
	toMotorposidx[9] = 25;
	toMotorvelidx[9] = nQ + 25;

    GearN << 25, 25, 16, 16, 50, 25, 25, 16, 16, 50;
    //from cassie controller
    Kp << 300, 300, 300, 300, 30, 300, 300, 300, 300, 30;
    Kd << 4, 4, 12, 12, 3, 4, 4, 12, 12, 3;
    for (int i = 0; i < nU; i++)
    {
    	Kp(i,0) /= GearN(i,0);
    	Kd(i,0) /= GearN(i,0);
    }

    TorqueBounds << -4.5, 4.5, -4.5, 4.5, -12.2, 12.2, -12.2, 12.2, -0.9, 0.9,
    		-4.5, 4.5, -4.5, 4.5, -12.2, 12.2, -12.2, 12.2, -0.9, 0.9;

    cout << "Torque Bounds: \n" << TorqueBounds << endl;

//    usleep(10e6);

}

void Cassie::GetInitU(vector<ActionMatrix>* init_u)
{
	ActionMatrix null_u = ActionMatrix::Zero();

	/*standing
	for (int i = 0; i < nU; i++)
		null_u(i,0) = 0.0;

	null_u(0,0) = 0.1139;
	null_u(5,0) = 0.1139;
	null_u(2,0) = -0.3325;
	null_u(7,0) = -0.3325;
	null_u(3,0) = 2.484;
	null_u(8,0) = 2.484;
	null_u(4,0) = -0.084;
	null_u(9,0) = -0.084;

	for (int i =0 ; i < N; i++)
	{
		(*init_u).push_back(null_u);
	}
	*/

	ifstream file ("cassie_walking_seed_u.csv");
	string value;
	int c = 0;
	int line = 0;
	while (file.good())
	{
		line++;
		if (c == 9)
			getline(file, value, '\n');
		else
			getline(file, value, ',');
		null_u(c,0) = atof(value.c_str());

		if (++c == 10)
		{
			c = 0;
			(*init_u).push_back(null_u);
			if ((*init_u).size() == N)
				return;
		}
	}
}

void Cassie::GetStateActionCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, Matrix<double, nU, 1> mU, double* l, bool bTerminalState)
{
	*l = 0;
	double cp_offset = 0.0;
	double cg_offsetXY = 0.0;
	double cg_offsetZ = 0.0;
	GetExtVars(m, dmain, &cp_offset, &cg_offsetXY, &cg_offsetZ);

//	*l += m_dCostCoeff_foot*pow(cp_offset,2.0);
//	*l += m_dCostCoeff_torsoH*pow(cg_offsetXY,2.0);
//	*l += m_dCostCoeff_torsoV*pow(cg_offsetZ,2.0);
	*l += m_dCostCoeff_xd*(pow(m_dTargetVel_mps - mX(nQ,0),2.0) + pow(mX(nQ+1,0),2.0));

	*l += m_dCostCoeff_foot*smooth_abs(0.0 - cp_offset, m_dAlpha_foot);
	*l += m_dCostCoeff_torsoH*smooth_abs(0.0 - cg_offsetXY, m_dAlpha_torsoH);
	*l += m_dCostCoeff_torsoV*smooth_abs(0.0 - cg_offsetZ, m_dAlpha_torsoV);
	*l += m_dCostCoeff_q*smooth_abs(((double)dmain->qpos[4]+0.01288), m_dAlpha_q);
	*l += m_dCostCoeff_q*smooth_abs((double)dmain->qpos[5], m_dAlpha_q);
	*l += m_dCostCoeff_q*smooth_abs((double)dmain->qpos[6], m_dAlpha_q);

//	*l += m_dCostCoeff_footYaw*pow(mX(7,0) + mX(12,0), 2.0);

//	*l += m_dCostCoeff_q*smooth_abs((double)dmain->qpos[6], m_dAlpha_q);
	if (bTerminalState)
	{
//		*l += 2.0*m_dCostCoeff_torsoV*smooth_abs(0.0 - cg_offsetZ, m_dAlpha_torsoV);
//		*l += m_dCostCoeff_xd*smooth_abs(0.0 - mX(nQ,0), m_dAlpha_xd);
//		*l += m_dCostCoeff_xd*smooth_abs(0.0 - mX(nQ+1,0), m_dAlpha_xd);
//		*l += m_dCostCoeff_xd*smooth_abs(0.0 - mX(nQ+2,0), m_dAlpha_xd);
//		*l += 100.0*m_dCostCoeff_torsoV*smooth_abs(0.0 - cg_offsetZ, m_dAlpha_torsoV);
//		*l += 1.0*m_dCostCoeff_foot*smooth_abs(0.0 - cp_offset, m_dAlpha_foot);
		return;
	}
	for (int k = 0; k < nU; k++)
	{
		*l += m_dCostCoeff_u[k]*(GearN(k,0)*GearN(k,0)*dmain->ctrl[k]*dmain->ctrl[k]);
	}

}

void Cassie::UpdateCostDerivatives(const mjModel* m, const mjData* dmain, Cost* C, bool bTerminal)
{

	//dCMdx needs to be calculated before calling this

	C->lu = Matrix<double, nU, 1>::Zero();
	C->luu = Matrix<double, nU, nU>::Zero();
	C->lux = Matrix<double, nU, nX>::Zero();

	double cp_offset = 0.0;
	double cg_offsetXY = 0.0;
	double cg_offsetZ = 0.0;
	GetExtVars(m, dmain, &cp_offset, &cg_offsetXY, &cg_offsetZ);

//	C->lx = -m_dCostCoeff_foot*2.0*(0.0 - cp_offset)*dCPdx;
//	C->lx += -m_dCostCoeff_torsoH*2.0*(0.0 - cg_offsetXY)*dCMxydx;
//	C->lx += -m_dCostCoeff_torsoV*2.0*(0.0 - cg_offsetZ)*dCMzdx;

//
//	C->lxx = m_dCostCoeff_foot*2.0*dCPdx*dCPdx.transpose();
//	C->lxx += m_dCostCoeff_torsoH*2.0*dCMxydx*dCMxydx.transpose();
//	C->lxx += m_dCostCoeff_torsoV*2.0*dCMzdx*dCMzdx.transpose();
//	C->lxx(nQ,nQ) += 2.0*m_dCostCoeff_xd;
//	C->lxx(nQ+1,nQ+1) += 2.0*m_dCostCoeff_xd;

	C->lx = m_dCostCoeff_foot*smooth_abs_deriv((cp_offset), m_dAlpha_foot)*dCPdx;
	C->lx += m_dCostCoeff_torsoH*smooth_abs_deriv((cg_offsetXY), m_dAlpha_torsoH)*dCMxydx;
	C->lx += m_dCostCoeff_torsoV*smooth_abs_deriv((cg_offsetZ), m_dAlpha_torsoV)*dCMzdx;
	C->lx(3,0) += m_dCostCoeff_q*smooth_abs_deriv((double)dmain->qpos[4], m_dAlpha_q);
	C->lx(4,0) += m_dCostCoeff_q*smooth_abs_deriv(((double)dmain->qpos[5]+0.01288), m_dAlpha_q);
	C->lx(5,0) += m_dCostCoeff_q*smooth_abs_deriv(((double)dmain->qpos[6]), m_dAlpha_q);

	C->lx(nQ,0) += -m_dCostCoeff_xd*2.0*(m_dTargetVel_mps - dmain->qvel[0]);
	C->lx(nQ+1,0) += -m_dCostCoeff_xd*2.0*(0.0 - dmain->qvel[1]);
//	C->lx(nQ+2,0) += -m_dCostCoeff_xd*2.0*(0.0 - dmain->qvel[2]);

//	C->lx(7,0) += 2.0*m_dCostCoeff_footYaw*(dmain->qpos[8]+dmain->qpos[22]);
//	C->lx(12,0) += 2.0*m_dCostCoeff_footYaw*(dmain->qpos[8]+dmain->qpos[22]);


//	C->lx(nQ,0) += 2.0*m_dCostCoeff_xd*dmain->qvel[0];
//	C->lx(nQ+1,0) += 2.0*m_dCostCoeff_xd*dmain->qvel[1];
//	C->lx(nQ+2,0) += 2.0*m_dCostCoeff_xd*dmain->qvel[2];

//	C->lx(5,0) += m_dCostCoeff_q*smooth_abs_deriv((double)dmain->qpos[6], m_dAlpha_q);
	if (bTerminal)
	{
//		C->lx += 1.0*m_dCostCoeff_foot*smooth_abs_deriv((0.0 - cp_offset), m_dAlpha_foot)*dCPdx;
//		C->lx += 2.0*m_dCostCoeff_torsoV*smooth_abs_deriv(cg_offsetZ, m_dAlpha_torsoV)*dCMzdx;
//		C->lx(nQ,0) += m_dCostCoeff_xd*smooth_abs_deriv((0.0 - dmain->qvel[0]), m_dAlpha_xd);
//		C->lx(nQ+1,0) += m_dCostCoeff_xd*smooth_abs_deriv((0.0 - dmain->qvel[1]), m_dAlpha_xd);
//		C->lx(nQ+2,0) += m_dCostCoeff_xd*smooth_abs_deriv((0.0 - dmain->qvel[2]), m_dAlpha_xd);
	}
	C->lxx = m_dCostCoeff_foot*smooth_abs_second_deriv((cp_offset), m_dAlpha_foot)*dCPdx*dCPdx.transpose();
	C->lxx += m_dCostCoeff_torsoH*smooth_abs_second_deriv((cg_offsetXY), m_dAlpha_torsoH)*dCMxydx*dCMxydx.transpose();
	C->lxx += m_dCostCoeff_torsoV*smooth_abs_second_deriv((cg_offsetZ), m_dAlpha_torsoV)*dCMzdx*dCMzdx.transpose();
	C->lxx(3,3) += m_dCostCoeff_q*smooth_abs_second_deriv((double)dmain->qpos[4], m_dAlpha_q);
	C->lxx(4,4) += m_dCostCoeff_q*smooth_abs_second_deriv(((double)dmain->qpos[5] + 0.01288), m_dAlpha_q);
	C->lxx(5,5) += m_dCostCoeff_q*smooth_abs_second_deriv((double)dmain->qpos[6], m_dAlpha_q);

	C->lxx(nQ,nQ) += m_dCostCoeff_xd;
	C->lxx(nQ+1,nQ+1) += m_dCostCoeff_xd;
//	C->lxx(nQ+2,nQ+2) += m_dCostCoeff_xd;

//	C->lxx(7,7) += 2.0*m_dCostCoeff_footYaw;
//	C->lxx(12,12) += 2.0*m_dCostCoeff_footYaw;

	if (bTerminal)
	{
//		C->lxx += 100.0*m_dCostCoeff_foot*smooth_abs_second_deriv((0.0 - cp_offset), m_dAlpha_foot)*dCPdx*dCPdx.transpose();
//		C->lxx += 2.0*m_dCostCoeff_torsoV*smooth_abs_second_deriv(cg_offsetZ, m_dAlpha_torsoV)*dCMzdx*dCMzdx.transpose();
//		C->lxx(nQ,nQ) += m_dCostCoeff_xd*smooth_abs_second_deriv((0.0 - dmain->qvel[0]), m_dAlpha_xd);
//		C->lxx(nQ+1,nQ+1) += m_dCostCoeff_xd*smooth_abs_second_deriv((0.0 - dmain->qvel[1]), m_dAlpha_xd);
//		C->lxx(nQ+2,nQ+2) += m_dCostCoeff_xd*smooth_abs_second_deriv((0.0 - dmain->qvel[2]), m_dAlpha_xd);
		return;
	}

	for (int i = 0; i < nU; i++)
	{
		C->lu(i,0) = m_dCostCoeff_u[i]*GearN(i,0)*dmain->ctrl[i];
		C->luu(i,i) = m_dCostCoeff_u[i]*GearN(i,0);
	}
}

void Cassie::IntegrateQuaternion(const mjData* dmain, Matrix<double, nX, nX>* A, int start_pos_idx, int start_w_idx)
{
	//quaternion derivatives
	double qw = dmain->qpos[start_pos_idx];
	double qx = dmain->qpos[start_pos_idx+1];
	double qy = dmain->qpos[start_pos_idx+2];
	double qz = dmain->qpos[start_pos_idx+3];
	double wx = dmain->qvel[start_w_idx];
	double wy = dmain->qvel[start_w_idx+1];
	double wz = dmain->qvel[start_w_idx+2];

//    //qwdot
//    (*A)(3,nQ+3) = -0.5*qx; // wx
//    (*A)(3,nQ+4) = -0.5*qy; // wy
//    (*A)(3,nQ+5) = -0.5*qz; // wz
	//qxdot
	(*A)(start_w_idx,nQ+start_w_idx) = 0.5*qw; // wx
	(*A)(start_w_idx,nQ+start_w_idx+1) = 0.5*qz; // wy
	(*A)(start_w_idx,nQ+start_w_idx+2) = -0.5*qy; // wz
	//qydot
	(*A)(start_w_idx+1,nQ+start_w_idx) = -0.5*qz; // wx
	(*A)(start_w_idx+1,nQ+start_w_idx+1) = 0.5*qw; // wy
	(*A)(start_w_idx+1,nQ+start_w_idx+2) = 0.5*qx; // wz
	//qzdot
	(*A)(start_w_idx+2,nQ+start_w_idx) = 0.5*qy; // wx
	(*A)(start_w_idx+2,nQ+start_w_idx+1) = -0.5*qx; // wy
	(*A)(start_w_idx+2,nQ+start_w_idx+2) = 0.5*qw; // wz

//    //qwdot
//    (*A)(3,3) = 0.0;	 // qw
//	(*A)(3,4) = -0.5*wx; // qx
//	(*A)(3,5) = -0.5*wy; // qy
//	(*A)(3,6) = -0.5*wz; // qz
	//qxdot
//	(*A)(3,3) = 0.5*wx;	 // qw
	(*A)(start_w_idx,start_w_idx) = 0.0; // qx
	(*A)(start_w_idx,start_w_idx+1) = -0.5*wz; // qy
	(*A)(start_w_idx,start_w_idx+2) = 0.5*wy; // qz
	//qydot
//	(*A)(4,3) = 0.5*wy;	 // qw
	(*A)(start_w_idx+1,start_w_idx) = 0.5*wz; // qx
	(*A)(start_w_idx+1,start_w_idx+1) = 0.0; // qy
	(*A)(start_w_idx+1,start_w_idx+2) = -0.5*wx; // qz
	//qzdot
//	(*A)(5,3) = 0.5*wz;	 // qw
	(*A)(start_w_idx+2,start_w_idx) = -0.5*wy; // qx
	(*A)(start_w_idx+2,start_w_idx+1) = 0.5*wx; // qy
	(*A)(start_w_idx+2,start_w_idx+2) = 0.0; // qz
}

void Cassie::IntegrateJacobians(const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B)
{
    //translation derivatives
    for (int i = 0; i < nQ; i++)
    	(*A)(i,nQ+i) = 1.0;

    IntegrateQuaternion(dmain, A, 3, 3);
    IntegrateQuaternion(dmain, A, 17, 16);
    IntegrateQuaternion(dmain, A, 31, 29);
//
//    for (int i = 6; i < nQ; i++)
//    	(*A)(i, nQ+i) = 1.0;

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

	double cp_offset = 0.0;
	double cg_offsetxy = 0.0;
	double cg_offsetz = 0.0;
	GetExtVars(mj_Model, dmain, &cp_offset, &cg_offsetxy, &cg_offsetz);

	//printf("CP offset: %f\t%f\t%f\n", cp_offset, cg_offsetxy, cg_offsetz);

	cout << "dCPdx matrix" << "\n" <<  dCPdx << "\n";
	cout << "dCMxydx matrix" << "\n" <<  dCMxydx << "\n";
	cout << "dCMzdx matrix" << "\n" <<  dCMzdx << "\n";

	m_Vis->Draw((mjData*)dmain);
#endif

}


void Cassie::GetStateFromMujoco(const mjData* d, Matrix<double, nX, 1>* mX)
{
	for (int i = 0; i < nQ; i++)
		(*mX)(i,0) = d->qpos[toMJposidx[i]];
	for (int i = 0; i < nQd; i++)
		(*mX)(nQ+i,0) = d->qvel[toMJvelidx[i]];
}

void Cassie::DrawTargetTrajectory(vector<StateMatrix> xnew)
{
	vector<Vector3d> pts;
	for (int i = 0; i < N; i++)
	{
		Vector3d pt = Vector3d::Zero();
		pt(0) = xnew[i](0,0);
		pt(1) = xnew[i](1,0);
		pt(2) = xnew[i](2,0);
		pts.push_back(pt);
	}
	m_Vis->SetTrajPoints(pts);
}

void Cassie::Jacobian(const mjModel* m, const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B, Cost* C, bool bTerminalState)
{
    int nv = m->nv;

    mjData* d = mj_makeData(m);
    // allocate stack space for result at center
    mjtNum* deriv = (mjtNum*) mju_malloc(sizeof(mjtNum)*m->nv*m->nv);
    mjMARKSTACK
    mjtNum* center = mj_stackAlloc(d, nv);
    mjtNum* warmstart = mj_stackAlloc(d, nv);

    Matrix<double, 32, 32> dpos = Matrix<double, 32, 32>::Zero();
    Matrix<double, 32, 32> dvel = Matrix<double, 32, 32>::Zero();
    Matrix<double, nX, nU> dacc = Matrix<double, nX, nU>::Zero();

    // prepare static schedule: range of derivative columns to be computed by this thread
    int isforward = 1;
    int istart = 0;
    int iend = m->nv;

    // copy state and control from dmain to thread-specific d
    d->time = dmain->time;
    mju_copy(d->qpos, dmain->qpos, m->nq);
    mju_copy(d->qvel, dmain->qvel, m->nv);
    mju_copy(d->qacc, dmain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dmain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dmain->ctrl, m->nu);

    // run full computation at center point (usually faster than copying dmain)
    if( isforward )
    {
        mj_forward(m, d);

        // extra solver iterations to improve warmstart (qacc) at center point
        for( int rep=1; rep<3; rep++ )
            mj_forwardSkip(m, d, mjSTAGE_VEL, 1);
    }
    else
        mj_inverse(m, d);

    // select output from forward or inverse dynamics
    mjtNum* output = (isforward ? d->qacc : d->qfrc_inverse);

//    double orig_cp_offset = 0.0;
//    double orig_cg_offsetXY = 0.0;
//    double orig_cg_offsetZ = 0.0;
//    GetExtVars(m, d, &orig_cp_offset, &orig_cg_offsetXY, &orig_cg_offsetZ);

    // save output for center point and warmstart (needed in forward only)
    mju_copy(center, output, nv);
    mju_copy(warmstart, d->qacc_warmstart, nv);

    Matrix<double, 32, 1> dCPdq = Matrix<double, 32, 1>::Zero();
    Matrix<double, 32, 1> dCMxydq = Matrix<double, 32, 1>::Zero();
    Matrix<double, 32, 1> dCMzdq = Matrix<double, 32, 1>::Zero();

    // select target vector and original vector for force or acceleration derivative
    mjtNum* target = d->ctrl;
    const mjtNum* original = dmain->ctrl;

    // finite-difference over force or acceleration: skip = mjSTAGE_VEL
    for( int i=0; i<nU; i++ )
    {
        // perturb selected target
        target[i] -= m_deps;

        // evaluate dynamics, with center warmstart
        if( isforward )
        {
            mju_copy(d->qacc_warmstart, warmstart, m->nv);
            mj_forwardSkip(m, d, mjSTAGE_VEL, 1);
        }
        else
            mj_inverseSkip(m, d, mjSTAGE_VEL, 1);

        // undo perturbation
        target[i] = original[i];
        mju_copy(center, output, nv);

        // perturb selected target
		target[i] += m_deps;

		// evaluate dynamics, with center warmstart
		if( isforward )
		{
			mju_copy(d->qacc_warmstart, warmstart, m->nv);
			mj_forwardSkip(m, d, mjSTAGE_VEL, 1);
		}
		else
			mj_inverseSkip(m, d, mjSTAGE_VEL, 1);

		target[i] = original[i];

        // compute column i of derivative 2
        for( int j=0; j<nQd; j++ )
            (*B)(j+nQ,i) = (output[toMJvelidx[j]] - center[toMJvelidx[j]])/(2.0*m_deps);
    }

    // finite-difference over velocity: skip = mjSTAGE_POS
    for( int i=0; i<m->nv; i++ )
    {
        // perturb velocity
        d->qvel[i] -= m_deps;

        // evaluate dynamics, with center warmstart
        if( isforward )
        {
            mju_copy(d->qacc_warmstart, warmstart, m->nv);
            mj_forwardSkip(m, d, mjSTAGE_POS, 1);
        }
        else
            mj_inverseSkip(m, d, mjSTAGE_POS, 1);

        // undo perturbation
        d->qvel[i] = dmain->qvel[i];
        mju_copy(center, output, nv);

        // perturb velocity
		d->qvel[i] += m_deps;

		// evaluate dynamics, with center warmstart
		if( isforward )
		{
			mju_copy(d->qacc_warmstart, warmstart, m->nv);
			mj_forwardSkip(m, d, mjSTAGE_POS, 1);
		}
		else
			mj_inverseSkip(m, d, mjSTAGE_POS, 1);

		// undo perturbation
		d->qvel[i] = dmain->qvel[i];

		for (int j = 0; j < m->nv; j++)
			dvel(j,i) = (output[j] - center[j])/(2.0*m_deps);

        // compute column i of derivative 1
//        for( int j=0; j<nQd; j++ )
//            (*A)(nQ+j,nQ+i) = (output[toMJvelidx[j]] - center[toMJvelidx[j]])/(2.0*m_deps);
    }

    // finite-difference over position: skip = mjSTAGE_NONE
    for( int i=istart; i<iend; i++ )
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

        int freeoffset = 0;
        if (m->jnt_type[jid]==mjJNT_FREE)
        	freeoffset = i;

        // apply quaternion or simple perturbation
        if( quatadr>=0 )
        {
            mjtNum angvel[3] = {0,0,0};
            angvel[dofpos] = -m_deps;
            mju_quatIntegrate(d->qpos+quatadr, angvel, 1);
        }
        else
        	d->qpos[m->jnt_qposadr[jid] + i - m->jnt_dofadr[jid]] -= m_deps;

        // evaluate dynamics, with center warmstart
        if( isforward )
        {
            mju_copy(d->qacc_warmstart, warmstart, m->nv);
            mj_forwardSkip(m, d, mjSTAGE_NONE, 1);
        }
        else
            mj_inverseSkip(m, d, mjSTAGE_NONE, 1);

		double cp_offsetn = 0.0;
		double cg_offsetXYn = 0.0;
		double cg_offsetZn = 0.0;
		GetExtVars(m, d, &cp_offsetn, &cg_offsetXYn, &cg_offsetZn);

//		m_Vis->Draw(d);


        // undo perturbation
        mju_copy(d->qpos, dmain->qpos, m->nq);
        mju_copy(center, output, nv);

        // apply quaternion or simple perturbation
		if( quatadr>=0 )
		{
			mjtNum angvel[3] = {0,0,0};
			angvel[dofpos] = m_deps;
			mju_quatIntegrate(d->qpos+quatadr, angvel, 1);
		}
		else
			d->qpos[m->jnt_qposadr[jid] + i - m->jnt_dofadr[jid]] += m_deps;

		// evaluate dynamics, with center warmstart
		if( isforward )
		{
			mju_copy(d->qacc_warmstart, warmstart, m->nv);
			mj_forwardSkip(m, d, mjSTAGE_NONE, 1);
		}
		else
			mj_inverseSkip(m, d, mjSTAGE_NONE, 1);

//		m_Vis->Draw(d);

		double cp_offsetp = 0.0;
		double cg_offsetXYp = 0.0;
		double cg_offsetZp = 0.0;
		GetExtVars(m, d, &cp_offsetp, &cg_offsetXYp, &cg_offsetZp);

		dCPdq(i,0) = (cp_offsetp - cp_offsetn)/(2.0*m_deps);
		dCMxydq(i,0) = (cg_offsetXYp - cg_offsetXYn)/(2.0*m_deps);
		dCMzdq(i,0) = (cg_offsetZp - cg_offsetZn)/(2.0*m_deps);

//		printf("%d\t%d\t%d\t%f\t%f\n", i,m->jnt_qposadr[jid],freeoffset,cg_offsetZp, cg_offsetZn);

		// undo perturbation
		mju_copy(d->qpos, dmain->qpos, m->nq);

        // compute column i of derivative 0
        for( int j=0; j<nv; j++ )
        {
        	deriv[j*nv + i] = (output[j] - center[j])/(2.0*m_deps);
        	dpos(j,i) = (output[j] - center[j])/(2.0*m_deps);
        }
    }

    for (int i = 0; i < nQd; i++)
    {
    	for (int j = 0; j < nQd; j++)
    	{
    		(*A)(i+nQ, j) = dpos(toMJvelidx[i], toMJvelidx[j]);
    		(*A)(i+nQ, j+nQ) = dvel(toMJvelidx[i], toMJvelidx[j]);
    	}
    	dCPdx(i, 0) = dCPdq(toMJvelidx[i], 0);
		dCMxydx(i, 0) = dCMxydq(toMJvelidx[i], 0);
		dCMzdx(i, 0) = dCMzdq(toMJvelidx[i], 0);
    }

    UpdateCostDerivatives(m, dmain, C, bTerminalState);

//    cout << "dpos:\n" << dpos << endl;
//    cout << "dvel:\n" << dvel << endl;
//
//    printf("deriv:\n");
//    mju_printMat(deriv, m->nv, m->nv);


//    cout << "A:\n" << (*A) << endl;
//    cout << "B:\n" << (*B) << endl;

    mjFREESTACK

    mj_deleteData(d);
}

//void Cassie::Jacobian(const mjModel* m, const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B, Cost* C, bool bTerminalState)
//{
//
//	/*
//	 * This function performs finite differences about the current state x and action u.
//	 *
//	 * A (12 x 12)
//	 * B (12 x 1)
//	 *
//	 * The output space u is actually reference joint positions (10) so the resulting actuator forces
//	 * will be obtained by applying a PD controller about each of the joints.  To differentiate x_dot
//	 * w.r.t control inputs u we will simulate the PD controller to find the equivalent motor torque
//	 * then perform finite differences.  Quaternion differentiation is performed following mujoco
//	 * derivative example.
//	 *
//	 * In the example derivative.cpp code they also allowed a couple iterations to improve qacc estimate
//	 * which they termed warmup iterations.  The mujoco states passed to this function were collected
//	 * during the rollout so the qacc value should already be "warmed up".
//	 *
//	 */
//
//	mjData* d = mj_makeData(m);
//	mjMARKSTACK
//	mjtNum* center = mj_stackAlloc(d, m->nv);
//    mjtNum* warmstart = mj_stackAlloc(d, m->nv);
//
//	d->time = dmain->time;
//	mju_copy(d->qpos, dmain->qpos, m->nq);
//	mju_copy(d->qvel, dmain->qvel, m->nv);
//	mju_copy(d->qacc, dmain->qacc, m->nv);
//	mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, m->nv);
//	mju_copy(d->userdata, dmain->userdata, m->nuserdata);
//	mju_copy(d->qfrc_applied, dmain->qfrc_applied, m->nv);
//	mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
//	mju_copy(d->ctrl, dmain->ctrl, m->nu);
//
//	mj_forward(m, d);
//
//	//warmup
//	for( int rep=1; rep<3; rep++ )
//	   mj_forwardSkip(m, d, mjSTAGE_VEL, 1);
//
//	// select output from forward dynamics
//    mjtNum* output = d->qacc;
//    double orig_cp_offset = 0.0;
//    double orig_cg_offsetXY = 0.0;
//    double orig_cg_offsetZ = 0.0;
//    GetExtVars(m, d, &orig_cp_offset, &orig_cg_offsetXY, &orig_cg_offsetZ);
//
//    // save output for warmstart
//    mju_copy(center, output, m->nv);
//    mju_copy(warmstart, d->qacc_warmstart, m->nv);
//
//    dCPdx = Matrix<double, nX, 1>::Zero();
//    dCMxydx = Matrix<double, nX, 1>::Zero();
//    dCMzdx = Matrix<double, nX, 1>::Zero();
//
//
//    Matrix<double, nU, 1> m_Original;
//    for (int k = 0; k < nU; k++)
//    	m_Original(k,0) = d->userdata[k];
//
//   // cout << "Original" << m_Original.transpose() << endl;
//    Matrix<double, nU, 1> m_Target = m_Original;
//
//	// finite-difference over reference angle: skip = mjSTAGE_VEL
//	for( int i=0; i<nU; i++ )
//	{
//		m_Target(i,0) += m_deps;
//		UpdateControl(d, m_Target);
//
//		mju_copy(d->qacc_warmstart, warmstart, m->nv);
//		mj_forwardSkip(m, d, mjSTAGE_VEL, 1);
//
//		m_Target(i,0) -= m_deps;
//
//		for (int j = 0; j < nQd; j++)
//			(*B)(nQ+j,i) = (output[toMJvelidx[j]] - center[toMJvelidx[j]])/(m_deps);
//	}
//
//	mju_copy(d->ctrl, dmain->ctrl, m->nu);
//
//	// finite-difference over velocity: skip = mjSTAGE_POS
//	for( int i=0; i < nQd; i++ )
//	{
//		//first subtract eps
//		d->qvel[toMJvelidx[i]] += m_deps;
//		UpdateControl(d, m_Original);
//
//		mju_copy(d->qacc_warmstart, warmstart, m->nv);
//		mj_forwardSkip(m, d, mjSTAGE_POS, 1);
//
//		// undo perturbation
//		mju_copy(d->qvel, dmain->qvel, m->nv);
//
//		for (int j = 0; j < nQd; j++)
//		{
//			(*A)(nQ+j,nQ+i) = (output[toMJvelidx[j]] - center[toMJvelidx[j]])/(m_deps);
//		}
//	}
//	mju_copy(d->ctrl, dmain->ctrl, m->nu);
//
//	// finite-difference over position: skip = mjSTAGE_NONE
//    for( int i=0; i<nQ; i++ )
//    {
//    	int k = toMJposidx[i];
//    	// get joint id for this dof
//	    int jid = m->dof_jntid[k];
//
//		// get quaternion address and dof position within quaternion (-1: not in quaternion)
//		int quatadr = -1, dofpos = 0;
//		if( m->jnt_type[jid]==mjJNT_BALL )
//		{
//			quatadr = m->jnt_qposadr[jid];
//			dofpos = k - m->jnt_dofadr[jid];
//		}
//		else if( m->jnt_type[jid]==mjJNT_FREE && k>=m->jnt_dofadr[jid]+3 )
//		{
//			quatadr = m->jnt_qposadr[jid] + 3;
//			dofpos = k - m->jnt_dofadr[jid] - 3;
//		}
//
//		// apply quaternion or simple perturbation
//		if( quatadr>=0 )
//		{
//			mjtNum angvel[3] = {0,0,0};
//			angvel[dofpos] = m_deps;
//			mju_quatIntegrate(d->qpos+quatadr, angvel, 1);
//		}
//		else
//			d->qpos[k] += m_deps;
//
//		UpdateControl(d, m_Original);
//
//		// evaluate dynamics, with center warmstart
//		mju_copy(d->qacc_warmstart, warmstart, m->nv);
//		mj_forwardSkip(m, d, mjSTAGE_NONE, 1);
//
//		double cp_offset = 0.0;
//		double cg_offsetXY = 0.0;
//		double cg_offsetZ = 0.0;
//		GetExtVars(m, d, &cp_offset, &cg_offsetXY, &cg_offsetZ);
//
//		//printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\n", i, orig_cp_offset, cp_offset, orig_cg_offsetXY, cg_offsetXY, orig_cg_offsetZ, cg_offsetZ);
//
//		dCPdx(i,0) = (cp_offset - orig_cp_offset)/m_deps;
//		dCMxydx(i,0) = (cg_offsetXY - orig_cg_offsetXY)/m_deps;
//		dCMzdx(i,0) = (cg_offsetZ - orig_cg_offsetZ)/m_deps;
//
//		// undo perturbation
//		mju_copy(d->qpos, dmain->qpos, m->nq);
//
//        for (int j = 0; j < nQd; j++)
//		{
//			(*A)(nQ+j,i) = (output[toMJvelidx[j]] - center[toMJvelidx[j]])/(m_deps);
//		}
//    }
//    mju_copy(d->ctrl, dmain->ctrl, m->nu);
//
//    UpdateCostDerivatives(m, dmain, C, bTerminalState);
//
//    mjFREESTACK
//
//    mj_deleteData(d);
//
//}

void Cassie::GetExtVars(const mjModel* m, const mjData* d, double* cp_offset, double* cg_offsetXY, double* cg_offsetZ)
{
        Vector3d left_foot;
        Vector3d right_foot;
        Vector3d torso;
        Vector3d com;

        GetBodyPosition(d, m_nLeftFootBodyID, &left_foot);
        GetBodyPosition(d, m_nRightFootBodyID, &right_foot);
        GetBodyPosition(d, m_nTorsoBodyID, &torso);
        GetCOMPosition(m, d, &com);

//        printf("Left Foot: %f\t%f\t%f\n",left_foot(0), left_foot(1), left_foot(2));
//        printf("Right Foot: %f\t%f\t%f\n",right_foot(0), right_foot(1), right_foot(2));
//        printf("Torso: %f\t%f\t%f\n",torso(0), torso(1), torso(2));
//        printf("COM: %f\t%f\t%f\n",com(0), com(1), com(2));

        (*cp_offset) = 0.01154 - sqrt(pow((left_foot(0)+right_foot(0))/2 - com(0),2) + pow((left_foot(1)+right_foot(1))/2 - com(1),2));
        (*cg_offsetXY) = 0.06826 - sqrt(pow(torso(0) - com(0),2) + pow(torso(1) - com(1),2));
        (*cg_offsetZ) = 1.013 - (torso(2)-(left_foot(2)+right_foot(2))/2.0);

}


void Cassie::GetBodyPosition(const mjData* dmain, int body_id, Vector3d* res)
{
	for (int i = 0; i < 3; i++)
		(*res)(i) = dmain->xipos[body_id*3 + i];
}

void Cassie::GetCOMPosition(const mjModel* m, const mjData* dmain, Vector3d* com)
{
	double dTotalMass = 0;

	for (int i = 0; i < 3; i++)
		(*com)(i) = 0.0;

	for (int i = 1; i < m->nbody; i++)
	{
		dTotalMass += m->body_mass[i];
		(*com)(0) += dmain->xipos[i*3]*m->body_mass[i];
		(*com)(1) += dmain->xipos[i*3 + 1]*m->body_mass[i];
		(*com)(2) += dmain->xipos[i*3 + 2]*m->body_mass[i];
	}

	for (int i = 0; i < 3; i++)
		(*com)(i) /= dTotalMass;
}

void Cassie::UpdateControl(mjData* d, Matrix<double, nU, 1> targ_q)
{
	Matrix<double, nX, 1> cur_x;
	GetStateFromMujoco(d, &cur_x);
	for (int i = 0; i < nU; i++)
	{
//		int pos_idx = toMotorposidx[i];
//		int vel_idx = toMotorvelidx[i];
		d->ctrl[i] = targ_q(i,0);//Kp(i,0)*(targ_q(i,0) - cur_x(pos_idx,0)) + Kd(i,0)*(0.0 - cur_x(vel_idx,0)); //joint PD controller
	}
}

void Cassie::LogClassSpecificVars(const mjModel* m, const mjData* dmain, ofstream* outfile)
{
	double cp_offset, cg_offsetXY, cg_offsetZ;
	GetExtVars(m, dmain, &cp_offset, &cg_offsetXY, &cg_offsetZ);

	(*outfile) << cp_offset << "," << cg_offsetXY << "," << cg_offsetZ << endl;

}


