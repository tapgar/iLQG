/*
 * SLIP.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#include "SLIP.h"
#include <iostream>

using namespace Eigen;
using namespace std;

SLIP::SLIP(bool save_vid) {
	mujoco_initialized = false;
	Init();
	m_Vis = new Visualizer(mj_Model, save_vid, "SLIP");
}

SLIP::~SLIP() {
	// TODO Auto-generated destructor stub
}

void SLIP::Init() {

	// Activate mujoco and load the model if this is the first instance
	if (!mujoco_initialized) {
		mj_activate("/home/tapgar/.mujoco/mjkey.txt");
		char error[1000] = "Could not load binary model";
		mj_Model = mj_loadXML("slip.xml", 0, error, 1000);
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

	mj_forward(mj_Model, mj_Data);

	//35 qpos states... get rid of qws
	for (int i = 0; i < nQ; i++)
	{
		toMJposidx[i] = i;
		toMJvelidx[i] = i;
	}


    TorqueBounds << -50.0, 50, -300, 300, -50, 50, -300, 300;

    cout << "Torque Bounds: \n" << TorqueBounds << endl;

//    usleep(10e6);

}

void SLIP::GetInitU(vector<ActionMatrix>* init_u)
{
	ActionMatrix null_u = ActionMatrix::Zero();

	for (int i = 0; i < nU; i++)
		null_u(i,0) = 0.0;

//	null_u(1,0) = -50;
//	null_u(3,0) = -50;

	for (int i =0 ; i < N; i++)
	{
		(*init_u).push_back(null_u);
	}
}

void SLIP::GetStateActionCost(const mjModel* m, const mjData* dmain, Matrix<double, nX, 1> mX, Matrix<double, nU, 1> mU, double* l, bool bTerminalState)
{
	*l = 0;
	double cp_offset = 0.0;
	double cg_offsetXY = 0.0;
	double cg_offsetZ = 0.0;
	GetExtVars(m, dmain, &cp_offset, &cg_offsetXY, &cg_offsetZ);

	*l += m_dCostCoeff_xd*(pow(m_dTargetVel_mps - mX(nQ,0),2.0));// + pow(mX(nQ+1,0),2.0));

	*l += m_dCostCoeff_foot*smooth_abs(0.0 - cp_offset, m_dAlpha_foot);
	*l += m_dCostCoeff_torsoH*smooth_abs(0.0 - cg_offsetXY, m_dAlpha_torsoH);
	*l += m_dCostCoeff_torsoV*smooth_abs(0.0 - cg_offsetZ, m_dAlpha_torsoV);

	for (int k = 0; k < nU; k++)
	{
		*l += m_dCostCoeff_u*dmain->ctrl[k]*dmain->ctrl[k];
	}

}

void SLIP::UpdateCostDerivatives(const mjModel* m, const mjData* dmain, Cost* C, bool bTerminal)
{

	//dCMdx needs to be calculated before calling this

	C->lu = Matrix<double, nU, 1>::Zero();
	C->luu = Matrix<double, nU, nU>::Zero();
	C->lux = Matrix<double, nU, nX>::Zero();

	double cp_offset = 0.0;
	double cg_offsetXY = 0.0;
	double cg_offsetZ = 0.0;
	GetExtVars(m, dmain, &cp_offset, &cg_offsetXY, &cg_offsetZ);

	C->lx = m_dCostCoeff_foot*smooth_abs_deriv((cp_offset), m_dAlpha_foot)*dCPdx;
	C->lx += m_dCostCoeff_torsoH*smooth_abs_deriv((cg_offsetXY), m_dAlpha_torsoH)*dCMxydx;
	C->lx += m_dCostCoeff_torsoV*smooth_abs_deriv((cg_offsetZ), m_dAlpha_torsoV)*dCMzdx;

	C->lx(nQ,0) += -m_dCostCoeff_xd*2.0*(m_dTargetVel_mps - dmain->qvel[0]);
//	C->lx(nQ+1,0) += -m_dCostCoeff_xd*2.0*(0.0 - dmain->qvel[1]);

	C->lxx = m_dCostCoeff_foot*smooth_abs_second_deriv((cp_offset), m_dAlpha_foot)*dCPdx*dCPdx.transpose();
	C->lxx += m_dCostCoeff_torsoH*smooth_abs_second_deriv((cg_offsetXY), m_dAlpha_torsoH)*dCMxydx*dCMxydx.transpose();
	C->lxx += m_dCostCoeff_torsoV*smooth_abs_second_deriv((cg_offsetZ), m_dAlpha_torsoV)*dCMzdx*dCMzdx.transpose();

	C->lxx(nQ,nQ) += m_dCostCoeff_xd;
//	C->lxx(nQ+1,nQ+1) += m_dCostCoeff_xd;


	for (int i = 0; i < nU; i++)
	{
		C->lu(i,0) = m_dCostCoeff_u*dmain->ctrl[i];
		C->luu(i,i) = m_dCostCoeff_u;
	}
}

void SLIP::IntegrateJacobians(const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B)
{
    //translation derivatives
    for (int i = 0; i < nQ; i++)
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


void SLIP::GetStateFromMujoco(const mjData* d, Matrix<double, nX, 1>* mX)
{
	for (int i = 0; i < nQ; i++)
		(*mX)(i,0) = d->qpos[toMJposidx[i]];
	for (int i = 0; i < nQd; i++)
		(*mX)(nQ+i,0) = d->qvel[toMJvelidx[i]];
}

void SLIP::DrawTargetTrajectory(vector<StateMatrix> xnew)
{
	vector<Vector3d> pts;
	for (int i = 0; i < N; i++)
	{
		Vector3d pt = Vector3d::Zero();
		pt(0) = xnew[i](0,0);
		pt(1) = 0.0;
		pt(2) = xnew[i](1,0);
		pts.push_back(pt);
	}
	m_Vis->SetTrajPoints(pts);
}

void SLIP::Jacobian(const mjModel* m, const mjData* dmain, Matrix<double, nX, nX>* A, Matrix<double, nX, nU>* B, Cost* C, bool bTerminalState)
{
    int nv = m->nv;

    mjData* d = mj_makeData(m);
    // allocate stack space for result at center
    mjtNum* deriv = (mjtNum*) mju_malloc(sizeof(mjtNum)*m->nv*m->nv);
    mjMARKSTACK
    mjtNum* center = mj_stackAlloc(d, nv);
    mjtNum* warmstart = mj_stackAlloc(d, nv);

    Matrix<double, nQ, nQ> dpos = Matrix<double, nQ, nQ>::Zero();
    Matrix<double, nQ, nQ> dvel = Matrix<double, nQ, nQ>::Zero();

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

    Matrix<double, nQ, 1> dCPdq = Matrix<double, nQ, 1>::Zero();
    Matrix<double, nQ, 1> dCMxydq = Matrix<double, nQ, 1>::Zero();
    Matrix<double, nQ, 1> dCMzdq = Matrix<double, nQ, 1>::Zero();

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

    mjFREESTACK

    mj_deleteData(d);
}

void SLIP::GetExtVars(const mjModel* m, const mjData* d, double* cp_offset, double* cg_offsetXY, double* cg_offsetZ)
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

        (*cp_offset) = (left_foot(0)+right_foot(0))/2 - com(0);
        (*cg_offsetXY) = torso(0) - com(0);
        (*cg_offsetZ) = 0.7 -( torso(2)-(left_foot(2)+right_foot(2))/2.0);
}


void SLIP::GetBodyPosition(const mjData* dmain, int body_id, Vector3d* res)
{
	for (int i = 0; i < 3; i++)
		(*res)(i) = dmain->xipos[body_id*3 + i];
}

void SLIP::GetCOMPosition(const mjModel* m, const mjData* dmain, Vector3d* com)
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

void SLIP::UpdateControl(mjData* d, Matrix<double, nU, 1> targ_q)
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

void SLIP::LogClassSpecificVars(const mjModel* m, const mjData* dmain, ofstream* outfile)
{
	double cp_offset, cg_offsetXY, cg_offsetZ;
	GetExtVars(m, dmain, &cp_offset, &cg_offsetXY, &cg_offsetZ);

	(*outfile) << cp_offset << "," << cg_offsetXY << "," << cg_offsetZ << endl;

}


