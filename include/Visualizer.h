/*
 * Visualizer.h
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#ifndef VISUALIZER_H_
#define VISUALIZER_H_

#include "mujoco.h"
#include "glfw3.h"
#include "stdio.h"
#include <vector>
#include "Common_Structs.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Visualizer {
public:
	Visualizer(mjModel *m, bool save_vid, const char* win_title);
	virtual ~Visualizer();

	int Init(bool save_video, const char* win_title);
	void Close();
	void Scroll(double xoffset, double yoffset);
	void Mouse_Button(int button, int act, int mods);
	void Mouse_Move(double xpos, double ypos);
	void Keyboard(int key, int scancode, int act, int mods);

	void Draw(mjData* data);
	void DrawWithPoints(mjData* data, vector<Vector3d> pts);

	void GetDisturbanceForce(double* Fx, double* Fy, double* Fz);

	void SetTrajPoints(vector<Vector3d> pts)
	{
		traj_pts = pts;
		m_bTrajInit = true;
	}

private:

	mjModel *mj_Model;

	GLFWwindow* m_Window;
	mjvCamera mj_Cam;
	mjvOption mj_Opt;
	mjvScene mj_Scn;
	mjrContext mj_Con;
	unsigned char* m_image_rgb;
	float* m_image_depth;

	int m_Width;
	int m_Height;
	bool m_bSaveVideo;

	FILE *fp;

	bool button_left = false;
	bool button_middle = false;
	bool button_right =  false;
	double cursor_lastx = 0;
	double cursor_lasty = 0;

	bool bOKtoComplete;
	bool bWaitForUserFeedback;


	vector<Vector3d> traj_pts;
	bool m_bTrajInit;

};

#endif /* VISUALIZER_H_ */
