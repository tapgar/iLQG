/*
 * Visualizer.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#include "Visualizer.h"

using namespace std;
using namespace Eigen;

Visualizer::Visualizer(mjModel *m, bool save_vid, const char* win_title) {
	mj_Model = m;
	m_Width = 1200;
	m_Height = 900;

	m_bTrajInit = false;

	Vector3d temp = Vector3d::Zero();
	for (int i = 0; i < N; i++)
	{
		traj_pts.push_back(temp);
	}

	Init(save_vid, win_title);

	bWaitForUserFeedback = PAUSE_VIS;
	bOKtoComplete = false;
}

Visualizer::~Visualizer() {
	// TODO Auto-generated destructor stub
}

static void window_close_callback(GLFWwindow* window)
{
	((Visualizer*)(glfwGetWindowUserPointer(window)))->Close();
}
static void scroll(GLFWwindow* window, double xoffset, double yoffset)
{
	((Visualizer*)(glfwGetWindowUserPointer(window)))->Scroll(xoffset, yoffset);
}
static void mouse_move(GLFWwindow* window, double xpos, double ypos)
{
	((Visualizer*)(glfwGetWindowUserPointer(window)))->Mouse_Move(xpos, ypos);
}
static void mouse_button(GLFWwindow* window, int button, int act, int mods)
{
	((Visualizer*)(glfwGetWindowUserPointer(window)))->Mouse_Button(button, act, mods);
}

static void keyboard(GLFWwindow* window, int key, int scancode, int act, int mods)
{
	((Visualizer*)(glfwGetWindowUserPointer(window)))->Keyboard(key, scancode, act, mods);
}


int Visualizer::Init(bool save_video, const char* win_title) {

	if (!glfwInit()) {
		mju_error("Could not initialize GLFW");
		return 1;
	}

	// Create window
	m_Window = glfwCreateWindow(m_Width, m_Height, win_title, NULL, NULL);
	glfwMakeContextCurrent(m_Window);
	glfwSwapInterval(1);

//	// Set up mujoco visualization objects
//	mj_Cam.lookat[0] = mj_Model->stat.center[0];
//	mj_Cam.lookat[1] = mj_Model->stat.center[1];
//	mj_Cam.lookat[2] = 1.0 + mj_Model->stat.center[2];
////
////	mj_Cam.distance = 0.5 * mj_Model->stat.extent;
////	mjv_moveCamera(mj_Model, mjMOUSE_ROTATE_H, 0.5, 0.0, &mj_Scn, &mj_Cam);
//	mj_Cam.type = mjCAMERA_FREE;
//
////	mj_Cam.type = mjCAMERA_TRACKING;
////	mj_Cam.trackbodyid = 1;
//	mj_Cam.distance = 0.3 * mj_Model->stat.extent;
//	mjv_moveCamera(mj_Model, mjMOUSE_ROTATE_H, 0.75, 0.0, &mj_Scn, &mj_Cam);

	mj_Cam.type = mjCAMERA_FIXED;
	mj_Cam.fixedcamid = 0;
	mjv_defaultOption(&mj_Opt);
	mj_Opt.flags[10] = 1;
	mj_Opt.flags[11] = 1;
	mj_Opt.flags[8] = 1;
//	mj_Opt.flags[15] = 1;
//	mj_Opt.flags[2] = 1;
//	mj_Opt.flags[3] = 1;

	//int present = glfwJoystickPresent(GLFW_JOYSTICK_1);
	//printf("Joystick present %d\n", present);

	mjr_defaultContext(&mj_Con);
	mjv_makeScene(&mj_Scn, 1E5);

	mjr_makeContext(mj_Model, &mj_Con, mjFONTSCALE_100);

	// Set callback for user-initiated window close events
	glfwSetWindowUserPointer(m_Window, this);
	glfwSetWindowCloseCallback(m_Window, window_close_callback);
	glfwSetCursorPosCallback(m_Window, mouse_move);
	glfwSetMouseButtonCallback(m_Window, mouse_button);
	glfwSetScrollCallback(m_Window, scroll);
	glfwSetKeyCallback(m_Window, keyboard);

	if (save_video)
	{
		m_image_rgb = (unsigned char*)malloc(3*m_Width*m_Height);
		m_image_depth = (float*)malloc(sizeof(float)*m_Width*m_Height);

		// create output rgb file
		fp = fopen("out/temp.out", "wb");
		if( !fp )
			mju_error("Could not open rgbfile for writing");
	}

	m_bSaveVideo = save_video;

	return 0;
}

void Visualizer::GetDisturbanceForce(double* Fx, double* Fy, double* Fz)
{
	int count;
	const float* axes = glfwGetJoystickAxes(GLFW_JOYSTICK_1, &count);

	(*Fx) = -axes[2]*200.0;
	(*Fy) = -axes[0]*200.0;
	(*Fz) = -axes[3]*200.0;
}

void Visualizer::DrawWithPoints(mjData* data, vector<Vector3d> pts)
{
	// Return early if window is closed
	if (!m_Window)
		return;

	// Set up for rendering
	glfwMakeContextCurrent(m_Window);
	mjrRect viewport = {0, 0, 0, 0};
	glfwGetFramebufferSize(m_Window, &viewport.width, &viewport.height);

	// Render scene
	mjv_updateScene(mj_Model, data, &mj_Opt, NULL, &mj_Cam, mjCAT_ALL, &mj_Scn);

	for (int i = 0; i < pts.size(); i++)
	{
		mjtNum pos[3];
		for (int j = 0; j < 3; j++)
			pos[j] = pts[i](j);
		mjtNum size[3] = {0.05,0.05,0.05};
		mjv_initGeom(&(mj_Scn.geoms[mj_Scn.ngeom]), mjGEOM_SPHERE, size, pos, NULL, NULL );
		mj_Scn.ngeom++;
	}
	mjv_addGeoms(mj_Model, data, &mj_Opt, NULL, mjCAT_DECOR, &mj_Scn);
	mjr_render(viewport, &mj_Scn, &mj_Con);

	// Show updated scene
	glfwSwapBuffers(m_Window);
	glfwPollEvents();

	if (m_bSaveVideo)
	{
		mjr_readPixels(m_image_rgb, m_image_depth, viewport, &mj_Con);
		// insert subsampled depth image in lower-left corner of rgb image
		const int NS = 3;           // depth image sub-sampling
		for( int r=0; r<m_Height; r+=NS )
		   for( int c=0; c<m_Width; c+=NS )
		   {
			  int adr = (r/NS)*m_Width + c/NS;
			  m_image_rgb[3*adr] = m_image_rgb[3*adr+1] = m_image_rgb[3*adr+2] = (unsigned char)((1.0f-m_image_depth[r*m_Width+c])*255.0f);
		   }

		 // write rgb image to file
		 fwrite(m_image_rgb, 3, m_Width*m_Height, fp);
	}

	if (bWaitForUserFeedback)
	{
		while (!bOKtoComplete)
			glfwPollEvents();
		bOKtoComplete = false;
	}
}

void Visualizer::Draw(mjData* data) {
	// Return early if window is closed
	if (!m_Window)
		return;

	// Set up for rendering
	glfwMakeContextCurrent(m_Window);
	mjrRect viewport = {0, 0, 0, 0};
	glfwGetFramebufferSize(m_Window, &viewport.width, &viewport.height);

	// Render scene
//	printf("Before Update Geoms: %d\n", mj_Scn.ngeom);
	mjv_updateScene(mj_Model, data, &mj_Opt, NULL, &mj_Cam, mjCAT_ALL, &mj_Scn);
//	printf("After Update Geoms: %d\n", mj_Scn.ngeom);

	if (m_bTrajInit)
	{

		for (int i = 0; i < N; i++)
		{
	//		traj_pts.push_back(&(mj_Scn.geoms[mj_Scn.ngeom+i]));
			//printf("Vis: %f, %f\n", (traj_pts[i])[0], (traj_pts[i])[2]);
			mjtNum pos[3];
			for (int j = 0; j < 3; j++)
				pos[j] = traj_pts[i](j);
			mjtNum size[3] = {0.05,0.05,0.05};
			mjv_initGeom(&(mj_Scn.geoms[mj_Scn.ngeom]), mjGEOM_SPHERE, size, pos, NULL, NULL );
			mj_Scn.ngeom++;
		}
//		printf("After Init Geoms: %d\n", mj_Scn.ngeom);
		mjv_addGeoms(mj_Model, data, &mj_Opt, NULL, mjCAT_DECOR, &mj_Scn);
//		printf("After Add Geoms: %d\n", mj_Scn.ngeom);
	}
//	mjtNum pos[3];
//	for (int j = 0; j < 3; j++)
//		pos[j] = data->xipos[4*3 + j];
//
//	mjtNum size[3] = {0.1,0.05,0.05};
//	mjv_initGeom(&(mj_Scn.geoms[mj_Scn.ngeom]), mjGEOM_SPHERE, size, pos, NULL, NULL );
//	mj_Scn.ngeom++;
//	mjv_addGeoms(mj_Model, data, &mj_Opt, NULL, mjCAT_DECOR, &mj_Scn);


	mjr_render(viewport, &mj_Scn, &mj_Con);


	// Show updated scene
	glfwSwapBuffers(m_Window);
	glfwPollEvents();


//	for (int i = 0; i < count; i++)
//		printf("Axis: %d\t%f\n", i, axes[i]);

	if (m_bSaveVideo)
	{
		mjr_readPixels(m_image_rgb, m_image_depth, viewport, &mj_Con);
		// insert subsampled depth image in lower-left corner of rgb image
		const int NS = 3;           // depth image sub-sampling
		for( int r=0; r<m_Height; r+=NS )
		   for( int c=0; c<m_Width; c+=NS )
		   {
			  int adr = (r/NS)*m_Width + c/NS;
			  m_image_rgb[3*adr] = m_image_rgb[3*adr+1] = m_image_rgb[3*adr+2] = (unsigned char)((1.0f-m_image_depth[r*m_Width+c])*255.0f);
		   }

		 // write rgb image to file
		 fwrite(m_image_rgb, 3, m_Width*m_Height, fp);
	}

	if (bWaitForUserFeedback)
	{
		while (!bOKtoComplete)
			glfwPollEvents();
		bOKtoComplete = false;
	}
}

// mouse button
void Visualizer::Mouse_Button(int button, int act, int mods)
{
    // update button state
    button_left =   (glfwGetMouseButton(m_Window, GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS);
    button_middle = (glfwGetMouseButton(m_Window, GLFW_MOUSE_BUTTON_MIDDLE)==GLFW_PRESS);
    button_right =  (glfwGetMouseButton(m_Window, GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS);

    // update mouse position
    glfwGetCursorPos(m_Window, &cursor_lastx, &cursor_lasty);
}



// mouse move
void Visualizer::Mouse_Move(double xpos, double ypos)
{
    // no buttons down: nothing to do
    if( !button_left && !button_middle && !button_right )
        return;

    // compute mouse displacement, save
    double dx = xpos - cursor_lastx;
    double dy = ypos - cursor_lasty;
    cursor_lastx = xpos;
    cursor_lasty = ypos;

    // get current window size
    int width, height;
    glfwGetWindowSize(m_Window, &width, &height);

    // get shift key state
    bool mod_shift = (glfwGetKey(m_Window, GLFW_KEY_LEFT_SHIFT)==GLFW_PRESS ||
                      glfwGetKey(m_Window, GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS);

    // determine action based on mouse button
    mjtMouse action;
    if( button_right )
        action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
    else if( button_left )
        action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
    else
        action = mjMOUSE_ZOOM;

    mjv_moveCamera(mj_Model, action, dx/height, dy/height, &mj_Scn, &mj_Cam);
}

void Visualizer::Scroll(double xoffset, double yoffset)
{
    // scroll: emulate vertical mouse motion = 5% of window height
    mjv_moveCamera(mj_Model, mjMOUSE_ZOOM, 0, -0.05*yoffset, &mj_Scn, &mj_Cam);
    printf("scroll callback %f, %f\n", xoffset, yoffset);
}

// keyboard
void Visualizer::Keyboard(int key, int scancode, int act, int mods)
{
    // do not act on release
    if( act==GLFW_RELEASE )
        return;

    bOKtoComplete = true;
}

void Visualizer::Close() {
    // Free mujoco objects
    mjv_freeScene(&mj_Scn);
    mjr_freeContext(&mj_Con);

    // Close window
    glfwDestroyWindow(m_Window);
    m_Window = NULL;
}
