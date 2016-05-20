//////////////////////////////////////////////////////////////////////////////////
//
// CAMERA.H
//
// Written by Shawn Simonson for Quadrion Engine 11/2005
//
// Provides facility for basic free-roaming camera object.
// Cameras support either Orthographic or Perspective projection. You will
// also find facility for frustum plane extraction herein.
//
//////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "qmath.h"



/*
CCameraObject

Camera object which maintains the view state.
*/

class CCameraObject
{
	public:

		CCameraObject();
		~CCameraObject();


		/*
		SetCamera

		Used for initial setting of the camera state.
		Sets camera position, look at position, and up vector

		"camPos" - position of camera in world space
		"lookAt" - position of the look at point in world space
		"up" - worldspace up vector (relative)
		*/
		void        SetCamera(vec3f camPos, vec3f lookAt, vec3f up);

		/*
		SetPerspective

		Sets camera parameters for generating perspective projection

		"fovy" - field of view (vertical) in degrees
		"aspect" - aspect ratio of the frustum
		"nearClip" - near clip plane distance from camera
		"farClip" - far clip plane distance from camera
		*/
		void        SetPerspective(float fovy, float aspect, float nearClip, float farClip);


		void        Update(int mouseX, int mouseY, int midScreenX, int midScreenY, float sensitivity);
		void		GetViewMatrix(mat4 out);
		void		GetProjectionMatrix(mat4 out);
//		void		Render();

		/*
		MoveX

		The MoveForward/MoveBack functions move forwards and backwards along the view vector
		The MoveRight/MoveLeft functions move left and right along the strafe vector

		"fMul" - multiplier for sensitivity (default is 0.1F)
		*/
		void        MoveForward(float fMul = 0.1F);
		void        MoveBack(float fMul = 0.1F);
		void        MoveRight(float fMul = 0.1F);
		void        MoveLeft(float fMul = 0.1F);


		/*
		GetFarClip/GetNearClip

		Get far clip distance/Get near clip distance
		*/
		const inline float      GetFarClip() const { return m_zFar; }
		const inline float      GetNearClip() const { return m_zNear; }

		vec4f GetVecFromScreenspace(uint32_t x, uint32_t y, uint32_t winWidth, uint32_t winHeight);

		vec3f GetAdjustedUpVec() { return m_adjustedUpVec; }

		/*
		GetFOV

		Get current cam field of view (degrees)
		*/
		const inline float      GetFOV() const { return m_curFOV; }

		const inline vec3f		GetPosition() const { return m_curPosition; }
		void					SetPosition(vec3f pos);

		const vec3f				GetDirection()
		{
			return m_curLookAt - m_curPosition;
		}

		const vec3f				GetLookAtDirection()
		{
			return m_curLookAt;
		}


	private:

		vec3f        m_curPosition;      // Current camera worldspace position
		vec3f        m_curLookAt;        // Current camera look at position
		vec3f        m_curUpVec;         // Current worldspace up vector (relative)

		vec3f		m_adjustedUpVec;	

		float       m_curFOV;           // Current camera field of view (degrees)
		float       m_curAspect;        // Current camera frustum aspect ratio
		float       m_zNear;            // Near clip distance
		float       m_zFar;             // Far clip distance
		float		m_theta;
		float		m_phi;

		float PI2 = PI2_F - 0.001f;
};

