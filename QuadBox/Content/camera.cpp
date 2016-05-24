#include "pch.h"

#include <sstream>
#include <algorithm>
#include "camera.h"

CCameraObject::CCameraObject()
{
	m_curFOV = 0.0F;
	m_curAspect = 0.0F;
	m_zNear = 0.0F;
	m_zFar = 0.0F;

	m_theta = 0;
	m_phi = 0;
}

CCameraObject::~CCameraObject()
{
	m_curFOV = 0.0F;
	m_curAspect = 0.0F;
	m_zNear = 0.0F;
	m_zFar = 0.0F;
}


void CCameraObject::SetCamera(vec3f camPos, vec3f lookAt, vec3f up)
{
	m_curPosition = camPos;
	m_curLookAt = lookAt;
	m_curUpVec = up;

	vec3f direction = m_curLookAt - m_curPosition;
	direction.normalize();

	m_phi = asin(direction.y);

	auto d = max(abs(direction.z), abs(direction.x));

	if (d == abs(direction.z))
	{
		m_theta = acos(direction.z / cos(m_phi));
	}
	else
	{
		m_theta = asin(direction.x / cos(m_phi));
	}
}

/*
rotate
Internally used helper for CCameraObject::RotateByMouse
*/
static void rotate(float* camPos, float* lookAtPos, float* upVec, float xDir)
{
	quat qRot, qView, qNewView;
	quat_rotate(qRot, deg_2_rad(xDir), upVec);

	qView[0] = lookAtPos[0] - camPos[0];
	qView[1] = lookAtPos[1] - camPos[1];
	qView[2] = lookAtPos[2] - camPos[2];
	qView[3] = 0.0f;

	quat_mul(qRot, qView, qNewView);
	quat_conjugate(qRot);
	quat_mul(qNewView, qRot, qNewView);

	lookAtPos[0] = camPos[0] + qNewView[0];
	lookAtPos[1] = camPos[1] + qNewView[1];
	lookAtPos[2] = camPos[2] + qNewView[2];
}


void CCameraObject::Update(int mouseX, int mouseY, int midScreenX, int midScreenY, float sensitivity)
{
	if (mouseX == 0 && mouseY == 0)
		return;

	//	Sensitivity should be something...
	if (sensitivity == 0)
		sensitivity = 0.001f;

	m_theta -= mouseX * sensitivity;
	m_phi += mouseY * sensitivity;

	if (m_phi > PI2)
		m_phi = PI2;

	if (m_phi < -PI2)
		m_phi = -PI2;

	float x = sin(m_theta) * cos(m_phi);
	float y = sin(m_phi);
	float z = cos(m_theta) * cos(m_phi);

	vec3f direction(x, y, z);
	direction.normalize();

	m_curLookAt = m_curPosition + direction;
}

/*
void CCameraObject::Render()
{
	vec3 pos;
	vec3 lookAt;
	vec3 upVec;

	pos[0] = m_curPosition.x;
	pos[1] = m_curPosition.y;
	pos[2] = m_curPosition.z;

	lookAt[0] = m_curLookAt.x;
	lookAt[1] = m_curLookAt.y;
	lookAt[2] = m_curLookAt.z;

	upVec[0] = m_curUpVec.x;
	upVec[1] = m_curUpVec.y;
	upVec[2] = m_curUpVec.z;

	mat4 la;
	mat4_load_lookat_rh(pos, lookAt, upVec, la);
	m_adjustedUpVec.set(la[1], la[5], la[9]);

	setMatrix(MATRIX_VIEW, la);


	mat4 p;
	mat4_load_perspective_rh(m_curFOV, m_curAspect, m_zNear, m_zFar, p);

	setMatrix(MATRIX_PROJECTION, p);
}
*/


void CCameraObject::GetViewMatrix(mat4 out)
{
	vec3 pos;
	vec3 lookAt;
	vec3 upVec;

	pos[0] = m_curPosition.x;
	pos[1] = m_curPosition.y;
	pos[2] = m_curPosition.z;

	lookAt[0] = m_curLookAt.x;
	lookAt[1] = m_curLookAt.y;
	lookAt[2] = m_curLookAt.z;

	upVec[0] = m_curUpVec.x;
	upVec[1] = m_curUpVec.y;
	upVec[2] = m_curUpVec.z;

	mat4_load_lookat_rh(pos, lookAt, upVec, out);
}

void CCameraObject::GetProjectionMatrix(mat4 out)
{
	mat4_load_perspective_rh(m_curFOV, m_curAspect, m_zNear, m_zFar, out);
}

void CCameraObject::SetPerspective(float fovy, float aspect, float nearClip, float farClip)
{
	m_curFOV = fovy;
	m_curAspect = aspect;
	m_zNear = nearClip;
	m_zFar = farClip;
}

void CCameraObject::MoveForward(float fMul)
{
	vec3f direction = m_curLookAt - m_curPosition;
	direction.normalize();
	direction *= fMul;
	m_curPosition += direction;
	m_curLookAt += direction;
}

void CCameraObject::MoveBack(float fMul)
{
	vec3f direction = m_curLookAt - m_curPosition;
	direction.normalize();
	direction *= fMul;
	m_curPosition -= direction;
	m_curLookAt -= direction;
}

void CCameraObject::MoveRight(float fMul)
{
	vec3f up = m_curUpVec;
	vec3f direction = m_curLookAt - m_curPosition;
	direction.normalize();

	vec3f right = direction.crossProd(up);
	right.normalize();
	right *= fMul;

	m_curPosition += right;
	m_curLookAt += right;
}

void CCameraObject::MoveLeft(float fMul)
{
	vec3f up = m_curUpVec;
	vec3f direction = m_curLookAt - m_curPosition;
	direction.normalize();

	vec3f left = up.crossProd(direction);
	left.normalize();
	left *= fMul;

	m_curPosition += left;
	m_curLookAt += left;
}

void CCameraObject::SetPosition(vec3f pos)
{
	vec3f direction = m_curLookAt - m_curPosition;
	direction.normalize();
	m_curPosition = pos;
	m_curLookAt = m_curPosition + direction;
}

vec4f CCameraObject::GetVecFromScreenspace(uint32_t x, uint32_t y, uint32_t winWidth, uint32_t winHeight)
{
/*
	float percX = x / (float)winWidth;
	float percY = y / (float)winHeight;

	vec3f farPCenter(0.0f, 0.0f, -m_zFar);
	vec3f vpUp(0.0f, 1.0f, 0.0f);
	float farPHeight = 2.0f * tanf(m_curFOV / 2.0f) * m_zFar;
	float farPWidth = farPHeight * m_curAspect;

	float interpX = (-farPWidth * 0.5f) + (farPWidth * percX);
	float interpY = (farPHeight * 0.5f) - (farPHeight * percY);

	vec3f z(interpX, interpY, -m_zFar);
	z.normalize();
	vec4f z4(z.x, z.y, z.z, 0.0f);
	mat4f V = getMatrix(MATRIX_VIEW);
	V.invert();
	vec4f res = V * z;

	return res;
*/
	return vec4f(1.0, 1.0, 1.0, 1.0);
}