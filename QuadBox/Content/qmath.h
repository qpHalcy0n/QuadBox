#ifndef __qikmath_h_
#define __qikmath_h_

#include <cstdint>
#include <math.h>
#include <string.h>
#include <stdlib.h>
//#include "qengine_export.h"

//#define max(x, y) (((x) > (y)) ? (x) : (y))
//#define min(x, y) (((x) < (y)) ? (x) : (y))


typedef float vec2[2];
typedef float vec3[3];
typedef float vec4[4];

typedef float mat3[12];
typedef float mat4[16];

typedef float quat[4];

const double PI		= 3.141592653589793238463;
const float  PI_F	= 3.14159265358979f;
const float  PI2_F	= 1.570796326794895f;





static float deg_2_rad(const float deg)
{
	return PI_F * deg / 180.0F;
}

static float rad_2_deg(const float rad)
{
	return 180.0F / PI_F * rad;
}

static void swap(float* l, float* r)
{
	float tmp = *l;
	*l = *r;
	*r = tmp;
}

static float FastISqrt(const float& a)
{
	float ha = 0.5F * a;
	float ret;
	int i = *(int*)&a;
	i = 0x5F3759D5 - (i >> 1);
	ret = *(float*)&i;
	ret = ret * (1.5F - ha * ret * ret);
	return ret;
}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
//-=-=-=-=-=-=-=-= VECTOR3 ROUTINES -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//


static void vec3_clear(vec3 l)
{
	l[0] = 0.0F;
	l[1] = 0.0F;
	l[2] = 0.0F;
}


static void vec3_set(vec3 l, float x, float y, float z)
{
	l[0] = x;
	l[1] = y;
	l[2] = z;
}

static void vec4_set(vec4 l, float x, float y, float z, float w)
{
	l[0] = x;
	l[1] = y;
	l[2] = z;
	l[3] = w;
}


// 3 component piecewise sub //
static void vec3_sub(const vec3 l, const vec3 r, vec3 ret)
{
	ret[0] = l[0] - r[0];
	ret[1] = l[1] - r[1];
	ret[2] = l[2] - r[2];
}

// 3 component piecewise sub and copy to lhs //
static void vec3_sub_replace(vec3 l, const vec3 r)
{
	l[0] -= r[0];
	l[1] -= r[1];
	l[2] -= r[2];
}

// 3 component piecewise add //
static void vec3_add(const vec3 l, const vec3 r, vec3 ret)
{
	ret[0] = l[0] + r[0];
	ret[1] = l[1] + r[1];
	ret[2] = l[2] + r[2];
}

// 3 component piecewise add and copy to lhs //
static void vec3_add_replace(vec3 l, const vec3 r)
{
	l[0] += r[0];
	l[1] += r[1];
	l[2] += r[2];
}

// Vector dot product //
static float vec3_dot(const vec3 l, const vec3 r)
{
	return l[0] * r[0] + l[1] * r[1] + l[2] * r[2];
}

// Vector cross product //
static void vec3_cross(const vec3 l, const vec3 r, vec3 ret)
{
	ret[0] = l[1] * r[2] - l[2] * r[1];
	ret[1] = l[2] * r[0] - l[0] * r[2];
	ret[2] = l[0] * r[1] - l[1] * r[0];
}

// Vector length as a scalar //
static float vec3_length(const vec3 v)
{
	return sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// Vector normalization (length unrecoverable) //
static void vec3_normalize(vec3 v)
{
	float iLen = 1.0F / vec3_length(v);
	v[0] *= iLen;
	v[1] *= iLen;
	v[2] *= iLen;
}


// Vector scale (scalar) //
static void vec3_scale(vec3 v, const float s)
{
	v[0] *= s;
	v[1] *= s;
	v[2] *= s;
}


// Vector copy (piecewise) //
static void vec3_copy(vec3 l, const vec3 r)
{
	l[0] = r[0];
	l[1] = r[1];
	l[2] = r[2];
}

static void vec4_copy(vec4 l, const vec4 r)
{
	l[0] = r[0];
	l[1] = r[1];
	l[2] = r[2];
	l[3] = r[3];
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= MATRIX4 ROUTINES -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//


static void mat4_load_identity(mat4 m)
{
	m[0] = 1.0F;	m[4] = 0.0F;	m[8] = 0.0F;		m[12] = 0.0F;
	m[1] = 0.0F;	m[5] = 1.0F;	m[9] = 0.0F;		m[13] = 0.0F;
	m[2] = 0.0F;	m[6] = 0.0F;	m[10] = 1.0F;		m[14] = 0.0F;
	m[3] = 0.0F;    m[7] = 0.0F;	m[11] = 0.0F;       m[15] = 1.0F;
}

static void mat4_load_scale(mat4 m, float sx, float sy, float sz)
{
	m[0] = sx;		m[4] = 0.0F;	m[8] = 0.0F;		m[12] = 0.0F;
	m[1] = 0.0F;	m[5] = sy;		m[9] = 0.0F;		m[13] = 0.0F;
	m[2] = 0.0F;	m[6] = 0.0F;	m[10] = sz;			m[14] = 0.0F;
	m[3] = 0.0F;    m[7] = 0.0F;	m[11] = 0.0F;       m[15] = 1.0F;
}

static void mat4_load_translation(mat4 m, float tx, float ty, float tz)
{
	m[0] = 1.0F;	m[4] = 0.0F;	m[8] = 0.0F;		m[12] = tx;
	m[1] = 0.0F;	m[5] = 1.0F;	m[9] = 0.0F;		m[13] = ty;
	m[2] = 0.0F;	m[6] = 0.0F;	m[10] = 1.0F;		m[14] = tz;
	m[3] = 0.0F;    m[7] = 0.0F;	m[11] = 0.0F;       m[15] = 1.0F;
}

static void mat4_load_x_rot(mat4 m, float ang)
{
	m[0] = 1.0f;	m[4] = 0.0f;		m[8] = 0.0F;		m[12] = 0.0F;
	m[1] = 0.0f;	m[5] = cosf(ang);	m[9] = -sinf(ang);	m[13] = 0.0F;
	m[2] = 0.0F;	m[6] = sinf(ang);	m[10] = cosf(ang);	m[14] = 0.0F;
	m[3] = 0.0F;	m[7] = 0.0F;		m[11] = 0.0F;		m[15] = 1.0F;
}

static void mat4_load_y_rot(mat4 m, float ang)
{
	m[0] = cosf(ang);	m[4] = 0.0f;	m[8] = sinf(ang);	m[12] = 0.0F;
	m[1] = 0.0f;		m[5] = 1.0f;	m[9] = 0.0F;		m[13] = 0.0F;
	m[2] = -sinf(ang);	m[6] = 0.0F;	m[10] = cosf(ang);	m[14] = 0.0F;
	m[3] = 0.0F;	    m[7] = 0.0F;	m[11] = 0.0F;		m[15] = 1.0F;	
}

static void mat4_load_z_rot(mat4 m, float ang)
{
	m[0] = cosf(ang);	m[4] = -sinf(ang);	m[8] = 0.0F;	m[12] = 0.0F;
	m[1] = sinf(ang);	m[5] = cosf(ang);	m[9] = 0.0F;	m[13] = 0.0F;
	m[2] = 0.0F;		m[6] = 0.0F;		m[10] = 1.0F;	m[14] = 0.0F;
	m[3] = 0.0F;	    m[7] = 0.0F;		m[11] = 0.0F;	m[15] = 1.0F;
}

static void mat4_transpose(mat4 m)
{
	swap(&m[1], &m[4]);
	swap(&m[2], &m[8]);
	swap(&m[3], &m[12]);
	swap(&m[6], &m[9]);
	swap(&m[7], &m[13]);
	swap(&m[11], &m[14]);
}

// Brute force 4x4 matrix mul //
static void mat4_mul_mat4(const mat4 l, const mat4 r, mat4 output)
{
	output[0] = (l[0] * r[0]) + (l[4] * r[1]) + (l[8] * r[2]) + (l[12] * r[3]);
	output[4] = (l[0] * r[4]) + (l[4] * r[5]) + (l[8] * r[6]) + (l[12] * r[7]);
	output[8] = (l[0] * r[8]) + (l[4] * r[9]) + (l[8] * r[10]) + (l[12] * r[11]);
	output[12] = (l[0] * r[12]) + (l[4] * r[13]) + (l[8] * r[14]) + (l[12] * r[15]);

	output[1] = (l[1] * r[0]) + (l[5] * r[1]) + (l[9] * r[2]) + (l[13] * r[3]);
	output[5] = (l[1] * r[4]) + (l[5] * r[5]) + (l[9] * r[6]) + (l[13] * r[7]);
	output[9] = (l[1] * r[8]) + (l[5] * r[9]) + (l[9] * r[10]) + (l[13] * r[11]);
	output[13] = (l[1] * r[12]) + (l[5] * r[13]) + (l[9] * r[14]) + (l[13] * r[15]);

	output[2] = (l[2] * r[0]) + (l[6] * r[1]) + (l[10] * r[2]) + (l[14] * r[3]);
	output[6] = (l[2] * r[4]) + (l[6] * r[5]) + (l[10] * r[6]) + (l[14] * r[7]);
	output[10] = (l[2] * r[8]) + (l[6] * r[9]) + (l[10] * r[10]) + (l[14] * r[11]);
	output[14] = (l[2] * r[12]) + (l[6] * r[13]) + (l[10] * r[14]) + (l[14] * r[15]);

	output[3] = (l[3] * r[0]) + (l[7] * r[1]) + (l[11] * r[2]) + (l[15] * r[3]);
	output[7] = (l[3] * r[4]) + (l[7] * r[5]) + (l[11] * r[6]) + (l[15] * r[7]);
	output[11] = (l[3] * r[8]) + (l[7] * r[9]) + (l[11] * r[10]) + (l[15] * r[11]);
	output[15] = (l[3] * r[12]) + (l[7] * r[13]) + (l[11] * r[14]) + (l[15] * r[15]);
}

// Tested Correct - SMS 05/24/15
static float mat4_det(const mat4 m)
{
	float d = (m[0] * m[5] * m[10] * m[15]) + (m[0] * m[9] * m[14] * m[7]) + (m[0] * m[13] * m[6] * m[11])
		      +(m[4] * m[1] * m[14] * m[11]) + (m[4] * m[9] * m[2] * m[15]) + (m[4] * m[13] * m[10] * m[3])
		      +(m[8] * m[1] * m[6] * m[15]) + (m[8] * m[5] * m[14] * m[3]) + (m[8] * m[13] * m[2] * m[7])
			  +(m[12] * m[1] * m[10] * m[7]) + (m[12] * m[5] * m[2] * m[11]) + (m[12] * m[9] * m[6] * m[3])
			  -(m[0] * m[5] * m[14] * m[11]) - (m[0] * m[9] * m[6] * m[15]) - (m[0] * m[13] * m[10] * m[7])
			  -(m[4] * m[1] * m[10] * m[15]) - (m[4] * m[9] * m[14] * m[3]) - (m[4] * m[13] * m[2] * m[11])
			  -(m[8] * m[1] * m[14] * m[7]) - (m[8] * m[5] * m[2] * m[15]) - (m[8] * m[13] * m[6] * m[3])
			  -(m[12] * m[1] * m[6] * m[11]) - (m[12] * m[5] * m[10] * m[3]) - (m[12] * m[9] * m[2] * m[7]);

	return d;
}


// Tested correct - SMS 05/24/15
static void mat4_invert(const mat4 m, mat4 out)
{
	float d = mat4_det(m);
	// Check for 0 determinant //
	if(fabs(d - 0.0f) < 0.00001f || d == 0.0f)
		return;

	float i = 1.0f / d;
	out[0] = i * ((m[5]*m[10]*m[15])+(m[9]*m[14]*m[7])+(m[13]*m[6]*m[11])-(m[5]*m[14]*m[11])-(m[9]*m[6]*m[15])-(m[13]*m[10]*m[7]));
	out[4] = i * ((m[4]*m[14]*m[11])+(m[8]*m[6]*m[15])+(m[12]*m[10]*m[7])-(m[4]*m[10]*m[15])-(m[8]*m[14]*m[7])-(m[12]*m[6]*m[11]));
	out[8] = i * ((m[4]*m[9]*m[15])+(m[8]*m[13]*m[7])+(m[12]*m[5]*m[11])-(m[4]*m[13]*m[11])-(m[8]*m[5]*m[15])-(m[12]*m[9]*m[7]));
	out[12] = i * ((m[4]*m[13]*m[10])+(m[8]*m[5]*m[14])+(m[12]*m[9]*m[6])-(m[4]*m[9]*m[14])-(m[8]*m[13]*m[6])-(m[12]*m[5]*m[10]));
	out[1] = i * ((m[1]*m[14]*m[11])+(m[9]*m[2]*m[15])+(m[13]*m[10]*m[3])-(m[1]*m[10]*m[15])-(m[9]*m[14]*m[3])-(m[13]*m[2]*m[11]));
	out[5] = i * ((m[0]*m[10]*m[15])+(m[8]*m[14]*m[3])+(m[12]*m[2]*m[11])-(m[0]*m[14]*m[11])-(m[8]*m[2]*m[15])-(m[12]*m[10]*m[3]));
	out[9] = i * ((m[0]*m[13]*m[11])+(m[8]*m[1]*m[15])+(m[12]*m[9]*m[3])-(m[0]*m[9]*m[15])-(m[8]*m[13]*m[3])-(m[12]*m[1]*m[11]));
	out[13] = i * ((m[0]*m[9]*m[14])+(m[8]*m[13]*m[2])+(m[12]*m[1]*m[10])-(m[0]*m[13]*m[10])-(m[8]*m[1]*m[14])-(m[3]*m[9]*m[2]));
	out[2] = i * ((m[1]*m[6]*m[15])+(m[5]*m[14]*m[3])+(m[13]*m[2]*m[7])-(m[1]*m[14]*m[7])-(m[5]*m[2]*m[15])-(m[13]*m[6]*m[3]));
	out[6] = i * ((m[0]*m[14]*m[7])+(m[4]*m[2]*m[15])+(m[12]*m[6]*m[3])-(m[0]*m[6]*m[15])-(m[4]*m[14]*m[3])-(m[12]*m[2]*m[7]));
	out[10] = i * ((m[0]*m[5]*m[15])+(m[4]*m[13]*m[3])+(m[12]*m[1]*m[7])-(m[0]*m[13]*m[7])-(m[4]*m[1]*m[15])-(m[12]*m[5]*m[3]));
	out[14] = i * ((m[0]*m[13]*m[6])+(m[4]*m[1]*m[14])+(m[12]*m[5]*m[2])-(m[0]*m[5]*m[14])-(m[4]*m[13]*m[2])-(m[12]*m[1]*m[6]));
	out[3] = i * ((m[1]*m[10]*m[7])+(m[5]*m[2]*m[11])+(m[9]*m[6]*m[3])-(m[1]*m[6]*m[11])-(m[5]*m[10]*m[3])-(m[9]*m[2]*m[7]));
	out[7] = i * ((m[0]*m[6]*m[11])+(m[4]*m[10]*m[3])+(m[8]*m[2]*m[7])-(m[0]*m[10]*m[7])-(m[4]*m[2]*m[11])-(m[8]*m[6]*m[3]));
	out[11] = i * ((m[0]*m[9]*m[7])+(m[4]*m[1]*m[11])+(m[8]*m[5]*m[3])-(m[0]*m[5]*m[11])-(m[4]*m[9]*m[3])-(m[8]*m[1]*m[7]));
	out[15] = i * ((m[0]*m[5]*m[10])+(m[4]*m[9]*m[2])+(m[8]*m[1]*m[6])-(m[0]*m[9]*m[6])-(m[4]*m[1]*m[10])-(m[8]*m[5]*m[2]));
}

static void mat4_mul_vec4SSE(const mat4 l, const vec4 r, vec4 output)
{
	float* r0 = (float*)&l[0];

#ifdef WIN32
	//__asm
	//{
	//	mov			esi, r					// Move address of r into esi
	//		mov 		edi, output				// Move address of out into edi

	//		mov			edx, r0
	//		movups		xmm4, [edx]				// Copy first col into xmm4    (16b)  from edx
	//		movups		xmm5, [edx + 0x10]		// Copy second col into xmm5   (16b)  from edx
	//		movups		xmm6, [edx + 0x20]		// Copy third col into xmm6    (16b)  from edx
	//		movups		xmm7, [edx + 0x30]		// Copy fourth col into xmm7   (16b)  from edx

	//		movups		xmm0, [esi]			// load vec into xmm0 from esi
	//		xorps		xmm2, xmm2			// Clear xmm2 to 0

	//		movups		xmm1, xmm0			// load vector into xmm1 from xmm0
	//		shufps		xmm1, xmm1, 0x00	// No shuf for x
	//		mulps		xmm1, xmm4			// Mul col0 by vec
	//		addps		xmm2, xmm1			// Add result to xmm2

	//		// Same for 2nd column //
	//		movups		xmm1, xmm0
	//		shufps		xmm1, xmm1, 0x55
	//		mulps		xmm1, xmm5
	//		addps		xmm2, xmm1

	//		// 3rd column //
	//		movups		xmm1, xmm0
	//		shufps		xmm1, xmm1, 0xAA
	//		mulps		xmm1, xmm6
	//		addps		xmm2, xmm1

	//		// 4th column //
	//		movups		xmm1, xmm0
	//		shufps		xmm1, xmm1, 0xFF
	//		mulps		xmm1, xmm7
	//		addps		xmm2, xmm1

	//		// Broadcast result to edi (out) //
	//		movups[edi], xmm2
	//}

	// ATT Syntax //
#else
	/*
	asm(".intel_syntax noprefix");
	asm("mov    esi, r");
	asm("mov    edi, out");
	asm("mov    edx, r0");
	asm("movups xmm4, [edx]");
	asm("movups xmm5, [edx + 0x10]");
	asm("movups xmm6, [edx + 0x20]");
	asm("movups xmm7, [edx + 0x30]");
	asm("movups xmm0, [esi]");
	asm("xorps  xmm2, xmm2");
	asm("movups xmm1, xmm0");
	asm("shufps xmm1, xmm1, 0x00");
	asm("mulps  xmm1, xmm4");
	asm("addps  xmm2, xmm1");
	asm("movups xmm1, xmm0");
	asm("shufps xmm1, xmm1, 0x55");
	asm("mulps  xmm1, xmm5");
	asm("addps  xmm2, xmm1");
	asm("movups xmm1, xmm0");
	asm("shufps xmm1, xmm1, 0xAA");
	asm("mulps  xmm1, xmm6");
	asm("addps  xmm2, xmm1");
	asm("movups xmm1, xmm0");
	asm("shufps xmm1, xmm1, 0xFF");
	asm("mulps  xmm1, xmm7");
	asm("addps  xmm2, xmm1");
	asm("movups [edi], xmm2");
	asm(".att_syntax noprefix");
	*/
#endif

}


// Brute force 4x4 matrix to vector mul //
static void mat4_mul_vec4(const mat4 l, const vec4 r, vec4 output)
{
	//	if(CPUSSE2)
	//	{
	//		mat4_mul_vec4SSE(l, r, out);
	//		return;
	//	}

	output[0] = l[0] * r[0] + l[4] * r[1] + l[8] * r[2] + l[12] * r[3];
	output[1] = l[1] * r[0] + l[5] * r[1] + l[9] * r[2] + l[13] * r[3];
	output[2] = l[2] * r[0] + l[6] * r[1] + l[10] * r[2] + l[14] * r[3];
	output[3] = l[3] * r[0] + l[7] * r[1] + l[11] * r[2] + l[15] * r[3];
}



static void mat4_load_lookat_rh(vec3 pos, vec3 lookAt, vec3 up, mat4 output)
{
	// Find normalized view vector //
	vec3 z, x, y;
	vec3_sub(lookAt, pos, z);
	vec3_normalize(z);

	// Construct x axis from Z X UP //
	vec3_cross(z, up, x);
	vec3_normalize(x);

	// Cam y then becomes already normalized X X Z //
	vec3_cross(x, z, y);

	z[0] = -z[0];
	z[1] = -z[1];
	z[2] = -z[2];
	vec3 ne = { -pos[0], -pos[1], -pos[2] };

	output[0] = x[0];        output[4] = x[1];    output[8] = x[2];        output[12] = vec3_dot(ne, x);
	output[1] = y[0];        output[5] = y[1];    output[9] = y[2];        output[13] = vec3_dot(ne, y);
	output[2] = z[0];        output[6] = z[1];    output[10] = z[2];       output[14] = vec3_dot(ne, z);
	output[3] = 0.0F;        output[7] = 0.0F;    output[11] = 0.0F;       output[15] = 1.0F;
}



static void mat4_load_perspective_rh(float fovy, float aspect, float nearClip, float farClip, mat4 output)
{

	float iTan = 1.0F / tanf(fovy * 0.5F);
	float dNF = (nearClip - farClip);

	output[0] = iTan / aspect;     output[4] = 0.0F;      output[8] = 0.0F;
	output[1] = 0.0F;              output[5] = iTan;      output[9] = 0.0F;
	output[2] = 0.0F;              output[6] = 0.0F;      output[10] = (farClip + nearClip) / dNF;
	output[3] = 0.0F;              output[7] = 0.0F;      output[11] = -1.0F;

	output[12] = 0.0F;
	output[13] = 0.0F;
	output[14] = (2.0F* farClip * nearClip) / dNF;
	output[15] = 0.0F;

}


static void mat4_copy(mat4 to, mat4 from)
{
	memcpy(to, from, sizeof(float) * 16);
}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= QUATERNION ROUTINES -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//



// Quaternion rotation //
static void quat_rotate(quat q, const float rad, vec3 axis)
{
	if ((axis[0] != 0 && axis[0] != 1) || (axis[1] != 0 && axis[1] != 1) || (axis[2] != 0 && axis[2] != 1))
		vec3_normalize(axis);

	float s = (float)sinf(rad * 0.5F);

	vec3 tm = { axis[0] * s, axis[1] * s, axis[2] * s };
	float w = (float)cosf(rad * 0.5F);
	float len = 1.0F / (float)sqrtf(tm[0] * tm[0] + tm[1] * tm[1] + tm[2] * tm[2] + w * w);
	tm[0] *= len;
	tm[1] *= len;
	tm[2] *= len;

	memcpy(q, tm, sizeof(float) * 3);
	q[3] = w;

}

// Quaternion multiplication //
static void quat_mul(quat l, quat r, quat output)
{
	output[0] = l[3] * r[0] + l[0] * r[3] + l[1] * r[2] - l[2] * r[1];
	output[1] = l[3] * r[1] - l[0] * r[2] + l[1] * r[3] + l[2] * r[0];
	output[2] = l[3] * r[2] + l[0] * r[1] - l[1] * r[0] + l[2] * r[3];
	output[3] = l[3] * r[3] - l[0] * r[0] - l[1] * r[1] - l[2] * r[2];
}

// Quaternion conjugation //
static void quat_conjugate(quat q)
{
	q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; q[3] = q[3];
}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
//				VECTOR CLASS DEFINITIONS					  //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

struct vec2f
{
	float x, y;

	vec2f(){}
	~vec2f(){}

	vec2f(const float* v) { x = v[0]; y = v[1]; }
	vec2f(const vec2f& v) { x = v.x; y = v.y; }
	vec2f(const float& nx, const float& ny) { x = nx; y = ny; }

	const inline void set(const float* v) { x = v[0]; y = v[1]; }
	const inline void set(const vec2f& v) { x = v.x; y = v.y; }
	const inline void set(const float& nx, const float& ny) { x = nx; y = ny; }

	const inline void operator*= (const float& s) { x *= s; y *= s; }
	const inline void operator/= (const float& s) { x /= s; y /= s; }
	const inline void operator+= (const vec2f& v) { x += v.x; y += v.y; }
	const inline void operator+= (const float* v) { x += v[0]; y += v[1]; }
	const inline void operator-= (const vec2f& v) { x -= v.x; y -= v.y; }
	const inline void operator-= (const float* v) { x -= v[0]; y -= v[1]; }
	const inline void operator=  (const vec2f& v) { x = v.x; y = v.y; }
	const inline void operator=  (const float* v) { x = v[0]; y = v[1]; }
	const inline bool operator== (const vec2f& v) { return ((x == v.x) && (y == v.y)); }
	const inline bool operator== (const float* v) { return ((x == v[0]) && (y == v[1])); }
	const inline bool operator!= (const vec2f& v) { return !(*this == v); }
	const inline bool operator!= (const float* v) { return !(*this == v); }


	float getLength()
	{
		return 1.0f / FastISqrt(x * x + y * y);
	}

	const inline float dotProd(const vec2f& v)
	{
		return (x * v.x) + (y * v.y);
	}

	const inline vec2f getPerpVector()
	{
		return vec2f(-y, x);
	}

	const inline void normalize()
	{
		float iLen = 1.0F / getLength();
		x *= iLen;
		y *= iLen;
	}
};

static vec2f operator+ (const vec2f& l, const vec2f& r)
{
	return vec2f(l.x + r.x, l.y + r.y);
}

static vec2f operator* (const vec2f& l, const float& r)
{
	return vec2f(l.x * r, l.y * r);
}


struct vec2i
{

	vec2i(){}
	~vec2i(){}

	vec2i(const int nx, const int ny) { x = nx; y = ny; }


	int32_t		x, y;
};

struct vec3f
{
	float x, y, z;

	vec3f(){}
	~vec3f(){}

	vec3f(const float* v) { x = v[0]; y = v[1]; z = v[2]; }
	vec3f(const vec3f& v) { x = v.x; y = v.y; z = v.z; }
	vec3f(const float& nx, const float& ny, const float& nz) { x = nx; y = ny; z = nz; }

	const inline void set(const float* v) { x = v[0]; y = v[1]; z = v[2]; }
	const inline void set(const vec3f& v) { x = v.x; y = v.y; z = v.z; }
	const inline void set(const float& nx, const float& ny, const float& nz) { x = nx; y = ny; z = nz; }

	const inline vec3f operator- () { return vec3f(-x, -y, -z); }
	const inline void operator*= (const float& s) { x *= s; y *= s; z *= s; }
	const inline void operator*= (const vec3f& v) { x *= v.x; y *= v.y; z *= v.z; }
	const inline void operator/= (const float& s) { float inv = 1.0F / s; x *= inv; y *= inv; z *= inv; }
	const inline void operator+= (const vec3f& v) { x += v.x; y += v.y; z += v.z; }
	const inline void operator+= (const float* v) { x += v[0]; y += v[1]; z += v[2]; }
	const inline void operator-= (const vec3f& v) { x -= v.x; y -= v.y; z -= v.z; }
	const inline void operator-= (const float* v) { x -= v[0]; y -= v[1]; z -= v[2]; }
	const inline void operator=  (const vec3f& v) { x = v.x; y = v.y; z = v.z; }
	const inline void operator=  (const float* v) { x = v[0]; y = v[1]; z = v[2]; }
	const inline bool operator== (const vec3f& v) { return ((x == v.x) && (y == v.y) && (z == v.z)); }
	const inline bool operator== (const float* v) { return ((x == v[0]) && (y == v[1]) && (z == v[2])); }
	const inline bool operator!= (const vec3f& v) { return !(*this == v); }
	const inline bool operator!= (const float* v) { return !(*this == v); }

	const inline vec3f lerp(const vec3f& op, const float& dv)
	{
		return vec3f(x + dv * (op.x - x), y + dv * (op.y - y), z + dv * (op.z - z));
	}

	float getLength()
	{
		return sqrt(x*x + y*y + z*z);
			//FastISqrt(x * x + y * y + z * z);
	}

	const inline float dotProd(const vec3f& v)
	{
		return ((x * v.x) + (y * v.y) + (z * v.z));
	}

	const inline void normalize()
	{
		float length = getLength();

		x /= length;
		y /= length;
		z /= length;

		//float mag = 1.0F / getLength(); x *= mag; y *= mag; z *= mag;
	}

	const inline vec3f crossProd(const vec3f& v)
	{
		return vec3f((y * v.z - z * v.y), (z * v.x - x * v.z), (x * v.y - y * v.x));
	}

	const inline vec3f crossProd(const float* v)
	{
		return vec3f((y * v[2] - z * v[1]), (z * v[0] - x * v[2]), (x * v[1] - y * v[0]));
	}

	const inline float tripleScalar(const vec3f& v, const vec3f& w)
	{
		return ((y * v.z - z * v.y) * w.x + (z * v.x - x * v.z) * w.y + (x * v.y - y * v.x) * w.z);
	}

	const inline void	assignIfLess(const vec3f& v)
	{
		x = (v.x < x) ? v.x : x;
		y = (v.y < y) ? v.y : y;
		z = (v.z < z) ? v.z : z;
	}

	const inline void	assignIfGreater(const vec3f& v)
	{
		x = (v.x > x) ? v.x : x;
		y = (v.y > y) ? v.y : y;
		z = (v.z > z) ? v.z : z;
	}
};

static vec3f operator- (const vec3f& l, const vec3f& r)
{
	return vec3f(l.x - r.x, l.y - r.y, l.z - r.z);
}

static vec3f operator+ (const vec3f& l, const vec3f& r)
{
	return vec3f(l.x + r.x, l.y + r.y, l.z + r.z);
}

static vec3f operator* (const vec3f& l, const float& r)
{
	return vec3f(l.x * r, l.y * r, l.z * r);
}



struct mat4f
{
	float m[16];

	mat4f() {}
	~mat4f() {}

	mat4f(const float* in) { memcpy(m, in, sizeof(float) * 16); }
	mat4f(const mat4f& in) { memcpy(m, in.m, sizeof(float) * 16); }

	const inline void set(const float* in) { memcpy(m, in, sizeof(float) * 16); }
	const inline void set(const mat4f& in) { memcpy(m, in.m, sizeof(float) * 16); }

	const inline void operator= (const float* in) { memcpy(m, in, sizeof(float) * 16); }
	const inline void operator= (const mat4f& in) { memcpy(m, in.m, sizeof(float) * 16); }

	const inline void operator*= (const mat4f& r) { mat4 res; mat4_mul_mat4(m, r.m, res); mat4_copy(m, res); }

	float getDeterminant() { return mat4_det(m); }
	mat4f getInverse() { mat4 inv; mat4_invert(m, inv); mat4f ret(inv); return ret; }

	float& operator[] (int i)
	{
		if (i > 15)
			return m[0];

		return m[i];
	}

	void transpose() { mat4_transpose(m); }
	void invert() { mat4 inv; mat4_invert(m, inv); mat4_copy(m, inv); }
	void loadIdentity() { mat4_load_identity(m); }
	void loadScale(float sx, float sy, float sz) { mat4_load_scale(m, sx, sy, sz); }
	void loadTranslation(float tx, float ty, float tz) { mat4_load_translation(m, tx, ty, tz); }

	void loadXRotation(float ang) { mat4_load_x_rot(m, ang); }
	void loadYRotation(float ang) { mat4_load_y_rot(m, ang); }
	void loadZRotation(float ang) { mat4_load_z_rot(m, ang); }
};


struct vec4f
{
	float x, y, z, w;

	vec4f() {}
	~vec4f(){}

	vec4f(const float* v) { x = v[0]; y = v[1]; z = v[2]; w = v[3]; }
	vec4f(const vec4f& v) { x = v.x; y = v.y; z = v.z; w = v.w; }
	vec4f(const vec3f& v) { x = v.x; y = v.y; z = v.z; w = 0.0f; }
	vec4f(const float& nx, const float& ny, const float& nz, const float& nw) { x = nx; y = ny; z = nz; w = nw; }

	const inline void set(const float* v) { x = v[0]; y = v[1]; z = v[2]; w = v[3]; }
	const inline void set(const vec4f& v) { x = v.x; y = v.y; z = v.z; w = v.w; }
	const inline void set(const float& nx, const float& ny, const float& nz, const float& nw) { x = nx; y = ny; z = nz; w = nw; }

	const inline void operator=  (const vec4f& v) { x = v.x; y = v.y; z = v.z; w = v.w; }
	const inline void operator=  (const vec3f& v) { x = v.x; y = v.y; z = v.z; w = 0.0f; }
	const inline void operator=  (const float* v) { x = v[0]; y = v[1]; z = v[2]; w = v[3]; }


	const inline void operator*= (const float& v) { x *= v; y *= v; z *= v; w *= w; }
	const inline void operator*= (const mat4f& m)
	{
		float operand[4] = { x, y, z, w };
		float out[4];
		mat4_mul_vec4SSE(m.m, operand, out);
		x = out[0];
		y = out[1];
		z = out[2];
		w = out[3];
	}

	const inline vec3f crossProd(const vec3f& v)
	{
		return vec3f((y * v.z - z * v.y), (z * v.x - x * v.z), (x * v.y - y * v.x));
	}

	const inline vec3f crossProd(const vec4f& v)
	{
		return vec3f((y * v.z - z * v.y), (z * v.x - x * v.z), (x * v.y - y * v.x));
	}

	const inline float dotProd(const vec3f& v)
	{
		return ((x * v.x) + (y * v.y) + (z * v.z));
	}


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



static vec4f operator* (const mat4f& l, const vec4f& r)
{
	vec4f out;
	vec4 oper, res;
	vec4_set(oper, r.x, r.y, r.z, r.w);
	mat4_mul_vec4(l.m, oper, res);
	out = res;

	return out;
}

static mat4f operator* (const mat4f& l, const mat4f& r)
{
	mat4 ret;
	mat4_mul_mat4(l.m, r.m, ret);

	mat4f F(ret);
	return F;
}

//////////////////////////////////////////////////////////////////////////////////////////

struct rectf
{
	float l, r, b, t;
};

struct texture_rect
{
	float leftU;
	float topV;
	float rightU;
	float bottomV;
};


static void QMATH_GET_TEXTURE_COORDINATES(const vec2f& srcDims, const rectf* src, const vec2f& destDims, const rectf* dest, texture_rect* pCoords)
{
	float tU, tV;
	pCoords->leftU = 0.0f;
	pCoords->topV = 1.0f;
	pCoords->rightU = 1.0f;
	pCoords->bottomV = 0.0f;

	if(src)
	{
		tU = 1.0f / srcDims.x;
		tV = 1.0f / srcDims.y;
		
		pCoords->leftU += src->l * tU;
		pCoords->topV -= src->t * tV;
		pCoords->rightU -= (srcDims.x - src->r) * tU;
		pCoords->bottomV += (srcDims.y - src->b) * tV;	
	}

	if(dest)
	{
		tU = 1.0f / destDims.x;
		tV = 1.0f / destDims.y;

		pCoords->leftU -= dest->l * tU;
//		pCoords->topV += dest->t * tV;
		pCoords->topV += tV;
		pCoords->rightU += (destDims.x - dest->r) * tU;
//		pCoords->bottomV -= (destDims.y - dest->b) * tV;
		pCoords->bottomV -= tV;

		if(src)
		{
			pCoords->bottomV = 1.0f - pCoords->bottomV;
			pCoords->topV = 1.0f - pCoords->topV;
		}
	}
}

static void QMATH_GET_TEXTURE_RECT(const vec2f& dims, rectf* rect)
{
	if(!rect)
		return;

	rect->l = 0.0f;
	rect->t = dims.y;
	rect->r = dims.x;
	rect->b = 0.0f;		
}

static void QMATH_INFLATE_RECT(rectf* rect, float dx, float dy)
{
	rect->r += dx;
	rect->l -= dx;
	rect->b -= dy;
	rect->t += dy;
}


static void QMATH_GET_SAMPLE4X4_OFFSETS(const int& bbWidth, const int& bbHeight, vec2f* avSampleOffsets)
{
	if(!avSampleOffsets)
		return;
	
	float tu = 1.0F / bbWidth;
	float tv = 1.0F / bbHeight;
	
	int index = 0;
	for(int y = 0; y < 4; ++y)
	{
		for(int x = 0; x < 4; ++x)
		{
			avSampleOffsets[index].x = (x - 1.5F) * tu;
			avSampleOffsets[index].y = (y - 1.5F) * tv;
			++index;
		}
	}	
}


static float QMATH_GET_GAUSSIAN_DIST(const float& x, const float& y, const float& rho)
{
	float g = 1.0f / sqrtf(2.0f * PI_F * rho * rho);
	g *= expf(-(x * x + y * y) / (2.0f * rho * rho));
	return g;
}


static void QMATH_GET_GAUSSIAN5X5_OFFSETS(unsigned int w, unsigned int h, vec2f* offsets, vec4f* weights, float mul)
{
	float tu = 1.0f / (float)w;
	float tv = 1.0f / (float)h;

	vec4f white(1.0f, 1.0f, 1.0f, 1.0f);
	float totalWeight = 0.0f;
	int idx = 0;

	for(int x = -2; x <= 2; ++x)
	{
		for(int y = -2; y <= 2; ++y)
		{
			if(abs(x) + abs(y) > 2)
				continue;

			offsets[idx].x = x * tu;
			offsets[idx].y = y * tv;
			weights[idx].x = white.x * QMATH_GET_GAUSSIAN_DIST((float)x, (float)y, 1.0f);
			weights[idx].y = white.y * QMATH_GET_GAUSSIAN_DIST((float)x, (float)y, 1.0f);
			weights[idx].z = white.z * QMATH_GET_GAUSSIAN_DIST((float)x, (float)y, 1.0f);
			weights[idx].w = white.w * QMATH_GET_GAUSSIAN_DIST((float)x, (float)y, 1.0f);
			totalWeight += weights[idx].x;
			++idx;
		}
	}

	for(int i = 0; i < idx; ++i)
	{
		weights[i].x /= totalWeight; weights[i].x *= mul;
		weights[i].y /= totalWeight; weights[i].y *= mul;
		weights[i].z /= totalWeight; weights[i].z *= mul;
		weights[i].w /= totalWeight; weights[i].w *= mul;
	}
}

static void QMATH_GET_SAMPLE2X2_OFFSETS(int w, int h, vec2f* avSampleOffsets)
{
	if(!avSampleOffsets)
		return;

	float tu = 1.0f / w;
	float tv = 1.0f / h;

	int idx = 0;
	for(int y = 0; y < 2; ++y)
	{
		for(int x = 0; x < 2; ++x)
		{
			avSampleOffsets[idx].x = (x - 0.5f) * tu;
			avSampleOffsets[idx].y = (y - 0.5f) * tv;
			++idx;
		}
	}
}


static void QMATH_GET_BLOOM_OFFSETS(int size, float texCoordOffset[15], vec4f* colorWeight, float dev, float mul)
{
	int i = 0;
	float tu = 1.0f / size;
	
	float weight = mul * QMATH_GET_GAUSSIAN_DIST(0.0f, 0.0f, dev);
	colorWeight[0].x = weight;
	colorWeight[0].y = weight;
	colorWeight[0].z = weight;
	colorWeight[0].w = 1.0f;
	
	texCoordOffset[0] = 0.0f;
	
	for(i = 1; i < 8; ++i)
	{
		weight = mul * QMATH_GET_GAUSSIAN_DIST((float)i, 0.0f, dev);
		texCoordOffset[i] = i * tu;
		colorWeight[i].x = weight;
		colorWeight[i].y = weight;
		colorWeight[i].z = weight;
		colorWeight[i].w = 1.0f;
	}	

	for(i = 8; i < 15; ++i)
	{
		colorWeight[i].set(colorWeight[i - 7]);
		texCoordOffset[i] = -(texCoordOffset[i - 7]);
	}
}


static void QMATH_CREATE_VERTEX_NORMALS(const vec3f* verts, const uint32_t& nVerts, const uint32_t* polys, const uint32_t& nPolys, float* norms)
{
	unsigned int* tmp;
	vec3f a, b, c, tmpNorm;
	vec3f tmpList[64];
	int counter;
	
	// iterate through every vertex //
	for(unsigned int j = 0; j < nVerts; ++j)
	{
		counter = 0;
		
		// check every poly for this vert //
		for(unsigned int v = 0; v < nPolys; ++v)
		{
			tmp = (unsigned int*)&polys[v * 3];
			if(*tmp != j)
			{
				++tmp;
				if(*tmp != j)
				{
					++tmp;
					if(*tmp != j)
						continue;
				}
			}
			
			// if we're here, this vert is in the tri, so fetch verts //
			unsigned int i1, i2, i3;
			i1 = polys[v * 3];
			i2 = polys[v * 3 + 1];
			i3 = polys[v * 3 + 2];
			
			a.x = verts[i1].x; 
			a.y = verts[i1].y;
			a.z = verts[i1].z;
			b.x = verts[i2].x;
			b.y = verts[i2].y;
			b.z = verts[i2].z;
			c.x = verts[i3].x;
			c.y = verts[i3].y;
			c.z = verts[i3].z;
			
			vec3f ab, ac;
			ac = a - c;
			ab = a - b;
//			tmpNorm = QMATH_VEC3F_CROSSPROD(ab, ac);
			tmpNorm = ac.crossProd(ab); //QMATH_VEC3F_CROSSPROD(ac, ab);
			tmpNorm.normalize();
//			QMATH_VEC3F_NORMALIZE(tmpNorm);

			// check through the buffer list for this vert and look for dupes //
			for(int z = 0; z < counter; ++z)
			{
				if(tmpNorm == tmpList[z])
					goto end;
			}
			
			// add to buffer list //
			tmpList[counter].x = tmpNorm.x;
			tmpList[counter].y = tmpNorm.y;
			tmpList[counter].z = tmpNorm.z;
			
			++counter;
			end:;
		}
		
		tmpNorm.x = 0.0F;
		tmpNorm.y = 0.0F;
		tmpNorm.z = 0.0F;
		for(int q = 0; q < counter; ++q)
		{
			tmpNorm.x += tmpList[q].x;
			tmpNorm.y += tmpList[q].y;
			tmpNorm.z += tmpList[q].z;
		}
		
		// copy to norm list //
//		QMATH_VEC3F_NORMALIZE(tmpNorm);
		tmpNorm.normalize();
		norms[j * 3] = tmpNorm.x;
		norms[j * 3 + 1] = tmpNorm.y;
		norms[j * 3 + 2] = tmpNorm.z;
	}
}

struct tan_index
{
	vec3f tan;
	int32_t vertRef;
};


static void QMATH_CREATE_TANGENT_SPACE(const vec3f* verts, const uint32_t& nVerts, const uint32_t* polys, const uint32_t& nPolys,
									   const vec2f* texcoords, const vec3f* norms, vec3f* tangent)
{
	tan_index* tan1 = static_cast<tan_index*>(calloc(nPolys * 3, sizeof(tan_index)));
	
	if(!tan1)
		return;
	
	for(unsigned int i = 0; i < nPolys; ++i)
	{
		unsigned int i1 = polys[i * 3];
		unsigned int i2 = polys[i * 3 + 1];
		unsigned int i3 = polys[i * 3 + 2];
		
		tan1[i * 3].vertRef = i1;
		tan1[i * 3 + 1].vertRef = i2;
		tan1[i * 3 + 2].vertRef = i3;
	
		vec3f a, b, c;
		a = verts[i1];
		b = verts[i2];
		c = verts[i3];
//		QMATH_VEC3F_COPY(a, verts[i1]);
//		QMATH_VEC3F_COPY(b, verts[i2]);
//		QMATH_VEC3F_COPY(c, verts[i3]);

		vec2f ta, tb, tc;
		ta = texcoords[i1];
		tb = texcoords[i2];
		tc = texcoords[i3];
//		QMATH_VEC2F_COPY(ta, texcoords[i1]);
//		QMATH_VEC2F_COPY(tb, texcoords[i2]);
//		QMATH_VEC2F_COPY(tc, texcoords[i3]);

		float x1 = b.x - a.x;
		float x2 = c.x - a.x;
		float y1 = b.y - a.y;
		float y2 = c.y - a.y;
		float z1 = b.z - a.z;
		float z2 = c.z - a.z;
		
		float s1 = tb.x - ta.x;
		float s2 = tc.x - ta.x;
		float t1 = tb.y - ta.y;
		float t2 = tc.y - ta.y;
		
		float r = 1.0F / (s1 * t2 - s2 * t1);
		vec3f sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r, (t2 * z1 - t1 * z2) * r);
		vec3f tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r, (s1 * z2 - s2 * z1) * r);

		tan1[i * 3].tan = tan1[i * 3].tan + sdir;
		tan1[i * 3 + 1].tan = tan1[i * 3 + 1].tan + sdir;
		tan1[i * 3 + 2].tan = tan1[i * 3 + 2].tan + sdir;
		
		vec3f n, t, tan;
		
		for(unsigned int j = 0; j < 3; ++j )
		{
//			QMATH_VEC3F_COPY(n, norms[tan1[i * 3 + j].vertRef]);
			n = norms[tan1[i * 3 + j].vertRef];
//			QMATH_VEC3F_COPY(t, tan1[i * 3 + j].tan);
			t = tan1[i * 3 + j].tan;
			
//			float dp3 = QMATH_VEC3F_DOTPROD(n, t);
			float dp3 = n.dotProd(t);
			n = n * dp3;
			tan = t - n;
			tan.normalize();
//			QMATH_VEC3F_NORMALIZE(tan);
//			QMATH_VEC3F_COPY(tangent[tan1[i * 3 + j].vertRef], tan);
			tangent[tan1[i * 3 + j].vertRef] = tan;
		}
	}
	
	if(tan1)
	{
		free(tan1);
		tan1 = NULL;
	}	
}

#endif
