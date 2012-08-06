#include "fixvector3d.h"

v3d v3d_add(v3d a, v3d b)
{
	a.x = fix16_add(a.x, b.x);
	a.y = fix16_add(a.y, b.y);
	a.z = fix16_add(a.z, b.z);
	return a;
}

v3d v3d_sub(v3d a, v3d b)
{
	a.x = fix16_sub(a.x, b.x);
	a.y = fix16_sub(a.y, b.y);
	a.z = fix16_sub(a.z, b.z);
	return a;
}

v3d v3d_mul(v3d a, fix16_t b)
{
	a.x = fix16_mul(a.x, b);
	a.y = fix16_mul(a.y, b);
	a.z = fix16_mul(a.z, b);
	return a;
}

v3d v3d_div(v3d a, fix16_t b)
{
	a.x = fix16_div(a.x, b);
	a.y = fix16_div(a.y, b);
	a.z = fix16_div(a.z, b);
	return a;
}

// Norm
fix16_t v3d_norm(v3d a)
{
	a.x >>= 8;
	a.y >>= 8;
	a.z >>= 8;
	
	return fix16_sqrt(v3d_dot(a, a)) * 256;
}

// Dot product
fix16_t v3d_dot(v3d a, v3d b)
{
	a.x = fix16_mul(a.x, b.x);
	a.y = fix16_mul(a.y, b.y);
	a.z = fix16_mul(a.z, b.z);
	return fix16_add(a.x, fix16_add(a.y, a.z));
}

// Cross product
v3d v3d_cross(v3d a, v3d b)
{
	v3d result;
	result.x = fix16_sub(fix16_mul(a.y, b.z), fix16_mul(a.z, b.y));
	result.y = fix16_sub(fix16_mul(a.z, b.x), fix16_mul(a.x, b.z));
	result.z = fix16_sub(fix16_mul(a.x, b.y), fix16_mul(a.y, b.x));
	return result;
}