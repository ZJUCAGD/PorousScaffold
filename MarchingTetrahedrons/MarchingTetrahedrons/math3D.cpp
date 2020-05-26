#include "math3D.h"

void normalize(Vector3D& v)
{
    float norm = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
	assert(norm > 0);
	v.x /= norm;
	v.y /= norm;
	v.z /= norm;
}
