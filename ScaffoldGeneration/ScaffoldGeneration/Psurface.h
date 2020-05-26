#ifndef PSURFACE_H
#define PSURFACE_H

#include "Isosurface.h"
#include "math3D.h"

class Psurface : public Isosurface
{

public:
	virtual float valueAt(float x, float y, float z) const
	{
		return (cosf(2 * PI*x) + cosf(2 * PI*y) + cosf(2 * PI*z));
	};
	virtual Vector3D gradientAt(float x, float y, float z) const
	{
		float gx = -sinf(2 * PI*x);
		float gy = -sinf(2 * PI*y);
		float gz = -sinf(2 * PI*z);

		Vector3D result = { gx, gy, gz };
		normalize(result);

		return result;
	};

};

#endif