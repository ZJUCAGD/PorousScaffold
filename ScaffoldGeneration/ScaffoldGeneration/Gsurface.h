#ifndef GSURFACE_H
#define GSURFACE_H

#include "Isosurface.h"
#include "math3D.h"

class Gsurface : public Isosurface
{

public:
    virtual float valueAt(float x, float y, float z) const
	{
		return (cosf(2 * PI*x) * sinf(2 * PI*y) + cosf(2 * PI*y) * sinf(2 * PI*z) + cosf(2 * PI*z) * sinf(2 * PI*x));
	};
    virtual Vector3D gradientAt(float x, float y, float z) const
	{
		float gx = cosf(2 * PI*x) * cosf(2 * PI*z) - sinf(2 * PI*x) * sinf(2 * PI*y);
		float gy = cosf(2 * PI*x) * cosf(2 * PI*y) - sinf(2 * PI*y) * sinf(2 * PI*z);
		float gz = cosf(2 * PI*y) * cosf(2 * PI*z) - sinf(2 * PI*x) * sinf(2 * PI*z);

		Vector3D result = { gx, gy, gz };
		normalize(result);

		return result;
	};

};

#endif
