#ifndef I2YSURFACE_H
#define I2YSURFACE_H

#include "Isosurface.h"
#include "math3D.h"

class I2Ysurface : public Isosurface
{
	
public:
	virtual float valueAt(float x, float y, float z) const
	{
		return (-2 * (sinf(4 * PI*x)*cosf(2 * PI*y)*sinf(2 * PI*z) + sinf(2 * PI*x)*sinf(4 * PI*y)*cosf(2 * PI*z) + cosf(2 * PI*x)*sinf(2 * PI*y)*sinf(4 * PI*z)) 
			+ cosf(4 * PI*x)*cosf(4 * PI*y) + cosf(4 * PI*y)*cosf(4 * PI*z) + cosf(4 * PI*x)*cosf(4 * PI*z));
	};
	virtual Vector3D gradientAt(float x, float y, float z) const
	{
		float gx = 2*cosf(4 * PI*x)*cosf(2 * PI*y)*sinf(2 * PI*z) + cosf(2 * PI*x)*sinf(4 * PI*y)*cosf(2 * PI*z) - sinf(2 * PI*x)*sinf(2 * PI*y)*sinf(4 * PI*z)
			+ sinf(4 * PI*x)*(cosf(4 * PI*y) + cosf(4 * PI*z));
		float gy = -sinf(4 * PI*x)*sinf(2 * PI*y)*sinf(2 * PI*z) + 2*sinf(2 * PI*x)*cosf(4 * PI*y)*cosf(2 * PI*z) + cosf(2 * PI*x)*cosf(2 * PI*y)*sinf(4 * PI*z)
			+ sinf(4 * PI*y)*(cosf(4 * PI*x) + cosf(4 * PI*z));
		float gz = sinf(4 * PI*x)*cosf(2 * PI*y)*cosf(2 * PI*z) - sinf(2 * PI*x)*sinf(4 * PI*y)*sinf(2 * PI*z) + 2*cosf(2 * PI*x)*sinf(2 * PI*y)*cosf(4 * PI*z)
			+ sinf(4 * PI*z)*(cosf(4 * PI*x) + cosf(4 * PI*y));

		Vector3D result = { gx, gy, gz };
		normalize(result);

		return result;
	};

};

#endif