#ifndef TubularPSURFACE_H
#define TubularPSURFACE_H

#include "Isosurface.h"
#include "math3D.h"

class TubularPsurface : public Isosurface
{
	
public:
	virtual float valueAt(float x, float y, float z) const
	{
		return (10 * (cosf(2 * PI*x) + cos(2 * PI*y) + cos(2 * PI*z)) - 5.1*(cosf(2 * PI*x)*cosf(2 * PI*y) + cosf(2 * PI*y)*cosf(2 * PI*z) + cosf(2 * PI*z)*cosf(2 * PI*x)) - 14.6);
	};
	virtual Vector3D gradientAt(float x, float y, float z) const
	{
		float gx = -10 * sinf(2 * PI*x) + 5.1*sinf(2 * PI*x)*(cosf(2 * PI*y) + cosf(2 * PI*z));
		float gy = -10 * sinf(2 * PI*y) + 5.1*sinf(2 * PI*y)*(cosf(2 * PI*x) + cosf(2 * PI*z));
		float gz = -10 * sinf(2 * PI*z) + 5.1*sinf(2 * PI*z)*(cosf(2 * PI*x) + cosf(2 * PI*y));

		Vector3D result = { gx, gy, gz };
		normalize(result);

		return result;
	};

};

#endif