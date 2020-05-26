#ifndef TubularGSURFACE_H
#define TubularGSURFACE_H

#include "Isosurface.h"
#include "math3D.h"

class TubularGsurface : public Isosurface
{
	
public:
	virtual float valueAt(float x, float y, float z) const
	{
		return (10 * (cosf(2 * PI*x)*sinf(2 * PI*y) + cosf(2 * PI*y)*sinf(2 * PI*z) + cosf(2 * PI*z)*sinf(2 * PI*x))
			- 0.5*(cosf(4 * PI*x)*cosf(4 * PI*y) + cosf(4 * PI*y)*cosf(4 * PI*z) + cosf(4 * PI*z)*cosf(4 * PI*x)) - 14);
	};
	virtual Vector3D gradientAt(float x, float y, float z) const
	{
		float gx = 10 * (-sinf(2 * PI*x) *sinf(2 * PI*y) + cosf(2 * PI*x)*cosf(2 * PI*z)) + sinf(4 * PI*x)*(cosf(4 * PI*y) + cosf(4 * PI*z));
		float gy = 10 * (-sinf(2 * PI*y) *sinf(2 * PI*z) + cosf(2 * PI*y)*cosf(2 * PI*x)) + sinf(4 * PI*y)*(cosf(4 * PI*y) + cosf(4 * PI*x));
		float gz = 10 * (-sinf(2 * PI*z) *sinf(2 * PI*x) + cosf(2 * PI*z)*cosf(2 * PI*y)) + sinf(4 * PI*z)*(cosf(4 * PI*x) + cosf(4 * PI*z));

		Vector3D result = { gx, gy, gz };
		normalize(result);

		return result;
	};

};

#endif