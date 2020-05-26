#ifndef LSURFACE_H
#define LSURFACE_H

#include "Isosurface.h"
#include "math3D.h"

class Lsurface : public Isosurface
{
	//0.5*(sin(2*X).*cos(Y).*sin(Z)+sin(2*Y).*cos(Z).*sin(X)+sin(2*Z).*cos(X).*sin(Y))-0.5*(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))+0.15
public:
	virtual float valueAt(float x, float y, float z) const
	{
		return (0.5*(sinf(4 * PI*x)*cosf(2 * PI*y)*sinf(2 * PI*z) + sinf(2 * PI*x)*sinf(4 * PI*y)*cosf(2 * PI*z) + cosf(2 * PI*x)*sinf(2 * PI*y)*sinf(4 * PI*z))
			- 0.5*(cosf(4 * PI*x)*cosf(4 * PI*y) + cosf(4 * PI*y)*cosf(4 * PI*z) + cosf(4 * PI*z)*cosf(4 * PI*x)) + 0.15);
	};
	virtual Vector3D gradientAt(float x, float y, float z) const
	{
		float gx = 0.5*(2 * cosf(4 * PI*x)*cosf(2 * PI*y)*sinf(2 * PI*z) + cosf(2 * PI*x)*sinf(4 * PI*y)*cosf(2 * PI*z) - sinf(2 * PI*x)*sinf(2 * PI*y)*sinf(4 * PI*z))
			+ sinf(4 * PI*x)*(cosf(4 * PI*y) + cosf(4 * PI*z));
		float gy = 0.5*(2 * cosf(4 * PI*y)*cosf(2 * PI*x)*sinf(2 * PI*z) + cosf(2 * PI*y)*sinf(4 * PI*z)*cosf(2 * PI*x) - sinf(2 * PI*y)*sinf(2 * PI*z)*sinf(4 * PI*x))
			+ sinf(4 * PI*y)*(cosf(4 * PI*x)+cosf(4 * PI*z));
		float gz = 0.5*(2 * cosf(4 * PI*z)*cosf(2 * PI*x)*sinf(2 * PI*y) + cosf(2 * PI*z)*sinf(4 * PI*x)*cosf(2 * PI*y) - sinf(2 * PI*z)*sinf(2 * PI*x)*sinf(4 * PI*y))
			+ sinf(4 * PI*z)*(cosf(4 * PI*x) +cosf(4 * PI*y));

		Vector3D result = { gx, gy, gz };
		normalize(result);

		return result;
	};

};

#endif