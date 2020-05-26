#ifndef ISOSURFACE_H
#define ISOSURFACE_H

#include "math3D.h"

class Isosurface
{

public:
    virtual float valueAt(float x, float y, float z) const = 0;

    // returns a normalized vector
	virtual Vector3D gradientAt(float x, float y, float z) const{
		const float epsilon = 0.0001;

		float dx = valueAt(x + epsilon, y, z) - valueAt(x - epsilon, y, z);
		float dy = valueAt(x, y + epsilon, z) - valueAt(x, y - epsilon, z);
		float dz = valueAt(x, y, z + epsilon) - valueAt(x, y, z - epsilon);

		Vector3D result = { dx, dy, dz };
		normalize(result);
		return result;
	};

};

#endif
