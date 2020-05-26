#ifndef MATH3D_H
#define MATH3D_H

#include <cmath>
#include <cassert>
#define PI  3.14159265358979323846
struct Triangle {
	int f1, f2, f3;
	bool isDegenerated = false;
};
class Point3D {
public:
	float x, y, z, value;
	int index = -1;
	Point3D()
	{
		x = 0;
		y = 0;
		z = 0;
		value = 0;
		index = -1;
	}
	Point3D(float _x, float _y, float _z,float _value=0,int _index=-1)
	{
		x = _x;
		y = _y;
		z = _z;
		value = _value;
		index = _index;
	}

	Point3D& operator/ (const double value)
	{
		this->x /= value;
		this->y /= value;
		this->z /= value;
		return *this;
	}
	Point3D& operator- (const Point3D& p)
	{
		this->x -= p.x;
		this->y -= p.y;
		this->z -= p.z;
		return *this;
	}
	Point3D& operator= (const Point3D& p)
	{
		this->x = p.x;
		this->y = p.y;
		this->z = p.z;
		this->value = p.value;
		this->index = p.index;
		return *this;
	}
};
class Vector3D {
	public:
		float x;
		float y;
		float z;
		bool isDuplicated = false;
		Vector3D()
		{
			x = 0;
			y = 0;
			z = 0;
		}
		Vector3D(float _x, float _y, float _z)
		{
			x = _x;
			y = _y;
			z = _z;
		}
		Vector3D operator+(Vector3D& value)
		{
			Vector3D result;
			result.x = x + value.x;
			result.y = y + value.y;
			result.z = z + value.z;
			return result;
		}
		Vector3D operator-(Vector3D& value)
		{
			Vector3D result;
			result.x = x - value.x;
			result.y = y - value.y;
			result.z = z - value.z;
			return result;
		}
		Vector3D operator*(double value)
		{
			Vector3D result;
			result.x = x * value;
			result.y = y * value;
			result.z = z * value;
			return result;
		}
		Vector3D operator+=(Vector3D& value)
		{
			x += value.x;
			y += value.y;
			z += value.z;
			return *this;
		}
		Vector3D operator/=(double value)
		{
			x /= value;
			y /= value;
			z /= value;
			return *this;
		}

		Vector3D& operator= (const Vector3D& value)
		{
			this->x = value.x;
			this->y = value.y;
			this->z = value.z;
			return *this;
		}


		inline bool operator<(const Vector3D& value) const
		{
			if (x < value.x)
			{
				return true;
			}
			else if (x == value.x)
			{
				if (y < value.y)
				{
					return true;
				}
				else if (y == value.y)
				{
					if (z < value.z)
					{
						return true;
					}
				}
			}
			return false;
		}
		bool operator== (Vector3D& value)
		{
			if (x == value.x &&y == value.y&&z == value.z)
			{
				return true;
			}
			return false;
		}
		bool operator!= (Vector3D& value)
		{
			if (x == value.x &&y == value.y&&z == value.z)
			{
				return false;
			}
			return true;
		}
		Vector3D cross(Vector3D& value)
		{
			Vector3D result;
			result.x = y*value.z - value.y*z;
			result.y = x*value.z - value.x*z;
			result.z = x*value.y - value.x*y;
			return result;
		}
		double dot(Vector3D& value)
		{
			return x*value.x + y*value.y + z*value.z;
		}
};


struct edge {
	int id1, id2;//端点在PointSet中对应的索引
	int index[2];//等值点在Vertex中的索引
};
void normalize(Vector3D& v);

#endif
