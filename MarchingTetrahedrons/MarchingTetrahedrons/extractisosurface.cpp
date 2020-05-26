#include "extractisosurface.h"
#include <ctime>
#include <windows.h>


void cutVertex(Vector3D& _point)
{//避免计算误差
	_point.x = 1.0*(int(_point.x * 1000000 + 0.5)) / 1000000;
	_point.y = 1.0*(int(_point.y * 1000000 + 0.5)) / 1000000;
	_point.z = 1.0*(int(_point.z * 1000000 + 0.5)) / 1000000;
}

//不考虑索引关系
#pragma region 

static inline void calInterVert(const Isosurface& surface, const Point3D& p1, const Point3D& p2, float isolevel, Vector3D& interP)
{

    float v1 = p1.value;
    float v2 = p2.value;

    float x, y, z;

    if (v2 == v1) 
	{
        x = (p1.x + p2.x) / 2.0f;
        y = (p1.y + p2.y) / 2.0f;
        z = (p1.z + p2.z) / 2.0f;
    } 
	else 
	{

        /*

         <----+-----+---+----->
              v1    |   v2
                 isolevel


         <----+-----+---+----->
              0     |   1
                  interp

         */


        // interp == 0: vert should be at p1
        // interp == 1: vert should be at p2
        float interp = (isolevel - v1) / (v2 - v1);
        float oneMinusInterp = 1 - interp;

        x = p1.x * oneMinusInterp + p2.x * interp;
        y = p1.y * oneMinusInterp + p2.y * interp;
        z = p1.z * oneMinusInterp + p2.z * interp;
    }
	interP.x = x;
	interP.y = y;
	interP.z = z;
	cutVertex(interP);
}

static inline void calInterPoint(const Isosurface& surface, const Point3D& p1, const Point3D& p2, float isolevel, Point3D& interP)
{

	float v1 = p1.value;
	float v2 = p2.value;

	float x, y, z;

	if (v2 == v1)
	{
		x = (p1.x + p2.x) / 2.0f;
		y = (p1.y + p2.y) / 2.0f;
		z = (p1.z + p2.z) / 2.0f;
	}
	else
	{

		/*

		<----+-----+---+----->
		v1    |   v2
		isolevel


		<----+-----+---+----->
		0     |   1
		interp

		*/


		// interp == 0: vert should be at p1
		// interp == 1: vert should be at p2
		float interp = (isolevel - v1) / (v2 - v1);
		float oneMinusInterp = 1 - interp;

		x = p1.x * oneMinusInterp + p2.x * interp;
		y = p1.y * oneMinusInterp + p2.y * interp;
		z = p1.z * oneMinusInterp + p2.z * interp;
	}
	interP.x = x;
	interP.y = y;
	interP.z = z;
	interP.value = isolevel;
}

static void extractFromTetrahedron(const Isosurface& surface, const Point3D p[4], float isolevel, vector<Vector3D>& Vertex, vector<Triangle>& Face)
{
	Vector3D intersectionP1, intersectionP2, intersectionP3, intersectionP4;
	Triangle face;
    /*

     Tetrahedron layout:

           0
           *
          /|
         / |
      3 *-----* 1
         \ | /
          \|/
           *
           2
     */

    unsigned char index = 0;
    for (int i = 0; i < 4; ++i)
        if (p[i].value < isolevel)
            index |= (1 << i);

    switch (index) {

        // we don't do anything if everyone is inside or outside
        case 0x00:
        case 0x0F:
            break;

        // only vert 0 is inside
        case 0x01:
			calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
			calInterVert(surface, p[0], p[3], isolevel, intersectionP2);
			calInterVert(surface, p[0], p[2], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}
            break;

        // only vert 1 is inside
        case 0x02:
			calInterVert(surface, p[1], p[0], isolevel, intersectionP1);
			calInterVert(surface, p[1], p[2], isolevel, intersectionP2);
			calInterVert(surface, p[1], p[3], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}
            break;

        // only vert 2 is inside
        case 0x04:
			calInterVert(surface, p[2], p[0], isolevel, intersectionP1);
			calInterVert(surface, p[2], p[3], isolevel, intersectionP2);
			calInterVert(surface, p[2], p[1], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}
            break;

        // only vert 3 is inside
        case 0x08:
			calInterVert(surface, p[3], p[1], isolevel, intersectionP1);
			calInterVert(surface, p[3], p[2], isolevel, intersectionP2);
			calInterVert(surface, p[3], p[0], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}
            break;

        // verts 0, 1 are inside
        case 0x03:
			calInterVert(surface, p[3], p[0], isolevel, intersectionP1);
			calInterVert(surface, p[2], p[0], isolevel, intersectionP2);
			calInterVert(surface, p[1], p[3], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}

			calInterVert(surface, p[2], p[1], isolevel, intersectionP4);
			if (intersectionP4 != intersectionP2 &&intersectionP4 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				Vertex.push_back(intersectionP4);
				face.f1 = Vertex.size() - 2;
				face.f2 = Vertex.size() - 3;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}
            break;

        // verts 0, 2 are inside
        case 0x05:
			calInterVert(surface, p[1], p[0], isolevel, intersectionP1);
			calInterVert(surface, p[3], p[0], isolevel, intersectionP2);
			calInterVert(surface, p[1], p[2], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}

			calInterVert(surface, p[2], p[3], isolevel, intersectionP4);
			if (intersectionP4 != intersectionP2 &&intersectionP4 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				Vertex.push_back(intersectionP4);
				face.f1 = Vertex.size() - 2;
				face.f2 = Vertex.size() - 3;
				face.f3 = Vertex.size()-1;
				Face.push_back(face);
			}
            break;

        // verts 0, 3 are inside
        case 0x09:
			calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
			calInterVert(surface, p[1], p[3], isolevel, intersectionP2);
			calInterVert(surface, p[0], p[2], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}

			calInterVert(surface, p[3], p[2], isolevel, intersectionP4);
			if (intersectionP4 != intersectionP2 &&intersectionP4 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				Vertex.push_back(intersectionP4);
				face.f1 = Vertex.size() - 2;
				face.f2 = Vertex.size() - 3;
				face.f3 = Vertex.size()-1;
				Face.push_back(face);
			}
            break;

        // verts 1, 2 are inside
        case 0x06:
			calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
			calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
			calInterVert(surface, p[1], p[3], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}

			calInterVert(surface, p[3], p[2], isolevel, intersectionP4);
			if (intersectionP4 != intersectionP2 &&intersectionP4 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				Vertex.push_back(intersectionP4);
				face.f1 = Vertex.size() - 2;
				face.f2 = Vertex.size() - 3;
				face.f3 = Vertex.size()-1;
				Face.push_back(face);
			}
            break;

        // verts 2, 3 are inside
        case 0x0C:
			calInterVert(surface, p[3], p[0], isolevel, intersectionP1);
			calInterVert(surface, p[1], p[3], isolevel, intersectionP2);
			calInterVert(surface, p[2], p[0], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}

			calInterVert(surface, p[2], p[1], isolevel, intersectionP4);
			if (intersectionP4 != intersectionP2 &&intersectionP4 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				Vertex.push_back(intersectionP4);
				face.f1 = Vertex.size() - 2;
				face.f2 = Vertex.size() - 3;
				face.f3 = Vertex.size()-1;
				Face.push_back(face);
			}
            break;

        // verts 1, 3 are inside
        case 0x0A:
			calInterVert(surface, p[1], p[0], isolevel, intersectionP1);
			calInterVert(surface, p[1], p[2], isolevel, intersectionP2);
			calInterVert(surface, p[3], p[0], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}

			calInterVert(surface, p[2], p[3], isolevel, intersectionP4);
			if (intersectionP4 != intersectionP2 &&intersectionP4 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				Vertex.push_back(intersectionP4);
				face.f1 = Vertex.size() - 2;
				face.f2 = Vertex.size() - 3;
				face.f3 = Vertex.size()-1;
				Face.push_back(face);
			}
            break;

        // verts 0, 1, 2 are inside
        case 0x07:
			calInterVert(surface, p[3], p[0], isolevel, intersectionP1);
			calInterVert(surface, p[3], p[2], isolevel, intersectionP2);
			calInterVert(surface, p[3], p[1], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}
            break;

        // verts 0, 1, 3 are inside
        case 0x0B:
			calInterVert(surface, p[2], p[1], isolevel, intersectionP1);
			calInterVert(surface, p[2], p[3], isolevel, intersectionP2);
			calInterVert(surface, p[2], p[0], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}
            break;

        // verts 0, 2, 3 are inside
        case 0x0D:
			calInterVert(surface, p[1], p[0], isolevel, intersectionP1);
			calInterVert(surface, p[1], p[3], isolevel, intersectionP2);
			calInterVert(surface, p[1], p[2], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}
            break;

        // verts 1, 2, 3 are inside
        case 0x0E:
			calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
			calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
			calInterVert(surface, p[0], p[3], isolevel, intersectionP3);
			if (intersectionP1 != intersectionP2 &&intersectionP1 != intersectionP3 &&intersectionP2 != intersectionP3)
			{
				Vertex.push_back(intersectionP1);
				Vertex.push_back(intersectionP2);
				Vertex.push_back(intersectionP3);
				face.f1 = Vertex.size() - 3;
				face.f2 = Vertex.size() - 2;
				face.f3 = Vertex.size() - 1;
				Face.push_back(face);
			}
            break;

        // what is this I don't even
        default:
            assert(false);
    }

}

static void extractTriangleOuter(const Isosurface& surface, const Point3D p[3], float isolevel, vector<Vector3D>& Vertex, vector<Triangle>& Face)
{
	Vector3D intersectionP1, intersectionP2;
	Triangle face;
	if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//0,1,2 outer
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//0,1 outer 2 inner
		calInterVert(surface, p[1], p[2], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(intersectionP1);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);

		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP1);
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//0,2 outer 1 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[1], p[2], isolevel, intersectionP2);
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP1);
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);

		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP2);
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0 outer 1,2 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP1);
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//1,2 outer 0 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);

		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//1 outer 0,2 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[1], p[2], isolevel, intersectionP2);
		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//2 outer 0,1 inner
		calInterVert(surface, p[1], p[2], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0,1,2 inner
	}

}

static void extractTriangleInner(const Isosurface& surface, const Point3D p[3], float isolevel, vector<Vector3D>& Vertex, vector<Triangle>& Face)
{
	Vector3D intersectionP1, intersectionP2;
	Triangle face;
	if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//0,1,2 outer
	}
	else if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//0,1 outer 2 inner
		calInterVert(surface, p[1], p[2], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//0,2 outer 1 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[1], p[2], isolevel, intersectionP2);
		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0 outer 1,2 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);

		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//1,2 outer 0 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP1);
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//1 outer 0,2 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[1], p[2], isolevel, intersectionP2);
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP1);
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);

		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP2);
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//2 outer 0,1 inner

		calInterVert(surface, p[1], p[2], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(intersectionP1);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);

		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP1);
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0,1,2 inner
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}

}

#pragma endregion



//考虑索引关系
#pragma region 

void addFace(vector<Triangle>& Face, int f1, int f2, int f3)
{
	Triangle face;
	if (f1 != f2&&f1 != f3&&f2 != f3){
		face.f1 = f1;
		face.f2 = f2;
		face.f3 = f3;
		Face.push_back(face);
	}
}

static inline void calInterVert2(const Isosurface& surface, const int edge_id, float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, vector<Vector3D>& Vertex, int& index, int sur_id)
{
	if (Edge[edge_id].index[sur_id] == -1)
	{
		float value1 = pointset[Edge[edge_id].id1].value;
		float value2 = pointset[Edge[edge_id].id2].value;

		float x, y, z;

		if (value2 == value1)
		{
			x = (pointset[Edge[edge_id].id1].x + pointset[Edge[edge_id].id2].x) / 2.0f;
			y = (pointset[Edge[edge_id].id1].y + pointset[Edge[edge_id].id2].y) / 2.0f;
			z = (pointset[Edge[edge_id].id1].z + pointset[Edge[edge_id].id2].z) / 2.0f;
		}
		else
		{
			float interp = (isolevel - value1) / (value2 - value1);
			interp = 1.0*(int(interp * 100000 + 0.5)) / 100000;
			float oneMinusInterp = 1 - interp;
			if (interp == 1)
			{
				if (pointset[Edge[edge_id].id2].index == -1)
				{
					x = pointset[Edge[edge_id].id2].x;
					y = pointset[Edge[edge_id].id2].y;
					z = pointset[Edge[edge_id].id2].z;
					Vertex.push_back(Vector3D(x, y, z));
					pointset[Edge[edge_id].id2].index = Vertex.size() - 1;
				}
				Edge[edge_id].index[sur_id] = pointset[Edge[edge_id].id2].index;
			}
			else if (interp == 0)
			{
				if (pointset[Edge[edge_id].id1].index == -1)
				{
					x = pointset[Edge[edge_id].id1].x;
					y = pointset[Edge[edge_id].id1].y;
					z = pointset[Edge[edge_id].id1].z;
					Vertex.push_back(Vector3D(x, y, z));
					pointset[Edge[edge_id].id1].index = Vertex.size() - 1;
				}
				Edge[edge_id].index[sur_id] = pointset[Edge[edge_id].id1].index;
			}
			else
			{
				x = pointset[Edge[edge_id].id1].x* oneMinusInterp + pointset[Edge[edge_id].id2].x* interp;
				y = pointset[Edge[edge_id].id1].y* oneMinusInterp + pointset[Edge[edge_id].id2].y* interp;
				z = pointset[Edge[edge_id].id1].z* oneMinusInterp + pointset[Edge[edge_id].id2].z* interp;
				Vertex.push_back(Vector3D(x, y, z));
				Edge[edge_id].index[sur_id] = Vertex.size() - 1;
			}
		}
		//cutVertex(p);
	}
	index = Edge[edge_id].index[sur_id];
}

static inline void getPointinsurface(const int id, vector<Point3D>& pointset, vector<Vector3D>& Vertex, int& index)
{
	if (pointset[id].index == -1)
	{
		Vertex.push_back(Vector3D(pointset[id].x, pointset[id].y, pointset[id].z));
		pointset[id].index = Vertex.size() - 1;
	}
	index = pointset[id].index;
}

static void extractTriangleFromTetrahedron(const Isosurface& surface, const Point3D p[4], const int edges[6], float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, vector<Vector3D>& Vertex, vector<Triangle>& Face,int sur_id)
{
	Triangle face;
	int index1 = -1, index2 = -1, index3 = -1, index4 = -1;
	unsigned char index = 0;
	for (int i = 0; i < 4; ++i)
	if (p[i].value < isolevel)
		index |= (1 << i);

	switch (index) {

		// we don't do anything if everyone is inside or outside
	case 0x00:
	case 0x0F:
		break;

		// only vert 0 is inside
	case 0x01:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[0]p[2]
		addFace(Face, index1, index2, index3);
		break;

		// only vert 1 is inside
	case 0x02:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[1]p[3]
		addFace(Face, index1, index2, index3);
		break;

		// only vert 2 is inside
	case 0x04:
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[2]p[3]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[1]p[2]
		addFace(Face, index1, index2, index3);
		break;

		// only vert 3 is inside
	case 0x08:
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[2]p[3]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[0]p[3]
		addFace(Face, index1, index2, index3);
		break;

		// verts 0, 1 are inside
	case 0x03:
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, Vertex, index4, sur_id);//p[1]p[2]
		addFace(Face, index1, index2, index3);
		addFace(Face, index3, index2, index4);
		break;

		// verts 0, 2 are inside
	case 0x05:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, Vertex, index4, sur_id);//p[2]p[3]
		addFace(Face, index1, index2, index3);
		addFace(Face, index3, index2, index4);
		break;

		// verts 0, 3 are inside
	case 0x09:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, Vertex, index4, sur_id);//p[2]p[3]
		addFace(Face, index1, index2, index3);
		addFace(Face, index3, index2, index4);
		break;

		// verts 1, 2 are inside
	case 0x06:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, Vertex, index4, sur_id);//p[2]p[3]
		addFace(Face, index1, index2, index3);
		addFace(Face, index3, index2, index4);
		break;

		// verts 2, 3 are inside
	case 0x0C:
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, Vertex, index4, sur_id);//p[1]p[2]
		addFace(Face, index1, index2, index3);
		addFace(Face, index3, index2, index4);
		break;

		// verts 1, 3 are inside
	case 0x0A:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, Vertex, index4, sur_id);//p[2]p[3]
		addFace(Face, index1, index2, index3);
		addFace(Face, index3, index2, index4);
		break;

		// verts 0, 1, 2 are inside
	case 0x07:
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[2]p[3]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[1]p[3]
		addFace(Face, index1, index2, index3);
		break;

		// verts 0, 1, 3 are inside
	case 0x0B:
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[2]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[0]p[2]
		addFace(Face, index1, index2, index3);
		break;

		// verts 0, 2, 3 are inside
	case 0x0D:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[1]p[2]
		addFace(Face, index1, index2, index3);
		break;

		// verts 1, 2, 3 are inside
	case 0x0E:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index3, sur_id);//p[0]p[3]
		addFace(Face, index1, index2, index3);
		break;

		// what is this I don't even
	default:
		assert(false);
	}
}

static void extractTriangleOuter2(const Isosurface& surface, const Point3D p[3], const int edges[3], const int Pindex[3], float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, vector<Vector3D>& Vertex, vector<Triangle>& Face,int sur_id)
{
	int index1 = -1, index2 = -1, index3 = -1, index4 = -1;
	Triangle face;
	if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//0,1,2 outer
		getPointinsurface(Pindex[0], pointset, Vertex, index1);
		getPointinsurface(Pindex[1], pointset, Vertex, index2);
		getPointinsurface(Pindex[2], pointset, Vertex, index3);
		addFace(Face, index1, index2, index3);
	}
	else if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//0,1 outer 2 inner
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[0], pointset, Vertex, index3);
		getPointinsurface(Pindex[1], pointset, Vertex, index4);
		addFace(Face, index3, index4, index1);
		addFace(Face, index3, index1, index2);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//0,2 outer 1 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[1]p[2]
		getPointinsurface(Pindex[0], pointset, Vertex, index3);
		getPointinsurface(Pindex[2], pointset, Vertex, index4);
		addFace(Face, index3, index1, index2);
		addFace(Face, index3, index2, index4);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0 outer 1,2 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[0], pointset, Vertex, index3);
		addFace(Face, index3, index1, index2);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//1,2 outer 0 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[1], pointset, Vertex, index3);
		getPointinsurface(Pindex[2], pointset, Vertex, index4);
		addFace(Face, index1, index3, index4);
		addFace(Face, index1, index4, index2);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//1 outer 0,2 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[1]p[2]
		getPointinsurface(Pindex[1], pointset, Vertex, index3);
		addFace(Face, index1, index3, index2);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//2 outer 0,1 inner
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[2], pointset, Vertex, index3);
		addFace(Face, index1, index3, index2);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0,1,2 inner
	}

}

static void extractTriangleInner2(const Isosurface& surface, const Point3D p[3], const int edges[3], const int Pindex[3], float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, vector<Vector3D>& Vertex, vector<Triangle>& Face,int sur_id)
{
	int index1 = -1, index2 = -1, index3 = -1, index4 = -1;
	Triangle face;
	if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//0,1,2 outer
	}
	else if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//0,1 outer 2 inner
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[2], pointset, Vertex, index3);
		addFace(Face, index1, index3, index2);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//0,2 outer 1 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[1]p[2]
		getPointinsurface(Pindex[1], pointset, Vertex, index3);
		addFace(Face, index1, index3, index2);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0 outer 1,2 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[1], pointset, Vertex, index3);
		getPointinsurface(Pindex[2], pointset, Vertex, index4);
		addFace(Face, index1, index3, index4);
		addFace(Face, index1, index4, index2);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//1,2 outer 0 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[0], pointset, Vertex, index3);
		addFace(Face, index3, index1, index2);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//1 outer 0,2 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[1]p[2]
		getPointinsurface(Pindex[0], pointset, Vertex, index3);
		getPointinsurface(Pindex[2], pointset, Vertex, index4);
		addFace(Face, index3, index1, index2);
		addFace(Face, index3, index2, index4);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//2 outer 0,1 inner
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, Vertex, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, Vertex, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[0], pointset, Vertex, index3);
		getPointinsurface(Pindex[1], pointset, Vertex, index4);
		addFace(Face, index3, index4, index1);
		addFace(Face, index3, index1, index2);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0,1,2 inner
		getPointinsurface(Pindex[0], pointset, Vertex, index1);
		getPointinsurface(Pindex[1], pointset, Vertex, index2);
		getPointinsurface(Pindex[2], pointset, Vertex, index3);
		addFace(Face, index1, index2, index3);
	}
}

static void extractTriangleFromSurface(const Isosurface& surface, const Point3D p[3], const int edges[3], const int Pindex[3], float maxisolevel, float minisolevel, vector<Point3D>& pointset, vector<edge>& Edge, vector<Vector3D>& Vertex, vector<Triangle>& Face)
{
	int index1 = -1, index2 = -1, index3 = -1, index4 = -1, index5 = -1;
	if (p[0].value >= maxisolevel && p[1].value >= maxisolevel && p[2].value >= maxisolevel)
	{//0,1,2 outer
	}
	else if (p[0].value >= maxisolevel && p[1].value >= maxisolevel && p[2].value < maxisolevel)
	{
		if (p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index1, 0);
			getPointinsurface(Pindex[2], pointset, Vertex, index2);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index3, 0);

			addFace(Face, index1, index2, index3);
		}
		else if (p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index1, 0);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index2, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index3, 1);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index4, 0);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
	}
	else if (p[0].value >= maxisolevel && p[1].value < maxisolevel && p[2].value >= maxisolevel)
	{
		if (p[1].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index1, 0);
			getPointinsurface(Pindex[1], pointset, Vertex, index2);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index3, 0);

			addFace(Face, index1, index2, index3);
		}
		else if (p[1].value < minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index1, 0);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index2, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index3, 1);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index4, 0);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
	}
	else if (p[0].value < maxisolevel && p[1].value >= maxisolevel && p[2].value >= maxisolevel)
	{
		if (p[0].value >= minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, Vertex, index1);
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index2, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index3, 0);

			addFace(Face, index1, index2, index3);
		}
		else if (p[0].value < minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index1, 1);
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index2, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index3, 0);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index4, 1);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
	}
	else if (p[0].value >= maxisolevel && p[1].value < maxisolevel && p[2].value < maxisolevel)
	{
		if (p[1].value >= minisolevel && p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index1, 0);
			getPointinsurface(Pindex[1], pointset, Vertex, index2);
			getPointinsurface(Pindex[2], pointset, Vertex, index3);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index4, 0);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
		else if (p[1].value >= minisolevel && p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index1, 0);
			getPointinsurface(Pindex[1], pointset, Vertex, index2);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index3, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index4, 1);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index5, 0);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
			addFace(Face, index1, index4, index5);
		}
		else if (p[1].value < minisolevel && p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index1, 0);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index2, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index3, 1);
			getPointinsurface(Pindex[2], pointset, Vertex, index4);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index5, 0);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
			addFace(Face, index1, index4, index5);
		}
		else if (p[1].value < minisolevel && p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index1, 0);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index2, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index3, 1);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index4, 0);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
	}
	else if (p[0].value < maxisolevel && p[1].value >= maxisolevel && p[2].value < maxisolevel)
	{
		if (p[0].value >= minisolevel && p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index1, 0);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index2, 0);
			getPointinsurface(Pindex[2], pointset, Vertex, index3);
			getPointinsurface(Pindex[0], pointset, Vertex, index4);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
		else if (p[0].value >= minisolevel && p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index1, 0);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index2, 0);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index3, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index4, 1);
			getPointinsurface(Pindex[0], pointset, Vertex, index5);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
			addFace(Face, index1, index4, index5);
		}
		else if (p[0].value < minisolevel && p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index1, 1);
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index2, 0);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index3, 0);
			getPointinsurface(Pindex[2], pointset, Vertex, index4);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index5, 1);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
			addFace(Face, index1, index4, index5);
		}
		else if (p[0].value < minisolevel && p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index1, 1);
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, Vertex, index2, 0);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index3, 0);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index4, 1);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
	}
	else if (p[0].value < maxisolevel && p[1].value < maxisolevel && p[2].value >= maxisolevel)
	{
		if (p[0].value >= minisolevel && p[1].value >= minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, Vertex, index1);
			getPointinsurface(Pindex[1], pointset, Vertex, index2);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index3, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index4, 0);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
		else if (p[0].value >= minisolevel && p[1].value < minisolevel)
		{
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index1, 0);
			getPointinsurface(Pindex[0], pointset, Vertex, index2);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index3, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index4, 1);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index5, 0);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
			addFace(Face, index1, index4, index5);
		}
		else if (p[0].value < minisolevel && p[1].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index1, 1);
			getPointinsurface(Pindex[1], pointset, Vertex, index2);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index3, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index4, 0);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index5, 1);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
			addFace(Face, index1, index4, index5);
		}
		else if (p[0].value < minisolevel && p[1].value < minisolevel)
		{
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index1, 1);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, Vertex, index2, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, Vertex, index3, 0);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index4, 1);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
	}
	else if (p[0].value < maxisolevel && p[1].value < maxisolevel && p[2].value < maxisolevel)
	{
		if (p[0].value > minisolevel && p[1].value > minisolevel && p[2].value > minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, Vertex, index1);
			getPointinsurface(Pindex[1], pointset, Vertex, index2);
			getPointinsurface(Pindex[2], pointset, Vertex, index3);

			addFace(Face, index1, index2, index3);
		}
		else if (p[0].value > minisolevel && p[1].value > minisolevel && p[2].value <= minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, Vertex, index1);
			getPointinsurface(Pindex[1], pointset, Vertex, index2);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index3, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index4, 1);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
		else if (p[0].value > minisolevel && p[1].value <= minisolevel && p[2].value > minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, Vertex, index1);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index2, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index3, 1);
			getPointinsurface(Pindex[2], pointset, Vertex, index4);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
		else if (p[0].value <= minisolevel && p[1].value > minisolevel && p[2].value > minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index1, 1);
			getPointinsurface(Pindex[1], pointset, Vertex, index2);
			getPointinsurface(Pindex[2], pointset, Vertex, index3);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index4, 1);

			addFace(Face, index1, index2, index3);
			addFace(Face, index1, index3, index4);
		}
		else if (p[0].value > minisolevel && p[1].value <= minisolevel && p[2].value <= minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, Vertex, index1);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index2, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index3, 1);

			addFace(Face, index1, index2, index3);
		}
		else if (p[0].value <= minisolevel && p[1].value > minisolevel && p[2].value <= minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, Vertex, index1, 1);
			getPointinsurface(Pindex[1], pointset, Vertex, index2);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index3, 1);

			addFace(Face, index1, index2, index3);
		}
		else if (p[0].value <= minisolevel && p[1].value <= minisolevel && p[2].value > minisolevel)
		{
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, Vertex, index1, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, Vertex, index2, 1);
			getPointinsurface(Pindex[2], pointset, Vertex, index3);

			addFace(Face, index1, index2, index3);
		}
		else if (p[0].value <= minisolevel && p[1].value <= minisolevel && p[2].value <= minisolevel)
		{
		}
	}
}

//提取大于minisolevel的三角网格
static void extractTriangleFromSurfaceMin(const Isosurface& surface, Point3D p[3], float isolevel, vector<Vector3D>& Vertex, vector<Triangle>& Face)
{
	Vector3D intersectionP1, intersectionP2;
	Triangle face;
	if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//0,1,2 outer
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//0,1 outer 2 inner
		calInterVert(surface, p[1], p[2], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(intersectionP1);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);

		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP1);
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//0,2 outer 1 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[1], p[2], isolevel, intersectionP2);
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP1);
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);

		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP2);
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0 outer 1,2 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(Vector3D(p[0].x, p[0].y, p[0].z));
		Vertex.push_back(intersectionP1);
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//1,2 outer 0 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);

		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//1 outer 0,2 inner
		calInterVert(surface, p[0], p[1], isolevel, intersectionP1);
		calInterVert(surface, p[1], p[2], isolevel, intersectionP2);
		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[1].x, p[1].y, p[1].z));
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//2 outer 0,1 inner
		calInterVert(surface, p[1], p[2], isolevel, intersectionP1);
		calInterVert(surface, p[0], p[2], isolevel, intersectionP2);
		Vertex.push_back(intersectionP1);
		Vertex.push_back(Vector3D(p[2].x, p[2].y, p[2].z));
		Vertex.push_back(intersectionP2);
		face.f1 = Vertex.size() - 3;
		face.f2 = Vertex.size() - 2;
		face.f3 = Vertex.size() - 1;
		Face.push_back(face);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0,1,2 inner
	}

}

//首先提取小于maxisolevel的三角网格，再输入extractTriangleFromSurfaceMin中去
static void extractTriangleFromSurfaceMax(const Isosurface& surface, const Point3D p[3], float maxisolevel, float minisolevel, vector<Vector3D>& Vertex, vector<Triangle>& Face)
{
	Point3D intersectionP1, intersectionP2;
	Point3D Tri[3];
	Triangle face;
	if (p[0].value >= maxisolevel && p[1].value >= maxisolevel && p[2].value >= maxisolevel)
	{//0,1,2 outer
	}
	else if (p[0].value >= maxisolevel && p[1].value >= maxisolevel && p[2].value < maxisolevel)
	{//0,1 outer 2 inner
		calInterPoint(surface, p[0], p[2], maxisolevel, intersectionP1);
		calInterPoint(surface, p[1], p[2], maxisolevel, intersectionP2);
		Tri[0] = intersectionP1;
		Tri[1] = intersectionP2;
		Tri[2] = p[2];
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);
	}
	else if (p[0].value >= maxisolevel && p[1].value < maxisolevel && p[2].value >= maxisolevel)
	{//0,2 outer 1 inner
		calInterPoint(surface, p[0], p[1], maxisolevel, intersectionP1);
		calInterPoint(surface, p[1], p[2], maxisolevel, intersectionP2);
		Tri[0] = intersectionP1;
		Tri[1] = p[1];
		Tri[2] = intersectionP2;
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);
	}
	else if (p[0].value < maxisolevel && p[1].value >= maxisolevel && p[2].value >= maxisolevel)
	{//1,2 outer 0 inner
		calInterPoint(surface, p[0], p[1], maxisolevel, intersectionP1);
		calInterPoint(surface, p[0], p[2], maxisolevel, intersectionP2);
		Tri[0] = p[0];
		Tri[1] = intersectionP1;
		Tri[2] = intersectionP2;
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);
	}
	else if (p[0].value >= maxisolevel && p[1].value < maxisolevel && p[2].value < maxisolevel)
	{//0 outer 1,2 inner
		calInterPoint(surface, p[0], p[1], maxisolevel, intersectionP1);
		calInterPoint(surface, p[0], p[2], maxisolevel, intersectionP2);
		Tri[0] = intersectionP1;
		Tri[1] = p[1];
		Tri[2] = p[2];
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);

		Tri[0] = intersectionP1;
		Tri[1] = p[2];
		Tri[2] = intersectionP2;
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);
	}
	else if (p[0].value < maxisolevel && p[1].value >= maxisolevel && p[2].value < maxisolevel)
	{//1 outer 0,2 inner
		calInterPoint(surface, p[0], p[1], maxisolevel, intersectionP1);
		calInterPoint(surface, p[1], p[2], maxisolevel, intersectionP2);
		Tri[0] = p[0];
		Tri[1] = intersectionP1;
		Tri[2] = intersectionP2;
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);

		Tri[0] = p[0];
		Tri[1] = intersectionP2;
		Tri[2] = p[2];
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);
	}
	else if (p[0].value < maxisolevel && p[1].value < maxisolevel && p[2].value >= maxisolevel)
	{//2 outer 0,1 inner
		calInterPoint(surface, p[1], p[2], maxisolevel, intersectionP1);
		calInterPoint(surface, p[0], p[2], maxisolevel, intersectionP2);
		Tri[0] = p[0];
		Tri[1] = p[1];
		Tri[2] = intersectionP1;
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);

		Tri[0] = p[0];
		Tri[1] = intersectionP1;
		Tri[2] = intersectionP2;
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);
	}
	else if (p[0].value < maxisolevel && p[1].value < maxisolevel && p[2].value < maxisolevel)
	{//0,1,2 inner
		Tri[0] = p[0];
		Tri[1] = p[1];
		Tri[2] = p[2];
		extractTriangleFromSurfaceMin(surface, Tri, minisolevel, Vertex, Face);
	}

}

#pragma endregion


void save_obj(string filename, vector<Vector3D>& Vertex, vector<Triangle>& Face,int type)
{
	ofstream fout(filename);
	for (int i = 0; i < Vertex.size(); i++)
	{
		if (!Vertex[i].isDuplicated)
		{
			fout << "v " << Vertex[i].x << " " << Vertex[i].y << " " << Vertex[i].z << "\n";
		}
	}
	if (type == 0)
	{
		for (int i = 0; i < Face.size(); i++)
		{
			if (!Face[i].isDegenerated)
			{
				fout << "f " << Face[i].f1 + 1 << " " << Face[i].f2 + 1 << " " << Face[i].f3 + 1 << "\n";
			}
		}
	}
	else
	{
		for (int i = 0; i < Face.size(); i++)
		{
			if (!Face[i].isDegenerated)
			{
				fout << "f " << Face[i].f1 + 1 << " " << Face[i].f3 + 1 << " " << Face[i].f2 + 1 << "\n";
			}
		}
	}
	fout.close();
}

void save_stl(string filename, vector<Vector3D>& Vertex, vector<Triangle>& Face, vector<Vector3D> FaceNormal, int type)
{
	ofstream fout(filename);

	if (type == 0)
	{
		for (int i = 0; i < Face.size(); i++)
		{
			fout << "  facet normal " << FaceNormal[i].x << " " << FaceNormal[i].y << " " << FaceNormal[i].z << "\n";
			fout << "    outer loop" << "\n";
			fout << "      vertex  " << Vertex[Face[i].f1].x << " " << Vertex[Face[i].f1].y << " " << Vertex[Face[i].f1].z << "\n";
			fout << "      vertex  " << Vertex[Face[i].f2].x << " " << Vertex[Face[i].f2].y << " " << Vertex[Face[i].f2].z << "\n";
			fout << "      vertex  " << Vertex[Face[i].f3].x << " " << Vertex[Face[i].f3].y << " " << Vertex[Face[i].f3].z << "\n";
			fout << "    endloop" << "\n";
			fout << "  endfacet" << "\n";
		}
	}
	else
	{
		for (int i = 0; i < Face.size(); i++)
		{
			fout << "  facet normal " << FaceNormal[i].x << " " << FaceNormal[i].y << " " << FaceNormal[i].z << "\n";
			fout << "    outer loop" << "\n";
			fout << "      vertex  " << Vertex[Face[i].f1].x << " " << Vertex[Face[i].f1].y << " " << Vertex[Face[i].f1].z << "\n";
			fout << "      vertex  " << Vertex[Face[i].f3].x << " " << Vertex[Face[i].f3].y << " " << Vertex[Face[i].f3].z << "\n";
			fout << "      vertex  " << Vertex[Face[i].f2].x << " " << Vertex[Face[i].f2].y << " " << Vertex[Face[i].f2].z << "\n";
			fout << "    endloop" << "\n";
			fout << "  endfacet" << "\n";
		}
	}
	fout.close();
}


void constructEdge(Array3D<Point3D>& PointGrid, size_t resolution, vector<edge>& Edge)
{
	for (size_t k = 0; k <= resolution; ++k)
	{
		for (size_t j = 0; j <= resolution; ++j)
		{
			for (size_t i = 0; i < resolution; ++i)
			{
				edge e;
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.index[0] = -1; e.index[1] = -1;
				Edge.push_back(e);
			}
		}
	}//index: 0 - res*(res+1)*(res+1)-1
	for (size_t k = 0; k <= resolution; ++k)
	{
		for (size_t i = 0; i <= resolution; ++i)
		{
			for (size_t j = 0; j < resolution; ++j)
			{
				edge e;
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k;
				e.index[0] = -1; e.index[1] = -1;
				Edge.push_back(e);
			}
		}
	}//index: res*(res+1)*(res+1) - 2*res*(res+1)*(res+1)-1
	for (size_t i = 0; i <= resolution; ++i)
	{
		for (size_t j = 0; j <= resolution; ++j)
		{
			for (size_t k = 0; k < resolution; ++k)
			{
				edge e;
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1;
				e.index[0] = -1; e.index[1] = -1;
				Edge.push_back(e);
			}
		}
	}//index: 2*res*(res+1)*(res+1) - 3*res*(res+1)*(res+1)-1
	for (size_t k = 0; k <= resolution; ++k)//ij
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			for (size_t i = 0; i < resolution; ++i)
			{
				edge e;
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k;
				e.index[0] = -1; e.index[1] = -1;
				Edge.push_back(e);
			}
		}
	}//index: 3*res*(res+1)*(res+1) - 3*res*(res+1)*(res+1)+res*res*(res+1)-1
	for (size_t j = 0; j <= resolution; ++j)//ik
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			for (size_t i = 0; i < resolution; ++i)
			{
				edge e;
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k+1;
				e.id2 = (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.index[0] = -1; e.index[1] = -1;
				Edge.push_back(e);
			}
		}
	}//index: 3*res*(res+1)*(res+1)+res*res*(res+1) - 3*res*(res+1)*(res+1)+2*res*res*(res+1)-1
	for (size_t i = 0; i <= resolution; ++i)//jk
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			for (size_t j = 0; j < resolution; ++j)
			{
				edge e;
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1;
				e.index[0] = -1; e.index[1] = -1;
				Edge.push_back(e);
			}
		}
	}//index: 3*res*(res+1)*(res+1)+2*res*res*(res+1) - 3*res*(res+1)*(res+1)+3*res*res*(res+1)-1
	for (size_t k = 0; k < resolution; ++k)
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			for (size_t i = 0; i < resolution; ++i)
			{
				edge e;
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + (k + 1);
				e.index[0] = -1; e.index[1] = -1;
				Edge.push_back(e);
			}
		}
	}//index: 3*res*(res+1)*(res+1)+3*res*res*(res+1) - 3*res*(res+1)*(res+1)+3*res*res*(res+1)+res*res*res-1
}

void extract_isosurface(const Isosurface& surface, float isolevel, size_t resolution, int type, vector<Vector3D>& Vertex, vector<Triangle>& Face, vector<edge>& Edge, vector<Point3D>& pointset, Array3D<Point3D>& PointGrid,int sur_id) // resolution indicates # of cubes
{
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			for (size_t k = 0; k < resolution; ++k)
			{
				/*
				 Coordinates:
				 z
				 |
				 |___ y
				 /
				 /
				 x

				 Cube layout:
	    		   4-------7
		    	  /|      /|
				 / |     / |
				 5-------6  |
				 |  0----|--3
				 | /     | /
				 |/      |/
				 1-------2

				 Tetrahedrons are:
				 0, 7, 3, 2
				 0, 7, 2, 6
				 0, 4, 6, 7
				 0, 6, 1, 2
				 0, 6, 1, 4
				 5, 6, 1, 4
				 */

				const Point3D v[8] = {
					{ PointGrid.get(i, j, k) },
					{ PointGrid.get(i + 1, j, k) },
					{ PointGrid.get(i + 1, j + 1, k) },
					{ PointGrid.get(i, j + 1, k) },
					{ PointGrid.get(i, j, k + 1) },
					{ PointGrid.get(i + 1, j, k + 1) },
					{ PointGrid.get(i + 1, j + 1, k + 1) },
					{ PointGrid.get(i, j + 1, k + 1) }
				};
				const Point3D tetrahedra[6][4] = {
					{ v[0], v[7], v[3], v[2] },
					{ v[0], v[7], v[2], v[6] },
					{ v[0], v[4], v[7], v[6] },
					{ v[0], v[1], v[6], v[2] },
					{ v[0], v[4], v[6], v[1] },
					{ v[5], v[1], v[6], v[4] }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tet_edge[6][6] = {
					{ 3 * len1 + 2 * len2 + i*resolution*resolution + k*resolution + j, //v[0]v[7]
					len1 + k*(resolution + 1)*resolution + i*resolution + j, //v[0]v[3]
					3 * len1 + k*resolution*resolution + j*resolution + i, //v[0]v[2]
					2 * len1 + i*(resolution + 1)*resolution + (j + 1)*resolution + k, //v[7]v[3]
					3 * len1 + len2 + (j + 1)*resolution*resolution + k*resolution + i,//v[7]v[2]
					k*(resolution + 1)*resolution + (j + 1)*resolution + i//v[3]v[2]
					},
					{ 3 * len1 + 2 * len2 + i*resolution*resolution + k*resolution + j,//v[0]v[7]
					3 * len1 + k*resolution*resolution + j*resolution + i,//v[0]v[2]
					3 * len1 + 3 * len2 + k*resolution*resolution + j*resolution + i,//v[0]v[6]
					3 * len1 + len2 + (j + 1)*resolution*resolution + k*resolution + i,//v[7]v[2]
					(k + 1)*(resolution + 1)*resolution + (j + 1)*resolution + i,//v[7]v[6]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + (j + 1)*resolution + k //v[2]v[6]
					},
					{ 2 * len1 + i*(resolution + 1)*resolution + j*resolution + k, //v[0]v[4]
					3 * len1 + 2 * len2 + i*resolution*resolution + k*resolution + j,//v[0]v[7]
					3 * len1 + 3 * len2 + k*resolution*resolution + j*resolution + i,//v[0]v[6]
					len1 + (k + 1)*(resolution + 1)*resolution + i*resolution + j, //v[4]v[7]
					3 * len1 + (k + 1)*resolution*resolution + j*resolution + i, //v[4]v[6]
					(k + 1)*(resolution + 1)*resolution + (j + 1)*resolution + i//v[7]v[6]
					},
					{ k*(resolution + 1)*resolution + j*resolution + i,//v[0]v[1]
					3 * len1 + 3 * len2 + k*resolution*resolution + j*resolution + i,//v[0]v[6]
					3 * len1 + k*resolution*resolution + j*resolution + i,//v[0]v[2]
					3 * len1 + 2 * len2 + (i + 1)*resolution*resolution + k*resolution + j, //v[1]v[6]
					len1 + k*(resolution + 1)*resolution + (i + 1)*resolution + j, //v[1]v[2]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + (j + 1)*resolution + k //v[6]v[2]
					},
					{ 2 * len1 + i*(resolution + 1)*resolution + j*resolution + k, //v[0]v[4]
					3 * len1 + 3 * len2 + k*resolution*resolution + j*resolution + i,//v[0]v[6]
					k*(resolution + 1)*resolution + j*resolution + i,//v[0]v[1]
					3 * len1 + (k + 1)*resolution*resolution + j*resolution + i, //v[4]v[6]
					3 * len1 + len2 + j*resolution*resolution + k*resolution + i,//v[4]v[1]
					3 * len1 + 2 * len2 + (i + 1)*resolution*resolution + k*resolution + j //v[1]v[6]
					},
					{ 2 * len1 + (i + 1)*(resolution + 1)*resolution + j*resolution + k, //v[5]v[1]
					len1 + (k + 1)*(resolution + 1)*resolution + (i + 1)*resolution + j, //v[5]v[6]
					(k + 1)*(resolution + 1)*resolution + j*resolution + i,//v[5]v[4]
					3 * len1 + 2 * len2 + (i + 1)*resolution*resolution + k*resolution + j, //v[1]v[6]
					3 * len1 + len2 + j*resolution*resolution + k*resolution + i,//v[1]v[4]
					3 * len1 + (k + 1)*resolution*resolution + j*resolution + i, //v[6]v[4]
					}
				};
				for (int t = 0; t < 6; ++t)
				{
					extractTriangleFromTetrahedron(surface, tetrahedra[t], tet_edge[t], isolevel, pointset, Edge, Vertex, Face, sur_id);
				}
			}
		}
	}
	////提取孔洞模型，面片法向翻转
	//if (type == 1)
	//{
	//	int index_temp;
	//	for (int i = 0; i < Face.size(); i++)
	//	{
	//		index_temp = Face[i].f2;
	//		Face[i].f2 = Face[i].f3;
	//		Face[i].f3 = index_temp;
	//	}
	//}
}

void close_isosurface_network(const Isosurface& surface, float isolevel, size_t resolution, int type, vector<Vector3D>& Vertex, vector<Triangle>& Face, vector<edge>& Edge, vector<Point3D>& pointset, Array3D<Point3D>& PointGrid,int sur_id)
{
	//提取六张表面生成封闭模型
	float x1, y1, z1, x2, y2, z2;
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			const Point3D v[8] = {
				{ PointGrid.get(i, j, 0) },
				{ PointGrid.get(i + 1, j, 0) },
				{ PointGrid.get(i + 1, j + 1, 0) },
				{ PointGrid.get(i, j + 1, 0) },
				{ PointGrid.get(i, j, resolution) },
				{ PointGrid.get(i + 1, j, resolution) },
				{ PointGrid.get(i + 1, j + 1, resolution) },
				{ PointGrid.get(i, j + 1, resolution) },
			};
			if (type == 0)
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[2], v[1] },
					{ v[0], v[3], v[2] },
					{ v[4], v[5], v[6] },
					{ v[4], v[6], v[7] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ 3 * len1 + j*resolution + i,//v[0]v[2]
					j*resolution + i,//v[0]v[1]
					len1 + (i + 1)*resolution + j//v[2]v[1]
					},
					{ len1 + i*resolution + j, //v[0]v[3]
					3 * len1 + j*resolution + i,//v[0]v[2]
					(j + 1)*resolution + i//v[3]v[2]
					},
					{ resolution*(resolution + 1)*resolution + j*resolution + i,//v[4]v[5]
					3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					len1 + resolution*(resolution + 1)*resolution + (i + 1)*resolution + j//v[5]v[6]
					},
					{ 3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					len1 + resolution*(resolution + 1)*resolution + i*resolution + j, //v[4]v[7]
					resolution*(resolution + 1)*resolution + (j + 1)*resolution + i//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleInner2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, Vertex, Face, sur_id);
				}
			}
			else
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[1], v[2] },
					{ v[0], v[2], v[3] },
					{ v[4], v[6], v[5] },
					{ v[4], v[7], v[6] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ j*resolution + i,//v[0]v[1]
					3 * len1 + j*resolution + i,//v[0]v[2]
					len1 + (i + 1)*resolution + j//v[2]v[1]
					},
					{ 3 * len1 + j*resolution + i,//v[0]v[2]
					len1 + i*resolution + j, //v[0]v[3]
					(j + 1)*resolution + i//v[3]v[2]
					},
					{ 3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					resolution*(resolution + 1)*resolution + j*resolution + i,//v[4]v[5]
					len1 + resolution*(resolution + 1)*resolution + (i + 1)*resolution + j//v[5]v[6]
					},
					{ len1 + resolution*(resolution + 1)*resolution + i*resolution + j, //v[4]v[7]
					3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					resolution*(resolution + 1)*resolution + (j + 1)*resolution + i//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleOuter2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, Vertex, Face, sur_id);
				}
			}
		}
	}
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			const Point3D v[8] = {
				{ PointGrid.get(i, 0, k) },
				{ PointGrid.get(i + 1, 0, k) },
				{ PointGrid.get(i + 1, 0, k + 1) },
				{ PointGrid.get(i, 0, k + 1) },
				{ PointGrid.get(i, resolution, k) },
				{ PointGrid.get(i + 1, resolution, k) },
				{ PointGrid.get(i + 1, resolution, k + 1) },
				{ PointGrid.get(i, resolution, k + 1) },
			};
			if (type == 0)
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[1], v[3] },
					{ v[1], v[2], v[3] },
					{ v[4], v[7], v[5] },
					{ v[5], v[7], v[6] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1 },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + k + 1, i*(resolution + 1)*(resolution + 1) + k + 1 },
					{ i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ k*(resolution + 1)*resolution + i,//v[0]v[1]
					2 * len1 + i*(resolution + 1)*resolution + k,//v[0]v[3]
					3 * len1 + len2 + k*resolution + i//v[1]v[3]
					},
					{ 2 * len1 + (i + 1)*(resolution + 1)*resolution + k, //v[1]v[2]
					3 * len1 + len2 + k*resolution + i,//v[1]v[3]
					(k + 1)*(resolution + 1)*resolution + i//v[2]v[3]
					},
					{ 2 * len1 + i*(resolution + 1)*resolution + resolution*resolution + k,//v[4]v[7]
					k*(resolution + 1)*resolution + resolution*resolution + i,//v[4]v[5]
					3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i//v[7]v[5]
					},
					{ 3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i,//v[5]v[7]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + resolution*resolution + k, //v[5]v[6]
					(k + 1)*(resolution + 1)*resolution + resolution*resolution + i//v[7]v[6]
					}
				};
				for (int t = 0; t < 4; ++t)
				//for (int t = 0; t < 2; ++t)
				{
					extractTriangleInner2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, Vertex, Face, sur_id);
				}
			}
			else
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[3], v[1] },
					{ v[1], v[3], v[2] },
					{ v[4], v[5], v[7] },
					{ v[5], v[6], v[7] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + k },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + k + 1 },
					{ i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ 2 * len1 + i*(resolution + 1)*resolution + k,//v[0]v[3]
					k*(resolution + 1)*resolution + i,//v[0]v[1]
					3 * len1 + len2 + k*resolution + i//v[1]v[3]
					},
					{ 3 * len1 + len2 + k*resolution + i,//v[1]v[3]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + k, //v[1]v[2]
					(k + 1)*(resolution + 1)*resolution + i//v[2]v[3]
					},
					{ k*(resolution + 1)*resolution + resolution*resolution + i,//v[4]v[5]
					2 * len1 + i*(resolution + 1)*resolution + resolution*resolution + k,//v[4]v[7]
					3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i//v[7]v[5]
					},
					{ 2 * len1 + (i + 1)*(resolution + 1)*resolution + resolution*resolution + k, //v[5]v[6]
					3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i,//v[5]v[7]
					(k + 1)*(resolution + 1)*resolution + resolution*resolution + i//v[7]v[6]
					}
				};
				for (int t = 0; t < 4; ++t)
				//for (int t = 0; t < 2; ++t)
				{
					extractTriangleOuter2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, Vertex, Face, sur_id);
				}
			}
		}
	}
	for (size_t j = 0; j < resolution; ++j)
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			const Point3D v[8] = {
				{ PointGrid.get(0, j, k) },
				{ PointGrid.get(0, j + 1, k) },
				{ PointGrid.get(0, j + 1, k + 1) },
				{ PointGrid.get(0, j, k + 1) },
				{ PointGrid.get(resolution, j, k) },
				{ PointGrid.get(resolution, j + 1, k) },
				{ PointGrid.get(resolution, j + 1, k + 1) },
				{ PointGrid.get(resolution, j, k + 1) },
			};
			if (type == 0)
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[2], v[1] },
					{ v[0], v[3], v[2] },
					{ v[4], v[5], v[6] },
					{ v[4], v[6], v[7] },
				};
				const int Pindex[4][3] = {
					{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1, (j + 1)*(resolution + 1) + k },
					{ j*(resolution + 1) + k, j*(resolution + 1) + k + 1, (j + 1)*(resolution + 1) + k + 1 },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1 },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ 3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					len1 + k*(resolution + 1)*resolution + j,//v[0]v[1]
					2 * len1 + (j + 1)*resolution + k//v[2]v[1]
					},
					{ 2 * len1 + j*resolution + k, //v[0]v[3]
					3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					len1 + (k + 1)*(resolution + 1)*resolution + j//v[3]v[2]
					},
					{ len1 + k*(resolution + 1)*resolution + resolution*resolution + j,//v[4]v[5]
					3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					2 * len1 + resolution*(resolution + 1)*resolution + (j + 1)*resolution + k//v[5]v[6]
					},
					{ 3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					2 * len1 + resolution*(resolution + 1)*resolution + j*resolution + k, //v[4]v[7]
					len1 + (k + 1)*(resolution + 1)*resolution + resolution*resolution + j//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleInner2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, Vertex, Face, sur_id);
				}
			}
			else
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[1], v[2] },
					{ v[0], v[2], v[3] },
					{ v[4], v[6], v[5] },
					{ v[4], v[7], v[6] },
				};
				const int Pindex[4][3] = {
					{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1 },
					{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1, j*(resolution + 1) + k + 1 },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ len1 + k*(resolution + 1)*resolution + j,//v[0]v[1]
					3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					2 * len1 + (j + 1)*resolution + k//v[2]v[1]
					},
					{ 3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					2 * len1 + j*resolution + k, //v[0]v[3]
					len1 + (k + 1)*(resolution + 1)*resolution + j//v[3]v[2]
					},
					{ 3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					len1 + k*(resolution + 1)*resolution + resolution*resolution + j,//v[4]v[5]
					2 * len1 + resolution*(resolution + 1)*resolution + (j + 1)*resolution + k//v[5]v[6]
					},
					{ 2 * len1 + resolution*(resolution + 1)*resolution + j*resolution + k, //v[4]v[7]
					3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					len1 + (k + 1)*(resolution + 1)*resolution + resolution*resolution + j//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleOuter2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, Vertex, Face, sur_id);
				}
			}
		}
	}
}


//根据边索引来求取等值点，判断当前边是否已经插入等值点，插入后直接调用已知点及其对应点索引
void close_isosurface_sheet(const Isosurface& surface, float maxisolevel, float minisolevel, size_t resolution, vector<Vector3D>& Vertex, vector<Triangle>& Face, vector<edge>& Edge, vector<Point3D>& pointset, Array3D<Point3D>& PointGrid)
{
	//提取六张表面生成封闭模型
	float x1, y1, z1, x2, y2, z2;
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			const Point3D v[8] = {
				{ PointGrid.get(i, j, 0) },
				{ PointGrid.get(i + 1, j, 0) },
				{ PointGrid.get(i + 1, j + 1, 0) },
				{ PointGrid.get(i, j + 1, 0) },
				{ PointGrid.get(i, j, resolution) },
				{ PointGrid.get(i + 1, j, resolution) },
				{ PointGrid.get(i + 1, j + 1, resolution) },
				{ PointGrid.get(i, j + 1, resolution) },
			};
			const Point3D triangle[4][3] = {
				{ v[0], v[2], v[1] },
				{ v[0], v[3], v[2] },
				{ v[4], v[5], v[6] },
				{ v[4], v[6], v[7] },
			};
			const int Pindex[4][3] = {
				{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) },
				{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
				{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution },
				{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution }
			};
			int len1 = resolution*(resolution + 1)*(resolution + 1);
			int len2 = resolution*resolution*(resolution + 1);
			const int tri_edge[4][3] = {
				{ 3 * len1 + j*resolution + i,//v[0]v[2]
				j*resolution + i,//v[0]v[1]
				len1 + (i + 1)*resolution + j//v[2]v[1]
				},
				{ len1 + i*resolution + j, //v[0]v[3]
				3 * len1 + j*resolution + i,//v[0]v[2]
				(j + 1)*resolution + i//v[3]v[2]
				},
				{ resolution*(resolution + 1)*resolution + j*resolution + i,//v[4]v[5]
				3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
				len1 + resolution*(resolution + 1)*resolution + (i + 1)*resolution + j//v[5]v[6]
				},
				{ 3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
				len1 + resolution*(resolution + 1)*resolution + i*resolution + j, //v[4]v[7]
				resolution*(resolution + 1)*resolution + (j + 1)*resolution + i//v[6]v[7]
				}
			};
			for (int t = 0; t < 4; ++t)
			{
				extractTriangleFromSurface(surface, triangle[t], tri_edge[t], Pindex[t], maxisolevel, minisolevel, pointset, Edge, Vertex, Face);
			}
		}
	}
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			const Point3D v[8] = {
				{ PointGrid.get(i, 0, k) },
				{ PointGrid.get(i + 1, 0, k) },
				{ PointGrid.get(i + 1, 0, k + 1) },
				{ PointGrid.get(i, 0, k + 1) },
				{ PointGrid.get(i, resolution, k) },
				{ PointGrid.get(i + 1, resolution, k) },
				{ PointGrid.get(i + 1, resolution, k + 1) },
				{ PointGrid.get(i, resolution, k + 1) },
			};
			const Point3D triangle[4][3] = {
				{ v[0], v[1], v[3] },
				{ v[1], v[2], v[3] },
				{ v[4], v[7], v[5] },
				{ v[5], v[7], v[6] },
			};
			const int Pindex[4][3] = {
				{ i*(resolution + 1)*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1 },
				{ (i + 1)*(resolution + 1)*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + k + 1, i*(resolution + 1)*(resolution + 1) + k + 1 },
				{ i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k },
				{ (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 }
			};
			int len1 = resolution*(resolution + 1)*(resolution + 1);
			int len2 = resolution*resolution*(resolution + 1);
			const int tri_edge[4][3] = {
				{ k*(resolution + 1)*resolution + i,//v[0]v[1]
				2 * len1 + i*(resolution + 1)*resolution + k,//v[0]v[3]
				3 * len1 + len2 + k*resolution + i//v[1]v[3]
				},
				{ 2 * len1 + (i + 1)*(resolution + 1)*resolution + k, //v[1]v[2]
				3 * len1 + len2 + k*resolution + i,//v[1]v[3]
				(k + 1)*(resolution + 1)*resolution + i//v[2]v[3]
				},
				{ 2 * len1 + i*(resolution + 1)*resolution + resolution*resolution + k,//v[4]v[7]
				k*(resolution + 1)*resolution + resolution*resolution + i,//v[4]v[5]
				3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i//v[7]v[5]
				},
				{ 3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i,//v[5]v[7]
				2 * len1 + (i + 1)*(resolution + 1)*resolution + resolution*resolution + k, //v[5]v[6]
				(k + 1)*(resolution + 1)*resolution + resolution*resolution + i//v[7]v[6]
				}
			};
			for (int t = 0; t < 4; ++t)
			//for (int t = 0; t < 2; ++t)
			{
				extractTriangleFromSurface(surface, triangle[t], tri_edge[t], Pindex[t], maxisolevel, minisolevel, pointset, Edge, Vertex, Face);
			}
		}
	}
	for (size_t j = 0; j < resolution; ++j)
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			const Point3D v[8] = {
				{ PointGrid.get(0, j, k) },
				{ PointGrid.get(0, j + 1, k) },
				{ PointGrid.get(0, j + 1, k + 1) },
				{ PointGrid.get(0, j, k + 1) },
				{ PointGrid.get(resolution, j, k) },
				{ PointGrid.get(resolution, j + 1, k) },
				{ PointGrid.get(resolution, j + 1, k + 1) },
				{ PointGrid.get(resolution, j, k + 1) },
			};
			const Point3D triangle[4][3] = {
				{ v[0], v[2], v[1] },
				{ v[0], v[3], v[2] },
				{ v[4], v[5], v[6] },
				{ v[4], v[6], v[7] },
			};
			const int Pindex[4][3] = {
				{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1, (j + 1)*(resolution + 1) + k },
				{ j*(resolution + 1) + k, j*(resolution + 1) + k + 1, (j + 1)*(resolution + 1) + k + 1 },
				{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1 },
				{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1 }
			};
			int len1 = resolution*(resolution + 1)*(resolution + 1);
			int len2 = resolution*resolution*(resolution + 1);
			const int tri_edge[4][3] = {
				{ 3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
				len1 + k*(resolution + 1)*resolution + j,//v[0]v[1]
				2 * len1 + (j + 1)*resolution + k//v[2]v[1]
				},
				{ 2 * len1 + j*resolution + k, //v[0]v[3]
				3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
				len1 + (k + 1)*(resolution + 1)*resolution + j//v[3]v[2]
				},
				{ len1 + k*(resolution + 1)*resolution + resolution*resolution + j,//v[4]v[5]
				3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
				2 * len1 + resolution*(resolution + 1)*resolution + (j + 1)*resolution + k//v[5]v[6]
				},
				{ 3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
				2 * len1 + resolution*(resolution + 1)*resolution + j*resolution + k, //v[4]v[7]
				len1 + (k + 1)*(resolution + 1)*resolution + resolution*resolution + j//v[6]v[7]
				}
			};
			for (int t = 0; t < 4; ++t)
			{
				extractTriangleFromSurface(surface, triangle[t], tri_edge[t], Pindex[t], maxisolevel, minisolevel, pointset, Edge, Vertex, Face);
			}
		}
	}
}

void solid_construction_Bspline(const Isosurface& surface, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, float isolevel, size_t resolution)
{
	int type = 0;
	vector<Vector3D> Vertex;
	vector<Triangle> Face;
	vector<edge> Edge;
	vector<Point3D> pointset;
	size_t pointRes = resolution + 1; // indicates the # of points per side
	Array3D<Point3D> PointGrid(pointRes, pointRes, pointRes);

	//double ***CValue;
	//ifstream cvaluefile("E:\\project\\huchuanfeng\\SemanticControlGrid\\src\\GeometryQTGUI\\PorosityDistribution_tooth.txt");
	//CValue = new double **[pointRes];
	//for (int i = 0; i < pointRes; i++)
	//{
	//	CValue[i] = new double *[pointRes];
	//	for (int j = 0; j < pointRes; j++)
	//	{
	//		CValue[i][j] = new double[pointRes];
	//	}
	//}
	//for (int i = 0; i < pointRes; i++)
	//{
	//	for (int j = 0; j < pointRes; j++)
	//	{
	//		for (int k = 0; k < pointRes; k++)
	//		{
	//			cvaluefile >> CValue[i][j][k];
	//		}
	//	}
	//}

	float xrange = xMax - xMin;
	float yrange = yMax - yMin;
	float zrange = zMax - zMin;
	
	for (size_t i = 0; i <= resolution; ++i) {
		float x = (float)i / resolution * xrange + xMin;
		for (size_t j = 0; j <= resolution; ++j) {
			float y = (float)j / resolution * yrange + yMin;
			for (size_t k = 0; k <= resolution; ++k) {
				float z = (float)k / resolution * zrange + zMin;
				float value = surface.valueAt(8*x, 6*y, 6*z);
				PointGrid.set(i, j, k, Point3D(x, y, z, value));
				pointset.push_back(Point3D(x, y, z, value));
			}
		}
	}
	constructEdge(PointGrid, resolution, Edge);

	//clock_t time_start1, time_end1;//计算时间
	//time_start1 = clock();
	extract_isosurface(surface, isolevel, resolution, type, Vertex, Face, Edge, pointset, PointGrid, 0);
	close_isosurface_network(surface, isolevel, resolution, type, Vertex, Face, Edge, pointset, PointGrid,0);
	//time_end1 = clock();
	//double duration1 = (static_cast<double>(time_end1 - time_start1));//返回毫秒

	//extract_isosurface(surface, isolevel - 0.3, resolution, type, Vertex, Face, Edge, pointset, PointGrid, 1);
	//////only for constructing IWP_sheet_solid 
	////int index_temp;
	////for (int i = 0; i < Face.size(); i++)
	////{
	////	index_temp = Face[i].f2;
	////	Face[i].f2 = Face[i].f3;
	////	Face[i].f3 = index_temp;
	////}
	//close_isosurface_sheet(surface, isolevel, isolevel - 0.3, resolution, Vertex, Face, Edge, pointset, PointGrid);

	Vector3D normal;
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < Face.size(); i++)
	{
		normal = surface.gradientAt(Vertex[Face[i].f1].x, Vertex[Face[i].f1].y, Vertex[Face[i].f1].z);
		glNormal3f(normal.x, normal.y, normal.z);
		glVertex3f(Vertex[Face[i].f1].x, Vertex[Face[i].f1].y, Vertex[Face[i].f1].z);
		normal = surface.gradientAt(Vertex[Face[i].f2].x, Vertex[Face[i].f2].y, Vertex[Face[i].f2].z);
		glNormal3f(normal.x, normal.y, normal.z);
		glVertex3f(Vertex[Face[i].f2].x, Vertex[Face[i].f2].y, Vertex[Face[i].f2].z);
		normal = surface.gradientAt(Vertex[Face[i].f3].x, Vertex[Face[i].f3].y, Vertex[Face[i].f3].z);
		glNormal3f(normal.x, normal.y, normal.z);
		glVertex3f(Vertex[Face[i].f3].x, Vertex[Face[i].f3].y, Vertex[Face[i].f3].z);
	}
	glEnd();


	save_obj("C:\\Users\\Admin\\Desktop\\Punit.obj",Vertex, Face, type);
}

void solid_construction_Bspline_save(const Isosurface& surface, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, float isolevel, size_t resolution,string s)
{
	int type = 1;
	vector<Vector3D> Vertex;
	vector<Triangle> Face;
	vector<edge> Edge;
	vector<Point3D> pointset;
	size_t pointRes = resolution + 1; // indicates the # of points per side
	Array3D<Point3D> PointGrid(pointRes, pointRes, pointRes);

	float xrange = xMax - xMin;
	float yrange = yMax - yMin;
	float zrange = zMax - zMin;

	for (size_t i = 0; i <= resolution; ++i) {
		float x = (float)i / resolution * xrange + xMin;
		for (size_t j = 0; j <= resolution; ++j) {
			float y = (float)j / resolution * yrange + yMin;
			for (size_t k = 0; k <= resolution; ++k) {
				float z = (float)k / resolution * zrange + zMin;
				float value = surface.valueAt(x, y, z);
				PointGrid.set(i, j, k, Point3D(x, y, z, value));
				pointset.push_back(Point3D(x, y, z, value));
			}
		}
	}
	constructEdge(PointGrid, resolution, Edge);

	extract_isosurface(surface, isolevel, resolution, type, Vertex, Face, Edge, pointset, PointGrid, 0);
	close_isosurface_network(surface, isolevel, resolution, type, Vertex, Face, Edge, pointset, PointGrid, 0);

	//extract_isosurface(surface, isolevel - 0.3, resolution, type, Vertex, Face, Edge, pointset, PointGrid, 1);
	//////only for constructing IWP_sheet_solid 
	////int index_temp;
	////for (int i = 0; i < Face.size(); i++)
	////{
	////	index_temp = Face[i].f2;
	////	Face[i].f2 = Face[i].f3;
	////	Face[i].f3 = index_temp;
	////}
	//close_isosurface_sheet2(surface, isolevel, isolevel - 0.3, resolution, Vertex, Face, Edge, pointset, PointGrid);
	////close_isosurface_sheet(surface, isolevel, isolevel - 0.4, resolution, Vertex, Face, PointGrid);

	save_obj(s, Vertex, Face, type);
}

void solid_construction_Tspline(const Isosurface& surface, float isolevel, size_t resolution)
{
	int type = 0;
	vector<Vector3D> Vertex;
	vector<Triangle> Face;
	vector<edge> Edge;
	vector<Point3D> pointset;
	size_t pointRes = resolution + 1; // indicates the # of points per side
	Array3D<Point3D> PointGrid(pointRes, pointRes, pointRes);

	//double ***CValue;
	//CValue = new double **[pointRes];
	//for (int i = 0; i < pointRes; i++)
	//{
	//	CValue[i] = new double *[pointRes];
	//	for (int j = 0; j < pointRes; j++)
	//	{
	//		CValue[i][j] = new double[pointRes];
	//	}
	//}
	//ifstream cvaluefile("E:\\project\\huchuanfeng\\SemanticControlGrid\\src\\GeometryQTGUI\\PorosityDistribution_tooth.txt");
	//for (int i = 0; i < pointRes; i++)
	//{
	//	for (int j = 0; j < pointRes; j++)
	//	{
	//		for (int k = 0; k < pointRes; k++)
	//		{
	//			cvaluefile >> CValue[i][j][k];
	//		}
	//	}
	//}

	Point3D ***Mesh;
	Mesh = new Point3D **[pointRes];
	for (int i = 0; i < pointRes; i++)
	{
		Mesh[i] = new Point3D *[pointRes];
		for (int j = 0; j < pointRes; j++)
		{
			Mesh[i][j] = new Point3D[pointRes];
		}
	}

	double max_z = -999;
	double min_z = 999;

	ifstream meshfile("E:\\project\\huchuanfeng\\SemanticControlGrid\\src\\GeometryQTGUI\\tsolidmesh_tooth.txt");
	for (int i = 0; i < pointRes; i++)
	{
		for (int j = 0; j < pointRes; j++)
		{
			for (int k = 0; k < pointRes; k++)
			{
				meshfile >> Mesh[i][j][k].x >> Mesh[i][j][k].y >> Mesh[i][j][k].z;
				max_z = (max_z>Mesh[i][j][k].z) ? max_z : Mesh[i][j][k].z;
				min_z = (min_z<Mesh[i][j][k].z) ? min_z : Mesh[i][j][k].z;
			}
		}
	}


	for (size_t i = 0; i <= resolution; ++i) {
		for (size_t j = 0; j <= resolution; ++j) {
			for (size_t k = 0; k <= resolution; ++k) {
				float x = Mesh[i][j][k].x, y = Mesh[i][j][k].y, z = Mesh[i][j][k].z;
				//float value = surface.valueAt(6 * x, 6 * y, 6 * z) + 0.9*((float)k / resolution - 0.7);
				double por = 0.3*(1 - (z - min_z) / (max_z - min_z)) + 0.7*(z - min_z) / (max_z - min_z);
				float value = surface.valueAt(6 * x, 6 * y, 6 * z) - (3.4644*por - 1.7322);
				//float value = surface.valueAt(15 * x, 15 * y, 15 * z) - CValue[i][j][k];
				//float value = surface.valueAt(10 * x, 10 * y, 10 * z) - 1.4*z + 0.8;
				//float value = surface.valueAt(6 * x, 6 * y, 6 * z) - 0.6*(CValue[i][j][k] - 7.5) / 15;
				PointGrid.set(i, j, k, Point3D(x, y, z, value));
				pointset.push_back(Point3D(x, y, z, value));
			}
		}
	}
	constructEdge(PointGrid, resolution, Edge);
	extract_isosurface(surface, isolevel, resolution, type, Vertex, Face, Edge, pointset, PointGrid, 0);
	//close_isosurface_network(surface, isolevel, resolution, type, Vertex, Face, Edge, pointset, PointGrid);

	extract_isosurface(surface, isolevel - 0.2, resolution, type, Vertex, Face, Edge, pointset, PointGrid, 1);

	////only for constructing IWP_sheet_solid 
	//int index_temp;
	//for (int i = 0; i < Face.size(); i++)
	//{
	//	index_temp = Face[i].f2;
	//	Face[i].f2 = Face[i].f3;
	//	Face[i].f3 = index_temp;
	//}

	close_isosurface_sheet(surface, isolevel, isolevel - 0.2, resolution, Vertex, Face, Edge, pointset, PointGrid);

	//Vector3D normal;
	//glBegin(GL_TRIANGLES);
	//for (int i = 0; i < Face.size(); i++)
	//{
	//	normal = surface.gradientAt(Vertex[Face[i].f1].x, Vertex[Face[i].f1].y, Vertex[Face[i].f1].z);
	//	glNormal3f(normal.x, normal.y, normal.z);
	//	glVertex3f(Vertex[Face[i].f1].x, Vertex[Face[i].f1].y, Vertex[Face[i].f1].z);
	//	normal = surface.gradientAt(Vertex[Face[i].f2].x, Vertex[Face[i].f2].y, Vertex[Face[i].f2].z);
	//	glNormal3f(normal.x, normal.y, normal.z);
	//	glVertex3f(Vertex[Face[i].f2].x, Vertex[Face[i].f2].y, Vertex[Face[i].f2].z);
	//	normal = surface.gradientAt(Vertex[Face[i].f3].x, Vertex[Face[i].f3].y, Vertex[Face[i].f3].z);
	//	glNormal3f(normal.x, normal.y, normal.z);
	//	glVertex3f(Vertex[Face[i].f3].x, Vertex[Face[i].f3].y, Vertex[Face[i].f3].z);
	//}
	//glEnd();


	save_obj("Tooth_TSporous100_new.obj", Vertex, Face, 1);
}

