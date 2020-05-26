#include <iostream>
#include <cstdlib>
#include <cassert>
#include <GL/glut.h>

#include "extractisosurface.h"
#include "porosityPereservation.h"
#include "curvature2threshold.h"

#include "Psurface.h"
#include "Gsurface.h"
#include "Dsurface.h"
#include "I_WPsurface.h"
#include "FRDsurface.h"
#include "Lsurface.h"
#include "TubularPsurface.h"
#include "TubularGsurface.h"
#include "I2Ysurface.h"

GLuint preDecimate(const Isosurface& surface,
                   float xMin, float xMax,
                   float yMin, float yMax,
                   float zMin, float zMax,
                   float isolevel,
                   size_t resolution)
{
    GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);

	solid_construction_Bspline(surface, xMin, xMax, yMin, yMax, zMin, zMax, isolevel, resolution);

	//solid_construction_Tspline(surface, isolevel, resolution);

	//solid_construction_Bspline_save(surface, xMin, xMax, yMin, yMax, zMin, zMax, isolevel, resolution,"unit.obj");

	int type_TPMS=0, type_structure=0;

	//porosity_pereservation(surface, isolevel, type_TPMS, type_structure);

	//curvature2threshold(type_TPMS, type_structure);

    glEndList();
    return list;
}

static float rotX;
static float rotY;

static struct { int x, y; } mouse = { -1, -1 };

GLuint list;

void render()
{
    glClearColor(0, 0.3, 0.6, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 1.0, 1.0, -1.0, 0.0 };
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glShadeModel (GL_SMOOTH);

    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    //glPushMatrix();
    //glRotatef(rotY, 1, 0, 0);
    //glRotatef(rotX, 0, 1, 0);
    //glutSolidCube(0.15);
    //glPopMatrix();

    glPushMatrix();
    glRotatef(rotY, 1, 0, 0);
    glRotatef(rotX, 0, 1, 0);

	//glScalef(0.1, 0.1, 0.1);
	glScalef(1, 1, 1);
	glTranslatef(-0.5, -0.5, -0.5);
    glCallList(list);
    glPopMatrix();

    glutSwapBuffers();
}

void changeSize(int w, int h)
{

    // Prevent a divide by zero, when window is too shor
    // (you cant make a window of zero width).
    if(h == 0)
        h = 1;
    float ratio = 1.0* w / h;

    // Use the Projection Matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // Set the viewport to be the entire window
    glViewport(0, 0, w, h);

    // Set the correct perspective.
    gluPerspective(45,ratio,1,1000);

    // Get Back to the Modelview
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Set up the camera
    gluLookAt(0, 0, -3, 0, 0, 0, 0, 1, 0);
}

void mouseDown(int button, int state, int x, int y)
{
    mouse.x = x;
    mouse.y = y;
}

void mouseDragged(int x, int y)
{

    float dx = x - mouse.x;
    float dy = y - mouse.y;

    rotX += dx / 10;
    rotY -= dy / 10;

    mouse.x = x;
    mouse.y = y;

    glutPostRedisplay();

}

void init()
{
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);
	glEnable(GL_CULL_FACE);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	int type = 1;


	if (type == 1){
		Psurface surface;
		list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	}
	else if (type == 2) {
		Gsurface surface;
		list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	}
	else if (type == 3) {
		Dsurface surface;
		list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	}
	else if (type == 4) {
		I_WPsurface surface;
		list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	}
	else if (type == 5) {
		FRDsurface surface;
		list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	}
	else if (type == 6) {
		Lsurface surface;
		list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	}
	else if (type == 7) {
		I2Ysurface surface;
		list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	}
	else if (type == 8) {
		TubularPsurface surface;
		list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	}
	else if (type == 9) {
		TubularGsurface surface;
		list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	}
	//Psurface surface;
	//list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	//Gsurface surface;
	//list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0.2, 100);
	//Dsurface surface;
	//list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	//I_WPsurface surface;
	//list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	//FRDsurface surface;
	//list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	//Lsurface surface;
	//list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	//I2Ysurface surface;
	//list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	//TubularPsurface surface;
	//list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);
	//TubularGsurface surface;
	//list = preDecimate(surface, 0, 1, 0, 1, 0, 1, 0, 100);

}

int main (int argc, char * argv[])
{

    assert(sizeof(char) == 1);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);

    glutInitWindowPosition(100, 100);
    glutInitWindowSize(800, 600);

    glutCreateWindow("Volume Rendering");

    init();

    glutDisplayFunc(render);
    glutReshapeFunc(changeSize);
    glutMouseFunc(mouseDown);
    glutMotionFunc(mouseDragged);

    glutMainLoop();
    return 0;
}
