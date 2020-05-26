#ifndef OPENGLWIDGET_H
#define OPENGLWIDGET_H
#include <QGLWidget>
#include <QMouseEvent>
#include <qgl.h>
#include <gl/gl.h>  
#include <gl/glu.h>  
#include <ScaffoldMesh.h>

class OpenGlWidget : public QGLWidget
{
	Q_OBJECT

public:
	explicit OpenGlWidget(QWidget *parent = 0);
	~OpenGlWidget();

	Vector3D center;
	bool _isreadfile;
	int paintType;
	ScaffoldMesh *Mesh;
	GLfloat diagonal_length;
	GLfloat cameraloc;

protected:
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
private:
	GLfloat rotationX;
	GLfloat rotationY;
	GLfloat rotationZ;
	QPoint lastPos;

	//vector<Vector3D> vertex_normal;
};

#endif // OPENGLWIDGET_H
