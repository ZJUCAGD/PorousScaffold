#include "OpenGlWidget.h"
#include <qdebug.h>

OpenGlWidget::OpenGlWidget(QWidget* parent)
: QGLWidget(parent)
{
	_isreadfile = false;
	setFormat(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer));
	rotationX = 0.0;
	rotationY = 0.0;
	rotationZ = 0.0;
	cameraloc = 4.0;
	diagonal_length = 2.0;
	paintType = 3;//默认绘制模式为面片
	center.x = 0.0;
	center.y = 0.0;
	center.z = 0.0;
}

OpenGlWidget::~OpenGlWidget()
{
}

void OpenGlWidget::initializeGL()
{
	qglClearColor(Qt::darkBlue);
	//加光照信息
	GLfloat sun_light_position[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat sun_light_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat sun_light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat sun_light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

	glLightfv(GL_LIGHT0, GL_POSITION, sun_light_position); //指定第0号光源的位置   
	glLightfv(GL_LIGHT0, GL_AMBIENT, sun_light_ambient); //GL_AMBIENT表示各种光线照射到该材质上，  
	//经过很多次反射后最终遗留在环境中的光线强度（颜色）  
	glLightfv(GL_LIGHT0, GL_DIFFUSE, sun_light_diffuse); //漫反射后~~  
	glLightfv(GL_LIGHT0, GL_SPECULAR, sun_light_specular);//镜面反射后~~~  

	glEnable(GL_LIGHT0); //使用第0号光照  
	glEnable(GL_LIGHTING); //在后面的渲染中使用光照  
	glEnable(GL_DEPTH_TEST); //这句是启用深度测试，这样，在后面的物体会被挡着，例如房子后面有棵树，如果不启用深度测试，  
	//你先画了房子再画树，树会覆盖房子的；但启用深度测试后无论你怎么画，树一定在房子后面（被房子挡着）
}

void OpenGlWidget::resizeGL(int width, int height)
{
	if (height == 0)	height = 1;
	float ratio = 1.0* width / height;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// 设置正确的投影矩阵
	gluPerspective(45, ratio, 0.1, 1000);
	//下面是设置模型视图矩阵
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0f, 1.0f, 0.0f);
}
//绘制函数
void OpenGlWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, cameraloc, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
	//glTranslatef(center.x, center.y, center.z);
	glRotatef(rotationX, 1.0, 0.0, 0.0);
	glRotatef(rotationY, 0.0, 1.0, 0.0);
	glRotatef(rotationZ, 0.0, 0.0, 1.0);
	glTranslatef(-center.x, -center.y, -center.z);

	if (_isreadfile)
	{
		if (paintType == 0)//点模型
		{
			glPointSize(2);
			glBegin(GL_POINTS);
			for (int i = 0; i < Mesh->Vertex.size(); i++)
			{
				glNormal3d(Mesh->VertexNormal[i].x, Mesh->VertexNormal[i].y, Mesh->VertexNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[i].x, Mesh->Vertex[i].y, Mesh->Vertex[i].z);
			}
			glEnd();
		}
		if (paintType == 1)//线框模型
		{
			for (int i = 0; i < Mesh->Face.size(); i++)
			{
				glBegin(GL_LINE_LOOP);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f1].x, Mesh->Vertex[Mesh->Face[i].f1].y, Mesh->Vertex[Mesh->Face[i].f1].z);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f2].x, Mesh->Vertex[Mesh->Face[i].f2].y, Mesh->Vertex[Mesh->Face[i].f2].z);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f3].x, Mesh->Vertex[Mesh->Face[i].f3].y, Mesh->Vertex[Mesh->Face[i].f3].z);
				glEnd();
			}
		}
		else if (paintType == 2)//面片加线框
		{
			//绘制wire
			glDepthRange(0, 0.9999);
			glDisable(GL_LIGHTING);
			glLineWidth(1.0f);
			glColor3f(0.0f, 0.0f, 0.0f);
			for (int i = 0; i < Mesh->Face.size(); i++)
			{
				glBegin(GL_LINE_LOOP);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f1].x, Mesh->Vertex[Mesh->Face[i].f1].y, Mesh->Vertex[Mesh->Face[i].f1].z);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f2].x, Mesh->Vertex[Mesh->Face[i].f2].y, Mesh->Vertex[Mesh->Face[i].f2].z);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f3].x, Mesh->Vertex[Mesh->Face[i].f3].y, Mesh->Vertex[Mesh->Face[i].f3].z);
				glEnd();
			}
			glEnable(GL_LIGHTING);
			// 绘制flat
			glDepthRange(0, 1);
			glShadeModel(GL_FLAT);
			for (int i = 0; i < Mesh->Face.size(); i++)
			{
				glBegin(GL_TRIANGLES);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f1].x, Mesh->Vertex[Mesh->Face[i].f1].y, Mesh->Vertex[Mesh->Face[i].f1].z);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f2].x, Mesh->Vertex[Mesh->Face[i].f2].y, Mesh->Vertex[Mesh->Face[i].f2].z);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f3].x, Mesh->Vertex[Mesh->Face[i].f3].y, Mesh->Vertex[Mesh->Face[i].f3].z);
				glEnd();
			}
		}
		else if (paintType == 3)//面片模型
		{
			glShadeModel(GL_FLAT);
			for (int i = 0; i < Mesh->Face.size(); i++)
			{
				glBegin(GL_TRIANGLES);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f1].x, Mesh->Vertex[Mesh->Face[i].f1].y, Mesh->Vertex[Mesh->Face[i].f1].z);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f2].x, Mesh->Vertex[Mesh->Face[i].f2].y, Mesh->Vertex[Mesh->Face[i].f2].z);
				glNormal3d(Mesh->FaceNormal[i].x, Mesh->FaceNormal[i].y, Mesh->FaceNormal[i].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f3].x, Mesh->Vertex[Mesh->Face[i].f3].y, Mesh->Vertex[Mesh->Face[i].f3].z);
				glEnd();
			}
		}
		else if (paintType == 4)//光滑模型
		{
			glShadeModel(GL_SMOOTH);
			for (int i = 0; i < Mesh->Face.size(); i++)
			{
				glBegin(GL_TRIANGLES);
				glNormal3d(Mesh->VertexNormal[Mesh->Face[i].f1].x, Mesh->VertexNormal[Mesh->Face[i].f1].y, Mesh->VertexNormal[Mesh->Face[i].f1].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f1].x, Mesh->Vertex[Mesh->Face[i].f1].y, Mesh->Vertex[Mesh->Face[i].f1].z);
				glNormal3d(Mesh->VertexNormal[Mesh->Face[i].f2].x, Mesh->VertexNormal[Mesh->Face[i].f2].y, Mesh->VertexNormal[Mesh->Face[i].f2].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f2].x, Mesh->Vertex[Mesh->Face[i].f2].y, Mesh->Vertex[Mesh->Face[i].f2].z);
				glNormal3d(Mesh->VertexNormal[Mesh->Face[i].f3].x, Mesh->VertexNormal[Mesh->Face[i].f3].y, Mesh->VertexNormal[Mesh->Face[i].f3].z);//法向量
				glVertex3d(Mesh->Vertex[Mesh->Face[i].f3].x, Mesh->Vertex[Mesh->Face[i].f3].y, Mesh->Vertex[Mesh->Face[i].f3].z);
				glEnd();
			}
		}
	}
}
void OpenGlWidget::mousePressEvent(QMouseEvent *event)
{
	lastPos = event->pos();
}
void OpenGlWidget::mouseMoveEvent(QMouseEvent *event)
{
	GLfloat dx = GLfloat(event->x() - lastPos.x()) / width();
	GLfloat dy = GLfloat(event->y() - lastPos.y()) / height();
	if (event->buttons() & Qt::LeftButton){
		rotationX += 180 * dy;
		rotationY += 180 * dx;
		updateGL();
	}
	else if (event->buttons() & Qt::RightButton){
		rotationX += 180 * dy;
		rotationZ += 180 * dx;
		updateGL();
	}
	lastPos = event->pos();
}
void OpenGlWidget::wheelEvent(QWheelEvent *event)
{
	if (event->delta() > 0)
	{
		cameraloc = cameraloc + diagonal_length / 5;
		updateGL();
	}
	else if (event->delta() < 0)
	{
		cameraloc = cameraloc - diagonal_length / 5;
		updateGL();
	}
}