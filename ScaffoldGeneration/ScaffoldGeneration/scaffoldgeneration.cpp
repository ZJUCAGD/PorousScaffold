#include "scaffoldgeneration.h"

ScaffoldGeneration::ScaffoldGeneration(QWidget *parent)
	: QMainWindow(parent)
{
	//ui.setupUi(this);

	type_TPMS = 0;
	type_structure = 0;

	QDir dir;
	rsrcPath = dir.currentPath();//获得当前路径
	centralWidget = new QWidget;
	setCentralWidget(centralWidget);
	glwidget = new OpenGlWidget;
	QGridLayout *centralLayout = new QGridLayout;
	centralLayout->addWidget(glwidget, 0, 0);
	centralWidget->setLayout(centralLayout);
	setWindowTitle(tr("ScaffoldGeneration Demo"));

	createActions();
	createMenus();
	createToolBars();
	createStatusBar();
	//creatTPMSDialog();

	glwidget->Mesh = NULL;
}

ScaffoldGeneration::~ScaffoldGeneration()
{
	if (!glwidget->Mesh)
	{
		delete glwidget->Mesh;
		glwidget->Mesh = NULL;
	}

}

/*************************************************
函数功能：创建各个action
*************************************************/
void ScaffoldGeneration::createActions()
{
	fileopenAct = new QAction(QIcon(rsrcPath + "/images/fileopen.png"),"Open TDF file", this);
	fileopenAct->setToolTip("Open TDF file");//设置按钮提示
	fileopenAct->setStatusTip("Open TDF file");//设置状态栏提示
	connect(fileopenAct, SIGNAL(triggered()), this, SLOT(fileOpen()));//action触发时调用函数fileOpen()

	filesaveAct = new QAction(QIcon(rsrcPath + "/images/filesave.png"), "Save Mesh", this);
	filesaveAct->setToolTip("Save Mesh");//设置按钮提示
	filesaveAct->setStatusTip("Save Mesh");//设置状态栏提示
	connect(filesaveAct, SIGNAL(triggered()), this, SLOT(fileSave()));//action触发时调用函数fileSave()
	filesaveAct->setEnabled(false);


	PointsAct = new QAction(QIcon(rsrcPath + "/images/points.png"),"Points", this);
	PointsAct->setStatusTip("Points");//设置状态栏提示
	connect(PointsAct, SIGNAL(triggered()), this, SLOT(PointsPaint()));
	PointsAct->setEnabled(false);
	PointsAct->setCheckable(true);

	WireAct = new QAction(QIcon(rsrcPath + "/images/wire.png"), "Wireframe", this);
	WireAct->setStatusTip("Wireframe");//设置状态栏提示
	connect(WireAct, SIGNAL(triggered()), this, SLOT(WirePaint()));
	WireAct->setEnabled(false);
	WireAct->setCheckable(true);

	FlatlinesAct = new QAction(QIcon(rsrcPath + "/images/flatlines.png"), "Flat Lines", this);
	FlatlinesAct->setStatusTip("Flat Lines");//设置状态栏提示
	connect(FlatlinesAct, SIGNAL(triggered()), this, SLOT(FlatlinesPaint()));
	FlatlinesAct->setEnabled(false);
	FlatlinesAct->setCheckable(true);

	FlatAct = new QAction(QIcon(rsrcPath + "/images/flat.png"),"Flat", this);
	FlatAct->setStatusTip("Flat");//设置状态栏提示
	connect(FlatAct, SIGNAL(triggered()), this, SLOT(FlatPaint()));
	FlatAct->setEnabled(false);
	FlatAct->setCheckable(true);

	SmoothAct = new QAction(QIcon(rsrcPath + "/images/smooth.png"), "Smooth", this);
	SmoothAct->setStatusTip("Smooth");//设置状态栏提示
	connect(SmoothAct, SIGNAL(triggered()), this, SLOT(SmoothPaint()));
	SmoothAct->setEnabled(false);
	SmoothAct->setCheckable(true);
}
/*************************************************
函数功能：创建菜单
*************************************************/
void ScaffoldGeneration::createMenus()
{
	FileMenu = menuBar()->addMenu("File");
	FileMenu->addAction(fileopenAct);
	FileMenu->addAction(filesaveAct);
	RenderMenu = menuBar()->addMenu("Render");
	RenderMenu->addAction(PointsAct);
	RenderMenu->addAction(WireAct);
	RenderMenu->addAction(FlatlinesAct);
	RenderMenu->addAction(FlatAct);
	RenderMenu->addAction(SmoothAct);
}
/*************************************************
函数功能：创建工具栏
*************************************************/
void ScaffoldGeneration::createToolBars()
{
	fileToolBar = addToolBar(QStringLiteral("File"));
	fileToolBar->addAction(fileopenAct);
	fileToolBar->addSeparator();
	fileToolBar->addAction(filesaveAct);
	SeeToolBar = addToolBar(QStringLiteral("Render"));
	SeeToolBar->addAction(PointsAct);
	SeeToolBar->addSeparator();
	SeeToolBar->addAction(WireAct);
	SeeToolBar->addSeparator();
	SeeToolBar->addAction(FlatlinesAct);
	SeeToolBar->addSeparator();
	SeeToolBar->addAction(FlatAct);
	SeeToolBar->addSeparator();
	SeeToolBar->addAction(SmoothAct);

}
/*************************************************
函数功能：创建状态栏
*************************************************/
void ScaffoldGeneration::createStatusBar()
{
	statusBar()->showMessage("Ready");
	statusBar()->setStyleSheet("QWidget{ border: 1px solid #d3d3d3; background - color:#FFFFFF }");
}
/*************************************************
函数功能：创建弹窗选择TPMS类型
*************************************************/
//void ScaffoldGeneration::creatTPMSDialog()
//{
//	dlg = new QDialog(this);
//	dlg->setFixedSize(240, 180);
//
//	PradioButton = new QRadioButton("P Surface", dlg);
//	PradioButton->setGeometry(QRect(20, 10, 100, 20));
//	GradioButton = new QRadioButton("G Surface", dlg);
//	GradioButton->setGeometry(QRect(20, 40, 100, 20));
//	DradioButton = new QRadioButton("D Surface", dlg);
//	DradioButton->setGeometry(QRect(20, 70, 100, 20));
//	IWPradioButton = new QRadioButton("I-WP Surface", dlg);
//	IWPradioButton->setGeometry(QRect(20, 100, 100, 20));
//	TPMSGroup = new QButtonGroup(dlg);
//	TPMSGroup->addButton(PradioButton, 0);
//	TPMSGroup->addButton(GradioButton, 1);
//	TPMSGroup->addButton(DradioButton, 2);
//	TPMSGroup->addButton(IWPradioButton, 3);
//	PradioButton->setChecked(true);
//
//	RodradioButton = new QRadioButton("Rod Structure", dlg);
//	RodradioButton->setGeometry(QRect(120, 10, 120, 20));
//	PoreradioButton = new QRadioButton("Pore Structure", dlg);
//	PoreradioButton->setGeometry(QRect(120, 40, 120, 20));
//	SheetradioButton = new QRadioButton("Sheet Structure", dlg);
//	SheetradioButton->setGeometry(QRect(120, 70, 120, 20));
//	StructureGroup = new QButtonGroup(dlg);
//	StructureGroup->addButton(RodradioButton, 0);
//	StructureGroup->addButton(PoreradioButton, 1);
//	StructureGroup->addButton(SheetradioButton, 2);
//	RodradioButton->setChecked(true);
//
//	btn = new QPushButton("OK", dlg);
//	btn->setGeometry(QRect(80, 130, 80, 40));
//
//	connect(PradioButton, SIGNAL(clicked(bool)), this, SLOT(TPMS_type()));
//	connect(GradioButton, SIGNAL(clicked(bool)), this, SLOT(TPMS_type()));
//	connect(DradioButton, SIGNAL(clicked(bool)), this, SLOT(TPMS_type()));
//	connect(IWPradioButton, SIGNAL(clicked(bool)), this, SLOT(TPMS_type()));
//
//	connect(RodradioButton, SIGNAL(clicked(bool)), this, SLOT(structure_type()));
//	connect(PoreradioButton, SIGNAL(clicked(bool)), this, SLOT(structure_type()));
//	connect(SheetradioButton, SIGNAL(clicked(bool)), this, SLOT(structure_type()));
//
//	connect(btn, SIGNAL(clicked(bool)), this, SLOT(pushButton_clicked()));
//
//	dlg->setModal(true);
//}
/*************************************************
函数功能：用于打开文件，点击打开后被调用
*************************************************/
void ScaffoldGeneration::fileOpen()
{
	//出现一个打开文件的对话框，选择了文件以后可以得到文件的路径
	QString fileName = QFileDialog::getOpenFileName(this, "Open", QString(), "TDF file (*.tdf);;All files (*.*)");
	if (fileName.isEmpty()) {
		//QMessageBox::information(NULL, tr("Path"), tr("You didn't select any files."));
	}
	else
	{
		if (!glwidget->Mesh)
		{
			delete glwidget->Mesh;
			glwidget->Mesh = NULL;
		}
		glwidget->Mesh = new ScaffoldMesh(fileName.toLocal8Bit().replace("/", "\\").toStdString());
		type_TPMS = glwidget->Mesh->type_TPMS;
		type_structure = glwidget->Mesh->type_structure;


		statusBar()->showMessage("Busy");
		glwidget->Mesh->scaffold_generation();
		glwidget->_isreadfile = true;
		glwidget->center = glwidget->Mesh->center;
		glwidget->diagonal_length = glwidget->Mesh->diagonal_length;
		glwidget->cameraloc = 2.0*glwidget->Mesh->diagonal_length;
		glwidget->update();
		statusBar()->showMessage(" ");

		//fileopenAct->setEnabled(false);
		filesaveAct->setEnabled(true);
		PointsAct->setEnabled(true);
		WireAct->setEnabled(true);
		FlatlinesAct->setEnabled(true);
		FlatAct->setEnabled(true);
		SmoothAct->setEnabled(true);
		FlatAct->setChecked(true);
	}
}
void ScaffoldGeneration::fileSave()
{
	//QString fileName = QFileDialog::getSaveFileName(this, "Save", QString(), "STL file(*.stl)");
	//if (fileName.isEmpty())
	//{
	//}
	//else
	//{
	//	statusBar()->showMessage("Saving");
	//	ofstream ss(fileName.toLocal8Bit().replace("/", "\\").toStdString());
	//	ss << "solid STL generated by ScaffoldGeneration" << "\n";
	//	if ((type_TPMS != 3 && type_structure == 1) || (type_TPMS == 3 && type_structure != 1))
	//	{
	//		for (int i = 0; i < glwidget->Mesh->Face.size(); i++)
	//		{
	//			ss << "  facet normal " << glwidget->Mesh->FaceNormal[i].x << " " << glwidget->Mesh->FaceNormal[i].y << " " << glwidget->Mesh->FaceNormal[i].z << "\n";
	//			ss << "    outer loop" << "\n";
	//			ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].z << "\n";
	//			ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].z << "\n";
	//			ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].z << "\n";
	//			ss << "    endloop" << "\n";
	//			ss << "  endfacet" << "\n";
	//		}
	//	}
	//	else
	//	{
	//		for (int i = 0; i < glwidget->Mesh->Face.size(); i++)
	//		{
	//			ss << "  facet normal " << glwidget->Mesh->FaceNormal[i].x << " " << glwidget->Mesh->FaceNormal[i].y << " " << glwidget->Mesh->FaceNormal[i].z << "\n";
	//			ss << "    outer loop" << "\n";
	//			ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].z << "\n";
	//			ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].z << "\n";
	//			ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].z << "\n";
	//			ss << "    endloop" << "\n";
	//			ss << "  endfacet" << "\n";
	//		}
	//	}
	//	ss << "endsolid vcg" << "\n";
	//	ss.close();
	//	statusBar()->showMessage("Saved");
	//}

	QString fileName = QFileDialog::getSaveFileName(this, "Save", QString(), "STL file(*.stl);;OBJ file (*.obj)");
	if (fileName.isEmpty())
	{
	}
	else
	{
		statusBar()->showMessage("Saving");
		ofstream ss(fileName.toLocal8Bit().replace("/", "\\").toStdString());
		if (fileName[fileName.size() - 1] == 'l')
		{
			ss << "solid STL generated by ScaffoldGeneration" << "\n";
			if ((type_TPMS != 3 && type_structure == 1) || (type_TPMS == 3 && type_structure != 1))
			{
				for (int i = 0; i < glwidget->Mesh->Face.size(); i++)
				{
					ss << "  facet normal " << glwidget->Mesh->FaceNormal[i].x << " " << glwidget->Mesh->FaceNormal[i].y << " " << glwidget->Mesh->FaceNormal[i].z << "\n";
					ss << "    outer loop" << "\n";
					ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].z << "\n";
					ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].z << "\n";
					ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].z << "\n";
					ss << "    endloop" << "\n";
					ss << "  endfacet" << "\n";
				}
			}
			else
			{
				for (int i = 0; i < glwidget->Mesh->Face.size(); i++)
				{
					ss << "  facet normal " << glwidget->Mesh->FaceNormal[i].x << " " << glwidget->Mesh->FaceNormal[i].y << " " << glwidget->Mesh->FaceNormal[i].z << "\n";
					ss << "    outer loop" << "\n";
					ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f1].z << "\n";
					ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f3].z << "\n";
					ss << "      vertex  " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].x << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].y << " " << glwidget->Mesh->Vertex[glwidget->Mesh->Face[i].f2].z << "\n";
					ss << "    endloop" << "\n";
					ss << "  endfacet" << "\n";
				}
			}
			ss << "endsolid vcg" << "\n";
			ss.close();
		}
		else
		{
			ss << "# Vertices: " << glwidget->Mesh->Vertex.size() << " Faces: " << glwidget->Mesh->Face.size() << endl;
			for (int i = 0; i < glwidget->Mesh->Vertex.size(); i++)
			{
				ss << "v " << glwidget->Mesh->Vertex[i].x << ' ' << glwidget->Mesh->Vertex[i].y << ' ' << glwidget->Mesh->Vertex[i].z << endl;
			}
			if ((type_TPMS != 3 && type_structure == 1) || (type_TPMS == 3 && type_structure != 1))
			{
				for (int i = 0; i < glwidget->Mesh->Face.size(); i++)
				{
					ss << "f " << glwidget->Mesh->Face[i].f1 + 1 << ' ' << glwidget->Mesh->Face[i].f2 + 1 << ' ' << glwidget->Mesh->Face[i].f3 + 1 << endl;
				}
			}
			else
			{
				for (int i = 0; i < glwidget->Mesh->Face.size(); i++)
				{
					ss << "f " << glwidget->Mesh->Face[i].f1 + 1 << ' ' << glwidget->Mesh->Face[i].f3 + 1 << ' ' << glwidget->Mesh->Face[i].f2 + 1 << endl;
					//ss << "f " << glwidget->Mesh->Face[i].f1 + 1 << ' ' << glwidget->Mesh->Face[i].f2 + 1 << ' ' << glwidget->Mesh->Face[i].f3 + 1 << endl;
				}
			}
			ss.close();
		}
		statusBar()->showMessage("Saved");
	}
}
void ScaffoldGeneration::PointsPaint()
{
	PointsAct->setChecked(true);
	WireAct->setChecked(false);
	FlatAct->setChecked(false);
	FlatlinesAct->setChecked(false);
	SmoothAct->setChecked(false);
	glwidget->paintType = 0;
	glwidget->update();
}
void ScaffoldGeneration::WirePaint()
{
	PointsAct->setChecked(false);
	WireAct->setChecked(true);
	FlatAct->setChecked(false);
	FlatlinesAct->setChecked(false);
	SmoothAct->setChecked(false);
	glwidget->paintType = 1;
	glwidget->update();
}
void ScaffoldGeneration::FlatlinesPaint()
{
	PointsAct->setChecked(false);
	WireAct->setChecked(false);
	FlatAct->setChecked(false);
	FlatlinesAct->setChecked(true);
	SmoothAct->setChecked(false);
	glwidget->paintType = 2;
	glwidget->update();
}
void ScaffoldGeneration::FlatPaint()
{
	PointsAct->setChecked(false);
	WireAct->setChecked(false);
	FlatAct->setChecked(true);
	FlatlinesAct->setChecked(false);
	SmoothAct->setChecked(false);
	glwidget->paintType = 3;
	glwidget->update();
}
void ScaffoldGeneration::SmoothPaint()
{
	PointsAct->setChecked(false);
	WireAct->setChecked(false);
	FlatAct->setChecked(false);
	FlatlinesAct->setChecked(false);
	SmoothAct->setChecked(true);
	glwidget->paintType = 4;
	glwidget->update();
}

/*************************************************
函数功能：为用户提供弹窗选择TPMS单元和体结构类型
*************************************************/
//void ScaffoldGeneration::type_choose()
//{
//	dlg->show();
//}
//void ScaffoldGeneration::pushButton_clicked()
//{
//	dlg->close();
//
//	statusBar()->showMessage("Busy");
//	glwidget->Mesh->type_TPMS = type_TPMS;
//	glwidget->Mesh->type_structure = type_structure;
//	glwidget->Mesh->scaffold_generation();
//	glwidget->_isreadfile = true;
//	glwidget->center = glwidget->Mesh->center;
//	glwidget->diagonal_length = glwidget->Mesh->diagonal_length;
//	glwidget->cameraloc = 2.0*glwidget->Mesh->diagonal_length;
//	glwidget->update();
//	statusBar()->showMessage(" ");
//
//	//fileopenAct->setEnabled(false);
//	filesaveAct->setEnabled(true);
//	PointsAct->setEnabled(true);
//	WireAct->setEnabled(true);
//	FlatlinesAct->setEnabled(true);
//	FlatAct->setEnabled(true);
//	SmoothAct->setEnabled(true);
//	FlatAct->setChecked(true);
//}
//void ScaffoldGeneration::TPMS_type()
//{
//	type_TPMS = TPMSGroup->checkedId();
//}
//void ScaffoldGeneration::structure_type()
//{
//	type_structure = StructureGroup->checkedId();
//}