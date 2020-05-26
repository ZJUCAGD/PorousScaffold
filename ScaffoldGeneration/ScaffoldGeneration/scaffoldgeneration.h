#ifndef SCAFFOLDGENERATION_H
#define SCAFFOLDGENERATION_H

//#include <QtWidgets/QMainWindow>
//#include "ui_scaffoldgeneration.h"

#include <QtWidgets/QMainWindow>
#include <QGridLayout>
#include <QFileDialog>
#include <QMessageBox>
#include <QButtonGroup>
#include <QRadioButton>
#include "OpenGlWidget.h"
#include <qaction.h>
#include <qtoolbar.h>
#include <qmenu.h>
#include <qmenubar.h>
#include <qstatusbar.h>
#include <QInputDialog>
#include <QDebug>
#include <QPushButton>

class ScaffoldGeneration : public QMainWindow
{
	Q_OBJECT

public:
	ScaffoldGeneration(QWidget *parent = 0);
	~ScaffoldGeneration();
private slots://槽函数
	void fileOpen();//打开图像文件
	void fileSave();//打开图像文件
	void PointsPaint();
	void WirePaint();
	void FlatPaint();
	void FlatlinesPaint();
	void SmoothPaint();
	//void pushButton_clicked();
	//void TPMS_type();
	//void structure_type();

private:
	int type_TPMS;//TPMS类型（P,G,D,IWP）
	int type_structure;//结构类型（Rod,Pore,Sheet）
	QString rsrcPath;//当前路径
	OpenGlWidget *glwidget;
	QWidget *centralWidget;
	QGridLayout *centralLayout;
	QMenu *FileMenu;//文件菜单	
	QMenu *RenderMenu;//文件菜单
	QToolBar *fileToolBar;//工具栏
	QToolBar *SeeToolBar;

	QAction *fileopenAct;//打开网格文件动作
	QAction *filesaveAct;//保存网格文件动作

	QAction *PointsAct;
	QAction *WireAct;
	QAction *FlatAct;
	QAction *FlatlinesAct;
	QAction *SmoothAct;

	////弹窗设置
	//QDialog* dlg;
	//QPushButton *btn;

	//QButtonGroup *TPMSGroup;
	//QRadioButton *PradioButton;
	//QRadioButton *GradioButton;
	//QRadioButton *DradioButton;
	//QRadioButton *IWPradioButton;

	//QButtonGroup *StructureGroup;
	//QRadioButton *RodradioButton;
	//QRadioButton *PoreradioButton;
	//QRadioButton *SheetradioButton;

	void createActions();//创建动作
	void createWidgets();//创建组件
	void createMenus();//创建菜单
	void createToolBars();//创建工具栏
	void createStatusBar();//创建状态栏
	//void creatTPMSDialog();//创建弹窗

	void type_choose();
};
#endif // SCAFFOLDGENERATION_H
