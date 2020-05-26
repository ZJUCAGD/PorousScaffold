/********************************************************************************
** Form generated from reading UI file 'scaffoldgeneration.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SCAFFOLDGENERATION_H
#define UI_SCAFFOLDGENERATION_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ScaffoldGenerationClass
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *ScaffoldGenerationClass)
    {
        if (ScaffoldGenerationClass->objectName().isEmpty())
            ScaffoldGenerationClass->setObjectName(QStringLiteral("ScaffoldGenerationClass"));
        ScaffoldGenerationClass->resize(600, 400);
        menuBar = new QMenuBar(ScaffoldGenerationClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        ScaffoldGenerationClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(ScaffoldGenerationClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        ScaffoldGenerationClass->addToolBar(mainToolBar);
        centralWidget = new QWidget(ScaffoldGenerationClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        ScaffoldGenerationClass->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(ScaffoldGenerationClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        ScaffoldGenerationClass->setStatusBar(statusBar);

        retranslateUi(ScaffoldGenerationClass);

        QMetaObject::connectSlotsByName(ScaffoldGenerationClass);
    } // setupUi

    void retranslateUi(QMainWindow *ScaffoldGenerationClass)
    {
        ScaffoldGenerationClass->setWindowTitle(QApplication::translate("ScaffoldGenerationClass", "ScaffoldGeneration", 0));
    } // retranslateUi

};

namespace Ui {
    class ScaffoldGenerationClass: public Ui_ScaffoldGenerationClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SCAFFOLDGENERATION_H
