/****************************************************************************
** Meta object code from reading C++ file 'scaffoldgeneration.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../scaffoldgeneration.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'scaffoldgeneration.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_ScaffoldGeneration_t {
    QByteArrayData data[9];
    char stringdata0[97];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ScaffoldGeneration_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ScaffoldGeneration_t qt_meta_stringdata_ScaffoldGeneration = {
    {
QT_MOC_LITERAL(0, 0, 18), // "ScaffoldGeneration"
QT_MOC_LITERAL(1, 19, 8), // "fileOpen"
QT_MOC_LITERAL(2, 28, 0), // ""
QT_MOC_LITERAL(3, 29, 8), // "fileSave"
QT_MOC_LITERAL(4, 38, 11), // "PointsPaint"
QT_MOC_LITERAL(5, 50, 9), // "WirePaint"
QT_MOC_LITERAL(6, 60, 9), // "FlatPaint"
QT_MOC_LITERAL(7, 70, 14), // "FlatlinesPaint"
QT_MOC_LITERAL(8, 85, 11) // "SmoothPaint"

    },
    "ScaffoldGeneration\0fileOpen\0\0fileSave\0"
    "PointsPaint\0WirePaint\0FlatPaint\0"
    "FlatlinesPaint\0SmoothPaint"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ScaffoldGeneration[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   49,    2, 0x08 /* Private */,
       3,    0,   50,    2, 0x08 /* Private */,
       4,    0,   51,    2, 0x08 /* Private */,
       5,    0,   52,    2, 0x08 /* Private */,
       6,    0,   53,    2, 0x08 /* Private */,
       7,    0,   54,    2, 0x08 /* Private */,
       8,    0,   55,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void ScaffoldGeneration::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ScaffoldGeneration *_t = static_cast<ScaffoldGeneration *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->fileOpen(); break;
        case 1: _t->fileSave(); break;
        case 2: _t->PointsPaint(); break;
        case 3: _t->WirePaint(); break;
        case 4: _t->FlatPaint(); break;
        case 5: _t->FlatlinesPaint(); break;
        case 6: _t->SmoothPaint(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject ScaffoldGeneration::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_ScaffoldGeneration.data,
      qt_meta_data_ScaffoldGeneration,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *ScaffoldGeneration::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ScaffoldGeneration::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_ScaffoldGeneration.stringdata0))
        return static_cast<void*>(const_cast< ScaffoldGeneration*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int ScaffoldGeneration::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
