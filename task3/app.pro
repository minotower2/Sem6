QT += widgets

CONFIG += c++11
QMAKE_CXXFLAGS += -W -Wall -Werror -Wfloat-equal -Wunused -Wnon-virtual-dtor -isystem $$[QT_INSTALL_HEADERS]

OBJECTS_DIR = build/obj
MOC_DIR     = build/moc
DESTDIR     = build/bin

CONFIG(debug, debug|release) {
    QMAKE_CXXFLAGS += -g
}

CONFIG(release, debug|release) {
    DEFINES += NDEBUG
    QMAKE_CXXFLAGS += -O2
}

CORE_H = core/Approximation.h 
CORE_S = core/Approximation.cpp

GUI_H = gui/GraphWidget.h gui/3DRender.h gui/Data.h gui/InfoDisplay.h
GUI_S = gui/GraphWidget.cpp gui/3DRender.cpp gui/Data.cpp gui/InfoDisplay.cpp

MODEL_H = model/Model.h
MODEL_S = model/Model.cpp

INTERFACE_H = interface/MainWindow.h
INTERFACE_S = interface/MainWindow.cpp interface/main.cpp

RESOURCES_H = resources/Storage.h
RESOURCES_S = resources/Storage.cpp

UTILS_H = utils/Geometry.h utils/CommonDefs.h
UTILS_S = utils/Geometry.cpp

EXTRA_DEBUG = 

HEADERS += \
		$$CORE_H \
		$$GUI_H \
		$$MODEL_H \
		$$INTERFACE_H \
		$$RESOURCES_H \ 
		$$UTILS_H

SOURCES += \
		$$CORE_S\
		$$GUI_S \
		$$MODEL_S \
		$$INTERFACE_S \
		$$RESOURCES_S \ 	
		$$UTILS_S \
		$$EXTRA_DEBUG
