# File generated by kdevelop's qmake manager. 
# ------------------------------------------- 
# Subdir relative project main directory: ./plugins/Engine/DynamicEngine/SDECDynamicEngine
# Target is a library:  

LIBS += -lSerialization \
        -lSphere \
        -lMath \
        -lEngine \
        -lInteraction \
        -lMultiMethods \
        -lSDECContactPhysics \
        -lSDECContactGeometry \
        -lActionForce \
        -lActionMomentum \
        -lBody \
        -lSDECParameters \
        -lRigidBodyParameters \
        -lSDECLinkGeometry \
        -lSDECLinkPhysics \
        -rdynamic 
INCLUDEPATH += $(YADEINCLUDEPATH) 
MOC_DIR = $(YADECOMPILATIONPATH) 
UI_DIR = $(YADECOMPILATIONPATH) 
OBJECTS_DIR = $(YADECOMPILATIONPATH) 
QMAKE_LIBDIR = ../../../../toolboxes/Libraries/Serialization/$(YADEDYNLIBPATH) \
               ../../../../plugins/Body/GeometricalModel/Sphere/$(YADEDYNLIBPATH) \
               ../../../../toolboxes/Libraries/Math/$(YADEDYNLIBPATH) \
               ../../../../yade/Engine/$(YADEDYNLIBPATH) \
               ../../../../yade/Interaction/Interaction/$(YADEDYNLIBPATH) \
               ../../../../toolboxes/Libraries/MultiMethods/$(YADEDYNLIBPATH) \
               ../../../../plugins/Interaction/InteractionPhysics/SDECContactPhysics/$(YADEDYNLIBPATH) \
               ../../../../plugins/Interaction/InteractionGeometry/SDECContactGeometry/$(YADEDYNLIBPATH) \
               ../../../../plugins/Engine/Action/ActionForce/$(YADEDYNLIBPATH) \
               ../../../../plugins/Engine/Action/ActionMomentum/$(YADEDYNLIBPATH) \
               ../../../../yade/Body/Body/$(YADEDYNLIBPATH) \
               ../../../../plugins/Body/BodyPhysicalParameters/SDECParameters/$(YADEDYNLIBPATH) \
               ../../../../plugins/Body/BodyPhysicalParameters/RigidBodyParameters/$(YADEDYNLIBPATH) \
               ../../../../plugins/Interaction/InteractionGeometry/SDECLinkGeometry/$(YADEDYNLIBPATH) \
               ../../../../plugins/Interaction/InteractionPhysics/SDECLinkPhysics/$(YADEDYNLIBPATH) \
               $(YADEDYNLIBPATH) 
QMAKE_CXXFLAGS_RELEASE += -lpthread \
                          -pthread 
QMAKE_CXXFLAGS_DEBUG += -lpthread \
                        -pthread 
DESTDIR = $(YADEDYNLIBPATH) 
CONFIG += debug \
          warn_on \
          dll 
TEMPLATE = lib 
HEADERS += SDECDynamicEngine.hpp 
SOURCES += SDECDynamicEngine.cpp 
