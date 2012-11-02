#Path to project directory
PROJECT_DIR = YOUR_PATH_HERE

#Paths to project dependencies
LCIO_DIR = YOUR_PATH_HERE
MARLIN_DIR = YOUR_PATH_HERE

PROJECT_INCLUDE_DIR = $(PROJECT_DIR)/include/
PROJECT_SOURCE_DIR  = $(PROJECT_DIR)/src/
PROJECT_TEST_DIR  = $(PROJECT_DIR)/tests/
PROJECT_LIBRARY_DIR = $(PROJECT_DIR)/lib/
PROJECT_BINARY_DIR = $(PROJECT_DIR)/

INCLUDES  = -I$(PROJECT_INCLUDE_DIR)
INCLUDES += -I$(LCIO_DIR)/include/
INCLUDES += -I$(MARLIN_DIR)/include/
INCLUDES += -I$(shell $(ROOTSYS)/bin/root-config --incdir)

CC = g++
CFLAGS = -c -Wall -g -w -fPIC -O2
ifdef BUILD_32BIT_COMPATIBLE
    CFLAGS += -m32
endif

SOURCES  = $(wildcard $(PROJECT_SOURCE_DIR)/*.cc)
SOURCES += $(wildcard $(PROJECT_TEST_DIR)/*.cc)

OBJECTS = $(SOURCES:.cc=.o)
DEPENDS = $(OBJECTS:.o=.d)

LIBS  = -L$(LCIO_DIR)/lib -llcio
LIBS += -L$(MARLIN_DIR)/lib -lMarlin
LIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)
ifdef BUILD_32BIT_COMPATIBLE
    LIBS += -m32
endif

LDFLAGS = $(LIBS) -Wl,-rpath

LIBRARY = $(PROJECT_LIBRARY_DIR)/libPandoraAnalysis.so

all: $(SOURCES) $(OBJECTS) AnalysePerformance AnalysePerformanceFull Calibrate ReclusterMonitoring
	$(CC) $(OBJECTS) $(LIBS) -shared -o $(LIBRARY)

AnalysePerformance:
	$(CC) $(INCLUDES) $(LIBS) $(PROJECT_TEST_DIR)AnalysePerformance.cc -o $(PROJECT_BINARY_DIR)AnalysePerformance

AnalysePerformanceFull:
	$(CC) $(INCLUDES) $(LIBS) $(PROJECT_TEST_DIR)AnalysePerformanceFull.cc -o $(PROJECT_BINARY_DIR)AnalysePerformanceFull

Calibrate:
	$(CC) $(INCLUDES) $(LIBS) $(PROJECT_TEST_DIR)Calibrate.cc -o $(PROJECT_BINARY_DIR)Calibrate

ReclusterMonitoring:
	$(CC) $(INCLUDES) $(LIBS) $(PROJECT_TEST_DIR)ReclusterMonitoring.cc -o $(PROJECT_BINARY_DIR)ReclusterMonitoring

$(LIBRARY): $(OBJECTS)
	$(CC) $(LDFLAGS) -fPIC $(OBJECTS) -o $@

-include $(DEPENDS)

.cc.o:
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -MP -MMD -MT $*.o -MT $*.d -MF $*.d -o $*.o $*.cc

clean:
	rm -f $(OBJECTS)
	rm -f $(DEPENDS)
	rm -f $(LIBRARY)
	rm -f $(PROJECT_BINARY_DIR)/AnalysePerformance
	rm -f $(PROJECT_BINARY_DIR)/AnalysePerformanceFull
	rm -f $(PROJECT_BINARY_DIR)/Calibrate
	rm -f $(PROJECT_BINARY_DIR)/ReclusterMonitoring
