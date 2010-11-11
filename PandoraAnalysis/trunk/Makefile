#Path to project directory
PROJECT_DIR = YOUR_PATH_HERE

#Paths to project dependencies
LCIO_DIR = YOUR_PATH_HERE
MARLIN_DIR = YOUR_PATH_HERE

PROJECT_INCLUDE_DIR = $(PROJECT_DIR)/include/
PROJECT_SOURCE_DIR  = $(PROJECT_DIR)/src/
PROJECT_LIBRARY_DIR = $(PROJECT_DIR)/lib/
PROJECT_BINARY_DIR = $(PROJECT_DIR)/

INCLUDES  = -I$(PROJECT_INCLUDE_DIR)
INCLUDES += -I$(LCIO_DIR)/include/
INCLUDES += -I$(MARLIN_DIR)/include/
INCLUDES += -I$(shell $(ROOTSYS)/bin/root-config --incdir)

CC = gcc
CFLAGS = -c -Wall -g -w -fPIC -m32
CFLAGS += $(INCLUDES)
ifdef BUILD_32BIT_COMPATIBLE
    CFLAGS += -m32
endif

SOURCES = $(PROJECT_SOURCE_DIR)MCPfoMaker.cc $(PROJECT_SOURCE_DIR)MCTree.cc $(PROJECT_SOURCE_DIR)PfoAnalysis.cc

OBJECTS = $(SOURCES:.cc=.o)

LIBS  = -L$(LCIO_DIR)/lib -llcio
LIBS += -L$(MARLIN_DIR)/lib -lMarlin
LIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)
ifdef BUILD_32BIT_COMPATIBLE
    LIBS += -m32
endif

LDFLAGS = $(LIBS) -Wl,-rpath

LIBRARY = $(PROJECT_LIBRARY_DIR)/libPandoraAnalysis.so

all: $(SOURCES) $(OBJECTS) analysePerformance
	$(CC) $(OBJECTS) $(LIBS) -shared -o $(LIBRARY)

analysePerformance:
	@echo Creating binary: $(PROJECT_BINARY_DIR)analysePerformance
	$(CC) $(INCLUDES) $(LIBS) $(PROJECT_SOURCE_DIR)AnalysePerformance.cc -o $(PROJECT_BINARY_DIR)analysePerformance
	@echo Created binary: $(PROJECT_BINARY_DIR)analysePerformance

$(LIBRARY): $(OBJECTS)
	$(CC) $(LDFLAGS) -fPIC $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $(DEFINES) $< -o $@

clean:
	rm -f $(OBJECTS)
	rm -f $(LIBRARY)
	rm -f analysePerformance
