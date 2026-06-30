# Makefile for TylerSort project

# Directories
SRC_DIR  := src
OBJ_DIR  := obj
BIN_DIR  := bin
INC_DIR  := include
INSTALL_DIR := $(HOME)/.local/bin

# Compiler and flags
CXX       := g++
CASORT_LIB := ../CASort/lib/libCASort.so
BREW_PREFIX := $(shell brew --prefix 2>/dev/null)
BOOST_LIB := 
ifneq ($(BREW_PREFIX),)
BOOST_LIB := -L$(BREW_PREFIX)/lib
endif
CXXFLAGS  := `root-config --cflags` -I$(HOME)/.local/include -I$(INC_DIR) -fPIC 
LDFLAGS   := `root-config --libs` -L$(HOME)/.local/lib $(BOOST_LIB) -lCASort -lboost_program_options -Wl,-rpath,$(HOME)/.local/lib
DEBUGFLAGS := -g -O0

# Target executable name
TARGET := $(BIN_DIR)/TylerSort

# Source and object files
SOURCES  := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS  := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

# Default target
all: $(TARGET)

# Debug target
debug: CXXFLAGS += $(DEBUGFLAGS)
debug: clean all

# Ensure CASort is installed before compiling or linking TylerSort
$(CASORT_LIB):
	$(MAKE) -C ../CASort install

$(OBJECTS): $(CASORT_LIB)

# Link object files into the executable
$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Install/uninstall targets
install: $(TARGET)
	@mkdir -p $(INSTALL_DIR)
	install -m 755 $(TARGET) $(INSTALL_DIR)/$(notdir $(TARGET))

uninstall:
	rm -f $(INSTALL_DIR)/$(notdir $(TARGET))

# Clean target
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# A quick test target
test: $(TARGET)
	$(TARGET) --caldir=data/energy_calibration --gsfile=data/gain_correction/70Ge_default.cags data/trees/root_data_70Ge_run208.mvmelst.bin_tree.root data/hists/test/test.hists.root

.PHONY: all clean debug