# Directories
SRC := src
LIB := lib
OBJ := obj
BIN := bin
MODELS := test_models
MODEL_FILE_EXT := .txt

# Compiler and flags
CC := gcc
CFLAGS := 
COMPILE_FLAGS := -Wall  -Wextra -Wshadow -I$(SRC) -I$(LIB) -L$(LIB)
LIBS := -lc -lgsl -lgslcblas -lm -lrc

DEBUG ?= 0

ifeq ($(DEBUG),1)
  BUILD_FLAGS := $(CFLAGS) $(COMPILE_FLAGS) -g -O0
else
  BUILD_FLAGS := $(CFLAGS) $(COMPILE_FLAGS) -O2
endif

# Target
TARGET := zmax

# Source files
SRCS := $(shell find $(SRC) -type f -name "*.c")

# Object files
OBJS := $(patsubst %.c, $(OBJ)/%.o, $(SRCS))

# Default target
all: $(BIN)/$(TARGET)

# Link
$(BIN)/$(TARGET): $(OBJS)
	@mkdir -p $(BIN)
	$(CC) $(BUILD_FLAGS) $^ $(LIBS) -o $@

# Compile all .c files to .o files
$(OBJ)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(BUILD_FLAGS) -c $< -o $@

# Run release
run: all
	@if [ -n "$(ARGS)" ]; then \
		./$(BIN)/$(TARGET) $(MODELS)/$(ARGS)$(MODEL_FILE_EXT); \
	else \
		./$(BIN)/$(TARGET); \
	fi

# Call make with debug
debug:
	$(MAKE) DEBUG=1 all

# Run debugger
drun: DEBUG := 1
drun: DEBUG := 1
drun: all
	@if [ -n "$(ARGS)" ]; then \
		gdb --args ./$(BIN)/$(TARGET) $(MODELS)/$(ARGS)$(MODEL_FILE_EXT); \
	else \
		gdb ./$(BIN)/$(TARGET); \
	fi

# Clean build artifacts
clean:
	rm -rf $(OBJ) $(BIN)/$(TARGET)

.PHONY: all run debug drun clean
