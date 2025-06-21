# Compiler and flags
CC := gcc
CFLAGS := 
COMPILE_FLAGS := -Wall -Wextra -Wshadow
LIBS := -lc -lgsl -lgslcblas -lm

DEBUG ?= 0

ifeq ($(DEBUG),1)
  BUILD_FLAGS := $(CFLAGS) $(COMPILE_FLAGS) $(LIBS) -g -O0
else
  BUILD_FLAGS := $(CFLAGS) $(COMPILE_FLAGS) $(LIBS) -O2
endif

# Directories
SRC := src
OBJ := obj
BIN := bin
MODELS := test_models

# Target
TARGET := simplex

# Source files
SRCS := $(shell find $(SRC) -type f -name "*.c")

# Object files
OBJS := $(patsubst %.c, $(OBJ)/%.o, $(SRCS))

# Default target
all: $(BIN)/$(TARGET)

# Link
$(BIN)/$(TARGET): $(OBJS)
	@mkdir -p $(BIN)
	$(CC) $(BUILD_FLAGS) $^ -o $@

# Compile all .c files to .o files
$(OBJ)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(BUILD_FLAGS) -I$(SRC) -c $< -o $@

# Run release
run: all
	@if [ -n "$(ARGS)" ]; then \
		./$(BIN)/$(TARGET) $(MODELS)/$(ARGS); \
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
		gdb --args ./$(BIN)/$(TARGET) $(MODELS)/$(ARGS); \
	else \
		gdb ./$(BIN)/$(TARGET); \
	fi

# Clean build artifacts
clean:
	rm -rf $(OBJ) $(BIN)/$(TARGET)

.PHONY: all run debug drun clean