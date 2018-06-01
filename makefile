# define the C compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -std=c++14 -Wall -O2
DEBUG_FLAGS = -std=c++14 -g -Wall -O2

# define any directories containing header files other than /usr/include
INCLUDES = -I./include

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
STAT_LIBS_PATH = -Wl,-Bstatic -L$(LIB_DIR)
STAT_LIBS = -llatticetester
DYN_LIBS_PATH = -Wl,-Bdynamic -L/usr/local/lib 
DYN_LIBS = -lntl -lgmp

# A few directories we need to be aware of
SRC_DIR = ./src
OBJ_DIR = ./obj
LIB_DIR = ./lib
BIN_DIR = ./bin
PRO_DIR = ./progs


# define the C source files
SRCS = $(wildcard $(SRC_DIR)/*.cc)
PROGS_CC = $(wildcard $(PRO_DIR)/*.cc)

# define the C object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
PROGS_O = $(PROGS_CC:$(PRO_DIR)/%.cc=$(PRO_DIR)/%.o)

all: mkdir lib bin

lib: lib_objects
	rm -f $(LIB_DIR)/liblatticetester.a
	ar rcs $(LIB_DIR)/liblatticetester.a $(OBJS)

lib_objects: $(OBJS)

bin: lib progs_objects

progs_objects: $(PROGS_O)

$(OBJ_DIR)/%.o:$(SRC_DIR)/%.cc
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(PRO_DIR)/%.o:$(PRO_DIR)/%.cc
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
	$(CC) $@ $(STAT_LIBS_PATH) $(STAT_LIBS) $(DYN_LIBS_PATH) $(DYN_LIBS) -o $(BIN_DIR)/$(@:progs/%.o=%) 

#==============================================================================
clean: clean_objects clean_bin clean_lib

clean_objects:
	rm -rf $(OBJ_DIR)
	rm -f $(PRO_DIR)/*.o

clean_bin:
	rm -rf $(BIN_DIR)
	rm -rf build

clean_lib:
	rm -rf $(LIB_DIR)

mkdir:
	mkdir -p bin
	mkdir -p obj
	mkdir -p lib
