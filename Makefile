CC = g++
TARGET_EXEC ?= main

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src

SRCS := $(wildcard $(SRC_DIRS)/*.cpp)
OBJS := $(SRCS:$(SRC_DIRS)/%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := /home/mushf/pce/include
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -g -Wno-unused-result
LDFLAGS := -lgmpxx -lgmp -lboost_iostreams -lboost_system -lboost_filesystem

$(TARGET_EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

# c++ source
$(BUILD_DIR)/%.cpp.o: $(SRC_DIRS)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	$(RM) $(BUILD_DIR)/*
	$(RM) $(TARGET_EXEC)

-include $(DEPS)

MKDIR_P ?= mkdir -p