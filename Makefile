# Variables
CC = g++
SRC = obj_render.cpp
TARGET = obj_render 
CFLAGS = -Wall -Wextra -std=c++11
LIBS =  -lsfml-graphics -lsfml-window -lsfml-system

# Default target
all: $(TARGET)

# Compile the target

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LIBS)

# Clean up
clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: all clean
