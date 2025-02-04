CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -O3 -Wall -Wextra -fopenmp
# CXXFLAGS = -std=c++17 -Iinclude -fsanitize=address -g -O1 -Wall -Wextra -fopenmp

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)
TARGET = dynamont

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f src/*.o $(TARGET)
