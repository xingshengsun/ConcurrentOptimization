# This library is free and distributed under
# Mozilla Public License Version 2.0.


CXX:=g++
#CXX:=x86_64-w64-mingw32-g++
DEBUG_FLAG:= -g -03
RELESE_FLAG:= -O3 -s -DNDEBUG
CURRENT_FLAGS:= $(RELESE_FLAG)
CURRENT_FLAGS += -std=c++11 -pthread -I./src -Wall -Wconversion -Wfatal-errors -Wextra
SRC:=./src
BIN:=./bin

LIBS:= # empty
G_LIBS:= -lGL -lGLU -lglut -lGLEW -lSDL -lSDL2main -lSDL2

all:
	@echo "***********************************************"
	@echo Run one of the following commandline examples:
	@echo ""
	@echo make clean
	@echo make design
	@echo "***********************************************"

design:
	$(CXX) $(CURRENT_FLAGS) $(SRC)/main.cpp $(SRC)/initialize.cpp $(SRC)/readDynaOutput.cpp $(SRC)/generateMesh.cpp $(SRC)/transferVariable.cpp $(SRC)/runDyna.cpp -o $(BIN)/Simulation.exe $(LIBS)
	#$(CXX) $(CURRENT_FLAGS) $(SRC)/VariableTransfer.cpp $(SRC)/initialize.cpp -o $(BIN)/VariableTransfer.exe $(LIBS)

clean:
	rm $(BIN)/*.exe
