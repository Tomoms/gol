CC = mpic++
CPPFLAGS = -Ofast -DDEBUG -march=native -mtune=native -mavx2 -mfma -mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store -Wall -Wextra -Wnon-virtual-dtor -Wcast-align -Wunused -Woverloaded-virtual --pedantic -fopenmp -lboost_mpi -lboost_serialization -L/u/dssc/tfonda/boost/lib -lmimalloc -L/u/dssc/tfonda/mimalloc/out/release
OUT = out
SRC = src
OBJS = $(addprefix $(OUT)/, $(patsubst %.cpp, %.o, $(notdir $(wildcard src/*.cpp))))
DEPS = $(wildcard $(OBJS:%=%.d))
INCLUDE_DIRS = $(SRC)/include /u/dssc/tfonda/mimalloc/include /u/dssc/tfonda/boost/include
INCLUDES = $(INCLUDE_DIRS:%=-I%)
TARGET = gol

$(TARGET): $(OBJS)
	$(CC) $(CPPFLAGS) $^ -o $@ $(INCLUDES)

$(OUT)/%.o: $(SRC)/%.cpp
	$(CC) -MD -MP -MF "$@.d" -c $(CPPFLAGS) $< -o $@ $(INCLUDES)

include $(DEPS)

clean:
	rm -r out/
	rm $(TARGET)

$(shell mkdir $(OUT))
