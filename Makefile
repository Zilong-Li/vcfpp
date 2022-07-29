# detect OS architecture and add flags
CXX      = g++
# CXXFLAGS = -std=c++11 -g -Wall -fsanitize=address
CXXFLAGS = -std=c++11 -Wall -g
LDFLAGS += -lhts
LIBS     = -lz -lm -lbz2 -llzma -lcurl -lpthread
INC      = -I.

OBJS = test-main.o bcf-reader.o

BINS = bcf-reader.bin

.PHONY: all test clean

all: $(BINS) $(OBJS)

test: $(BINS) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< ${INC}

%.bin: %.o test-main.o
	${CXX} ${CXXFLAGS} -o $@ $< test-main.o $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o *.bin
