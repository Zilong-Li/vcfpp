# detect OS architecture and add flags
CXX      = g++
# CXXFLAGS = -std=c++11 -g -Wall -fsanitize=address
CXXFLAGS = -std=c++11 -Wall
LDFLAGS += -lhts
LIBS     = -lz -lm -lbz2 -llzma -lcurl -lpthread
INC      = -I. -I..

OBJS = test-main.o bcf-reader.o bcf-writer.o

BINS = bcf-reader.bin bcf-writer.bin

.PHONY: all test clean

all: $(BINS) $(OBJS)

test: $(BINS) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< ${INC}

%.bin: %.o test-main.o
	${CXX} ${CXXFLAGS} -o $@ $< test-main.o $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o *.bin

test:
	./bcf-reader.bin --success
