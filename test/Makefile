# detect OS architecture and add flags
CXX      = g++
# CXXFLAGS = -std=c++11 -g -Wall -fsanitize=address
CXXFLAGS = -std=c++11 -Wall
HTSINC   = /usr/local/include
HTSLIB   = /usr/local/lib
LDFLAGS  =  -L$(HTSLIB) -Wl,-rpath,$(HTSLIB)
LIBS     = -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread
INC      = -I. -I..

OBJS = test-main.o bcf-variant.o bcf-reader.o bcf-writer.o example.o

BINS = bcf-variant.bin bcf-reader.bin bcf-writer.bin example.bin

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
	./bcf-writer.bin --success
	./bcf-variant.bin --success
	./example.bin --success
