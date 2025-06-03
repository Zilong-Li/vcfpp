HTSINC   = /usr/local/include
HTSLIB   = /usr/local/lib
# detect OS architecture and add flags
CXX      = g++
CXXFLAGS = -std=c++11 -Wall -g -O2 -fsanitize=address
LDFLAGS  =  -L$(HTSLIB) -Wl,-rpath,$(HTSLIB)
LIBS     = -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread
INC      = -I. -I.. -I${HTSINC}

OBJS = test-main.o bcf-header.o bcf-variant.o bcf-reader.o bcf-writer.o read-stdin.o

BINS = bcf-header.bin bcf-variant.bin bcf-reader.bin bcf-writer.bin read-stdin

.PHONY: all test clean

all: $(BINS) $(OBJS)

test: $(BINS) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< ${INC}

%.bin: %.o test-main.o
	${CXX} ${CXXFLAGS} -o $@ $< test-main.o $(LDFLAGS) $(LIBS)

read-stdin: read-stdin.o 
	${CXX} ${CXXFLAGS} -o $@ $< $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o *.bin read-stdin

test:
	cat test-vcf-read.vcf | ./read-stdin
	./bcf-header.bin --success 2> /dev/null
	./bcf-reader.bin --success 2> /dev/null
	./bcf-writer.bin --success 2> /dev/null
	./bcf-variant.bin --success 2> /dev/null
