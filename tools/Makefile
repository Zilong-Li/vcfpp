HTSINC   = /usr/local/include
HTSLIB   = /usr/local/lib
VCFPP   = ../
CXX      = g++
CXXFLAGS = -std=c++11 -Wall -O3
INC      =  -I$(HTSINC) -I$(VCFPP)
LDFLAGS  =  -L$(HTSLIB) -Wl,-rpath,$(HTSLIB)
LIBS     = -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread
# OBJS     = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
BINS     = $(patsubst %.cpp, %, $(wildcard *.cpp))

.PHONY: all clean test

all: $(BINS)

%: %.cpp
	${CXX} ${CXXFLAGS} -o $@ $< ${INC} $(LIBS) $(LDFLAGS)

clean:
	rm -f $(BINS)

test:
	./vcf_addDS -i test.vcf.gz -o out.vcf.gz
	./vcf_addDS -i test.vcf.gz -s Ind3,Ind1 -o out.vcf.gz
	./vcf_addINFO -i test.vcf.gz -o out.vcf.gz
	./vcf_addINFO -i test.vcf.gz -s Ind0,Ind2 -o out.vcf.gz
