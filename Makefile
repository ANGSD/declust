#modied from htslib makefile
#FLAGS=-O3 -std=c++11
FLAGS=-ggdb -std=c++11

CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

all: decluster


# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif


-include $(OBJ:.o=.d)

ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR) $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp >$*.d

decluster: $(OBJ)
	$(CXX) $(FLAGS)  -o decluster *.o $(HTS_LIBDIR) -lz -llzma -lbz2 -lpthread -lcurl -lgsl 
else
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

decluster: $(OBJ)
	$(CXX) $(FLAGS)  -o decluster *.o -lz -llzma -lbz2 -lpthread -lcurl -lhts -lgsl
endif

clean:	
	rm  -f decluster *.o *.d


testbams := $(wildcard tests/test*bam)

test: $(testbams)
	for bam in $(testbams); do \
		./decluster -0 -w $${bam} -o $${bam%.bam}.test 2> $${bam%.bam}.log; \
		diff $${bam%.bam}.test.hist.txt $${bam%.bam}.expected.hist.txt; \
		bash -v -c "diff <(sed 1d $${bam%.bam}.test.dupstat.txt) <(sed 1d $${bam%.bam}.expected.dupstat.txt)"; \
	done

cleantest:
	rm -v tests/*.test* tests/*.log tests/*bamin.test* tests/*histin.test*
