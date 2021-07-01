#modied from htslib makefile
#FLAGS=-O3 -std=c++11
FLAGS=-ggdb -std=c++11

CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

all: superduper


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

superduper: $(OBJ)
	$(CXX) $(FLAGS)  -o superduper *.o $(HTS_LIBDIR) -lz -llzma -lbz2 -lpthread -lcurl -lgsl 
else
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

superduper: $(OBJ)
	$(CXX) $(FLAGS)  -o superduper *.o -lz -llzma -lbz2 -lpthread -lcurl -lhts -lgsl
endif

clean:	
	rm  -f superduper *.o *.d


testsams := $(wildcard tests/*sam)

test: $(testsams)
	for sam in $(testsams); do \
		samtools view -bS $${sam} > $${sam%.sam}.bam; \
		./superduper -w $${sam%.sam}.bam -o $${sam%.sam}.out 2> $${sam%.sam}.log; \
		diff $${sam%.sam}.out.hist.txt $${sam%.sam}.expected.hist.txt; \
	done

bam: $(testsams)
	for sam in $(testsams); do \
		samtools view -bS $${sam} > $${sam%.sam}.bam; \
	done

testclean:
	rm -v tests/*.bam tests/*.out* tests/*.log
