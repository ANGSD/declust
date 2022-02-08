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
testsmallbams := $(wildcard tests/small.bam)

test: $(testbams) $(testsmallbams)
	for bam in $(testbams); do \
		./decluster -0 -w $${bam} -o $${bam%.bam}.test 2> $${bam%.bam}.log; \
		diff $${bam%.bam}.test.hist.txt $${bam%.bam}.expected.hist.txt; \
		bash -c "diff <(sed 1d $${bam%.bam}.test.dupstat.txt) <(sed 1d $${bam%.bam}.expected.dupstat.txt)"; \
	done
	for bam in $(testsmallbams);do \
		./decluster $${bam} -r 42 -o $${bam%.bam}_bamin.test 2> $${bam%.bam}_bamin.log; \
		./decluster $${bam%.bam}_bamin.test.noClusterDuplicates.bam -r 42 -o $${bam%.bam}_ncld_bamin.test 2> $${bam%.bam}_ncld_bamin.log; \
		./decluster -H $${bam%.bam}_bamin.hist.txt -r 42 -o $${bam%.bam}_histin.test 2> $${bam%.bam}_histin.log; \
		diff $${bam%.bam}_bamin.test.table.txt $${bam%.bam}_bamin.table.expected.txt; \
		diff $${bam%.bam}_histin.test.table.txt $${bam%.bam}_histin.table.expected.txt; \
		diff $${bam%.bam}_ncld_bamin.test.table.txt $${bam%.bam}_ncld_bamin.table.expected.txt; \
		diff $${bam%.bam}_bamin.test.table_defect.txt $${bam%.bam}_bamin.table_defect.expected.txt; \
		diff $${bam%.bam}_histin.test.table_defect.txt $${bam%.bam}_histin.table_defect.expected.txt; \
		diff $${bam%.bam}_ncld_bamin.test.table_defect.txt $${bam%.bam}_ncld_bamin.table_defect.expected.txt; \
		diff $${bam%.bam}_bamin.test.table.txt $${bam%.bam}_histin.test.table.txt; \
		diff $${bam%.bam}_bamin.test.table.txt $${bam%.bam}_ncld_bamin.test.table.txt; \
		diff $${bam%.bam}_histin.test.table.txt $${bam%.bam}_ncld_bamin.test.table.txt; \
		diff $${bam%.bam}_bamin.test.table_defect.txt $${bam%.bam}_histin.test.table_defect.txt; \
		diff $${bam%.bam}_bamin.test.table_defect.txt $${bam%.bam}_ncld_bamin.test.table_defect.txt; \
		diff $${bam%.bam}_histin.test.table_defect.txt $${bam%.bam}_ncld_bamin.test.table_defect.txt; \
	done

	#addtestaux

cleantest:
	rm -v tests/*.test* tests/*.log tests/*bamin.test* tests/*histin.test*
