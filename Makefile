CXX	= g++


ROOTCFLAGS   	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)

#FORMAT: MC_2016, MC_2017, MC_2018, DATA_2016, DATA_2017, DATA_2018
#SUBFORMAT: MC_DUMMY, DATA_2016B, ..., DATA_2017B, ..., DATA_2018B, ...
#INPUT: TFILE,TCHAIN
#NANOAOD: NANOAODV7, NANOAODV9

CXXFLAGS	= -O -Wall -fPIC -D $(FORMAT) -D $(SUBFORMAT) -D $(NANOAOD) -D $(INPUT)
CXXFLAGS	+= $(ROOTCFLAGS)	

SOFLAGS		= -O -shared
SOFLAGS 	+= $(ROOTGLIBS)

#main:	Ana.cxx libReader.so libReader_dict.so
#	$(CXX) $(CXXFLAGS) $(ROOTGLIBS) libReader.so libReader_dict.so $< -o $@

#Task: copy mydict main clean
Task: clean mydict main clean

#copy:
#	cp -f src/Reader_$(FORMAT).h  src/Reader.h

mydict:
	@rootcint Reader_dict.cxx -c src/Reader.h
	@rootcint Processor_dict.cxx -c src/Processor.h

main:	Ana.cxx Global.o Reader.o Processor.o BTagCalibrationStandalone.o LumiMaskFilter.o Reader_dict.o Processor_dict.o Selector.o VH_selection.o
	#$(CXX) $(CXXFLAGS) $(ROOTGLIBS) $^ -o $@
	$(CXX) $(CXXFLAGS) $(ROOTGLIBS) $^ -o $@ yaml-cpp/build/libyaml-cpp.a

#main:	Ana.cxx Global.o Reader.o Processor.o Reader_dict.o Processor_dict.o Selector.o ZbSelection.o
#	$(CXX) $(CXXFLAGS) $(ROOTGLIBS) $^ -o $@

#lib%.so: %.o
#	$(CXX) $(SOFLAGS) -o $@ $^

#Ana.o: Ana.cxx
#	$(CXX) $(CXXFLAGS) -c $^ -o $@

#StdArg.o: StdArg.hpp
#	$(CXX) $(CXXFLAGS) -c $^ -o $@

%.o:	src/%.cxx
	$(CXX) $(CXXFLAGS) -c $^ -o $@

%.o:	src/%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

%.o:	src/%.cc
	#$(CXX) $(CXXFLAGS) -I yaml-cpp/include -lyaml-cpp -c $^ -o $@
	$(CXX) $(CXXFLAGS) -I yaml-cpp/include -c $^ -o $@

%.o:	%.cxx
	$(CXX) $(CXXFLAGS) -c $^ -o $@

%.o:	src/%.h
	$(CXX) $(CXXFLAGS) -c $^ -o $@

%.o:	%.hpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

Reader_dict.o: Reader_dict.cxx
	$(CXX) $(CXXFLAGS) -c $^ -o $@

Processor_dict.o: Processor_dict.cxx
	$(CXX) $(CXXFLAGS) -c $^ -o $@

#LumiMaskFilter1.o: LumiMaskFilter.cc
#	$(CXX) $(CXXFLAGS) -I yaml-cpp/include -lyaml-cpp -c $^ -o $@

clean:
	rm -rf main *.o *.so *_dict.cxx
