#daid 04 Nov 2011 14:18:58  

CXX=g++
COPT=-O2 -march=native -pipe -fomit-frame-pointer
ROOTFLAGS=$(shell root-config --cflags)
CXXFLAGS=$(COPT) -Wall -Werror -fPIC $(ROOTFLAGS)
SRC=$(wildcard ./*.C)
ROOTLIBS=$(shell root-config --libs)
LIBS=$(ROOTLIBS)
OBJS=$(SRC:%.C=%.o)
TGT=tb

$(TGT): $(OBJS)
	@echo "CXX $@ $<"
	@$(CXX) -o $@ $(LIBS) $(CXXFLAGS) $<
	@echo
	@echo "binary: $(shell pwd)/$(TGT)"

%.o: %.C 
	@echo "CXX $<"
	@$(CXX) -c $(CXXFLAGS) $<

clean: 
	@echo "cleaning..."
	@rm -f $(TGT) $(OBJS)
