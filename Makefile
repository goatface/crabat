# Copyright: 
# 2010, 2012, 2013 daid KAHL
#
# This file is part of crabat
#
# crabat is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# crabat is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with crabat.  If not, see <http://www.gnu.org/licenses/>.
#  
# In order to publish results obtained from crabat or any
# future derivative software, the authors must make reference
# to the paper describing it (in preparation): D. Kahl et al.
# ``Algorithms for Active Target Data Analysis'' 2013.  This is an 
# additional term to the license as stipulated in section 7b.
# Publication includes, but is not limited to: books, peer-
# reviewed journals, conference proceedings, annual reports, etc.
# 

# Makefile for crabat
# Dependencies:
# ROOT (http://root.cern.ch)
# KaliVeda (http://indra.in2p3.fr/KaliVedaDoc/)
# Doxygen (http://www.stack.nl/~dimitri/doxygen/)
CXX=g++
#COPT=-O0 -g # valgrind
#COPT=-O2 -march=native -pipe #debug
#COPT=-DDEBUG -I"" -O0 -g3 -Wall -c -fmessage-length=0 # another debug
#CDEBUG=-ansi -pendantic -g # some problems? 
COPT=-O2 -march=native -pipe -fomit-frame-pointer # run fastest
ROOTFLAGS=$(shell root-config --cflags)
CXXFLAGS=$(COPT) -Wall -Werror -fPIC $(ROOTFLAGS)
FC=gfortran
FFLAGS=-lstdc++
FSRC=$(wildcard ./enewzsub/*.f)
LD=g++ $(FFLAGS)
ROOTLIBS=$(shell root-config --libs)
#KVLIBS=-L$(KVROOT)/lib -lKVMultiDet -lrange
RANGELIBS=-L/usr/local/lib -lrange
LIBS=-lXpm -lX11 -lm $(ROOTLIBS) $(KVLIBS) $(RANGELIBS) -lstdc++
#INCLUDE=$(addprefix -I,$(KVROOT)/include)
COBJS=dictCal.o dictAnaly.o run.o
#FOBJS=$(FSRC:%.f=%.o)
#ARCHIVE=libenewzlib.a
DOCCMD=doxygen
DOCCFG=doc/Doxyfile
DOCS=doc/html/index.html
TGT=run

all: tb $(ARCHIVE) $(TGT) $(DOCS) 

tb: force_look 
	@printf "make tbrowser\n"
	@cd tbrowser; ${MAKE} ${MFLAGS}

$(ARCHIVE):  
	@ar r $(ARCHIVE)
	@ranlib $(ARCHIVE)

$(TGT): $(FOBJS) $(COBJS)
	@printf "LD $(TGT)\n\tcopy of any errors recorded in log/linker.log\n"
	$(LD) -O  $(COBJS) $(FOBJS) $(ARCHIVE) -o $@ $(LIBS) 2>&1 | tee log/linker.log 
	@printf "\nCompile successful!\nbinary: $(shell pwd)/$(TGT)\n\n"

$(DOCS): $(TGT) 
	@printf "Generating documentation...\n\tall output and errors directed to log/doxygen.log\n"
	@$(DOCCMD) $(DOCCFG)  &> log/doxygen.log 

force_look:
	@true

clean:	
	@echo "cleaning..."
	@rm -f $(TGT) $(ARCHIVE) dict*
	@rm -rf doc/html doc/latex
	@find . -maxdepth 2 -name "*.o" | xargs -I{} rm {}

%.o: %.cxx Analyzer.cxx Analyzer_config.cxx run.h
	@echo "CXX $< $@"
	@$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ 

%.o: %.f
	@echo "FC $@"
	@$(FC) -c $< -o $@ 

dictCal.cxx: Calibration.h linkdefCal.h Calibration.cxx
	@echo "ROOT $@"
	@rootcint -f $@ -c $(CXXFLAGS) $< linkdefCal.h

dictAnaly.cxx: Analyzer.h linkdefAnaly.h Analyzer.cxx Analyzer_config.cxx
	@echo "ROOT $@"
	@rootcint -f $@ -c $(CXXFLAGS) $< linkdefAnaly.h
