
bin_PROGRAMS = treemix threepop fourpop f4ratio
DISTCHECK_CONFIGURE_FLAGS=LDFLAGS=-L/opt/local/lib CXXFLAGS=-I/opt/local/include
treemix_SOURCES = TreeMix.cpp gzstream.cpp PopGraph.cpp nnls.cpp PhyloPop_params.cpp CountData.cpp CmdLine.cpp GraphState2.cpp

threepop_SOURCES = three_pop.cpp gzstream.cpp PopGraph.cpp nnls.cpp PhyloPop_params.cpp CountData.cpp CmdLine.cpp GraphState2.cpp
fourpop_SOURCES = four_pop.cpp gzstream.cpp PopGraph.cpp nnls.cpp PhyloPop_params.cpp CountData.cpp CmdLine.cpp GraphState2.cpp
f4ratio_SOURCES = f4_ratio.cpp gzstream.cpp PopGraph.cpp nnls.cpp PhyloPop_params.cpp CountData.cpp CmdLine.cpp GraphState2.cpp

noinst_HEADERS = gzstream.h PopGraph.h nnls.h Settings.hpp PhyloPop_params.h CountData.h CmdLine.h GraphState2.h plotting_funcs.R 
