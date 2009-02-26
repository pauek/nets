# 
#  Copyright (c) 2007, Pau Fernández
#

OBJECTS=postscript.o graph.o
CXXFLAGS+=-Wall -I.. -O2 -DNDEBUG 

all: lib 

lib: $(OBJECTS)
	ar r libgraphics.a $(OBJECTS) 

postscript.o: postscript.C postscript.H
graph.o: graph.C graph.H

clean:
	rm -rf *.o libgraphics.a
