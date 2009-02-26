# 
#  Copyright (c) 2007, Pau Fernández
#

CXXFLAGS+=-Wall -I.. -O2 # -DNDEBUG
LDFLAGS=-lstdc++ 

NOGRAPHICS=Gconvert Gstats Gdstats Gpnstats Gmanip \
	   Gtopovrlp Gbtwns Gcommunity Gmodmeas Ggen

GRAPHICS=Gmaslov Ghierclust Gtreeanlz 

all: utils
utils: $(NOGRAPHICS) $(GRAPHICS)

convert.o: convert.C io.H
hierclust.o: hierclust.C hierclust.H
gen.o: gen.C gen.H

$(NOGRAPHICS): G%: %.o
	gcc -o $@ $(subst G,,$@).o $(LDFLAGS)

$(GRAPHICS): G%: %.o
	gcc -o $@ $(subst G,,$@).o $(LDFLAGS) -L../graphics -lgraphics

clean:
	rm -rf *.o $(GRAPHICS) $(NOGRAPHICS)
