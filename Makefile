CC=g++ -g
TARGET=ana
OBJS=myevent.o ana.o
HEADERS=myevent.h
ROOTLIBS=`root-config --libs`

TTARGET=tracking
TOBJS=basictracking.o track.o
THEADERS=basictracking.h

.SUFFIXES :.o.cc
.SUFFIXES :.o.C

.cc.o:
	$(CC) $(FFLAGS) `root-config --cflags` -c $<

.C.o:
	$(CC) $(FFLAGS) `root-config --cflags` -c $<

ana:$(OBJS) $(HEADERS)
	$(CC) $(CPPFLAGS) -o $(TARGET) $(OBJS) $(ROOTLIBS)

track:$(TOBJS) $(THEADERS)
	$(CC) $(CPPFLAGS) -o $(TTARGET) $(TOBJS) $(ROOTLIBS) -lMinuit


clean:
	rm -rf *.o *.d $(TARGET) $(TTARGET)

