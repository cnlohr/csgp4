all : test

checkProg : checkProg.c csgp4.h
	gcc -g -Os -flto -o $@ $< -lm
	objdump -S $@ > $@.lst

checkProg.float : checkProg.c
	gcc -g -Os -flto -o $@ $< -lm -DCSGP4_USE_FLOAT=1

checkProgSimple : checkProgSimple.c
	gcc -g -Os -flto -o $@ $< -lm -DCSGP4_USE_FLOAT=1 -pedantic -Wall

trackonly : trackonly.c csgp4.h
	gcc -g -Os -flto -o $@ $< -lm

spacestations.txt : 
	wget "https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle" -O spacestations.txt

test : checkProg spacestations.txt trackonly checkProg.float checkProgSimple
	./checkProg spacestations.txt
	./checkProg.float spacestations.txt
	./checkProgSimple

clean :
	rm -rf *.o *~ checkProg trackonly checkProg.float checkProgSimple


