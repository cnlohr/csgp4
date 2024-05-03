all : test

checkProg : checkProg.c csgp4.h
	gcc -g -Os -flto -o $@ $< -lm
	objdump -S $@ > $@.lst

checkProg.float : checkProg.c
	gcc -g -Os -flto -o $@ $< -lm -DCSGP4_USE_FLOAT=1
	objdump -S $@ > $@.lst

checkProgSimple : checkProgSimple.c
	gcc -g -Os -flto -o $@ $< -lm -DCSGP4_USE_FLOAT=1 -pedantic -Wall
	objdump -S $@ > $@.lst

trackonly : trackonly.c csgp4.h
	gcc -g -Os -flto -o $@ $< -lm

spacestations.txt : 
	wget "https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle" -O spacestations.txt

test : checkProg spacestations.txt trackonly checkProg.float checkProgSimple
	./checkProg spacestations.txt
	./checkProg.float spacestations.txt
	./checkProgSimple
	size checkProg checkProg.float checkProgSimple

clean :
	rm -rf *.o *~ checkProg trackonly checkProg.float checkProgSimple


