all : test

checkProg : checkProg.c sattrack.h
	gcc -g -Os -flto -o $@ $< -lm
	objdump -S $@ > $@.lst

trackonly : trackonly.c
	gcc -g -Os -flto -o $@ $< -lm

spacestations.txt : 
	wget "https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle" -O spacestations.txt

test : checkProg spacestations.txt trackonly
	./checkProg spacestations.txt


clean :
	rm -rf *.o *~ checkProg trackonly


