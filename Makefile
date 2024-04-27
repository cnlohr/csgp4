all : test

checkProg : checkProg.c sattrack.h
	gcc -g -Os -o $@ $< -lm
	objdump -S $@ > $@.lst


spacestations.txt : 
	wget "https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle" -O spacestations.txt

test : checkProg spacestations.txt
	./checkProg spacestations.txt


clean :
	rm -rf *.o *~ checkProg


