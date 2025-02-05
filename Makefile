all : test

ifeq ($(GITHUB_TOKEN),)
define GH_ADDSTATUS
	true
endef
else
define GH_ADDSTATUS
	( (cat status.txt | tr '\n' '<br>' > content.txt); (echo '{"state":$(2),"context":$(1)}' | jq  --rawfile content content.txt '."description" |= $$content' > payload.txt); \
	curl -L \
		-X POST \
		-H "Accept: application/vnd.github+json" \
		-H "Authorization: Bearer ${GITHUB_TOKEN}" \
		-H "X-GitHub-Api-Version: 2022-11-28" \
		https://api.github.com/repos/${GITHUB_REPOSITORY}/statuses/${GITHUB_WORKFLOW_SHA} \
		--data-binary payload.txt; cat payload.txt)
endef
# Need context and state: (error, failure, pending, success)
# Can have target_url, description, also btw --silent --output /dev/null 
endif

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
	((./checkProg spacestations.txt | tee status.txt) && $(call GH_ADDSTATUS,"Double precision","success")) || $(call GH_ADDSTATUS,"Double precision","failure")
	((./checkProg.float spacestations.txt | tee status.txt) && $(call GH_ADDSTATUS,"Single precision","success")) || $(call GH_ADDSTATUS,"Single precision","failure")
	((./checkProgSimple | tee status.txt) && $(call GH_ADDSTATUS,"Simple Test","success")) || $(call GH_ADDSTATUS,"Simple Test","failure")
	size checkProg checkProg.float checkProgSimple

clean :
	rm -rf *.o *~ checkProg trackonly checkProg.float checkProgSimple


