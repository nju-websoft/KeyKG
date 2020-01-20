SUBDIR = ./bcPLL ./sspexp_codes

all:$(SUBDIR)
	@echo make end!

$(SUBDIR):ECHO
	$(MAKE) -C $@

ECHO:  
	@echo $(SUBDIRS)
	@echo begin compile

clean:
	rm -rf ./bin
