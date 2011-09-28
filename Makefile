BIN = shuffle 
VERSION := $(shell awk '/VERSION/ {print $$3}' src/globals.h| sed 's/\"\(.*\)\"/\1/')

all: $(BIN)

$(BIN): force_look
	cd src; $(MAKE) $(MFLAGS)
	mv src/$(BIN) .

clean :
	$(RM) $(BIN)
	cd src; $(MAKE) clean

tarball:
	git archive --format=tar --prefix=shuffle-$(VERSION)/ HEAD:src|bzip2 >shuffle-$(VERSION).tar.bz2

force_look:
	true
