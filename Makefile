BIN = shuffle 

all: $(BIN)

$(BIN): force_look
	cd src; $(MAKE) $(MFLAGS)
	mv src/$(BIN) .

clean :
	$(RM) $(BIN)
	cd src; $(MAKE) cleana

force_look:
	true