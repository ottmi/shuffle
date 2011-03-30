BIN = shuffle 

all: $(BIN)

$(BIN): force_look
	cd src; $(MAKE) $(MFLAGS)
	mv src/$(BIN) .

clean :
	$(RM) $(BIN)
	cd src; $(MAKE) clean

force_look:
	true
