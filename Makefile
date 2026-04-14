CRYSTAL   = crystal
FLAGS     = -Dpreview_mt -Dexecution_context
BINDIR    = bin

.PHONY: all minimap2 paftools spec clean

all: minimap2 paftools

minimap2:
	$(CRYSTAL) build $(FLAGS) src/main.cr -o $(BINDIR)/minimap2

paftools:
	$(CRYSTAL) build $(FLAGS) src/paftools_main.cr -o $(BINDIR)/paftools

spec:
	$(CRYSTAL) spec $(FLAGS)

clean:
	rm -f $(BINDIR)/minimap2 $(BINDIR)/paftools
