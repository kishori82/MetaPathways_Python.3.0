ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

all: 
	make -C bwa
	make -C FAST
	make -C prodigal
	make -C rpkm
	make -C trnascan

.PHOHY: clean
clean:
	make -C bwa clean
	make -C FAST clean
	make -C prodigal clean
	make -C rpkm clean
	make -C trnascan clean

.PHONY: install
install: bwa/bwa
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	install bwa/bwa $(DESTDIR)$(PREFIX)/bin/
#
	install FAST/fastal $(DESTDIR)$(PREFIX)/bin/
	install FAST/fastdb $(DESTDIR)$(PREFIX)/bin/
#
	install prodigal/prodigal $(DESTDIR)$(PREFIX)/bin/
#
	install rpkm/rpkm $(DESTDIR)$(PREFIX)/bin/
#
	install trnascan/trnascan-1.4 $(DESTDIR)$(PREFIX)/bin/

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/bwa
	rm -f $(DESTDIR)$(PREFIX)/bin/fastal  $(DESTDIR)$(PREFIX)/bin/fastdb
	rm -f $(DESTDIR)$(PREFIX)/bin/prodigal
	rm -f $(DESTDIR)$(PREFIX)/bin/rpkm
	rm -f $(DESTDIR)$(PREFIX)/bin/trnascan-1.4
