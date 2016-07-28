# The location of the expat directory
CC=gcc  
LEX=lex  
LEXFLAGS=-lfl
CFLAGS=-C

TARS=RAxML-7.1.0.tar hmmer-2.3.2.tar prodigal.v2_00.tar rpkm.tar trnascan.tar wise2.2.0.tar
SUBDIRS=RAxML-7.1.0 hmmer-2.3.2 prodigal.v2_00 rpkm_c trnascan wise2.2.0


buildall:
	$(MAKE) $(CFLAGS) trnascan 
	$(MAKE) $(CFLAGS) rpkm_c 
	$(MAKE) $(CFLAGS) wise2.2.0/src all
	$(MAKE) $(CFLAGS) RAxML-7.1.0 -f Makefile.gcc all
	$(MAKE) $(CFLAGS) prodigal.v2_00 all
	#$(MAKE) $(CFLAGS) hmmer-2.3.2

cleanall:
	$(MAKE) $(CFLAGS) trnascan clean
	$(MAKE) $(CFLAGS) rpkm_c clean
	$(MAKE) $(CFLAGS) wise2.2.0/src clean
	$(MAKE) $(CFLAGS) RAxML-7.1.0 -f Makefile.gcc clean
	$(MAKE) $(CFLAGS) prodigal.v2_00 clean
	#$(MAKE) $(CFLAGS) hmmer-2.3.2 clean

unzip: $(TARS)
	$(foreach tar, $(TARS), tar -xvf $(tar);)

cleanunzip:
	rm -rf $(SUBDIRS)
