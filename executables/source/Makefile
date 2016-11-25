# The location of the expat directory
CC=gcc  
LEX=lex  
LEXFLAGS=-lfl
CFLAGS=-C

TARS= prodigal.v2_00.tar rpkm.tar trnascan.tar LAST-Plus.tar bwa-0.7.7.tar
SUBDIRS= prodigal.v2_00 rpkm trnascan  LAST-Plus bwa-0.7.7




all: untar_folders clean_folders  build_folders install 


.PHONY: build_folders
build_folders:
	$(MAKE) $(CFLAGS) trnascan 
	$(MAKE) $(CFLAGS) rpkm 
	$(MAKE) $(CFLAGS) prodigal.v2_00 all
	$(MAKE) $(CFLAGS) LAST-Plus
	$(MAKE) $(CFLAGS) bwa-0.7.7

.PHONY: clean_folders
clean_folders:
	$(MAKE) $(CFLAGS) trnascan clean
	$(MAKE) $(CFLAGS) rpkm clean
	$(MAKE) $(CFLAGS) prodigal.v2_00 clean
	$(MAKE) $(CFLAGS) LAST-Plus clean
	$(MAKE) $(CFLAGS) bwa-0.7.7 clean

.PHONY: untar_folders
untar_folders: $(TARS)
	$(foreach tar, $(TARS), tar -xvf $(tar);)

tar_folders: $(SUBDIRS)
	$(foreach folder, $(SUBDIRS), tar -cvvf $(folder).tar $(folder);)

remove_folders:
	rm -rf $(SUBDIRS)

install:
	cp  trnascan/trnascan-1.4  ../${TARGET}/
	cp  LAST-Plus/lastal+  ../${TARGET}/ 
	cp  LAST-Plus/lastdb+  ../${TARGET}/ 
	cp  bwa-0.7.7/bwa  ../${TARGET}/ 
	cp  prodigal.v2_00/prodigal  ../${TARGET}/ 
	cp  rpkm/rpkm  ../${TARGET}/ 
