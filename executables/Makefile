# The location of the expat directory
CC=gcc  
LEX=lex  
LEXFLAGS=-lfl
CFLAGS=-C

NCBI_BLAST=ncbi-blast-2.5.0+-x64-linux.tar.gz
NCBI_BLAST_VER=ncbi-blast-2.5.0+
BINARY_FOLDER=executables/linux

BLASTP=$(BINARY_FOLDER)/blastp
RPKM=$(BINARY_FOLDER)/rpkm
BWA=$(BINARY_FOLDER)/bwa
TRNASCAN=$(BINARY_FOLDER)/trnascan-1.4
LASTPlus=$(BINARY_FOLDER)/lastal+
PRODIGAL=$(BINARY_FOLDER)/prodigal
MICROBE_CENSUS=microbe_census
METAPATHWAYS_DB=MetaPathways_DBs
METAPATHWAYS_DB_TAG=Metapathways_DBs_2016-04.tar.xz



all: $(BINARY_FOLDER) $(PRODIGAL)  $(LASTPlus)  $(BWA) $(TRNASCAN)  $(RPKM) $(BLASTP) $(MICROBE_CENSUS) $(METAPATHWAYS_DB)

$(TRNASCAN):  
	$(MAKE) $(CFLAGS) executables/source/trnascan 
	mv executables/source/trnascan/trnascan-1.4 $(BINARY_FOLDER)/

$(RPKM):  
	$(MAKE) $(CFLAGS) executables/source/rpkm 
	mv executables/source/rpkm/rpkm $(BINARY_FOLDER)/

$(BWA):  
	$(MAKE) $(CFLAGS) executables/source/bwa 
	mv executables/source/bwa/bwa $(BINARY_FOLDER)/

$(PRODIGAL):  
	$(MAKE) $(CFLAGS) executables/source/prodigal 
	mv executables/source/prodigal/prodigal $(BINARY_FOLDER)/

$(LASTPlus):  
	$(MAKE) $(CFLAGS) executables/source/LAST-Plus
	mv executables/source/LAST-Plus/lastal+ $(BINARY_FOLDER)/
	mv executables/source/LAST-Plus/lastdb+ $(BINARY_FOLDER)/

$(BLASTP): $(NCBI_BLAST) 
	@echo -n "Extracting the binaries for BLAST...." 
	tar --extract --file=$(NCBI_BLAST)  $(NCBI_BLAST_VER)/bin
	mv $(NCBI_BLAST_VER)/bin/blastx  executables/linux/
	mv $(NCBI_BLAST_VER)/bin/blastp  executables/linux/
	mv $(NCBI_BLAST_VER)/bin/blastn  executables/linux/
	mv $(NCBI_BLAST_VER)/bin/makeblastdb  executables/linux/
	rm -rf  $(NCBI_BLAST_VER)
	rm -rf  $(NCBI_BLAST)
	@echo "done" 

$(NCBI_BLAST):
	@echo -n "Downloading BLAST from NCBI website...." 
	wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$(NCBI_BLAST)
	@echo "done" 

.ONESHELL:
$(MICROBE_CENSUS):
	@echo  "Installing MicrobeCensus...." 
	cd executables/source/MicrobeCensus/
	sudo python setup.py install
	@echo "Testing installation, may take a few minutes...." 
	cd tests; python test_microbe_census.py

$(METAPATHWAYS_DB_TAR):
	@echo  "Fetching the databases...." 
	aws s3 cp s3://wbfogdog/a2ac7fc4db0bfae6c05ca12a5818792d/Metapathways_DBs_2016-04.tar.xz .

$(METAPATHWAYS_DB): $(METAPATHWAYS_DB_TAR)
	@echo  "Unzipping the database...." 
	tar -xvf Metapathways_DBs_2016-04.tar.xz 
  

$(BINARY_FOLDER): 
	@if [ ! -d $(BINARY_FOLDER) ]; then mkdir $(BINARY_FOLDER); fi



build_folders:
	$(MAKE) $(CFLAGS) trnascan 
	$(MAKE) $(CFLAGS) rpkm 
	$(MAKE) $(CFLAGS) prodigal.v2_00 all
	$(MAKE) $(CFLAGS) LAST-Plus
	$(MAKE) $(CFLAGS) bwa-0.7.7

clean_folders:
	$(MAKE) $(CFLAGS) trnascan clean
	$(MAKE) $(CFLAGS) rpkm clean
	$(MAKE) $(CFLAGS) prodigal.v2_00 clean
	$(MAKE) $(CFLAGS) LAST-Plus clean
	$(MAKE) $(CFLAGS) bwa-0.7.7 clean

untar_folders: $(TARS)
	$(foreach tar, $(TARS), tar -xvf $(tar);)

tar_folders: $(SUBDIRS)
	$(foreach folder, $(SUBDIRS), tar -cvvf $(folder).tar $(folder);)

remove_folders:
	rm -rf $(SUBDIRS)

install_macosx:
	cp  trnascan/trnascan-1.4  ../macosx/
	cp  LAST-Plus/lastal+  ../macosx/ 
	cp  LAST-Plus/lastdb+  ../macosx/ 
	cp  bwa-0.7.7/bwa  ../macosx/ 
	cp  prodigal.v2_00/prodigal  ../macosx/ 
	cp  rpkm/rpkm  ../macosx/ 

remove_macosx:
	rm -rf  ../macosx/trnascan-1.4 
	rm -rf ../macosx/lastal+  
	rm -rf ../macosx/lastdb+  
	rm -rf ../macosx/bwa  
	rm -rf ../macosx/prodigal
	rm -rf ../macosx/rpkm 

install_redhat:
	cp  trnascan/trnascan-1.4  ../redhat/
	cp  LAST-Plus/lastal+  ../redhat/ 
	cp  LAST-Plus/lastdb+  ../redhat/ 
	cp  bwa-0.7.7/bwa  ../redhat/ 
	cp  prodigal.v2_00/prodigal  ../redhat/ 
	cp  rpkm/rpkm  ../redhat/ 

remove_redhat:
	rm -rf  ../redhat/trnascan-1.4 
	rm -rf ../redhat/lastal+  
	rm -rf ../redhat/lastdb+  
	rm -rf ../redhat/bwa  
	rm -rf ../redhat/prodigal
	rm -rf ../redhat/rpkm 

install_ubuntu:
	cp  trnascan/trnascan-1.4  ../ubuntu/
	cp  LAST-Plus/lastal+  ../ubuntu/ 
	cp  LAST-Plus/lastdb+  ../ubuntu/ 
	cp  bwa-0.7.7/bwa  ../ubuntu/ 
	cp  prodigal.v2_00/prodigal  ../ubuntu/ 
	cp  rpkm/rpkm  ../ubuntu/ 

remove_ubuntu:
	rm -rf  ../ubuntu/trnascan-1.4 
	rm -rf ../ubuntu/lastal+  
	rm -rf ../ubuntu/lastdb+  
	rm -rf ../ubuntu/bwa  
	rm -rf ../ubuntu/prodigal
	rm -rf ../ubuntu/rpkm 
