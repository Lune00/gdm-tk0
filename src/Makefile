# root makefile. Delegate to source subdirs.

PACKAGE=.

SUBDIRS	= algo struct struct/Bodies2D struct/Inter2D struct/CohesionLaw analyse systems systems_A processing 

all:
	@@for p in $(SUBDIRS); do \
		echo '============= building ' $(PACKAGE)/$$p; \
		make -C $(PACKAGE)/$$p --no-print-directory  all; \
	done

clean:
	@@for p in $(SUBDIRS); do \
		echo '============= cleaning ' $(PACKAGE)/$$p; \
		make -C $(PACKAGE)/$$p --no-print-directory clean;\
	done
	
	
depend:
		@@for p in $(SUBDIRS); do \
			echo '============= depend ' $(PACKAGE)/$$p; \
			make -C $(PACKAGE)/$$p --no-print-directory depend;\
		done
		
cleanlib:
	@@echo '============= Clean lib ';\
	cd ../lib/;\
	rm *.o;\
	cd ../src/;
		
		
