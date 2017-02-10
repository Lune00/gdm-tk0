# root makefile. Delegate to source subdirs.

PACKAGE=.

SUBDIRS	= disk2cluster disk2polyg his2mgp spl2cin 

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
		

		
