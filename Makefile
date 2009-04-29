DIRS = bddcmumps bddcpp lobpcg_driver mumps_adaptive multilevel

##################################################################
# Targets
##################################################################

all: 
	@ \
        for i in ${DIRS}; \
        do \
          echo "Making $$i ..."; \
          (cd $$i && $(MAKE)); \
          echo ""; \
        done

clean: 
	@ \
        for i in ${DIRS}; \
        do \
          if [ -d $$i ]; \
          then \
            echo "Cleaning $$i ..."; \
            (cd $$i &&  $(MAKE) $@); \
          fi; \
        done

