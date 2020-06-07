SUBDIRS	= src/surftomo src/Joint src/gravity 
.PHONY: subdirs $(SUBDIRS)
#$(SUBDIRS):
#	$(MAKE) -C $@

subdirs: $(SUBDIRS)
	@for dir in ${SUBDIRS}; do \
		${MAKE} -C $$dir; \
	done

clean:
	@for dir in ${SUBDIRS}; do \
		${MAKE} -C $$dir clean; \
	done
	rm src/*.o