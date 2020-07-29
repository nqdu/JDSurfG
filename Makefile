SUBDIRS	= src/surftomo src/gravity src/JointTomo
.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)
	@for dir in ${SUBDIRS}; do \
		${MAKE} -C $$dir; \
	done

clean:
	@for dir in ${SUBDIRS}; do \
		${MAKE} -C $$dir clean; \
	done