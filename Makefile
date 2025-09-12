.PHONY: all clean single double test

SINGLE_DIRS := src/Single/Recip src/Single/Sqrt src/Single/RSqrt
DOUBLE_DIRS := src/Double/Recip src/Double/Sqrt src/Double/RSqrt

# Default: build and run tests for both
all: test

single:
	@for dir in $(SINGLE_DIRS); do \
		echo "===> Building $$dir"; \
		$(MAKE) -C $$dir all || exit $$?; \
	done
	
double: single
	@for dir in $(DOUBLE_DIRS); do \
		echo "===> Building $$dir"; \
		$(MAKE) -C $$dir all || exit $$?; \
	done

test: single double

clean:
	@for dir in $(SINGLE_DIRS) $(DOUBLE_DIRS); do \
		echo "===> Cleaning $$dir"; \
		$(MAKE) -C $$dir clean || exit $$?; \
	done
