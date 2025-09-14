CC      = clang
CFLAGS  = -DRN -DSINGLE -DTABLE -DFRACTTABLE=25 -DFRACBITS=27 -g
LDFLAGS = -lm
AR      = ar
ARFLAGS = rcs

LIBOBJS = dtof_nearest.o dtof_zero.o dtof_pinf.o dtof_ninf.o rsqrt.o float_utils.o
LIB     = librsqrt_sp.a

TESTOBJS = test.o

TARGET  = test

.PHONY: all lib test run clean

all: lib test run clean

# Build static library
lib: $(LIB)

$(LIB): $(LIBOBJS)
	$(AR) $(ARFLAGS) $@ $^

test: $(LIBOBJS) $(TESTOBJS)
	$(CC) -o $(TARGET) $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Run local test
run: test
	./$(TARGET)

getdate:
	echo "date $$(date "+%Y-%m-%d")" > test.log

clean:
	rm -f $(LIBOBJS) $(TESTOBJS) $(LIB) $(TARGET)