CC      = gcc
CFLAGS  = -DRN -DSINGLE -DTABLE -DFTABLE=makePosRecipTable \
          -DFLINEARTABLE=makelinearTable -DFQTABLE=makeqTable -DFLIE=file_name \
          -DFRACTTABLE=25 -DFRACBITS=31 -g
LDFLAGS = -lm
AR      = ar
ARFLAGS = rcs

LIBOBJS = dtof_nearest.o dtof_zero.o dtof_pinf.o dtof_ninf.o recip.o float_utils.o
LIB     = librecip_sp.a

# Objects used only for Single/Recip's own testing
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