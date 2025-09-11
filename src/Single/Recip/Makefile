CC      = clang
CFLAGS  = -DRN -DSINGLE -DTABLE -DFTABLE=makePosRecipTable \
          -DFLINEARTABLE=makelinearTable -DFQTABLE=makeqTable -DFLIE=file_name \
          -DFRACTTABLE=25 -DFRACBITS=31 -g
LDFLAGS = -lm

OBJS = dtof_nearest.o dtof_zero.o dtof_pinf.o dtof_ninf.o test.o recip.o float_utils.o
TARGET  = test

.PHONY: clean

all: makedtof output clean

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

output: $(TARGET)
	./$(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

makedtof:
	$(CC) -c -frounding-math dtof_nearest.c dtof_zero.c dtof_pinf.c dtof_ninf.c

getdate:
	echo "date $$(date "+%Y-%m-%d")" > test.log

clean:
	rm -f $(OBJS) *.o a