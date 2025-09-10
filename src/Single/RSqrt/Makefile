CC=clang
CFLAGS=-DRN -DSINGLE -DTABLE -DFRACTTABLE=25 -DFRACBITS=27 -lm -g
OBJS=dtof_nearest.o dtof_zero.o dtof_pinf.o dtof_ninf.o

.PHONY: clean

all: makedtof output clean

output: test
	./a

test: test.c
	$(CC) -g -o a rsqrt.c test.c float_utils.c $(CFLAGS) -DGENQTABLE $(OBJS)

makedtof:
	$(CC) -c -frounding-math dtof_nearest.c dtof_zero.c dtof_pinf.c dtof_ninf.c

getdate:
	echo "date $$(date "+%Y-%m-%d")" > test.log

clean:
	rm -f $(OBJS) *.o a