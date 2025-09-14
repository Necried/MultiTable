CC      = clang
CFLAGS  = -Wall -Wextra -g \
           -DTABLE -DTESTCOUNT=2 -DFRACTTABLE=25 -DFRACBITS=58 \
					 -DDEBUG \
           -I../../Single/Sqrt
LDFLAGS = -L../../Single/Sqrt -lsqrt_sp -lm

TARGET  = test
OBJS    = float_utils.o sqrt_dp.o test.o

# Path to dependency library
DEPS_DIR = ../../Single/Sqrt
LIB      = $(DEPS_DIR)/libsqrt_sp.a

# Default target
all: $(TARGET) output clean

$(TARGET): $(OBJS) $(LIB)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

$(LIB):
	$(MAKE) -C $(DEPS_DIR) lib

# Object file rules
sqrt_dp.o: sqrt_dp.c sqrt_dp.h $(DEPS_DIR)/sqrt.h $(DEPS_DIR)/float_utils.h
	$(CC) $(CFLAGS) -c $< -o $@

test.o: test.c sqrt_dp.h
	$(CC) $(CFLAGS) -c $< -o $@

# Run program
output: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
	$(MAKE) -C $(DEPS_DIR) clean