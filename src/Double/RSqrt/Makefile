CC      = clang
CFLAGS  = -Wall -Wextra -g \
           -DTABLE -DTESTCOUNT=2 -DFRACTTABLE=25 -DFRACBITS=58 \
					 -DDEBUG \
           -I../../Single/RSqrt
LDFLAGS = -L../../Single/RSqrt -lrsqrt_sp -lm

TARGET  = test
OBJS    = float_utils.o rsqrt_dp.o test.o

# Path to dependency library
DEPS_DIR = ../../Single/RSqrt
LIB      = $(DEPS_DIR)/librsqrt_sp.a

# Default target
all: $(TARGET) output clean

$(TARGET): $(OBJS) $(LIB)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

$(LIB):
	$(MAKE) -C $(DEPS_DIR) lib

# Object file rules
rsqrt_dp.o: rsqrt_dp.c rsqrt_dp.h $(DEPS_DIR)/rsqrt.h $(DEPS_DIR)/float_utils.h
	$(CC) $(CFLAGS) -c $< -o $@

test.o: test.c rsqrt_dp.h
	$(CC) $(CFLAGS) -c $< -o $@

# Run program
output: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
	$(MAKE) -C $(DEPS_DIR) clean