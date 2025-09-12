CC      = clang
CFLAGS  = -Wall -Wextra -g \
           -DTABLE -DTESTCOUNT=8388608 -DFRACTTABLE=25 -DFRACBITS=58 \
					 -DDEBUG \
           -I../../Single/Recip
LDFLAGS = -L../../Single/Recip -lrecip_sp -lm

TARGET  = test
OBJS    = float_utils.o recip_dp.o test.o

# Path to dependency library
DEPS_DIR = ../../Single/Recip
LIB      = $(DEPS_DIR)/librecip_sp.a

# Default target
all: $(TARGET) output clean

$(TARGET): $(OBJS) $(LIB)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

$(LIB):
	$(MAKE) -C $(DEPS_DIR) lib

# Object file rules
recip_dp.o: recip_dp.c recip_dp.h $(DEPS_DIR)/recip.h $(DEPS_DIR)/float_utils.h
	$(CC) $(CFLAGS) -c $< -o $@

test.o: test.c recip_dp.h
	$(CC) $(CFLAGS) -c $< -o $@

# Run program
output: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
	$(MAKE) -C $(DEPS_DIR) clean