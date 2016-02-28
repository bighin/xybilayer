TARGET = xypcc
LIBS = 
CC = gcc-mp-5
CFLAGS = -O2 -Wall -ansi -fopenmp -std=c11
LDFLAGS = -lncurses

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c libprogressbar/*.c))
HEADERS = $(wildcard *.h libprogressbar/*.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o libprogressbar/*.o
	-rm -f $(TARGET)
