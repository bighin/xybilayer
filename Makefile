TARGET = xypcc
LIBS = 
CC = gcc-mp-12
#CC = clang-mp-3.8
CFLAGS = -O2 -Wall -fopenmp -std=c11
#CFLAGS = -g -pg
LDFLAGS = -lncurses

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c libprogressbar/*.c thr_pool/*.c))
HEADERS = $(wildcard *.h libprogressbar/*.h thr_pool/*.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o libprogressbar/*.o thr_pool/*.o
	-rm -f $(TARGET)
