CC=gcc
CPP=g++
CFLAGS=-g -Wall -Werror -Wextra -ansi -pedantic
CPPFLAGS=-g -Wall -Werror -Wextra -ansi -pedantic
LDFLAGS=
LIBS=-lm
SOURCES=$(wildcard *.c)
SOURCESCPP=$(wildcard *.cpp)
OBJECTS=$(SOURCES:%.c=%.o) $(SOURCESCPP:%.cpp=%.o)
TARGET=exercise9

.PHONY=all debug clean

all: $(SOURCES) $(TARGET)

debug: CC += -DDEBUG -DPRINT_DEBUG -g
debug: CPP += -DDEBUG -DPRINT_DEBUG -g
debug: $(SOURCES) $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CPP) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
