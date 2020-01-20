#CXXFLAGS = -std=c++11 -Wl,--no-as-needed -pthread -w -O3
CXXFLAGS = -std=c++11 -w -O3
OUTPUT=../bin/
PROGRAMS=$(wildcard *.cpp *.h)

all:../bin/bcPLL

../bin/bcPLL:$(PROGRAMS) 
	g++ $(CXXFLAGS) -o $@ $^

clean:
	rm
