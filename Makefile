CC = g++
CFLAGS = -lm -lgsl -lgslcblas

.PHONY: all
all: main.exe

.PHONY: debug
debug: main-g.exe

main.exe: main.cpp  update.cpp
	$(CC) -O3 -Wall  main.cpp  update.cpp -o main.exe $(CFLAGS) -Dfolder="\"results"\"

main-g.exe: main.cpp   update.cpp
	$(CC) -O3 -g -Wall  main.cpp  update.cpp -o main-g.exe $(CFLAGS) -Dfolder="\"results"\"


.PHONY: clean
clean:
	rm -rf *exe

.PHONY: install
install:
	mkdir -p ./results/trace
	cp *.cpp ./results
	cp LICENSE ./results
	cp *.h ./results

