all: find90

find90: find90.cpp
	g++ `root-config --libs --cflags` -Wall -Wextra -Werror find90.cpp -o find90 
