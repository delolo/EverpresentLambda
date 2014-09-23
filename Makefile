test: src/test.cpp
	g++ -o test.out src/test.cpp include/Tools.cpp

main: src/main.cpp
	g++ -o main.out src/main.cpp include/tools/random.cpp