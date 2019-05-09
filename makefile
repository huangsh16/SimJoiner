obj = $(wildcard *.cpp)

main.exe : $(obj)
	g++ -std=c++11 -O3 -o main.exe $(obj)