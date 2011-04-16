#For updating or adding to the library:
SOURCE = ../src/*.c ../src/*/*.c
GLYINC = ../inc
default: $(SOURCE)
	gcc -c -Wall -I$(GLYINC) $(SOURCE) -g
	ar rvs libglylib.a *.o
clean:	
	rm *.o
	rm *.a
