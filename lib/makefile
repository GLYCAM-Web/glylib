#For updating or adding to the library:
SOURCE = ../src/*/*.c
GLYINC = ../inc
default: $(SOURCE)
	gcc -Wall -I$(GLYINC) -c $(SOURCE) -g
	ar rvs libglylib.a *.o
	rm -f *.o
	g++ -Wall -I$(GLYINC) -c $(SOURCE) -g
	ar rvs libglylib++.a *.o
	rm -f *.o
justC: $(SOURCE)
	gcc -Wall -I$(GLYINC) -c $(SOURCE) -g
	ar rvs libglylib.a *.o
	rm -f *.o
justCpp: $(SOURCE)
	g++ -Wall -I$(GLYINC) -c $(SOURCE) -g
	ar rvs libglylib++.a *.o
	rm -f *.o
clean:	
	rm *.o
	rm *.a
