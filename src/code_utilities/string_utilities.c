/** \file string_utilities.c Functions for writing variable strings.

Contains utilities that produce numbers in certain formats or that remove
initial and trailing whitespace, etc.

*/
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

#include <mylib.h>

/************* prune_string_whitespace() ***************/
/** Function to prune leading and trailing whitespace from a string 
*/
const char *prune_string_whitespace(const char *I){
int 	i=0, ///< A counter
	len=0, ///< Length (meaning varies)
	LEN=0, ///< Length of incoming string
	start=0, ///< Start point for pruned string
	stop=0; ///< End point for pruned string
char 	*O; ///< Outgoing string

LEN=strlen(I);
start=strspn(I," \t\n\r\f\v"); 
//printf("LEN is %d and start is %d\n",LEN,start);
if(start>=LEN){ // if the string is all whitespace
	O=(char*)calloc(1,sizeof(char));
	O[0]='\0';
	return O;}

stop=LEN-1; ///< Begin at the last character
//printf("stop is %d\n",stop);
i=-1;
while(i!=0){
	switch(I[stop]){
		case ' ':
		case '\t':
		case '\n':
		case '\r':
		case '\f':
		case '\v':
			stop --;
			break;
		default:
			i=0;
			break;
		}
	if(stop<0){stop=i=0;}
	}

if(stop<start){
	printf("In prune_string_whitespace, stop is less than start.  This is a coding bug.\nExiting.\n");
	exit(1);}

len=stop-start+1;
O=(char*)calloc(len+1,sizeof(char));
for(i=0;i<len;i++){O[i]=I[i+start];}
O[len]='\0';

return O;
}

/** Functions that produce variable-format output for numbers and strings.

There must be some already-written way to do this, but I suspect I can write
the functions before I can figure that out.  So, these functions will generate
strings when one does not know ahead of time the output format.

For example:

Assume you know you need to write out a double variable in exponential 
notation, but that the program will determine the proper number of 
decimal places and character width.  These functions will help.

So, where normally, you might write:

- fprintf("%18.8e",x); 
	- that's a right justified exponential, 18 chars wide with 8 decimal places

You can do, instead:

- char j,*s;
- int w,d;

- w = figure_out_width();
- d = figure_out_number_decimal_places();
- j = figure_out_justification(); // 'L' or 'l' for left, 'R' or 'r' for right

- s = strdup(get_exp_string(x,j,w,d));

- fprintf("%s",s);

Normally, this is a lot of bother.  But, at times, it can be useful.  

I'm lazy, so the maximum number of digits that can be printed is about 145.
For float values, that number is divided into pre- and post-decimal, which makes
for some limits.  Hopefully, these limits are generous enough for most applications.
This is hard-coded in the number->string functions below.  There are much cooler 
ways to do this, and anyone who wants to code them is very welcome to do so.
*/

/************* get_char_string() ***************/
const char * get_char_string(const char *x, char j, int w){
char *s=NULL,*t=NULL;
int i,k;

k=w-strlen(x); // difference between outgoing and incoming lengths
/// If the incoming string (x) is longer than outgoing size (w),
/// 	write something to stderr and print what fits
if(k<0){ fprintf(stderr,"\nget_char_string: incoming string is longer than requested outgoing size.\n\t  Doing the best we can.\n\n"); }
s=(char*)calloc(w+1,sizeof(char)); ///< Allocate memory for new string
// Since clear allocating, we don't need to worry about the final '\0', but might anyway

switch (j) {
	case 'L':
	case 'l':
		if(k<=0){for(i=0;i<w;i++){s[i]=x[i];}} ///< If incoming string is long
		else{
			t=(char*)calloc(k+1,sizeof(char));
			for(i=0;i<k;i++){t[i]=' ';} // space padding
			t[k]='\0'; // for safety
			if((strlen(t)+strlen(x))!=w){mywhine("get_char_string: coding bug in s/t/w string handling, left.");}
			strcat(s,x);
			strcat(s,t);
			}
		break;
	case 'R':
	case 'r': 
		if(k<=0){for(i=w-1;i>=0;i--){s[i]=x[i-k];}} ///< If incoming string is long
		else{
			t=(char*)calloc(k+1,sizeof(char));
			for(i=0;i<k;i++){t[i]=' ';} // space padding
			t[k]='\0'; // for safety
			if((strlen(t)+strlen(x))!=w){mywhine("get_char_string: coding bug in s/t/w string handling, right.");}
			strcat(s,t);
			strcat(s,x);
			}
		break;
	default:
		mywhine("get_char_string: unrecognized justification specification");
		break;
	}
if(t!=NULL){free(t);}
s[w]='\0'; // for safety
return s;
}

/************* get_float_string() ***************/
const char * get_float_string(double x, char j, int w, int d){
char *s=NULL,*t=NULL,*xs,num[151],cdig[2];
int i=0,k=0,tst=0,len=0,start=0,end=0,digit=0,checkpt=0,first_non_nine=0;

if(w>150){mywhine("In get_float_string: Cannot have total width greater than 150 characters.");}
if(d>74){mywhine("In get_float_string: Cannot have number of decimal places greater than 74.");}
cdig[0]=cdig[1]='\0';

if(x>=0){
	tst=(int)(log10(x))+1;
/*printf("int(log10(x)) is %d and tst is %d\n",(int)(log10(x)),tst);*/
	}
if(x<0){
	tst=(int)(log10(-x))+1;
/*printf("int(log10(-x)) is %d and tst is %d\n",(int)(log10(-x)),tst);*/
	}
/* set for case where fabs(x)<1 */
if(fabs(x)<1.0){tst=1;} 
if(x==0){tst=1;}
k=w-tst; // difference between outgoing length and length of integer part of float value
/// If the incoming string (x) is longer than outgoing size (w),
/// 	write something to stderr and print what fits
if(k<0){
	fprintf(stderr,"\nget_float_string: WARNING!!!  Integer portion of incoming value exceeds requested outgoing width.\n");
	fprintf(stderr,"\tReturning non-numeric string. \n\n");
	s=(char*)calloc((w+1),sizeof(char)); ///< Allocate memory for new string
	for(i=0;i<w;i++){s[i]='*';}
	return s;
	}
/** Otherwise, print a really wide version of the incoming number to a really wide string
  Here, the string is 150 chars wide, and 74 of those chars are decimal digits.
  Another 74 are allocated for the integer part of the float representation.
  That leaves two spaces, one for a decimal point and another for '-' if needed. 
  At present, this function will not return a '+' sign for positive values. */
sprintf(num,"%150.74f",x); ///< 150 wide, minus +/- designation and decimal point and half of that for the fraction
len=tst+d+2; ///< tst=integer part ; d=fraction part ; 2 = space for . & sign (\0 comes later)
if(d==0){len=tst+1;} ///< tst=integer part ; 1 = space for sign and . (\0 comes later)
if(w<len){ //< print a statement if we have to give back fewer decimal places to fit the width, but do it anyway
	fprintf(stdout,"\nget_float_string: Float width including requested decimal places exceeds requested outgoing width.\n\tRounding.  Hope that's ok.\n\n");
	len=w;}

start=75-tst-1; ///< the starting point in the big number 

if(start<0){fprintf(stderr,"get_float_string: WARNING!! Incoming number too large to round up.\n");
	fprintf(stderr,"\tReturning non-numeric string.\n\n");
	s=(char*)calloc((w+1),sizeof(char)); ///< Allocate memory for new string
	for(i=0;i<w;i++){s[i]='*';}
	return s; // return the bad s and stop bothering
	}
//if(start<0){mywhine("get_float_string: Problem in first assignment of start value.");}
end=start+len-1; ///< the ending point in the big number
if(end>=150){mywhine("get_float_string: Problem in first assignment of end value.");}

checkpt=end+1; ///< check the integer where the null character will finally go
//printf("checkpt is >>%d<<\n",checkpt);
if((checkpt<150)&&(num[checkpt]=='.')){checkpt++;} ///< if there is a decimal point there, move over by one, 
		///< but only if that won't put us over 150 chars (shouldn't ever happen, really) 
//printf("end is >>%d<< num[end] is >>%c<<\n",end,num[end]);
//printf("checkpt is >>%d<< num[checkpt] is >>%c<<\n",checkpt,num[checkpt]);
xs=(char*)calloc((len+1),sizeof(char)); ///< Allocate memory for new string


cdig[0]=num[checkpt];
digit=atoi(cdig); 
/*printf("num[checkpt] (init) is %c, num[end] is %c\n",num[checkpt],num[end]);*/
/// See if we have a "9X' situation where we have to round up
if((num[end]=='9')&&(digit>4)){
		// Check to the left until a non-nine digit is found
		for(first_non_nine=checkpt-1;first_non_nine>=start;first_non_nine--){
			if((num[first_non_nine]!='9')&&(num[first_non_nine]!='.')){break;}
			}
//printf("first_non_nine is >>%d<< num[first_non_nine] is >>%c<<\n",first_non_nine,num[first_non_nine]);
//printf("first_non_nine+1 is >>%d<< num[first_non_nine+1] is >>%c<<\n",first_non_nine+1,num[first_non_nine+1]);
		if((num[first_non_nine]=='.')||(num[first_non_nine]=='9')){mywhine("get_float_string: Problem with coding in first non-nine.");}
		if(first_non_nine==start){// If the first space:
//printf("first non-nine is at start\n");
			xs[0]=num[start];///< Print the first character (space or '-') 
			xs[1]='1';///< Print a 1 in the next space
			/// For the rest of the string-1, for every 9 add a zero (one space over).  Copy any decimals as-is
			for(i=1;i<len;i++){
//printf("i is %d ; num[start+i] is %c\n",i,num[start+i]);
				switch (num[start+i]) {
					case '9': xs[i+1]='0';
//printf("\t xs[1+i] is %c\n",xs[i+1]);
						break;
					case '.': xs[i+1]='.';
//printf("\t xs[1+i] is %c\n",xs[i+1]);
						break;
					default: mywhine("get_float_string: Problem with coding in all nines switch.");
					}
				}
			xs[len]='\0';///< Ensure null termination
//printf("strlen(xs) is %d ; xs is >>%s<<\n",strlen(xs),xs);
			}
		else{
//printf("first non-nine is not at start\n");
			for(i=start;i<first_non_nine;i++){xs[i-start]=num[i];} ///< Copy as usual up to "first_not_nine"
			/// Increment the "first_not_nine" digit 
			cdig[0]=num[first_non_nine];
			digit=atoi(cdig); 
//printf("digit is %d ; cdig is %s\n",digit,cdig);
			if((digit<0)||(digit>8)){mywhine("get_float_string: Digit outside range of 0-8 but shouldn't be.");}
			digit++; 
			sprintf(cdig,"%1d",digit);
			xs[first_non_nine-start]=cdig[0]; 
//printf("digit is now %d ; cdig is now %s ; xs is >>%s<<\n",digit,cdig,xs);
			/// Change other nines to zeros
			for(i=(first_non_nine+1);i<(len+start);i++){
//printf("i is %d ; num[i] is %c\n",i,num[i]);
				switch (num[i]) {
					case '9': xs[i-start]='0';
//printf("\t xs[1+i] is %c\n",xs[i+1]);
						break;
					case '.': xs[i-start]='.';
//printf("\t xs[1+i] is %c\n",xs[i+1]);
						break;
					default: mywhine("get_float_string: Problem with coding in all nines switch.");
					}
				}
			xs[len]='\0';///< Ensure null termination
			}
	}
else{ ///< if that situation isn't present
/*printf("num[checkpt] (else) is %c\n",num[checkpt]);*/
	switch(num[checkpt]){
		case '9':	
		case '8':
		case '7':
		case '6':
		case '5':	// round up the last number
			for(i=0;i<len-1;i++){xs[i]=num[start+i];} ///< Copy as usual up to the end-1 space
			/// Increment the "end-1" digit 
			cdig[0]=num[start+len-1];
//printf("cdig is >>%s<<\n",cdig);
			digit=atoi(cdig); 
//printf("digit is >>%d<<\n",digit);
			if((digit<0)||(digit>8)){mywhine("get_float_string: Digit outside range of 0-8 but shouldn't be.");}
		digit++; 
//printf("digit is now >>%d<<\n",digit);
			sprintf(cdig,"%1d",digit);
//printf("cdig is now >>%s<<\n",cdig);
			xs[len-1]=cdig[0]; 
			xs[len]='\0'; ///< Ensure null termination
			break;
		default: // numbers 4 and lower and also the null terminator
			for(i=0;i<len;i++){xs[i]=num[start+i];} ///< Copy as usual up to the end space
			xs[len]='\0'; ///< Ensure termination.
			break;
		}
	}
	
// Since clear allocating, we don't need to worry about the final '\0', but might anyway
s=(char*)calloc(w,sizeof(char)); ///< Allocate memory for new string

k=w-len; ///< Calculated the padding needed (if any)
/*printf("len is %d and strlen(xs) is %d [w=%d ; d=%d]\n",len,strlen(xs),w,d);*/
if(k<0){mywhine("get_float_string: w less than len, but it shouldn't be.");}
switch (j) {
	case 'L':
	case 'l':
		t=(char*)calloc(k+1,sizeof(char));
		for(i=0;i<k;i++){t[i]=' ';} // space padding
		t[k]='\0'; // for safety
		if((strlen(t)+strlen(xs))!=w){mywhine("get_float_string: coding bug in s/t/w string handling, left.");}
		strcat(s,xs);
		strcat(s,t);
		break;
	case 'R':
	case 'r': 
		t=(char*)calloc(k+1,sizeof(char));
		for(i=0;i<k;i++){t[i]=' ';} // space padding
		t[k]='\0'; // for safety
		if((strlen(t)+strlen(xs))!=w){mywhine("get_float_string: coding bug in s/t/w string handling, right.");}
		strcat(s,t);
		strcat(s,xs);
		break;
	default:
		mywhine("get_float_string: unrecognized justification specification");
		break;
	}
free(t);
free(xs);
s[w]='\0'; // for safety
return s;
}


//START HERE -- rewrite this for exponential notation
/* */
/************* get_exp_string() ***************/
const char * get_exp_string(double x, char j, int w, int d){

printf("\n\nThe function get_exp_string hasn't been written yet.  It's very close, though!\n");
printf("Perhaps you might consider writing it?  The basic template is in place... \n");
printf("Essentially, you just need to alter the get_float_string function.\n");
printf("Unfortunately, however, at the time, this program will have to exit.\n");
printf("Bye!\n\n");
exit(0);

return "hello, world";
}
/*
const char * get_exp_string(double x, char j, int w, int d){
char *s,*t,*xs,*num[151],cdig[2];
int i=0,k=0,tst=0,len=0,start=0,digit=0;

if(w>150){mywhine("In get_float_string: Cannot have total width greater than 150 characters.");}
if(d>74){mywhine("In get_float_string: Cannot have number of decimal places greater than 74.");}
cdig[0]=cdig[1]='\0';

tst=log10(x)+1;
k=w-tst; // difference between outgoing length and length of integer part of float value
/// If the incoming string (x) is longer than outgoing size (w),
/// 	write something to stderr and print what fits
if(k<0){
	fprintf(stderr,"\nget_float_string: WARNING!!!  Integer portion of incoming value exceeds requested outgoing width.\n");
	fprintf(stderr,"\tReturning non-numeric string. \n\n");
	s=(char*)calloc((w+1),sizeof(char)); ///< Allocate memory for new string
	for(i=0;i<w;i++){s[i]='*';}
	return s;
	}
*/
/** Otherwise, print a really wide version of the incoming number to a really wide string
  Here, the string is 150 chars wide, and 74 of those chars are decimal digits.
  Another 74 are allocated for the integer part of the float representation.
  That leaves two spaces, one for a decimal point and another for '-' if needed. 
  At present, this function will not return a '+' sign for positive values. */
/*
sprintf(num,"%150.74f",x): ///< 150 wide, minus +/- designation and decimal point and half of that for the fraction
len=tst+d+3; ///< tst=integer part ; d=fraction part ; 3 = space for \0, . & sign
if(d==0){len=tst+2;} ///< tst=integer part ; 2 = space for \0 and sign.
if(w<len){ //< print a statement if we have to give back fewer decimal places to fit the width, but do it anyway
	fprintf(stdout,"\nget_float_string: Float width including requested decimal places exceeds requested outgoing width.\n\tRounding.  Hope that's ok.\n\n");
	len=w;}

start=75-tst-1; ///< the starting point in the big number 
if(start<0){fprintf(stderr,"get_float_string: WARNING!! Incoming number too large to round up.\n");
	fprintf(stderr,"\tReturning non-numeric string.\n\n");
	s=(char*)calloc((w+1),sizeof(char)); ///< Allocate memory for new string
	for(i=0;i<w;i++){s[i]='*';}
	return s; // return the bad s and stop bothering
	}
//if(start<0){mywhine("get_float_string: Problem in first assignment of start value.");}
end=start+len+1; ///< the ending point in the big number
if(end>=150){mywhine("get_float_string: Problem in first assignment of end value.");}

checkpt=end+1; ///< check the integer where the null character will finally go
if((checkpt<150)&&(num[checkpt]=='.')){checkpt++;} ///< if there is a decimal point there, move over by one, 
		///< but only if that won't put us over 150 chars (shouldn't ever happen, really) 
xs=(char*)calloc((len+1),sizeof(char)); ///< Allocate memory for new string
switch(num[checkpt]){
	case '9':	// do stuff for having found a nine
		// Check to the left until a non-nine digit is found
		for(first_non_nine=checkpt-1;first_non_nine>=start;first_non_nine--){if((num[first_non_nine]!='9')&&(num[first_non_nine]!='.')){break;}}
		if((num[first_non_nine]=='.')||(num[first_non_nine]=='9')){mywhine("get_float_string: Problem with coding in first non-nine.");}
		if(first_non_nine==start){// If the first space:
			xs[0]=num[start];///< Print the first character (space or '-') 
			xs[1]='1';///< Print a 1 in the next space
			/// For the rest of the string-1, for every 9 add a zero (one space over).  Copy any decimals as-is
			for(i=1;i<len;i++){
				switch (num[start+i]) {
					'9': xs[i+1]=0;
						break;
					'.': xs[i+1]='.';
						break;
					default: mywhine("get_float_string: Problem with coding in all nines switch.");
					}
				}
			xs[len]='\0';///< Ensure null termination
			}
		else{
			for(i=0;i<first_non_nine;i++){xs[i]=num[start+i];} ///< Copy as usual up to "first_not_nine"
			/// Increment the "first_not_nine" digit 
			cdig[0]=num[first_non_nine];
			digit=atoi(cdig); 
			if((digit<0)||(digit>8)){mywhine("get_float_string: Digit outside range of 0-8 but shouldn't be.");}
			digit++; 
			sprintf(cdig,"%1d",digit);
			xs[first_non_nine]=cdig[0]; 
			/// Change other nines to zeros
			for(i=(first_non_nine+1);i<len;i++){
				switch (num[start+i]) {
					'9': xs[i]=0;
						break;
					'.': xs[i]='.';
						break;
					default: mywhine("get_float_string: Problem with coding in all nines switch.");
					}
				}
			xs[len]='\0';///< Ensure null termination
			}
		break;
	case '8':
	case '7':
	case '6':
	case '5':	// round up the last number
		for(i=0;i<len-1;i++){xs[i]=num[start+i];} ///< Copy as usual up to the end-1 space
		/// Increment the "end-1" digit 
		cdig[0]=num[len-1];
		digit=atoi(cdig); 
		if((digit<0)||(digit>8)){mywhine("get_float_string: Digit outside range of 0-8 but shouldn't be.");}
		digit++; 
		sprintf(cdig,"%1d",digit);
		xs[len-1]=cdig[0]; 
		xs[len]='\0'; ///< Ensure null termination
		break;
	default: // numbers 4 and lower and also the null terminator
		for(i=0;i<len;i++){xs[i]=num[start+i];} ///< Copy as usual up to the end space
		xs[len]='\0'; ///< Ensure termination.
		break;
	}

// Since clear allocating, we don't need to worry about the final '\0', but might anyway
s=(char*)calloc(w,sizeof(char)); ///< Allocate memory for new string

k=w-len; ///< Calculated the padding needed (if any)
if(k<0){mywhine("get_float_string: w less than len, but it shouldn't be.");}
switch (j) {
	case 'L':
	case 'l':
		t=(char*)calloc(k+1,sizeof(char));
		for(i=0;i<k;i++){t[i]=' ';} // space padding
		t[k]='\0'; // for safety
		if((strlen(t)+strlen(x))!=w){mywhine("get_float_string: coding bug in s/t/w string handling, left.");}
		strcat(s,x);
		strcat(s,t);
		break;
	case 'R':
	case 'r': 
		t=(char*)calloc(k+1,sizeof(char));
		for(i=0;i<k;i++){t[i]=' ';} // space padding
		t[k]='\0'; // for safety
		if((strlen(t)+strlen(x))!=w){mywhine("get_float_string: coding bug in s/t/w string handling, right.");}
		strcat(s,t);
		strcat(s,x);
		break;
	default:
		mywhine(stderr,"get_float_string: unrecognized justification specification");
		break;
	}

s[w]='\0'; // for safety
return s;
}
*/
