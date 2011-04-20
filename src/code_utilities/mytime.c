// Function written by B. Lachele Foley, 2007
#include <mylib.h>
//#include "../inc/mylib.h"

// Usage:
//	const char * mylocaltimestring,mytime(format)
//		-for example:
//  	strcpy(mylocaltimestring,mytime("format"));
// where:
// 	format is one of:
//		[keyword]	[example]
// 		YYYYMMDD 	20070321
//		YYYYMMDDhhmmss	20070321152110
//		pretty		Wed, 21 March, 2007
//		longpretty	15:12:10, Wed, 21 March, 2007

const char * mytime(const char* timehow){
char timestring[1001],month[101],day[101];
 char* tm;
int DATE;
struct tm TM;
time_t LT;

// Get date
LT = time(NULL);
TM = *localtime(&LT);
switch(TM.tm_mon){
	case 0:		strcpy(month,"January");
			break; 
	case 1:		strcpy(month,"February");
			break; 
	case 2:		strcpy(month,"March");
			break; 
	case 3:		strcpy(month,"April");
			break; 
	case 4:		strcpy(month,"May");
			break; 
	case 5:		strcpy(month,"June");
			break; 
	case 6:		strcpy(month,"July");
			break; 
	case 7:		strcpy(month,"August");
			break; 
	case 8:		strcpy(month,"September");
			break; 
	case 9:		strcpy(month,"October");
			break; 
	case 10:	strcpy(month,"November");
			break; 
	case 11:	strcpy(month,"December");
			break; 
	default:	mywhine("Error assigning month in mytime.c");
	}
switch(TM.tm_wday){
	case 0:		strcpy(day,"Monday");
			break;
	case 1:		strcpy(day,"Tuesday");
			break;
	case 2:		strcpy(day,"Wednesday");
			break;
	case 3:		strcpy(day,"Thursday");
			break;
	case 4:		strcpy(day,"Friday");
			break;
	case 5:		strcpy(day,"Saturday");
			break;
	case 6:		strcpy(day,"Sunday");
			break;
	default:	mywhine("Error assigning week day in mytime.c");
	}
if(strcmp(timehow,"YYYYMMDD")==0){
//	YYYYMMDD 	20070321
	DATE=19000000 + TM.tm_year*10000 + (TM.tm_mon+1)*100 + TM.tm_mday;
	sprintf(timestring,"%d",DATE);
	}
if(strcmp(timehow,"YYYYMMDDhhmmss")==0){
//	YYYYMMDDhhmmss	20070321152110
	DATE=190000 + TM.tm_year*100 + (TM.tm_mon+1);
	sprintf(timestring,"%d",DATE);
	DATE=TM.tm_mday*1000000 + 10000*TM.tm_hour + 100*TM.tm_min + TM.tm_sec;
	sprintf(timestring,"%s%d",timestring,DATE);
	}
if(strcmp(timehow,"pretty")==0){
//	pretty		Wednesday, 21 March, 2007
	DATE=1900+TM.tm_year;
	sprintf(timestring,"%s, %d %s, %d",day,TM.tm_mday,month,DATE);
	}
if(strcmp(timehow,"longpretty")==0){
//	longpretty	15:12:10, Wednesday, 21 March, 2007
	DATE=1900+TM.tm_year;
	sprintf(timestring,"%02d:%02d:%02d, %s, %d %s, %d",TM.tm_hour,TM.tm_min,\
TM.tm_sec,day,TM.tm_mday,month,DATE);
	} 
 tm = timestring;
return tm;
}
