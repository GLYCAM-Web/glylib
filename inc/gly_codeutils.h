/** \file  gly_codeutils.h
\brief Header file that loads utilities to facilitate programming
*/
#if !defined(GLYLIB_CODING_UTILITIES)
#define GLYLIB_CODING_UTILITIES

/** \addtogroup CODE_UTILS
 * @{
 */
/******************** STRUCTURES ******************/

/** An array of characters for times when a 2D char array might be confusing */
typedef struct { 
	char *T;
} llcharset;


/******************** FUNCTIONS ******************/

/** Prints usage.*/
void  myusage(const char *USAGE); ///< Prints usage message and exits(1).

/** Prints error message and exits(1). 
Prints the following: 

"This program encountered the following error and will exit:

	whinetext " 
*/
void  mywhine(const char *whinetext); ///< Prints error message and exits(1).

/** Prints usage and an error message and exits(1). See mywhine for more information.
 Prints a usage statement after the whinetext*/
void  mywhineusage(const char *whinetext,const char *USAGE);///< Prints error message and usage and exits(1).

/**  Returns a the time in one of four popular formats.
*
* Usage:
*
*      const char * mylocaltimestring,mytime(format)
*
*              -for example:
*
*      strcpy(mylocaltimestring,mytime("format"));
*
* where:
*
*      format is one of:
*
*              [keyword]       [example]
*
*              YYYYMMDD        20070321
*
*              YYYYMMDDhhmmss  20070321152110
*
*              pretty          Wed, 21 March, 2007
*
*              longpretty      15:12:10, Wed, 21 March, 2007
*/
const char * mytime(const char*); ///< Gets the time in a nice format

/** Prunes leading and trailing whitespace from string I*/
const char *prune_string_whitespace(const char *I); 

/** Returns a character string of width w characters and justified as
per the character in j (l or L or r or R for left or right) based on the
string provided in x. */
const char * get_char_string(const char *x, char j, int w);

/** Returns a character string corresponding to a floating point 
representation of the number supplied by x.  The return string will have a 
width of w characters, d of which will be integers appearing to the right of 
the decimal point.   The number will be justified as per the character in j 
(l or L or r or R for left or right). */
const char * get_float_string(double x, char j, int w, int d);

/** Returns a character string corresponding to a exponential notation
representation of the number supplied by x.  The return string will have a 
width of w characters, d of which will be integers appearing to the right of 
the decimal point.   The number will be justified as per the character in j 
(l or L or r or R for left or right). */
const char * get_exp_string(double x, char j, int w, int d);

/** @}*/

#endif
