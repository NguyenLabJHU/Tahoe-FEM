

/********************************************************



	VAX-style Date and Time functions.

	

	These have become de-facto FORTRAN standards.

	

	Also includes a FORTRAN-90 style datetime function.

	

	IMT 23 Sep 97

	

	GES 01 Mar 99



********************************************************/







#include <time.h>	/* GES */



#include "f2c.h"



/*

	Returns a date string of the form DD-MMM-YY

	in the leftmost 9 characters of s

	

	Typical FORTRAN call syntax:

	

		character today*9

		call date( today )

*/



int date_( char* s, ftnlen len )

{

	/* GES 01Mar99 */

	/* IMT 05Aug99 bug fix:  s need not be null terminated and might be too 

	                         small to hold terminating null inserted by strftime */

	time_t 		t1;

	struct tm 	*now;

	char		temp[10];	/* allow room for terminating null */

	int			maxLen;

	int			i;

	

	time( &t1 );							/* System time */

	now = localtime( &t1 ); 				/* format as struct */						

	strftime( temp, 10, "%d-%b-%y", now ); 	/* use ANSI strftime to format */

	

	/* Copy over at most the first 9 characters */

	maxLen = (len < 9) ? len : 9;

	for ( i = 0; i < maxLen; i++ )

		s[i] = temp[i];



	return 1;

}







/*

	Returns a time string of the form HH:MM:SS	

	in the leftmost 8 characters of s



	Typical FORTRAN call syntax:

	

		character now*8

		call time( now )

*/



int time_( char *s, ftnlen len )

{

	/* GES 01Mar99 */

	/* IMT 05Aug99 bug fix:  s need not be null terminated and might be too 

	                         small to hold terminating null inserted by strftime */

	time_t 		t1;

	struct tm 	*now;

	char		temp[9];	/* allow room for terminating null */

	int			maxLen;

	int			i;

	

	time( &t1 );							/* System time */

	now = localtime( &t1 ); 				/* format as struct */

	strftime( temp, 9, "%I:%M:%S", now ); 	/* use ANSI strftime to format */

	

	/* Copy over at most the first 8 characters */

	maxLen = (len < 8) ? len : 8;

	for ( i = 0; i < maxLen; i++ )

		s[i] = temp[i];



	

	return 1;

}









/*

	Return the following as integer values:

	

		year: 	the year

		month:	the month

		day:	the day



	Typical FORTRAN call syntax:

	

		integer month, day, year

		call idate( month, day, year )

*/



int idate_( integer* month, integer* day, integer* year )

{

	/* GES 01Mar99 */

	time_t t1;

	struct tm *now;

	

	time( &t1 );				/* System time */

	now = localtime( &t1 ); 	/* format as struct */



	*month = now->tm_mon + 1;

	*day = now->tm_mday;

	*year = now->tm_year + 1900;

	

	return 1;

}









/*

	Return the following:

	

		sDate 		a string of the form YYYYMMDD	in the leftmost 8 characters

		sTime 		a string of the form HHMMSS.SSS	in the leftmost 10 characters

		zone 		a string of the form +HHMM in the leftmost 5 characters	

		        	( contains the time zone difference with respect to GMT)

		values[0] 	the year

		values[1]	the month

		values[2]	the day

		values[3]	the GMT time difference in minutes

		values[4]	the hour

		values[5]	the minutes

		values[6]	the seconds

		values[7]	the milliseconds (always zero)	



	Typical FORTRAN call syntax:

	

		character sdate*6, stime*10, szone*5

		dimension idatev(8)

		call date( sdate, stime, szone, idatev )

		

	NOTE:  this function assumes the arrays are large enough to accept the output

*/





int datetime_( char* sDate, char* sTime, char* zone, integer* values, 

				ftnlen dateLen, ftnlen timeLen, ftnlen zoneLen )

{

	time_t 		t1;

	struct tm 	*now;

	struct tm 	localNow;

	struct tm 	utcNow;

	integer 	zoneTime;

	char		temp[11];		/* allow room for terminating null */

	int			maxLen;

	int			i;

	

	time( &t1 );				/* System time */

	now = gmtime( &t1 );		/* format as struct in utc time */

	utcNow = *now;

	now = localtime( &t1 ); 	/* format as struct in local time */

	localNow = *now;



	/* Set up the values array */

	values[0] = localNow.tm_year + 1900;

	values[1] = localNow.tm_mon + 1;

	values[2] = localNow.tm_mday;

	

	/* Handle non-whole hour timezones correctly; produce value in minutes */

	values[3] = (localNow.tm_hour - utcNow.tm_hour) * 60;

	values[3] += localNow.tm_min - utcNow.tm_min;				

	

	values[4] = localNow.tm_hour;

	values[5] = localNow.tm_min;

	values[6] = localNow.tm_sec;

	values[7] = 0;

	

	

	/* Fill in the date string */

	strftime( temp, 11, "%Y%y%m%d", now ); 	/* use ANSI strftime to format */	

	/* Copy over at most the first 8 characters */

	maxLen = (dateLen < 8) ? dateLen : 8;

	for ( i = 0; i < maxLen; i++ )

		sDate[i] = temp[i];

	

	/* Fill in the time string, zeroing out the milliseconds */

	strftime( temp, 11, "%H%M%S", now );

	temp[6] = '.';

	temp[7] = temp[8] = temp[9] = '0';

	/* Copy over at most the first 11 characters */

	maxLen = (timeLen < 11) ? timeLen : 11;

	for ( i = 0; i < maxLen; i++ )

		sTime[i] = temp[i];

	

	/* Fill in the zone string */

	/* ...Get sign */

	if ( values[3] < 0 )

	{

		temp[0] = '-';

		zoneTime = -values[3];				// Force positive for following computations

	}

	else

	{

		temp[0] = '+';

		zoneTime = values[3];	

	}	

	/* ...Get hours */

	sprintf( temp+1, "%2d", zoneTime/60 );

	/* ...Get the minutes */

	sprintf( temp+3, "%2d", zoneTime % 60 );

	/* Copy over at most the first 5 characters */

	maxLen = (zoneLen < 5) ? zoneLen : 5;

	for ( i = 0; i < maxLen; i++ )

		zone[i] = temp[i];

	

	

	return 1;

}

