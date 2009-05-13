#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "adhoc_utils.h"

void get_timestamp(char ** stamp) {
	char * WEEK_DAYS[] = {"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};
	char * MONTHS[] = {
		"Jan", "Feb", "Mar", "Apr", 
		"May", "Jun", "Jul", "Aug", 
		"Sep", "Oct", "Nov", "Dec"
	};

	time_t t;
	struct tm *lt;

	t = time(NULL);
	lt = localtime(&t);

	asprintf(stamp, "%s %s %02d, %04d %02d:%02d:%02d", 
		WEEK_DAYS[lt->tm_wday],
		MONTHS[lt->tm_mon],
		lt->tm_mday,
		lt->tm_year+1900,
		lt->tm_hour,
		lt->tm_min,
		lt->tm_sec
	);
}
