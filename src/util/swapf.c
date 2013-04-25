/* bytes 1,2,3,4 switched to 4,3,2,1. */
#include <stdio.h>
#include <math.h>
main(ac,ap)
char *ap[];
int ac;
{
	char d[4],tmp;
	while(read(0,d,4)) {
	tmp=d[0];
	d[0]=d[3];
	d[3]=tmp;
	tmp=d[1];
	d[1]=d[2];
	d[2]=tmp;
	write(1,d,4);
	}
}
