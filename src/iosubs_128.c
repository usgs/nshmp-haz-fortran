/* io subroutines for hazard programs, last modified 10/00 */
/* these subroutines also used for opening binary files used in finite-diff
	programs and elsewhere. Mods of 2000 just for conformity with C standards */
/* the header structure is only relevant for the puthead and gethead routines */
#include <stdio.h>
#include <fcntl.h>
FILE *fp,*fopen();
int fd;
struct header
{
   char name[128][6];
   int period;
   int nlev;
   float xlev[20];
   float extra[10];
};
struct header *headr;
void openr_(name,len)
char name[];
int len;
{
       int i;
       char st[80];
/* modified to 80 chars to allow long directory path   */
       for(i=0; i<80; i++) st[i]= '\0';
       for(i=0; i<80 && name[i] != ' ' && name[i] != '\n' ; i++)
          st[i]=name[i];
       if((fp=fopen(st,"r"))==NULL) {
          fprintf(stderr,"cant open %s\n",st);
       exit(1);
       }
}
void openw_(name,len)
char name[];
int len;
{
       int i;
       char st[80];
       for(i=0; i<80; i++) st[i]= '\0';
       for(i=0; i<80 && name[i] != ' ' && name[i] != '\n' ; i++)
          st[i]=name[i];
       if((fp=fopen(st,"w+"))==NULL) {
          fprintf(stderr,"cant open %s\n",st);
       exit(1);
       }
}
void openwx_(fpx,name,len)
char name[];
int len;
int *fpx;
{
       int i;
       char st[80];
       for(i=0; i<80; i++) st[i]= '\0';
       for(i=0; i<80 && name[i] != ' ' && name[i] != '\n' ; i++)
          st[i]=name[i];
       if((*fpx=creat(st,0664)) < 0) {
          fprintf(stderr,"cant open %s\n",st);
       exit(1);
       }
       fprintf(stderr,"fpx %d\n",*fpx);
}
void openrx_(fpx,name,len)
/* open file descriptor for readonly (enables several open files for reading) */
char name[];
int len;
int *fpx;
{
       int i;
       char st[80];
       for(i=0; i<80; i++) st[i]= '\0';
       for(i=0; i<80 && name[i] != ' ' && name[i] != '\n' ; i++)
          st[i]=name[i];
       if((*fpx=open(st,O_RDONLY)) < 0) {
          fprintf(stderr,"cant open %s\n",st);
       exit(1);
       }
       fprintf(stderr,"fpx %d\n",*fpx);
}
void close_()
{
       if(fclose(fp)!= 0) {
          fprintf(stderr,"cant close \n");
          exit(1); }
}
void closex_()
{
       if(close(fd)!= 0) {
          fprintf(stderr,"cant close  %d\n",fd);
          exit(99); }
}
void getbuf_(buf,bufsiz,readn)
int *bufsiz,*readn;
short int *buf;
{
       *readn= fread(buf,2,*bufsiz,fp);
}
void getbuf2_(buf2,bufsiz,readn)
int *bufsiz,*readn;
float *buf2;
{
       *readn= fread(buf2,4,*bufsiz,fp);
}
void getbuf3_(buf2,bufsiz,readn,skip)
int *bufsiz, *readn;
int *skip;
float *buf2;
{
       fseek(fp, *skip, 0);
       *readn= fread(buf2,4,*bufsiz,fp);
/*
	fprintf(stderr,"file number,%d\n",fileno(fp));
 fprintf(stderr,"readn bufsiz,%i,%i,\n",*readn,*bufsiz);
	 fprintf(stderr,"after fread feof %d\n",feof(fp));
					*/
}
void putbuf2_(buf2,bufsiz,readn)
int *bufsiz,*readn;
float *buf2;
{
       *readn= fwrite(buf2,4,*bufsiz,fp);
}
void putbuf_(buf2,bufsiz,readn)
int *bufsiz,*readn;
short int *buf2;
{
       *readn= fwrite(buf2,2,*bufsiz,fp);
}
void putbufx_(fpx,buf2,bufsiz,readn)
unsigned int *bufsiz,*readn;
float *buf2;
int *fpx;
{
       *readn= write(*fpx,buf2,*bufsiz*4);
}
void getbufx_(fpx,buf2,bufsiz,readn)
unsigned int *bufsiz,*readn;
float *buf2;
int *fpx;
{
       *readn= read(*fpx,buf2,*bufsiz*4);
}
void putbufcx_(fpx,buf2,bufsiz,writen)
unsigned int *bufsiz,*writen;
char *buf2;
int *fpx;
{
       *writen= write(*fpx,buf2,*bufsiz);
}
void puthead_(fpx,headr,bufsiz,readn)
struct header *headr;
unsigned int *bufsiz,*readn;
int *fpx;
{
       *readn= write(*fpx,headr,*bufsiz);
}
void gethead_(headr,bufsiz,readn)
struct header *headr;
int *bufsiz,*readn;
{
      *readn= fread(headr,*bufsiz,1,fp);
}
void getheadx_(fpx,headr,bufsiz,readn)
struct header *headr;
int *bufsiz,*readn;
int *fpx;
{
      *readn= read(*fpx,headr,*bufsiz);
}
void getbufc_(buf,bufsiz,readn)
int *bufsiz,*readn;
short int *buf;
{
       *readn= fread(buf,1,*bufsiz,fp);
}
void putbufc_(buf2,bufsiz,readn)
int *bufsiz,*readn;
char *buf2;
{
       *readn= fwrite(buf2,1,*bufsiz,fp);
}
