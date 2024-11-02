/************************************************************************/
/*                                                                      */
/* CPDBE (CoPy DataBase Entries)                                        */
/*                                                                      */
/* SJF 000430 from LSDBE by SJF from RMDBE.C by JQ                      */
/*                                                                      */
/*                                                                      */
/************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#define FOREVER 1
#define TRUE                    1
#define FALSE                   0

char *progname;
char tubestr[500*4];                   /* ridiculously big number! */
char db_path[FILENAME_MAX];
char default_dbpath[]="/usr/dbgt";
/*char default_dbpath[]="./"; */

extern int optind; /* index of first non option command line arg */

typedef 
enum { None, CopyPedsAndVars, CopyTubesOff } Action;

Action action=None;

void usage();
void getoptions(int , char **);
FILE *output=NULL;

main(int argc, char **argv)
{

  char *runid;
  char *runid2;
  char hrc_name[FILENAME_MAX];
  char ped_dbfile[FILENAME_MAX];
  char n2_dbfile[FILENAME_MAX];
  char toff_dbfile[FILENAME_MAX];
  int runno;
  int date;
  int ntoff;
  int *toff;
  int i;
  int j;
  int nskip;
  int intdum;
  
  output=stdout;
  progname = basename(*argv);
  
  strcpy(db_path,default_dbpath);
  getoptions(argc,argv);

  if((argc-optind)!=3){        /* must give runid and date */
    usage();
    exit(EXIT_FAILURE);
  }

  runid=*(argv+optind);
  runid2=*(argv+optind+1);
  sscanf(*(argv+optind+2),"%d",&date);

  sscanf(runid,"%*c%*c%d",&runno);  /* assume run name is XX.......  */

  /* make name like "/usr/dbgt/hrc97a" so that only                   */
  /* .cpeds or .n2gains.... needs to be appended to make database name */
  strcpy(hrc_name,db_path);
  if(*(hrc_name+strlen(hrc_name)-1)!='/')
    strcat(hrc_name,"/");
  strcat(hrc_name,mk_hrc_name(date));

  /* Pedestals: */
  strcpy(ped_dbfile,hrc_name);
  strcat(ped_dbfile,".cpeds");

  /* N2 Gains: */
  strcpy(n2_dbfile,hrc_name);
  strcat(n2_dbfile,".n2gains");

  /* Tubes off: */
  strcpy(toff_dbfile,hrc_name);
  strcat(toff_dbfile,".ntubelist");

  switch(action)
    {
    case CopyTubesOff:
      if(chk_toff_(runid,&ntoff,toff_dbfile)!=0){
	exit(EXIT_FAILURE);
      }
      else{
	toff=(int *)malloc(ntoff*sizeof(float));
	read_toff_(runid,&ntoff,toff,toff_dbfile);
	write_toff_(runid2,&ntoff,toff,toff_dbfile);
      }
      break;

    case CopyPedsAndVars:
      if(chk_peds_(runid,&date,ped_dbfile)!=0){
	exit(EXIT_FAILURE);
      }
      else{
	float peds[10000];
	float vars[10000];
	float nums[3];
	char mode[100]; /* what's this for ?? */
	int nchan;
	int events;
	
	if(read_peds_(runid,&date,ped_dbfile,peds,vars,
		      &events,nums,mode,&nchan)!=0){
	  exit(EXIT_FAILURE);
	}
	else{
	  float *printme;
	  fflush(stdout);
	  fprintf(stdout,"Events: %d, Channels: %d, Mode: %s\n",
		  events,nchan,mode);

	  write_peds_(runid2,&date,ped_dbfile,peds,vars,
		      &events,nums,mode,&nchan);
	}
      }
      break;
    }
  
}

void getoptions(int argc, char **argv)
{
  extern char *optarg;
  int c;

  while(FOREVER){
    c = getopt(argc, argv, "pPoOh");
    switch(c){
    case 'p':
    case 'P':
      action=CopyPedsAndVars;
      break; 
    case 'o':
    case 'T':
      action=CopyTubesOff;
      break; 
    case 'h':
      usage();
      exit(EXIT_SUCCESS);
      break;
      
    case -1:
      return;
      
    default:
      printf("%s: ignoring unknown option %c \n",progname,c);
      break;
    }
  }
  
  return;
}


void usage()
{
  printf("Usage: %s [options] run_name run_name2 date\n",progname);
  printf("       where run_name is in form 'gtnnnn...' and date is an \n");
  printf("       integer such as 970101.\n");
  printf("\nBy default %s does nothing, ",progname);
  printf("at least one option must be selected,\n");
  printf("so, I guess, they are not really optional!\n");
  printf("Options are (%%s is string argument):\n");
  printf("   -P     copy peds and pedvars\n"); 
}

