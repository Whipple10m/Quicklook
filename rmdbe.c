/************************************************************************/
/*                                                                      */
/* RMDBE (ReMove DataBase Entries)                                      */
/*                                                                      */
/* JQ 970507                                                            */
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
int replace_tubes=FALSE;
int list_tubes=FALSE;
int del_n2gains=FALSE;
int del_peds=FALSE;
int del_toff=FALSE;
int more_toff=FALSE;

void usage();
void getoptions(int , char **);

main(int argc, char **argv)
{

  char *runid;
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
  
  progname = basename(*argv);
  
  strcpy(db_path,default_dbpath);
  getoptions(argc,argv);

  if((argc-optind)!=2){        /* must give runid and date */
    usage();
    exit(EXIT_FAILURE);
  }
  else if (!(replace_tubes || list_tubes || del_n2gains || del_peds 
             || del_toff || more_toff)){/* ie no options were given */
    usage();
    exit(EXIT_FAILURE);
  }
	   

  runid=*(argv+optind);
  sscanf(*(argv+optind+1),"%d",&date);

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


  if (del_peds)
    del_peds_(runid,&date,ped_dbfile);

  if (del_n2gains)
    del_n2gains_(&runno,&date,n2_dbfile); 
    
  if (del_toff)
    del_toff_(runid,toff_dbfile); 

  if (list_tubes){
    if(chk_toff_(runid,&ntoff,toff_dbfile)!=0){
      exit(EXIT_FAILURE);
    }
    else{
      toff=(int *)malloc(ntoff*sizeof(float));
      read_toff_(runid,&ntoff,toff,toff_dbfile);
      printf("\n");
      printf("%d tubes off for run %s in %s :\n",
	     ntoff,runid,toff_dbfile);
      for(i=0;i<ntoff;i+=10){
	for(j=0;(j<10)&&((i+j)<ntoff);j++)
	  printf("%5d",*(toff+i+j));
	printf("\n");
      }
    }
  }
 
  if (replace_tubes){
    if(chk_toff_(runid,&ntoff,toff_dbfile)!=0)
      printf("Entries not found for %s in %s, continuing....\n",
	     runid,toff_dbfile);
    else
      del_toff_(runid,toff_dbfile); 
    
    ntoff=0;
    toff=(int *)malloc(1*sizeof(float));  /* needed before reallocing */
    if(*(tubestr+strlen(tubestr)-1)!=',')
      strcat(tubestr,",");        /* hack to eliminate special case of */
                                  /* reading last entry of string      */ 

    while(sscanf(tubestr,"%d,%n",&intdum,&nskip)!=EOF){
      if((toff=realloc(toff,sizeof(int)*(ntoff+1)))==NULL){
	printf("Cannot allocate Memory, exiting...\n");
	exit(EXIT_FAILURE);
      }
      sscanf(tubestr,"%d,",toff+ntoff);
      strcpy(tubestr,tubestr+nskip);     /* chop off first entry each time */
      ntoff++;
    }
    printf("\n%d tubes to go off for run %s in %s :\n",
	   ntoff,runid,toff_dbfile);
    for(i=0;i<ntoff;i+=10){
      for(j=0;(j<10)&&((i+j)<ntoff);j++)
	printf("%5d",*(toff+i+j));
      printf("\n");
    }
    write_toff_(runid,&ntoff,toff,toff_dbfile); 
  }


  if (more_toff){
    if(chk_toff_(runid,&ntoff,toff_dbfile)!=0){
      ntoff=0;
      toff=(int *)malloc(1*sizeof(float));
      printf("Entries not found for %s in %s, continuing....\n",
	     runid,toff_dbfile);      
    }
    else{
      toff=(int *)malloc(ntoff*sizeof(float));
      read_toff_(runid,&ntoff,toff,toff_dbfile);
      del_toff_(runid,toff_dbfile); 
    }

    if(*(tubestr+strlen(tubestr)-1)!=',')
      strcat(tubestr,",");        /* hack to eliminate special case of */
                                  /* reading last entry of string      */ 
    while(sscanf(tubestr,"%d,%n",&intdum,&nskip)!=EOF){
      if((toff=realloc(toff,sizeof(int)*(ntoff+1)))==NULL){
	printf("Cannot allocate Memory, exiting...\n");
	exit(EXIT_FAILURE);
      }
      sscanf(tubestr,"%d,",toff+ntoff);
      strcpy(tubestr,tubestr+nskip);     /* chop off first entry each time */
      ntoff++;
    }
    
    printf("\n%d tubes to go off for run %s in %s :\n",
	   ntoff,runid,toff_dbfile);
    for(i=0;i<ntoff;i+=10){
      for(j=0;(j<10)&&((i+j)<ntoff);j++)
	printf("%5d",*(toff+i+j));
      printf("\n");
    }
    write_toff_(runid,&ntoff,toff,toff_dbfile); 
  }
}




void getoptions(int argc, char **argv)
{
  extern char *optarg;
  int c;

  while(FOREVER){
    c = getopt(argc, argv, "hd:lr:m:ntpa");
    switch(c){
    case 'a':
      del_n2gains=TRUE;
      del_peds=TRUE;
      del_toff=TRUE;
      break;
    case 'h':         
      usage();
      exit(EXIT_SUCCESS);
      break;
    case 'd':
      sscanf(optarg,"%s",db_path);
      break; 
    case 'r':
      sscanf(optarg,"%s",tubestr);
      replace_tubes=TRUE;
      break; 
    case 'l':
      list_tubes=TRUE;
      break;
    case 'n':
      del_n2gains=TRUE;
      break;
    case 'm':
      sscanf(optarg,"%s",tubestr);
      more_toff=TRUE;
      break;
    case 'p':
      del_peds=TRUE;
      break;
    case 't':
      del_toff=TRUE;
      break;
    case -1:
      return;
      break;
    default:
      printf("%s: ignoring unknown option %c \n",progname,c);
      break;
    }
  }
  
  return;
}


void usage()
{
  printf("Usage: %s [options] run_name date\n",progname);
  printf("       where run_name is in form 'gtnnnn...' and date is an \n");
  printf("       integer such as 970101.\n");
  printf("\nBy default %s does nothing, ",progname);
  printf("at least one option must be selected,\n");
  printf("so, I guess, they are not really optional!\n");
  printf("Options are (%%s is string argument):\n");
  printf("   -n     delete entries from n2gains database\n"); 
  printf("   -p     delete entries from pedestals/pedvars database\n"); 
  printf("   -t     delete entries from tubes off database\n"); 
  printf("   -a     delete all entries (ie n2, peds and tubes off)\n");
  printf("   -h     help (this message)\n"); 
  printf("   -d %%s  change path to database files [default is %s]\n",
	 default_dbpath);
  printf("   -l     list tubes off from database\n");
  printf("   -m %%s  turn off more tubes \n");
  printf("          (example: if tubes 1,2 are off then %s -m 3,4 gt1234 ",
	 progname);
  printf("970101\n");
  printf("          results in tubes 1,2,3,4 being turned off)\n" );
  printf("   -r %%s  replace tubes off in database\n");   
  printf("Options can be combined but some combinations may have adverse ");
  printf("effects\n");
}

