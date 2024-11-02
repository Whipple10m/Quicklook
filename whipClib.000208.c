/*
	whipClib.c
	961209
	J.Quinn and J.Buckley

	C libraries for Whipple data analysis.
*/

/* Major Modifications:                                             */
/*                                                                  */
/* 971202 - Error handling to include case of read but not write    */
/*          permission for database files included.                 */
/*        - All "return 2" (i.e. an error occurred) replaced with   */
/*          "exit(2)" since some of the calling programs do not     */
/*          have error handling.                                    */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

#define	TELE_10	1
#define TELE_11 2
#define	RETURN_ERROR	2
#define RETURN_SUCCESS	0

char *mk_hrc_name(int);
void make_hrc_name_(int *, char *);
char *rm_blanks(char *);
char *rootname(char *);
void num_chan_(int *, int *, int *, int *);
int read_peds_(char *, int *, char *, float *, float *,int *, 
	       float *, char *, int *); 
int chk_peds_(char *, int *, char *);
int write_peds_(char *, int *, char *, float *, float *,int *, 
		float *, char *, int *); 
int del_peds_(char *, int *, char *);
int read_n2gains_(int *, int *, char *, float *, float *,char *,
		  int *, int *); 
int chk_n2gains_(int*, int *, char *);
int write_n2gains_(int *, int *, char *, float *, float *,char *,
		  int *, int *); 
int del_n2gains_(int *, int *, char *);
int get_coords_(int *, int *, float *, float *, float *, int *) ;
void print_error(char *, char *) ;
void remove_pmt(int, double *, double *, int *) ;
int chk_toff_(char *, int *, char *);       
int read_toff_(char *, int *, int *, char *);       
int write_toff_(char *, int *, int *, char *);    
int del_toff_(char *, char *);    

/************************************************************************/
char *mk_hrc_name(int date)
{
  char name[6];
  char year_half[2];
  int year;
  
  date=date%1000000; /* Y2K complience or is it ?? */
  year=(int) (date+0.5)/10000;
  strcpy(year_half,((date-year*10000 < 600) ? "a" : "b"));
  sprintf(name,"hrc%2.2d%s",year,year_half);
  return name;
}

/************************************************************************/

void make_hrc_name_(int *date, char *name)
/* version of make_hrc_name which will work with FORTRAN */
{
  char year_half[2];
  int year;
  int Y2Kdate=*date;
  int YEARdate;
  if(Y2Kdate<800000)Y2Kdate+=1000000;
  YEARdate=Y2Kdate%1000000;

  year=(int) (YEARdate+0.5)/10000;
  strcpy(year_half,((YEARdate-year*10000 < 600) ? "a" : "b"));
  sprintf(name,"hrc%2.2d%s",year,year_half);
  return;
}
/*******************************************************************/
char *rm_blanks(char *string)
{
  char *new_string;
  int i, len_string;
  
  i=len_string=strlen(string);
  if (i>0){
    new_string=(char *)malloc(len_string*sizeof(char));
    strcpy(new_string,string);
    while((new_string[--i]==' ') && (i>=0))
      ;
    new_string[++i]='\0';
  }
  else new_string=string;
  
  return new_string;
}

/******************************************************************/
/* rootname:  (JQ 961111)                                         */
/* takes a character string (pointer) and returns the root part   */
/* (pointer). Example "/temp/name.tmp" -> "name"                  */
/* Use in same way as basename() and dirname()                    */
/* Modified 970224 so that only last extension is removed.        */
/* "/temp/name.ext1.ext2" -> "name.ext1"                          */
/******************************************************************/
     
char *rootname(char *instring)
{
  char *dot;
  char *newstring;

  if(strlen(instring)==NULL) 
    return (instring);
  newstring=(char *)malloc(sizeof(char)*strlen(instring));
  strcpy(newstring,basename(instring));

  if ((dot=strrchr(newstring,'.'))!= NULL)
    *dot='\0';
  
  return(newstring);
}

/******************************************************************/

void num_chan_(int *telescope, int *date, int *nadc, int *npmt)
{
  int Y2Kdate=*date;
  if(Y2Kdate<800000)Y2Kdate+=1000000;

  switch (*telescope){
  case TELE_10:
    if (Y2Kdate < 961203){
      *nadc=120;
      *npmt=109;
    }
    else if (Y2Kdate < 970901){
      *nadc=156;
      *npmt=151;
    }      
    else if (Y2Kdate < 990901){
      *nadc=336;
      *npmt=331;
    }      
    else {
      *nadc=492;
      *npmt=490;
    }
    break;
  case TELE_11:
    *nadc=120;
    *npmt=109;
    break;
  default:
    printf("num_chan_: %d - unknown telescope\n",*telescope);
    *nadc=0;
    *npmt=0;
    break;
  }

  return;
}

/*************************************************************************/
/* 	get_coords()
	961209
	JB

	Generates coordinates and neighbor list for a camera composed
	of *npmts 1" phototubes.   The algorithm assumes that the same
	number of PMTs is removed from each corner to make the total number
	work out, and that the missing tubes are not numbered - i.e. there
	are no empty channels. 
*/

void print_error(char *funname, char *estring)
{
        printf("%s ERROR: %s\n", funname, estring) ;
}

void remove_pmt(int ipmt, double *x, double *y, int *pn)
{
        register int    i ;

        for(i=ipmt;i<(*pn-1);i++) {
          x[i] = x[i+1] ;
          y[i] = y[i+1] ;
        }
        (*pn)-- ;
}

#include "GIIICamera.h"
int get_coords_GIII(int *telescope, int *Y2Kdate, 
		    float *wxdeg, float *wydeg, float *wradius,
		    int *nbr_list)
{
  int i;
  for(i=0;i<490;i++)
    {
      int j;
      wxdeg[i]=GIIIXcoord[i];
      wydeg[i]=GIIIYcoord[i];
      wradius[i]=GIIIRadius[i];
      for(j=0;j<7;j++)nbr_list[i*7+j]=GIIIneighbours[i][j];
    }
  return(RETURN_SUCCESS);
}

#define LSPACE  ((Y2Kdate<970901)?0.259:0.2436)
#define LSPACE2 (LSPACE*LSPACE)
#define TOLERANCE (LSPACE2*0.1)
#define RPMT    0.125

int get_coords_(int *telescope, int *date, 
		float *wxdeg, float *wydeg, float *wradius,
		int *nbr_list)
{
        register int    i,j,k,iring,jring,ipmt,iside ;
        double  dx0, dx_c60, dy_s60 ;
        double  sin60 ;
        double  cos60 ;
        double  dx[6], dy[6] ;
        double  drings ;
        double  x, y, dxp, dyp, dist2 ;
	double	*xpmt, *ypmt ;
        int     nrings, npmt, mpmt, ncorner ;
        int     ncorner_short, ncorner_long, nside ;
        int     lastpmt ;
        int     nsuggest1, nsuggest2 ;
        int     jmin,jmax,imin,imax ;
        int     inbr ;
        char    *pav ;
	char	*funname ;
	int	nchans ;
	int Y2Kdate;

	Y2Kdate=*date;
	if(Y2Kdate<800000)Y2Kdate+=1000000;

	if(Y2Kdate > 990901)return get_coords_GIII(telescope, &Y2Kdate, 
						   wxdeg, wydeg, wradius,
						   nbr_list);


/*	funname = (char *)malloc(FILENAME_MAX*sizeof(char));
	strcpy(funname,"get_coords") ; */

	funname="get_coords";
                                   
	printf("get_coords: telescope: %d, date: %d\n",*telescope, Y2Kdate) ;

	printf("%d %d %d %d\n",*telescope, Y2Kdate, npmt, nchans);
        num_chan_(telescope, &Y2Kdate, &nchans, &npmt) ;
	printf("%d %d %d %d\n",*telescope, Y2Kdate, npmt, nchans);


	if((*telescope == TELE_10) && (Y2Kdate < 930800)) {
          print_error(funname,
           "Tube coordinates not available for this utdate\n") ;
          return(RETURN_ERROR) ;
        }
	else if(*telescope == TELE_11) {
	  print_error(funname,
           "11m tube coordinates not yet available\n") ;
	  return(RETURN_ERROR) ;
        }
        drings = -0.5 + 0.5*sqrt((double)1.0+(4.0/3.0)*((double)npmt-1.0)) ;
        printf("drings = %lf\n", drings) ;
        nrings = (int)(drings + 0.5) ;
        if(fabs(drings - (double)nrings) > 1.0e-4) {
          printf("%d does not give an even number of rings\n", npmt) ;
          nrings = (int)(drings + 1.0) ;
          mpmt = 1 + 3*nrings*(nrings + 1) ;
          printf("mpmt:%d,npmt:%d\n",mpmt,npmt) ;
          ncorner = (mpmt-npmt)/6 ;
          if (((mpmt-npmt)%6 != 0)||(ncorner%2 == 0)) {
            printf("and number of corner tubes is not divisible by 6\n") ;
            printf("or is not an odd number\n") ;
            nsuggest1 = mpmt ;
            if(ncorner%2 == 0) nsuggest2 = mpmt - (ncorner+1) * 6 ;
            else nsuggest2 = mpmt - (ncorner) * 6 ;
            printf("try using %d or %d PMTs instead\n",
                   nsuggest1, nsuggest2) ;
            print_error(funname,"%s: bad number of tubes\n") ;
            return(RETURN_ERROR) ;
          }
        }
        else {
          nrings = (int)(drings + 0.5) ;
          mpmt = 1 + 3*nrings*(nrings + 1) ;    /* Check */
          ncorner = 0 ;
        }

        nside = nrings - ncorner ;
        ncorner_short = (ncorner - 1)/2 ;
        ncorner_long = (ncorner + 1)/2 ;

        printf("Number of PMTs in %d rings: %d\n", nrings, mpmt) ;
        printf("Number of PMTs to remove from ea. corner: %d\n",ncorner) ;
        printf("Number of PMTs on ea. side of outer ring: %d\n",nside) ;

        if((xpmt=(double *)calloc((size_t)(mpmt+1),sizeof(double)))==NULL){
          print_error(funname,
            "can not allocate memory for xpmt") ;
	  return(RETURN_ERROR) ;
	}
        if((ypmt=(double *)calloc((size_t)(mpmt+1),sizeof(double)))==NULL){
          print_error(funname,
            "can not allocate memory for ypmt") ;
	  return(RETURN_ERROR) ;
	}
        for(i=0;i<npmt*7;i++) { /* Initialize to the impossible tube */
          nbr_list[i] = npmt+1 ;
        }

        dx0 = LSPACE ;
        sin60 = sqrt(3.0)/2.0 ;
        cos60 = 0.5 ;
        dy_s60 = LSPACE * sin60 ;
        dx_c60 = LSPACE * cos60 ;
                                /* Initialize the relative      */
                                /* displacement vectors for     */
                                /* the six sides of the         */
                                /* hexagonal ring of pmts       */
        dx[0] = -1.0*dx_c60 ; dy[0] = -1.0*dy_s60 ;
        dx[1] = -1.0*dx0 ;    dy[1] = 0.0 ;
        dx[2] = -1.0*dx_c60 ; dy[2] = dy_s60 ;
        dx[3] = dx_c60 ;      dy[3] = dy_s60 ;
        dx[4] = dx0 ;         dy[4] = 0.0 ;
        dx[5] = dx_c60 ;      dy[5] = -1.0*dy_s60 ;

               xpmt[0] = ypmt[0] = 0.0 ;
        ipmt = 2 ;
                                /* For each ring of the camera  */
        for(iring=1; (iring<=nrings) && (ipmt<mpmt) ; iring++) {
          xpmt[--ipmt] = dx0*(double)iring ;
          ypmt[ipmt++] = 0.0 ;  /* For each side of the hex.    */
          for(iside=0; (iside<6) && (ipmt<mpmt); iside++) {
            for(i=0; (i<iring) && (ipmt<mpmt) ; i++) {
              xpmt[ipmt] = xpmt[ipmt-1] + dx[iside] ;
              ypmt[ipmt] = ypmt[ipmt-1] + dy[iside] ;
              ipmt++ ;
            } /* end for i */
          } /* end for iside */
        } /* end for iring */

        if(ncorner > 0) {
                                /* If the last ring is part.    */
                                /* filled, remove the corner    */
                                /* tubes, and renumber seq.     */
          lastpmt = 1 + 3*(nrings-1)*nrings ;
          ipmt = mpmt-1 ;
          do {
            for(i=0;i<ncorner_short;i++) remove_pmt(ipmt--,xpmt,ypmt,&mpmt) ;
            ipmt -= nside ;
            for(i=0;i<ncorner_long;i++) remove_pmt(ipmt--,xpmt,ypmt,&mpmt) ;
          } while(ipmt >= lastpmt) ;
        }
                        /* Build array of neighbor coordinates  */
                        /* Make a 1-d array for compat. with    */
                        /* fortran.  Only look in adjacent      */
                        /* rings to save time, and keep the     */
                        /* run time from being a function of    */
                        /* n^2                                  */
        nbr_list[0] = 2 ; nbr_list[1] = 3 ; nbr_list[2] = 4 ;
        nbr_list[3] = 5 ; nbr_list[4] = 6 ; nbr_list[5] = 7 ;
        for(iring=1; iring<=nrings ; iring++) {
          imin = 1 + 3*(iring-1)*iring ;
          imax = 1 + 3*iring*(iring+1) ;
          for(i=imin; (i<imax)&&(i<npmt); i++) {
            x = xpmt[i] ;
            y = ypmt[i] ;
            inbr = 0 ;
            for(jring=iring-1;jring<iring+2;jring++) {
              if((jring>=0)&&(jring<=nrings)) {
                if(jring==0) jmin = 0 ;
                else jmin = 1 + 3*(jring-1)*jring ;
                jmax = 1 + 3*jring*(jring+1) ;
                for(j=jmin; (j<jmax)&&(j<npmt); j++) {
                  if(j != i) {
                    dxp = xpmt[j] - x;
                    dyp = ypmt[j] - y ;
                    dist2 = dxp*dxp + dyp*dyp ;
                    if(dist2 < (LSPACE2 + TOLERANCE)) {
                      nbr_list[i*7+inbr++] = j+1 ;
                    } /* end if dist2 */
                  } /* end if j */
                } /* end for j */
              } /* end if jring */
            } /* end for jring */
          } /* end for i */
        } /* end for iring */
	for(i=0;i<npmt;i++) {
	  wxdeg[i] = xpmt[i] ;
	  wydeg[i] = ypmt[i] ;
	  wradius[i]=LSPACE;
	}
	return(RETURN_SUCCESS) ;
}

/***********************************************************************/
/* name,date,dbfile,peds,pedvars,event_cnt,nums[3],mode,nchan          */
/* returns: 0 success, 1 peds not found, 2 error (lock file or db not  */
/* found)  */
/***********************************************************************/
int read_peds_(char *run_id, int *date, char *db_file, float *peds, 
	       float *pedvars, int *event_cnt, float *nums, 
	       char *mode, int *n_chan)

{
  FILE *fp, *lock_fp;
  char line[80];
  int n_match;             /* number of matches for each info line of */
                           /* database file                           */
  int n_entries;
  int i,wait=0;
  char lock_file[]="/tmp/cpeds.lock";
  char sys_cmd[100];
  
  char tmp_run_id[FILENAME_MAX];
  char tmp_mode[8];
  int tmp_event_cnt;
  float tmp_nums[3];
  int tmp_date;
  int tmp_n_entries;
  
  int Y2Kdate=*date;
  if(Y2Kdate<800000)Y2Kdate+=1000000;

  strcpy(run_id,rm_blanks(run_id));
  strcpy(mode,rm_blanks(mode));
  strcpy(db_file,rm_blanks(db_file));

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("read_peds_: /tmp/cpeds.lock exists, trying again in 20 secs\n");
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("/tmp/cpeds.lock still exists after 2 mins, exiting....\n");
    exit(2);
  }

  if((fp = fopen(db_file, "r")) == NULL){
    printf("read_peds_ error: cannot find %s.\n",db_file);
    exit(2);
  }

  while(fgets(line, 80, fp) != NULL){
    n_match=sscanf(line,"%s %d %f %f %f%7c %d %d",tmp_run_id, 
		   &tmp_event_cnt,&tmp_nums[0],&tmp_nums[1],&tmp_nums[2],
		   tmp_mode,&tmp_date,&tmp_n_entries); 
    if (n_match<8)
      n_entries=120;
    else
      n_entries=tmp_n_entries;

    /* If run number the same as filename read peds */
    if((strcmp(tmp_run_id,run_id)==0)) {
      *event_cnt=tmp_event_cnt;
      *(tmp_mode+8)='\0';
      strcpy(mode,tmp_mode+4);
      for(i=0;i<3;i++)
	nums[i]=tmp_nums[i];
      *n_chan=n_entries;
      for(i=0; i< n_entries;i++)
	fscanf(fp,"%f",&peds[i]);
      for(i=0; i< n_entries;i++)
	fscanf(fp,"%f",&pedvars[i]);
      fclose(fp);
      printf("read_peds_: read entry for %s on %d from %s\n",
	     run_id,Y2Kdate,db_file);
      return 0;
    }else{
      for(i=0;i < 2*((n_entries+9)/10);i++)
	fgets(line,80,fp);
    }
  }

  fclose(fp);
  printf("read_peds_: cannot find entry for %s on %d in %s\n",
	   run_id,Y2Kdate,db_file);
  return 1;
}


/************************************************************************/
int chk_peds_(char *run_id, int *date, char *db_file)
/* returns: 0 success, 1 peds or database not found, 
           2 error (lock file)  */
{
  int tmp_event_cnt;
  int i,wait=0;
  int n_entries;
  int n_match;
  float tmp_nums[3];
  int tmp_n_entries, tmp_date;  
  FILE *fp,*lock_fp;
  char lock_file[]="/tmp/cpeds.lock";
  char sys_cmd[100];
  char line[80];
  char tmp_run_id[FILENAME_MAX];
  char tmp_mode[7];
  int Y2Kdate=*date;
  if(Y2Kdate<800000)Y2Kdate+=1000000;

  strcpy(run_id,rm_blanks(run_id));
  strcpy(db_file,rm_blanks(db_file));

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("chk_peds_: /tmp/cpeds.lock exists, trying again in 20 secs\n");
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("/tmp/cpeds.lock still exists after 2 mins, exiting....\n");
    exit(2);
  }

  if((fp = fopen(db_file, "r")) == NULL){
    printf("chk_peds_ : cannot find %s.\n",db_file);
    return 1;
  }

  while(fgets(line, 80, fp) != NULL){
    n_match=sscanf(line,"%s %d %f %f %f%7c %d %d",tmp_run_id,&tmp_event_cnt,
		   &tmp_nums[0],&tmp_nums[1],&tmp_nums[2],tmp_mode,
		   &tmp_date,&tmp_n_entries); 
    if (n_match<8){
      n_entries=120;
      tmp_date=0;
    }
    else
      n_entries=tmp_n_entries;

    /* If run number the same as filename read peds */
    if((strcmp(tmp_run_id,run_id)==0)){
      fclose(fp);
      printf("chk_peds_: entry exists for %s on %d in %s \n",
	     run_id,Y2Kdate,db_file);
      return 0;
    }else{
      for(i=0;i < 2*((n_entries+9)/10);i++)
	fgets(line,80,fp);
    }
  }

  fclose(fp);
  printf("chk_peds_: entry does not exist for %s on %d in %s\n",
	   run_id,Y2Kdate,db_file);
  return 1;
}


/************************************************************************/
int write_peds_(char *run_id, int *date, char *db_file, float *peds, 
		float *pedvars, int *event_cnt, float *nums, 
		char *mode, int *n_chan)
/* returns: 0 success, 2 error (lock file)  */

{
  int n_entries,i,j;
  int wait=0;
  char sys_cmd[100];
  FILE *fp,*lock_fp;
  char lock_file[]="/tmp/cpeds.lock";
  int Y2Kdate=*date;
  if(Y2Kdate<800000)Y2Kdate+=1000000;

  strcpy(run_id,rm_blanks(run_id));
  strcpy(mode,rm_blanks(mode));
  strcpy(db_file,rm_blanks(db_file));
  
  if (*n_chan<=0){
    printf("write_peds_ error: n_chan <=0, exiting....\n");
    exit(2);
  }

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("write_peds_: /tmp/cpeds.lock exists, trying again in 20 secs\n");
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("/tmp/cpeds.lock still exists after 2 mins, exiting....\n");
    exit(2);
  }
    

  if((lock_fp=fopen(lock_file,"w+"))==NULL){ /* cannot open lock-file */
    printf("write_peds_ error: cannot create %s, exiting....\n",lock_file);
    exit(2);
  }
  sprintf(sys_cmd,"/usr/bin/chmod +rwx %s",lock_file);
  system(sys_cmd);
  
  
  /*  n_entries=((int) (((float)*n_chan +9.0)/10.0))*10;;*/
  
  if ((fp=fopen(db_file,"r"))!=NULL){  /* database already exists */
    fclose(fp);
    if((fp=fopen(db_file,"a+"))==NULL){
      printf("write_peds_: error, can open %s for read but not write,\n",
	     db_file);
      printf("exiting...\n");
      fclose(lock_fp);  /* close and delete lock file */
      sprintf(sys_cmd,"rm -f %s",lock_file);
      system(sys_cmd);
      exit(2);
    }
  }
  else{                               /* doesn't exist so create it */
    if((fp=fopen(db_file,"a+"))==NULL){
      printf("write_peds_: error opening %s, exiting...\n");
      fclose(lock_fp);  /* close and delete lock file */
      sprintf(sys_cmd,"rm -f %s",lock_file);
      system(sys_cmd);
      exit(2);
    }
    sprintf(sys_cmd,"/usr/bin/chmod +rwx %s",db_file);
    system(sys_cmd);
  }

  fprintf(fp,"%7s%7i%7.3f%7.3f%7.3f%7s% 6i %6i\n",run_id,*event_cnt,
	  nums[0],nums[1],nums[2],mode,Y2Kdate,*n_chan);
  
  for(i=0; i < *n_chan;i+=10){
    for (j=0; (j < 10)&&((i+j)<*n_chan); j++)
      fprintf(fp,"%7.3f",*(peds+i+j));
    fprintf(fp,"\n");
  }
  for(i=0; i < *n_chan;i+=10){
    for (j=0;(j < 10)&&((i+j)<*n_chan); j++)
      fprintf(fp,"%7.3f",*(pedvars+i+j));
    fprintf(fp,"\n");
  }
  printf("write_peds_: wrote entry for %s on %d in %s \n",
	 run_id,Y2Kdate,db_file);
  
  fclose(fp);
  fclose(lock_fp);
  sprintf(sys_cmd,"rm -f %s",lock_file);
  system(sys_cmd);
  
  return 0;
}

/************************************************************************/
int del_peds_(char *run_id, int *date, char *db_file)
/* returns: 0 success, 1 cannot find entries, 2 error (lock file)  */
{
  FILE *fp, *lock_fp, *tmp_fp;
  char line[80];
  int n_match;             /* number of matches for each info line of */
                           /* database file                           */
  int n_entries;
  int i,wait=0;
  int found=0;
  char lock_file[]="/tmp/cpeds.lock";
  char tmp_file[]="/tmp/tmp.cpeds";
  char sys_cmd[100];

  char tmp_run_id[FILENAME_MAX];
  int tmp_event_cnt;
  float tmp_nums[3];
  int tmp_date;
  int tmp_n_entries;
  char tmp_mode[7];
  int Y2Kdate=*date;
  if(Y2Kdate<800000)Y2Kdate+=1000000;

  strcpy(run_id,rm_blanks(run_id));
  strcpy(db_file,rm_blanks(db_file));
  
  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("del_peds_: /tmp/cpeds.lock exists, trying again in 20 secs\n");
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("/tmp/cpeds.lock still exists after 2 mins, exiting....\n");
    exit(2);
  }
  
  /* open lock file */
  if((lock_fp=fopen(lock_file,"w+"))==NULL){ /* cannot open lock-file */
    printf("del_peds_ error: cannot create %s, exiting....\n",lock_file);
    exit(2);
  }
  sprintf(sys_cmd,"/usr/bin/chmod +rwx %s",lock_file);
  system(sys_cmd);
  
  
  /* open temp file */
  if((tmp_fp=fopen(tmp_file,"w+"))==NULL){
    printf("del_peds_ error: cannot create %s, exiting....\n",tmp_file);
    fclose(lock_fp);  /* close and delete lock file */
    sprintf(sys_cmd,"rm -f %s",lock_file);
    system(sys_cmd);
    exit(2);
  }

  if((fp = fopen(db_file, "r")) == NULL){
    printf("del_peds_ error: cannot open %s, exiting....\n",db_file);
    fclose(tmp_fp);    /* close tmp file */
    sprintf(sys_cmd,"rm -f %s",tmp_file);
    system(sys_cmd); 
    fclose(lock_fp);    /* close lock file */
    sprintf(sys_cmd,"rm -f %s",lock_file);
    system(sys_cmd); 
    exit(2);
  }
  
  while(fgets(line, 80, fp) != NULL){
    n_match=sscanf(line,"%s %d %f %f %f%7c %d %d",tmp_run_id, &tmp_event_cnt,
		   &tmp_nums[0],&tmp_nums[1],&tmp_nums[2],tmp_mode,
		   &tmp_date,&tmp_n_entries); 
    if (n_match<8){
      n_entries=120;
      tmp_date=0;
    }
    else
      n_entries=tmp_n_entries;
    
    /* If run number the same as filename read peds */
    if((strcmp(tmp_run_id,run_id)==0)){
      found=1;
      for(i=0;i < (int)(2.0*((float)n_entries/10.0)+0.9);i++)
	fgets(line,80,fp);
    }else{
      if(fprintf(tmp_fp,line) < 0){
	printf("Error writing %s, exiting...\n",tmp_file);
	fclose(tmp_fp);    /* close tmp file */
	sprintf(sys_cmd,"rm -f %s",tmp_file);
	system(sys_cmd); 
	fclose(lock_fp);    /* close lock file */
	sprintf(sys_cmd,"rm -f %s",lock_file);
	system(sys_cmd); 
	exit(2);
      }
      for(i=0;i < (int)(2.0*((float)n_entries/10.0)+0.9);i++){
	fgets(line,80,fp);
	if(fprintf(tmp_fp,line) < 0){
	  printf("Error writing %s, exiting...\n",tmp_file);
	  fclose(tmp_fp);    /* close tmp file */
	  sprintf(sys_cmd,"rm -f %s",tmp_file);
	  system(sys_cmd); 
	  fclose(lock_fp);    /* close lock file */
	  sprintf(sys_cmd,"rm -f %s",lock_file);
	  system(sys_cmd); 
	  exit(2);
	}
      }
    }
  }
  

  fclose(fp);
  fclose(tmp_fp);
  sprintf(sys_cmd,"cp %s %s~",db_file,db_file);
  system(sys_cmd); 
  sprintf(sys_cmd,"rm -f %s",db_file);
  system(sys_cmd); 
  sprintf(sys_cmd,"cp %s %s",tmp_file,db_file);
  system(sys_cmd); 
  sprintf(sys_cmd,"rm -f %s",tmp_file);
  system(sys_cmd); 

  fclose(lock_fp);
  sprintf(sys_cmd,"rm -f %s",lock_file);
  system(sys_cmd);

  if (found==1){
    printf("del_peds_: deleted entry for %s on %d from %s \n",
	   run_id,Y2Kdate,db_file);
    return 0;
  }else{
    printf("del_peds_: could not find entry for %s on %d in %s \n",
	   run_id,Y2Kdate,db_file);
    return 1;
  }
}

/************************************************************************/

int read_n2gains_(int *run_id, int *date, char *db_file, float *gains, 
	       float *gainvars, char *code, int *ped_id, int *n_chan)
/* returns: 0 success, 1 gains not found, 
           2 error (lock file or db not found)  */

{
  FILE *fp;
  FILE *lock_fp;
  char line[90];
  int n_match;             /* number of matches for each info line of */
                           /* database file                           */
  int n_entries;
  int i;
  int n_lines;
  int entries_per_line;    /* previously database was 16, from now on 10 */
  int wait=0;
  char lock_file[]="/tmp/gn2gains.lock";
  char sys_cmd[100];

  char tmp_code[4];
  int tmp_run_id;
  int tmp_ped_id;
  int tmp_date;
  int tmp_n_entries;
  
  int Y2Kdate=*date;
  if(Y2Kdate<800000)Y2Kdate+=1000000;

  strcpy(code,rm_blanks(code));
  strcpy(db_file,rm_blanks(db_file));

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("read_n2gains_: %s exists, trying again in 20 secs\n",lock_file);
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("%s still exists after 2 mins, exiting....\n",lock_file);
    exit(2);
  }

  if((fp = fopen(db_file, "r")) == NULL){
    printf("read_n2gains_ error: cannot find %s, exiting....\n",db_file);
    exit(2);
  }

  while(fgets(line, 90, fp) != NULL){
    n_match=sscanf(line,"%d%3c %d %d %d",&tmp_date,tmp_code,
		   &tmp_run_id,&tmp_ped_id,&tmp_n_entries); 
    if (n_match<5){
      n_entries=120;
      entries_per_line=16;
    }
    else{
      n_entries=tmp_n_entries;
      entries_per_line=10;
    }

    /* If run number the same as filename read gains */
    if(tmp_run_id==*run_id) {
      *ped_id=tmp_ped_id;
      *(tmp_code+3)='\0';
      strcpy(code,tmp_code+1);
      *n_chan=n_entries;
      for(i=0; i< n_entries;i++)
	fscanf(fp,"%f",&gains[i]);
      for(i=0; i< n_entries;i++)
	fscanf(fp,"%f",&gainvars[i]);
      fclose(fp);
      printf("read_n2gains_: read entry for %04i on %d from %s\n",
	     *run_id,Y2Kdate,db_file);
      return 0;
    }else{
      n_lines=2*((int) ((float)n_entries/(float)entries_per_line+0.9));
      for(i=0;i < n_lines;i++)
	fgets(line,90,fp);
    }
  }
  
  fclose(fp);
  printf("read_n2gains_: cannot find entry for %04i on %d in %s\n",
	   *run_id,Y2Kdate,db_file);
  return 1;
}


/************************************************************************/
int chk_n2gains_(int *run_id, int *date, char *db_file)
/* returns: 0 success, 1 gains or db not found, 
           2 error (lock file)  */
{
  int i;
  int wait=0;
  int n_entries;
  int n_match;
  int n_lines;
  int entries_per_line;    /* previously database was 16, from now on 10 */

  int tmp_n_entries;
  int tmp_date;  
  int tmp_run_id;
  int tmp_ped_id;
  char tmp_code[4];

  FILE *fp;
  FILE *lock_fp;
  char lock_file[]="/tmp/gn2gains.lock";
  char sys_cmd[100];
  char line[90];

  int Y2Kdate=*date;
  if(Y2Kdate<800000)Y2Kdate+=1000000;

  strcpy(db_file,rm_blanks(db_file));

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("chk_n2gains_: %s exists, trying again in 20 secs\n",lock_file);
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("%s still exists after 2 mins, exiting....\n",lock_file);
    exit(2);
  }

  if((fp = fopen(db_file, "r")) == NULL){
    printf("chk_n2gains_ : cannot find %s.\n",db_file);
    return 1;
  }

  while(fgets(line, 90, fp) != NULL){
    n_match=sscanf(line,"%d%3c %d %d %d",&tmp_date,tmp_code,
		   &tmp_run_id,&tmp_ped_id,&tmp_n_entries); 
    if (n_match<5){
      n_entries=120;
      entries_per_line=16;
    }
    else{
      n_entries=tmp_n_entries;
      entries_per_line=10;
    }


    /* If run number the same as filename read gains */
    if(tmp_run_id==*run_id){
      fclose(fp);
      printf("chk_n2gains_: entry exists for %04i on %d in %s \n",
	     *run_id,Y2Kdate,db_file);
      return 0;
    }else{
      n_lines=2*((int) ((float)n_entries/(float)entries_per_line+0.99));
      for(i=0;i < n_lines;i++)
	fgets(line,90,fp);
    }
  }

  fclose(fp);
  printf("chk_n2gains_: entry does not exist for %04i on %d in %s\n",
	   *run_id,Y2Kdate,db_file);
  return 1;
}


/************************************************************************/
int write_n2gains_(int *run_id, int *date, char *db_file, float *gains, 
		   float *gainvars, char *code, int *ped_id, int *n_chan)
  
/* returns: 0 success, 2 error (lock file)  */

{
  int n_entries,i,j;
  int wait=0;
  char sys_cmd[100];
  FILE *fp,*lock_fp;
  char lock_file[]="/tmp/gn2gains.lock";

  int Y2Kdate=*date;
  if(Y2Kdate<800000)Y2Kdate+=1000000;

  strcpy(code,rm_blanks(code));
  strcpy(db_file,rm_blanks(db_file));

  if (*n_chan<=0){
    printf("write_n2gains_ error: n_chan <=0, exiting....\n");
    exit(2);
  }

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("write_n2gains_: %s exists, ",lock_file);
    printf ("trying again in 20 secs\n");
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("%s still exists after 2 mins, exiting....\n",lock_file);
    exit(2);
  }
    
  if((lock_fp=fopen(lock_file,"w+"))==NULL){
    printf("write_n2gains_ error: cannot create %s, exiting...\n",lock_file);
    exit(2);
  }
  sprintf(sys_cmd,"/usr/bin/chmod +rwx %s",lock_file);
  system(sys_cmd);

  /*  n_entries=((int) (((float)*n_chan +9.0)/10.0))*10; */

  if ((fp=fopen(db_file,"r"))!=NULL){  /* database already exists */
    fclose(fp);
    if((fp=fopen(db_file,"a+"))==NULL){
      printf("write_n2gains_: error, can open %s for read but not write,",
	     db_file);
      printf(" exiting...\n");
      fclose(lock_fp);  /* close and delete lock file */
      sprintf(sys_cmd,"rm -f %s",lock_file);
      system(sys_cmd);
    }
  }
  else{                               /* doesn't exist so create it */
    if((fp=fopen(db_file,"a+"))==NULL){
      printf("write_n2gains_: error, cannot open %s, exiting.... \n",db_file);
      fclose(lock_fp);  /* close and delete lock file */
      sprintf(sys_cmd,"rm -f %s",lock_file);
      system(sys_cmd);
      exit(2);
    }
    else {
      sprintf(sys_cmd,"/usr/bin/chmod +rwx %s",db_file);
      system(sys_cmd);
    }
  }

  fprintf(fp,"%6i%3s %04i %04i %d\n",Y2Kdate,code,*run_id,*ped_id,*n_chan);
  
  for(i=0; i < *n_chan;i+=10){
    for (j=0; (j < 10)&&((i+j)<*n_chan); j++)
      fprintf(fp,"%7.3f",*(gains+i+j));
    fprintf(fp,"\n");
  }
  for(i=0; i < *n_chan;i+=10){
    for (j=0; (j < 10)&&((i+j)<*n_chan); j++)
      fprintf(fp,"%7.3f",*(gainvars+i+j));
    fprintf(fp,"\n");
  }
  printf("write_n2gains_: wrote entry for %04i on %d in %s \n",
	 *run_id,Y2Kdate,db_file);
  
  fclose(fp);
  fclose(lock_fp);
  sprintf(sys_cmd,"rm -f %s",lock_file);
  system(sys_cmd);
  
  return 0;
}

/************************************************************************/
int del_n2gains_(int *run_id, int *date, char *db_file)
/* returns: 0 success, 1 cannot find entries, 2 error (lock file)  */
{
  FILE *fp, *lock_fp, *tmp_fp;
  char line[90];
  int n_match;             /* number of matches for each info line of */
                           /* database file                           */
  int n_entries;
  int entries_per_line;
  int i;
  int n_lines;
  int wait=0;
  int found=0;
  char lock_file[]="/tmp/n2gains.lock";
  char tmp_file[]="/tmp/tmp.n2gains";
  char sys_cmd[100];

  int tmp_run_id;
  int tmp_date;
  int tmp_n_entries;
  int tmp_ped_id;
  char tmp_code[3];
  
  int Y2Kdate=*date;
  if(Y2Kdate<800000)Y2Kdate+=1000000;

  strcpy(db_file,rm_blanks(db_file));

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("del_n2gains_: %s exists, trying again in 20 secs\n",lock_file);
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("%s still exists after 2 mins, exiting....\n",lock_file);
    exit(2);
  }
  
  /* open lock file */
  if ((lock_fp=fopen(lock_file,"w+"))==NULL){
    printf("del_n2gains_ error: cannot create %s, exiting....\n",lock_file);
    exit(2);
  }
  sprintf(sys_cmd,"/usr/bin/chmod +rwx %s",lock_file);
  system(sys_cmd);
  
  /* open tmp file */
  if((tmp_fp=fopen(tmp_file,"w+"))==NULL){
    printf("del_n2gains_ error: cannot create %s, exiting....\n",tmp_file);
    exit(2);
  }

  if((fp = fopen(db_file, "r")) == NULL){
    printf("del_n2gains_ error: cannot find %s, exiting....\n",db_file);
    fclose(tmp_fp);    /* close tmp file */
    sprintf(sys_cmd,"rm -f %s",tmp_file);
    system(sys_cmd); 
    fclose(lock_fp);    /* close lock file */
    sprintf(sys_cmd,"rm -f %s",lock_file);
    system(sys_cmd); 
    exit(2);
  }

  while(fgets(line, 90, fp) != NULL){
    n_match=sscanf(line,"%d%3c %d %d %d",&tmp_date,tmp_code,
		   &tmp_run_id,&tmp_ped_id,&tmp_n_entries); 
    
    if (n_match<5){
      n_entries=120;
      entries_per_line=16;
    }
    else{
      n_entries=tmp_n_entries;
      entries_per_line=10;
    }

    n_lines=2*((int) ((float)n_entries/(float)entries_per_line+0.9));
    
    /* If run number the same as filename read gains */
    if(tmp_run_id==*run_id){
      found=1;
      for(i=0;i < n_lines;i++)
	fgets(line,90,fp);
    }else{
      if(fprintf(tmp_fp,line) < 0){
	printf("Error writing %s, exiting...\n",tmp_file);
	fclose(tmp_fp);    /* close tmp file */
	sprintf(sys_cmd,"rm -f %s",tmp_file);
	system(sys_cmd); 
	fclose(lock_fp);    /* close lock file */
	sprintf(sys_cmd,"rm -f %s",lock_file);
	system(sys_cmd); 
	exit(2);
      }
      for(i=0;i < n_lines;i++){
	fgets(line,90,fp);
	if(fprintf(tmp_fp,line) < 0){
	  printf("Error writing %s, exiting...\n",tmp_file);
	  fclose(tmp_fp);    /* close tmp file */
	  sprintf(sys_cmd,"rm -f %s",tmp_file);
	  system(sys_cmd); 
	  fclose(lock_fp);    /* close lock file */
	  sprintf(sys_cmd,"rm -f %s",lock_file);
	  system(sys_cmd); 
	  exit(2);
	}
      }
    }
  }

  fclose(fp);
  fclose(tmp_fp);
  
  sprintf(sys_cmd,"cp %s %s~",db_file,db_file);
  system(sys_cmd); 
  sprintf(sys_cmd,"rm -f %s",db_file);
  system(sys_cmd); 
  sprintf(sys_cmd,"cp %s %s",tmp_file,db_file);
  system(sys_cmd); 
  sprintf(sys_cmd,"rm -f %s",tmp_file);
  system(sys_cmd); 

  fclose(lock_fp);
  sprintf(sys_cmd,"rm -f %s",lock_file);
  system(sys_cmd);

  if (found==1){
    printf("del_n2gains_: deleted entry for %04i on %d from %s \n",
	   *run_id,Y2Kdate,db_file);
    return 0;
  }else{
    printf("del_n2gains_: could not find entry for %04i on %d in %s \n",
	   *run_id,Y2Kdate,db_file);
    return 1;
  }
}

/*************************************************************************/
/* int chk_toff_(char *runid, int *ntoff, char *dbfile)                  */
/* ntoff is ptr to number of tubes off, which is returned,               */
/* returns: 0 - success, 1 - entry or db not found,                      */
/*                                       2 - error (lock)                */
/*************************************************************************/

int chk_toff_(char *run_id, int *ntoff, char *db_file)       

{
  FILE *fp, *lock_fp;
  int i,wait=0;
  int tmp_ntoff;
  char line[2500];  /* (541*4)+100=2246 < 2500!! so all tubes can be off!*/ 
  char sys_cmd[100];
  char lock_file[]="/tmp/cpeds.lock";
  char tmp_run_id[FILENAME_MAX];
  
  strcpy(run_id,rm_blanks(run_id));
  strcpy(db_file,rm_blanks(db_file));

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("chk_toff_: %s exists, trying again in 20 secs\n",lock_file);
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("%s still exists after 2 mins, exiting....\n",lock_file);
    exit(2);
  }

  if((fp = fopen(db_file, "r")) == NULL){
    printf("chk_toff_ : cannot find %s\n",db_file);
    return 1;
  }


  while(fgets(line, 2500, fp) != NULL){
    sscanf(line,"%s %i",tmp_run_id,&tmp_ntoff);
    if (strcmp(tmp_run_id,run_id)==0){
      printf("chk_toff_: found entry for %s in %s\n",run_id,db_file);
      *ntoff=tmp_ntoff;
      fclose(fp);
      return 0;
    } 
  }
  printf("chk_toff_: entry does not exist for %s in %s\n",
	 run_id,db_file);
  fclose(fp);
  return 1;
}



/*************************************************************************/
/* int read_toff_(char *runid, int *ntoff, int *toff, char *dbfile)      */
/* ntoff is ptr to number of tubes off, toff is ptr to array of tubes off*/
/* returns: 0 - success, 1 - entry not found,                            */
/*                                       2 - error (lock or no database  */
/*************************************************************************/

int read_toff_(char *run_id, int *ntoff, int *toff, char *db_file)       

{
  FILE *fp, *lock_fp;
  int i,wait=0;
  int nskip;         /* number of characters to skip at start of each */
                     /* entry before actual tubes off are listed */
  int tube;
  int tmp_ntoff;
  char line[2500];  /* (541*4)+100=2246 < 2500!! so all tubes can be off!*/ 
  char sys_cmd[100];
  char lock_file[]="/tmp/toff.lock";
  char tmp_run_id[FILENAME_MAX];
  
  strcpy(run_id,rm_blanks(run_id));
  strcpy(db_file,rm_blanks(db_file));

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("read_toff_: %s exists, trying again in 20 secs\n",lock_file);
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("%s still exists after 2 mins, exiting....\n",lock_file);
    exit(2);
  }

  if((fp = fopen(db_file, "r")) == NULL){
    printf("read_toff_ error: cannot find %s, exiting....\n",db_file);
    exit(2);
  }


  while(fgets(line, 2500, fp) != NULL){
    sscanf(line,"%s %i",tmp_run_id,&tmp_ntoff);
    if (strcmp(tmp_run_id,run_id)==0){
      *ntoff=tmp_ntoff;
      sscanf(line,"%s %d%n",tmp_run_id,&tmp_ntoff,&nskip);
      for (i=0;i<*ntoff;i++){
	strcpy(line,line+nskip);
	sscanf(line,"%d",(toff+i));
	sscanf(line,"%d%n",&tube,&nskip);
      }
      printf("read_toff_: read entry for %s in %s\n",run_id,db_file);
      fclose(fp);
      return 0;
    } 
  }
  printf("read_toff_: entry does not exist for %s in %s\n",
	 run_id,db_file);
  fclose(fp);
  return 1;
}


/*************************************************************************/
/* int write_toff_(char *runid, int *ntoff, int *toff, char *dbfile)      */
/* ntoff is ptr to number of tubes off, toff is ptr to array of tubes off*/
/* returns: 0 - success, 2 - error (lock or no database)                 */
/*************************************************************************/

int write_toff_(char *run_id, int *ntoff, int *toff, char *db_file)       

{
  FILE *fp, *lock_fp;
  int i,wait=0;
  int nskip;         /* number of characters to skip at start of each */
                     /* entry before actual tubes off are listed */
  int tube;
  int tmp_ntoff;
  char line[2500];  /* (541*4)+100=2246 < 2500!! so all tubes can be off!*/ 
  char sys_cmd[100];
  char lock_file[]="/tmp/toff.lock";
  char tmp_run_id[FILENAME_MAX];
  
  strcpy(run_id,rm_blanks(run_id));
  strcpy(db_file,rm_blanks(db_file));

  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("write_toff_: %s exists, trying again in 20 secs\n",lock_file);
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("%s still exists after 2 mins, exiting....\n",lock_file);
    exit(2);
  }


  if((fp = fopen(db_file, "a")) == NULL){
    printf("write_toff_ error: cannot find %s, exiting....\n",db_file);
    exit(2);
  }


  if((lock_fp=fopen(lock_file,"w+"))==NULL){
    printf("write_toff_ error: cannot create %s, exiting...\n",lock_file);
    exit(2);
  }
  sprintf(sys_cmd,"/usr/bin/chmod +rwx %s",lock_file);
  system(sys_cmd);

  
  if ((fp=fopen(db_file,"r"))!=NULL){  /* database already exists */
    fclose(fp);
    if((fp=fopen(db_file,"a+"))==NULL){
      printf("write_toff_: error, can open %s for read but not write,",
	     db_file);
      printf(" exiting...\n");
      fclose(lock_fp);  /* close and delete lock file */
      sprintf(sys_cmd,"rm -f %s",lock_file);
      system(sys_cmd);
    }
  }
  else{                               /* doesn't exist so create it */
    if((fp=fopen(db_file,"a+"))==NULL){
      printf("write_toff_: error, cannot open %s, exiting.... \n",db_file);
      fclose(lock_fp);  /* close and delete lock file */
      sprintf(sys_cmd,"rm -f %s",lock_file);
      system(sys_cmd);
      exit(2);
    }
    else {
      sprintf(sys_cmd,"/usr/bin/chmod +rwx %s",db_file);
      system(sys_cmd);
    }
  }
  
  


  fprintf(fp,"%s %3d",run_id,*ntoff);
  for (i=0;i<*ntoff;i++)
    fprintf(fp," %3d",*(toff+i));
  fprintf(fp,"\n");
  
  printf("write_toff_: wrote entry for %s in %s\n",
	 run_id,db_file);
  fclose(fp);
  
  fclose(lock_fp);
  sprintf(sys_cmd,"rm -f %s",lock_file);
  system(sys_cmd);

  return 0;
}


int del_toff_(char *run_id, char *db_file)
/* returns: 0 success, 1 cannot find entries, 2 error (lock file)  */
{
  FILE *fp, *lock_fp, *tmp_fp;
  char line[2500];
  int i,wait=0;
  int found=0;
  char lock_file[]="/tmp/toff.lock";
  char tmp_file[]="/tmp/tmp.toff";
  char sys_cmd[100];

  char tmp_run_id[FILENAME_MAX];

  strcpy(run_id,rm_blanks(run_id));
  strcpy(db_file,rm_blanks(db_file));
  
  while (((lock_fp=fopen(lock_file,"r"))!=NULL) && (wait<120)){
    fclose(lock_fp);
    printf("del_toff_: %s exists, trying again in 20 secs\n",lock_file);
    sprintf(sys_cmd,"sleep %d",20);
    system(sys_cmd);
    wait+=20;
  }
  if(wait==120){
    printf("%s still exists after 2 mins, exiting....\n",lock_file);
    exit(2);
  }

  /* open lock file */
  if((lock_fp=fopen(lock_file,"w+"))==NULL){
    printf("del_toff_ error: cannot create %s, exiting...\n",lock_file);
    exit(2);
  }
  sprintf(sys_cmd,"/usr/bin/chmod +rwx %s",lock_file);
  system(sys_cmd);
  
  /* open temp file */
  if((tmp_fp=fopen(tmp_file,"w+"))==NULL){
    printf("del_toff_ error: cannot create %s, exiting...\n",tmp_file);
    fclose(lock_fp);    /* close lock file */
    sprintf(sys_cmd,"rm -f %s",lock_file);
    system(sys_cmd);     
    exit(2);
  }
  
  if((fp = fopen(db_file, "r")) == NULL){
    printf("del_toff_ error: cannot find %s, exiting....\n",db_file);
    fclose(tmp_fp);    /* close tmp file */
    sprintf(sys_cmd,"rm -f %s",tmp_file);
    system(sys_cmd); 
    fclose(lock_fp);    /* close lock file */
    sprintf(sys_cmd,"rm -f %s",lock_file);
    system(sys_cmd); 
    exit(2);
  }
  
  while(fgets(line, 2500, fp) != NULL){
    sscanf(line,"%s",tmp_run_id);
    /* If run number the same as filename read peds */
    if((strcmp(tmp_run_id,run_id)==0))
      found=1;
    else{
      if(fprintf(tmp_fp,line) < 0){
	printf("Error writing %s, exiting...\n",tmp_file);
	fclose(tmp_fp);    /* close tmp file */
	sprintf(sys_cmd,"rm -f %s",tmp_file);
	system(sys_cmd); 
	fclose(lock_fp);    /* close lock file */
	sprintf(sys_cmd,"rm -f %s",lock_file);
	system(sys_cmd); 
	exit(2);
      }
    }
  }
  
  

  fclose(fp);
  fclose(tmp_fp);
  
  sprintf(sys_cmd,"cp %s %s~",db_file,db_file);
  system(sys_cmd); 
  sprintf(sys_cmd,"rm -f %s",db_file);
  system(sys_cmd); 
  sprintf(sys_cmd,"cp %s %s",tmp_file,db_file);
  system(sys_cmd); 
  sprintf(sys_cmd,"rm -f %s",tmp_file);
  system(sys_cmd); 

  fclose(lock_fp);
  sprintf(sys_cmd,"rm -f %s",lock_file);
  system(sys_cmd);

  if (found==1){
    printf("del_toff_: deleted entry for %s from %s \n",
	   run_id,db_file);
    return 0;
  }else{
    printf("del_toff_: could not find entry for %s in %s \n",
	   run_id,db_file);
    return 1;
  }
}
