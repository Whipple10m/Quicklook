*---------------------------------------------------------------------
*--	GN2GAINS
*---------------------------------------------------------------------
*-- GV @ FLWO 	! Oct 3, 1988 !
*---------------------------------------------------------------------
*-- Produces nitrogen gains from files prepared with REDUCE1.
*-- All the nitrogen gains are stored in hrc.n2gains.
*---------------------------------------------------------------------
*-- Modified on Nov 8, 1988 by GV
*-- Now comments contained in the GRALP header are saved to disk to
*-- provide additional info in choosing the best nitrogen file for the
*-- night.
*-- mjl 18th Nov. 2 inch tubes are now given identical treatment.
*-- GV 22nd Nov, 1988: files whose runid is ge 603 use integer*2
*-- to store the events.
*-- Mod. by GV on 02/10/89: events which have at least one tube below
*-- threshold (=50dc) -after pedestals subtraction- are not included
*-- in the computation. This is done to get rid of Cherenkov events 
*-- which could contaminate the nitrogen data; it also gets rid of 
*-- weak nitrogen flashes.
*---------------------------------------------------------------------
c-- Modified to run on Viktor and to put info into current
c-- working directory. Also, must convert vax to IEEE reals using
c-- vaxflt.c and vaxdbl.c MP920505.
*---------------------------------------------------------------------
c-- Modified for calibration and run information in directory
c-- /vik/observer/quicklook/info    ADK 921105
c---------------------------------------------------------------------
c-- Automatic switching to proper file based on date.  Michael Punch
c-- method.                            ADK  931129
c---------------------------------------------------------------------
c-- Modified 940421 by JB to work with 11m data files
c-- i.e. reads reduce2 format files (i.e. fz2red format) and updates
c-- the 11m data base
*---------------------------------------------------------------------
*-- Modified by JB 940424: event threshold TEMPORARILY removed
*-- to allow preliminary analysis of 11m data with attenuated N2 pulser
*-- and events with bad channels
*---------------------------------------------------------------------
*-- Modified 940916 by JB to read new 10m and 11m data
*-- denoted by prefix 'gt' and 'ge', and modified to automatically
*-- change the protection on a new hrc.n2gains file to rwx for everyone.
*-- If the first two characters are other than 'ge', ten meter data is
*-- assumed.  This is a TEMPRORARY measure, and should be changed later.
*-- Also note that this program currently has no lock mechanism.  This 
*-- should probably be changed at some later time.
*---------------------------------------------------------------------
*-- Modified 940920 by JB.  Appropriate exit code is now returned.  
*---------------------------------------------------------------------
*-- Modified 941015 by JB.  Vax to IEEE conversion only for old data.
*---------------------------------------------------------------------
*---------------------------------------------------------------------
*-- Modified 961205 by JQ to handle variable no's of tubes. The reading
*-- and writing of the database is handled by a library of C routines
*---------------------------------------------------------------------
*-- MAC 970805 Modified to read in the path for the peds and tubelist 
*-- databases from dat.config so that we can write these values to any 
*-- directory we want.  
*----------------------------------------------------------------------
*-- MAC 970918 Increased maximum number of channels to 600
*------------------------------------------------------------------------
*-- MAC 971114 Changed the level for deciding that an event is not
*-- a N2 event to be based on a function of the fraction of the maximum
*-- number of tubes available (i.e., not off from gcpeds) rather than a 
*-- fixed number so that we can automatically correct for tube changes and
*-- also still get reasonable numbers when there are a lot of PMTs turned
*-- off.
*-----------------------------------------------------------------------
*-- PM 000718 : Minor modifications to formats and spacings
*-----------------------------------------------------------------------

	program n2g

	integer max_chan                   ! JQ
	parameter (max_chan=600)           ! as of 961202

	integer utdate,iunit,dunit,krunno,
     +  punit,npmt,event(max_chan),proc8,           !JQ
     +  gunit,cunit
        integer	nrevt,nonzeroevents
	integer*2 shortevent(max_chan)
	real peds(max_chan),pedvar(max_chan),revent(max_chan),mean, meanmean,
     +       gainaccum(max_chan),gainsquared(max_chan),gains(max_chan),
     +       gainvar(max_chan),thres,naccum(max_chan)          !JQ
	character n2file*8,source*2,id*4,n2*2,sn2*2
        character dbfile*80, comdev*20
C	character infdir*11   !This length must be exact (JQ: Why?-use lnblnk)
C	character eminfdir*11
	character infdir*80   !JQ
	character eminfdir*80
	integer*4 yy
	character yyc*2,hy*1
	logical short,exist
*-- JB
	character mode*3,run_id*4
        integer date,code,stdur,ut,st
        integer nevents,noff,tubesoff(max_chan)
        character srcname*20,skyq*6,comms*404,trig*6/' '/
        character*4 pdid
	integer*4 npdid
	real azimuth,elevation,ra,dec
	real*8 duration,mjd,frjd,time1,time2,time3
*	character*20 sthresh
	integer ithrp
        character*80 syscmd	! System command string

*-- MAC Stuff related to reading in dat.config file.

        integer*4 NPATHS
        parameter(NPATHS=4)

        character cfgfile*10              !Config file name
        parameter(cfgfile='dat.config')

        character path(NPATHS)*80,pthcal*80,pthntub*80,line*80
        equivalence (path(2),pthcal),(path(3),pthntub)

        logical eof

*-- JQ  C database functions  
	integer*4 read_peds
        integer*4 read_toff
	integer*4 chk_peds
	integer*4 read_n2gains
	integer*4 chk_n2gains
	integer*4 write_n2gains
	integer*4 del_n2gains

*-- JQ  other necessary variables:
	integer nadc
	integer lnblnk          ! this is a function
	integer*4 teleid
	integer*4 db_status     ! value returned by C database functions 
	integer*4 chi           ! channel i
	character*80 db_file,ntoff_file
	character*4 ped_mode    ! dummy - for reading peds
	integer*4 event_cnt     ! dummy - for reading peds
	real*4 av_val(3)        ! dummy - for reading peds
	integer*4 n_chan        ! dummy - for reading peds
*-- JQ
* New C function usage from fortran:
* i=read_n2gains(run_id,date,db_file,gains,gainvars,code,ped_id,nadc)
*   where  run_id     integer (4 digits max)
*          date       integer*4
*          db_file    character string
*          gains      real*4 array
*          gainvars   real*4 array
*          ped_id     integer*4
*          code       char*3
*          nadc      integer*4
*    the value returned (in i) is for error handling, generally 
*    0 means success, 1 means not found in database, 2 error 
*    (lock file or database file not found)
*
* i=write_n2gains(run_id,date,db_file,gains,gainvars,code,ped_id,nadc)
*   if database does not exist then it is created.     
*
* i=chk_n2gains(run_id,date,db_file)  
*   checks if gains&gainvars exist in database.
*
* i=del_n2gains(run_id,date,db_file)  
*   deletes (all occurrances of) gains&gainvars from database (note: 
*   a backup of the database before deletion of entries is made -
*   it has the same name as the database but with an appended ~)
* 	



*-- Initializations
*-- Mod by JB 940424
*	parameter (thres=50.0)
	data dunit/1/,iunit/6/,punit/2/,gunit/3/,
     +  n2/'n2'/,sn2/'N2'/, 
     +  cunit/1/,comdev/'nitrogen.coms'/,
     +  thres/50.0/,
     +  infdir/'/usr/dbgt/'/,
     +  eminfdir/'/usr/dbge/'/
*-- JQ          ,npmt/109/,

*-- MAC did this config file stuff.
*
* Before anything else, check for the configuration file.  Quit if
* it's not there.
*
      inquire(file=cfgfile,exist=exist)
      if (.not.exist) then
	 write(6,810)cfgfile
	 write(6,830)cfgfile
	 call exit(1)
      end if
*
* Open the configuration file
*
      open(10,file=cfgfile,status='old',readonly)
*
* The parameterize section of cfgfiledat.config is denoted by a line with
* "par" as the first three chracters.  Find this.
*
      line = '   '
      do while (line(1:3) .ne. 'par')
	 call getline(10,line,eof)
         if (eof) then 
           write(*,820)cfgfile
           write(*,830)cfgfile
           call exit(1)
         end if
      end do
*
* Parameterize section of cfgfile found.
*
* Lines beginning with a space or a # ar comment lines.  Skip these and
* also the path for the data file so that we can read the path for the 
* calibration files which is where the pedestals are.
*
* Loop until the proper number of non-comment lines are read
*
      i=1
      do while(i.le.NPATHS)
         call getline(10,line,eof)
         path(i) = line
* Insert a trailing slash in path names
         j=lnblnk(path(i))
         if ((path(i).ne.'default').and.
     >        (path(i)(j:j).ne.'/')) path(i)(j+1:j+1)='/'
*	 write(*,866)'paths: ',i,path(i)
	 i = i + 1
      end do

* Close the configuration file
      close(10)

* Print out the database paths
      write(*,860) pthcal
      write(*,861) pthntub

*-- Get date, file name and pedestals id
*	call intarg(1,utdate)
*	call getarg(2,n2file)
	call getarg(1,n2file)
*	call getarg(3,sthresh)
*	read(sthresh,'(f10.2)') thres
*       write(*,*) 'thres: ',thres
        thres = 50.0
	source=n2file(1:2)
	id=n2file(3:6)
	pdid=n2file(3:6)
	read(pdid,'(i4)') npdid
*-- If file name is like n2xxxx.r or ccxxxx.r where cc=gt,ge,or gs then go on
*-- this is TEMPORARY 
	if((source.eq.'n2').or.(source.eq.'gt').or.(source.eq.'ge').or.
     +   (source.eq.'gs'))then
	  open(dunit,file=n2file,form='unformatted',status='unknown')
*-- New read for red2 format


	  read (dunit)run_id,nevents,duration,stdur,mode,srcname,date,mjd,
     +     frjd,ra,dec,ut,st,azimuth,elevation,skyq,comms
	  read(run_id,'(i4)')krunno
*-- Hopefully the "date" read in will be >1,000,000 for years after Y2K
*-- so all date related logic should still work (i said hopefully) SJF 990906
	  utdate = date
	  if(utdate.eq.0)then
5	    write(6,'(''Enter UT date, and nitrogen file name'',
     +      '' : '',$)')
	    read(5,*)utdate,n2file
*
* 941015 JB -- changed from n2file.eq.'' to lnblnk(n2file).lt.2 to
* be compatible with new version of DEC fortran
*
	    if((utdate.eq.0).or.(lnblnk(n2file).lt.2))goto 5
	  endif

*         write(*,*)'Run ID     : ',run_id
*         write(*,*)'Source     : ',srcname
*         write(*,*)'Date       : ',date
*	  write(*,*)'MJD        : ',mjd
*         write(*,*)'No. Events : ',nevents
          
          if(utdate.lt.940800) then
	     call vaxdbl(duration)
          end if
*         write(*,*)'Duration   : ',duration
*         write(*,*)'ST Duration: ',stdur
*         write(*,*)'Comments   : '
*         write(*,*) comms

*	  read(dunit)knoev,krunl,krec,kelt,fname,krunno,ksmode,sname,kdate,
*    +    kra,kdec,kut,kst,kovflow,az,el,trig,observ,skyq,comms
*	  call vaxflt(az)
*	  call vaxflt(el)
*


* ** Figure out which pedestal file we want; then open it
*
        yy=int((float(date)+0.5)/10000.)
        write(yyc,'(i2.2)')mod(yy,100)
        if (utdate-yy*10000.lt.600) then
           hy='a'
        else
           hy='b'
        end if


*-- JQ 961205 - fitst check if gains are already in database:
*-- MAC 970805 - Change to use pthcal instead of infdir and eminfdir
	  if(source.eq.'ge') then
             if (pthcal.eq.'default') pthcal=eminfdir
             if (pthntub.eq.'default') pthntub=eminfdir
*MAC	     db_file=eminfdir(1:lnblnk(eminfdir))//'hrc'//yyc//hy//
*MAC     +	                                 '.n2gains'//char(0)
	  else
             if (pthcal.eq.'default') pthcal=infdir
             if (pthntub.eq.'default') pthntub=infdir
*MAC	     db_file=infdir(1:lnblnk(infdir))//'hrc'//yyc//hy//
*MAC     +                             '.n2gains'//char(0)
	  end if
          db_file=pthcal(1:lnblnk(pthcal))//'hrc'//yyc//hy//
     +         '.n2gains'//char(0)
          ntoff_file=pthntub(1:lnblnk(pthntub))//'hrc'//yyc//hy//
     +         '.ntubelist'//char(0)

	  db_stat=chk_n2gains(npdid,utdate,db_file)
	  if (db_stat.EQ.0) then
	     write(6,*)'Gains already exist, exiting...'
	     write(6,*)
	     call exit(0)
	  endif
*----
*-- Get pedestals:	  

	if(source.eq.'ge') then
* JB changed back to teleid=2 for 11m, teleid=1 for 10m
* Use pthcal instead of eminfdir
	   teleid=2
*MAC	   db_file=eminfdir(1:lnblnk(eminfdir))//'hrc'//yyc//hy//
*MAC     +                   '.cpeds'//char(0)
*
*          open(unit=punit,file=eminfdir//'hrc'//yyc//hy//'.cpeds',
*     +     status='unknown')
*           dbfile=eminfdir//'hrc'//yyc//hy//'.cpeds'
	else
	   teleid=1
*MAC	   db_file=infdir(1:lnblnk(infdir))//'hrc'//yyc//hy//
*MAC     +                   '.cpeds'//char(0)

*          open(unit=punit,file=infdir//'hrc'//yyc//hy//'.cpeds',
*     +     status='unknown')
*           dbfile=infdir//'hrc'//yyc//hy//'.cpeds'
	endif
        call num_chan(teleid,date,nadc,npmt)
        db_file=pthcal(1:lnblnk(pthcal))//'hrc'//yyc//hy//
     +       '.cpeds'//char(0)
        ntoff_file=pthcal(1:lnblnk(pthcal))//'hrc'//yyc//hy//
     +       '.ntubelist'//char(0)
	         ! char(0) because strings must be terminated by '\0' in C


*-- Get pedestals
* JB -- pdid gives the run number of the
* nitrogen file.  This means that unlike the old version of n2gains,
* we use injected pedistals calculated previously using
* the program cpeds from the nitrogen file itself,
* rather than using electronic pedistals from
* that night. 
*	call getpeds(punit,date,pdid,peds,pedvar,*25)
*--JQ	call getnpeds(punit,pdid,peds,pedvar,*25)
*	type *,'peds: ',peds
*--JQ	close(punit)

	n2file=n2file(1:lnblnk(n2file))
	db_status=read_peds(n2file//char(0),utdate,db_file,peds,pedvar,
     +                  event_cnt,av_val,ped_mode,nadc)
	if (db_status.EQ.1) then
	   goto 25
	else if (db_status.EQ.2) then
	   goto 25
	end if
	db_status=read_toff(n2file//char(0),noff,tubesoff,ntoff_file)


*-- Decide here if events are long or short
* JB -- convert from test run_id to integer run number
	if(source.eq.'ge') then
	  short=.TRUE.
	else
	  if((krunno.lt.603).and.(date.lt.881101))then
	    short=.FALSE.
	  else
	    short=.TRUE.
	  endif
	endif
*	write(*,*)'nevents: ',nevents

*-- Get mean of mean adc value of all tubes for all events - SJF June 99
	meanmean=0.0
	nonzeroevents=0

*-- Loop on the events
	  do 10 i=1,nevents
	    if(short)then
              if(date.lt.940800) then
	        read(dunit,end=15)code,time1,(shortevent(chi),chi=1,nadc) !JQ
*                 call vaxdbl(time1)
              else
		 read(dunit,end=15)code,time1,time2,time3,
     +	              (shortevent(chi),chi=1,nadc) !JQ
              endif

*	      if (i.lt.10) then
*                 write(*,*)i,':','code: ',code,
*     >                ' event: ',shortevent
*              end if
*-- Map "short" event into the "long" one (software compatibility).
	      do k=1,nadc         !JQ
		  event(k)=shortevent(k)
              end do
	    else
              if(date.lt.940800) then
	        read(dunit,end=15)code,time1,(event(chi),chi=1,nadc)  !JQ
*	        call vaxdbl(time1)
              else
	        read(dunit,end=15)code,time1,time2,time3,
     +                         (event(chi),chi=1,nadc)    !JQ
              endif
	    endif
*-- Process hypothetical nitrogen code 9's (it has happened before).
*	    write(*,*)'Raw ADC counts: ',i
*	    call show109(iunit,event) 

	    if(code.ne.8.and.code.ne.9)goto 10

*-- Skip saturated events. Subtract pedestals
*-- Skip Cherenkov and weak nitrogen events
	    ithrp = 0

	    do k=1,npmt  !JQ
*-- TEMPORARY Modification by JB 940424 to
*-- skip events where at least 10 adc values lie below threshold -
*-- a more forgiving criteria allowing for some bad channels
*-- as in the preliminary 11m data
		if(event(k).gt.1023)then
		  isat=isat+1
*	          type *,'pmt: ',k,' saturated'
		  goto 10
		endif
		revent(k)=event(k)-peds(k)
		
		if(revent(k).lt.thres)then
	          ithrp=ithrp+1
*	          type *,'pmt: ',k,' less than threshold'
*Change this to account for periods where lots of tubes are off
	          if (ithrp.gt.
     >                 (nint(0.25*float(nadc-noff))+noff)) then 
    	            ithr = ithr + 1
		    goto 10
	          endif
		endif
	     end do

*-- Get mean adc value of all tubes 
	    mean=0.0

*
* Modified 941017 by JB -- old version calculated mean of all 109 PMTs
* regardless of zeroed channels.  This might have led to a shift in the
* relative gains depending on how many PMTs were turned off, and could
* thus have resulted in a systematic error in the size parameter.  In
* this version, average only over nonzero channels making the assumption
* that for a nitrogen run revent(k) should never fluctuate down to zero.
*
            nrevt=0.0
	    do k=1,npmt		!JQ
	       if(revent(k).gt.thres) then ! If nonzero
		  mean=mean+revent(k)
		  nrevt = nrevt + 1.0
	       endif
	    end do
	    
*
* 941017 JB -- was: mean=mean/109.0
*
	    mean=mean/nrevt 
	    meanmean=meanmean+mean
	    nonzeroevents=nonzeroevents+1

*-- Gains are derived by accumulating the ratio of the mean adc value
*-- for an event to the ped-subtracted adc value for a given pmt
*-- Accumulate gains and gains squared tube by tube
	    do k=1,npmt  !JQ
*
* Modified 941017 by JB -- Only accumulate for events with nonzero
* adc values
*
              if(revent(k).gt.thres) then
	        gainaccum(k)=gainaccum(k)+mean/revent(k)
	        gainsquared(k)=gainsquared(k)+(mean/revent(k))**2
*
* Modified 941017 by JB -- Allow for the possiblility of tubes with
* no signal.  Only count events with nonzero adc values instead of all
* code-8 events
*
                naccum(k) = naccum(k) + 1.0
              endif
            end do
 	    proc8=proc8+1
10	  continue
	  close(dunit)

	  meanmean=meanmean/nonzeroevents;
	  write(*,*)'Mean ADC value ',meanmean

*-- Now get gains and their variances
15	  if(proc8.le.0)then
	    write(6,200)n2file
	    call exit(1)
	  else
	    do i=1,npmt
              if(naccum(i).gt.0.001) then
*
* Modified 941017 by JB -- was: .../proc8 now: .../naccum(i)
*
		gains(i)=gainaccum(i)/naccum(i)
		gainvar(i)=sqrt(gainsquared(i)/naccum(i)-gains(i)**2)
              end if
            end do
	  end if

*-- Trap tubes which could have bad HV (wild gain)
	  do i=1,npmt
	    if(gains(i).gt.5.0.or.gains(i).lt.0.1) then
              write(6,210)i
	      type *,' NG ** Gain is being set to zero to avoid'
	      type *,'       error in data base. Tube must be turned off!'
              gains(i) = 0.0	! Set gain to zero to effectively zero
              gainvar(i) = 1.0	! tube in case user is too dumb to do
				! this him/herself.  This is useful
				! for automated quicklook analysis
	    end if
          end do

*-- Update the gains file. The gains file is opened in append mode
*
* ** If file exists, open it for append; else, create it
*
* JB 
*
*JQ: this is now all taken care of by C database routines
*	if(source.eq.'ge') then
*             inquire(file=eminfdir//'hrc'//yyc//hy//'.n2gains',exist=exist)
*          if (exist) then
*            open(unit=gunit,file=eminfdir//'hrc'//
*     +       yyc//hy//'.n2gains',access='append',status='unknown')
*          else
*            open(unit=gunit,file=eminfdir//'hrc'//
*     +       yyc//hy//'.n2gains',status='new')
*            syscmd = 'chmod 777 '//eminfdir//'hrc'//yyc//hy//'.n2gains'
*            call system(syscmd)
*          end if
*	else
*        inquire(file=infdir//'hrc'//yyc//hy//'.n2gains',exist=exist)
*          if (exist) then
*            open(unit=gunit,file=infdir//'hrc'//
*     +       yyc//hy//'.n2gains',access='append',status='unknown')
*          else
*            open(unit=gunit,file=infdir//'hrc'//
*     +       yyc//hy//'.n2gains',status='new')
*            syscmd = 'chmod 777 '//infdir//'hrc'//yyc//hy//'.n2gains'
*            call system(syscmd)
*          end if
*	endif
*	write(gunit,300)date,sn2,id,npdid
*	write(gunit,400)gains
*	write(gunit,400)gainvar
*	close (gunit)

*MAC Use pthcal instead of infdir
	  if(source.eq.'ge') then
*MAC	   db_file=eminfdir(1:lnblnk(eminfdir))//'hrc'//yyc//hy//
*MAC     +                   '.n2gains'//char(0)
	   db_file=pthcal(1:lnblnk(pthcal))//'hrc'//yyc//hy//
     +                   '.n2gains'//char(0)
	else
*MAC	   db_file=infdir(1:lnblnk(infdir))//'hrc'//yyc//hy//
*MAC     +                   '.n2gains'//char(0)
	   db_file=pthcal(1:lnblnk(pthcal))//'hrc'//yyc//hy//
     +                   '.n2gains'//char(0)
	   
	end if


	db_stat=write_n2gains(npdid,utdate,db_file,gains,gainvar,
     +            sn2,npdid,nadc)
	if (db_stat.EQ.2) then
	   write(6,*)'Error in write_n2_gains, exiting...'
	   write(6,*)
	   call exit(1)
	end if

*-- Update the comment file
	open(cunit,file=comdev,access='append',status='unknown')
	write(cunit,600)n2file,utdate,proc8,isat,ithr
	write(cunit,611)comms
	close(cunit)

*-- Message to standard output
	write(6,601)dbfile(1:30),n2file,utdate,proc8,isat,ithr
	write(6,611)comms
	else
	  write(6,100)n2file
	endif
	goto 35

*-- Peds not found
25	write(6,500)n2file,pdid
        call exit(1)
35	call exit(0)

*-- Formats
100	format('  NG ** This is not a nitrogen file: ',a)
200	format('  NG ** No code 8s processed: ',a)
210	format('  NG ** Wild gain factor in channel ',i3)
300	format(i6,1x,a2,1x,a4,1x,i4)
400	format(16f5.2)
500	format('  NG ** Peds not found for given file: ',a,' id: ',a)
531	format('SOURCE: ',a,' FILE: ',a) 
600	format('  NG ** File: ',a,' Date: ',i6,/,
     +       '           Events processed: ',i6,' saturated: ',i6,
     +       ' sub-threshold: ',i6)
601     format('  NG ** File: ',a30,' Run: ',a6,' Date: ',i6,/,
     +         ' Events processed: ',i6,' saturated: ',i6,
     +       ' sub-threshold: ',i6)             
610	format(' Comments: ',/,5(80a,/),a4)
611     format(' Comments: ',/,5(80a,/),a4)
 810    format('  NG ** Configuration file ',a,' not found.')
 820    format('  NG ** Section "par" in ',a,' not found.')
 825    format('  NG ** Section "cut" in ',a,' not found.')
 830    format(  '  A list of standard config files may be found in ',
     &       /,'  /usr/sc/config.  Copy one to ',a,' in',
     &       /,'  your current working directory.',
     &       /,'Exiting...')
 840    format('  NG ** Unable to read enough lines from ',a)
 850    format('  NG ** Unable to read ',a)
 860    format('   NG ** Using ',a30,' for database path')
 861    format('   NG ** Using ',a30,' for tubelist path')

	end
