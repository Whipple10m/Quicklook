	program gcpeds
*
* ** Uses code1 and code2 events to accumulate the pedestals and
* ** variances for a run.  Rejects event if any tube has more than
* ** 50 counts.  Method of calculation is same as that used in the
* ** traditional 'pedestals' program which operated on code 8's in
* ** pedestals files taken at beginning and end of night.
*
* ** It has some other features:  it uses '/tmp/cpeds.lock' to
* **  lock out additional users when program is being used;  it
* **  finds the median value of pedestal variance for each bank
* **  of discriminators and writes that on first line of entry
* **  into pedestal file (av_val);  the pedestal file which is
* **  appended to is called infdir/hrc(yr/a or b).cpeds.  infdir
* **  is currently 14 chars long and is /usr/users/db.

* VC 940126 - based on preliminary version by ADK 
* 
* JB 940421 - Modified by JB to work with 10m and 11m data
* ** This version checks to see if the first two characters of
* ** the filename are 'em'.  This prefix is used to flag 11m data.
* ** If this prefix is detected, the database path is taken to
* ** be /usr/users/db11 rather than /usr/users/db.
* ** Note: This program should also be run on all n2 files prior
* ** to running n2gains, since these pedistals are more appropriate
* ** (and more readily available) than electronic pedistals.
*  
* JB 940917 - Modified to allow for a reorganization of 10m and
* ** 11m database directories following changeover to GDF format
* ** data.  Curently the prefix 'ge' means 11m data and results
* ** in data being stored in /home/dbge.  Otherwise 10m data
* ** is assumed, and the database in /home/dbgt is used.
* ** The program has also been changed to allow any user to create
* ** the first code-peds data base file.
*
* JB 941027 - Modified by JB to turn off phototubes with pedestal 
* ** variances which exceed some lower or upper limit.
*
* JB 950414 - Modified by JB to read new data format with three
* ** time variables: gpsutc, livetime and osctime, even though
* ** these are not used in this program
*
* MAC 950930 - Modified by MAC to avoid problem with tubes above
* ** channel 109 exceeding threshold
*
* JB 951017 - Modified so that gcpeds always returns with proper
* ** error code.
*
* JQ 961201 - Modified to use C functions for reading/writing
* database. The pedestal and pedvar array sizes been increased 
* to handle the extra channels which have been added to the camera
* and will be in the data stream from 961202. The C routines
* handle variable no's of channels and can read the old database
* entries.
*
* MAC 970805 Modified to read in the path for the peds and tubelist 
* databases from dat.config so that we can write these values to any 
* directory we want.  
*
* MAC 970918 Changed maximum number of tubes in declarations below to
* 600.
*
* MAC 971113 Change choice of 1/4, 1/2, 3/4 values of pedestals to
* be based on m_pmts minus the number of zero pedestals to reduce
* the effects of bright fields or periods when we have a lot of dead tubes.
*
* PM 000718 : Minor modifications to formats and spacings
*
	implicit none

	integer m_pmts,n_pmts,n_fail
*-- JQ	parameter (m_pmts=109, n_pmts=120)
	integer max_chans
	parameter (max_chans=600)
	
        integer*4 MAXOFF	! Maximum number of tubes which can
        parameter (MAXOFF=600)	! be turned off. Changed from 15 by JQ 970506
	real rsig
	integer i,j,k,l,m,istep,imin
	real av_val(3)
	logical exist
	character*11 infdir
	character*11 eminfdir
	character*9  oldinfdir
	integer yy
	character yyc*2,hy*1
	integer*2 pmts(max_chans)
	integer date,code,stdur,ut,st
	integer evnt_cnt,evnt_chr,nevents
	character source*20,skyq*6,comms*404,trig*6/' '/
	character*4 pdid
	real azimuth,elevation,ra,dec
	real*8 duration,mjd,frjd,time1,time2,time3
	integer*4 sum(max_chans),sqrsum(max_chans)
	real*4 avg(max_chans),wrms(max_chans),rms(max_chans),count
	character mode*3,run_id*4,run_name*6
*-- JB - string denoting telescope
	character telsrc*2
        character*80 fname
	character*80 syscmd	! system command string
	integer	pcutoff		! Cut off the tail of the ped dist
				! at this value.
        character*80 ftubelist	! Database file containing tubes to be 
                                ! turned off
        real*4 flowthresh,fupthresh ! Fraction of average pedestal variance
                                  ! giving the upper and lower limit for
                                  ! turning off tubes
	real*4 lowerlimit,upperlimit
        integer*4 ntubesoff(MAXOFF),nnt
        data flowthresh /0.6/	! Threshold to detect turned-off tubes
        data fupthresh /1.5/	! Threshold for stars
*-- SF	data pcutoff /200/	! TEMPORARY should be around 50
	data pcutoff /75/	! TEMPORARY should be around 50

	data infdir/'/usr/dbgt/'/
*-- JB
	data eminfdir/'/usr/dbge/'/
        data oldinfdir /'/usr/db/'/

*-- MAC Stuff related to reading in dat.config file.

        integer*4 NPATHS
        parameter(NPATHS=4)

        character cfgfile*10              !Config file name
        parameter(cfgfile='dat.config')

        character path(NPATHS)*80,pthcal*80,pthntub*80,line*80
        equivalence (path(2),pthcal),(path(3),pthntub)

        logical eof

*-- JQ  C database functions:
	integer*4 read_peds
        integer*4 write_peds
        integer*4 chk_peds
        integer*4 del_peds
	integer*4 chk_toff
	integer*4 write_toff

*-- JQ  other necessary variables:
	integer lnblnk
	integer*4 teleid
	integer*4 db_status     ! value returned by C database functions 
	integer*4 chi           ! channel i
	integer*4 ped_events
	character*80 db_file
*-- MAC Data variables
        logical first
        data first/.true./
*-- JQ
* New C function usage from fortran:
* i=read_peds(run_name,date,db_file,peds,pedvars,event_cnt,av_val,
*             mode, nadc)
*   where  run_name     character string (including gt/ge part)
*          date       integer*4
*          db_file    character string
*          peds       real*4 array
*          pedvars    real*4 array
*          event_cnt  integer*4
*          av_val     real*4 array (3 elements)
*          mode       char*3
*          nadc      integer*4
*    the value returned (in i) is for error handling, generally 
*    0 means success, 1 means not found in database, 2 error 
*    (lock file or database file not found)
*
* i=write_peds(run_name,date,db_file,peds,pedvars,event_cnt,av_val,
*              mode,nadc)
*   if database does not exist then it is created.     
*
* i=chk_peds(run_name,date,db_file)  
*   checks if peds&pedvars exist in database.
*
* i=del_peds(run_name,date,db_file)  
*   deletes (all occurrances of) peds&pedvars from database (note: 
*   a backup of the database before deletion of entries is made -
*   it has the same name as the database but with an appended ~)
*

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
      write(*,860) pthcal
      write(*,861) pthntub
      
*--
*-- Get name of file
	call getarg(1,run_name)
1       if (run_name .eq. ' ')  then
2	  write(6,11)
	  read(5,12)i,run_name
	  if (i .ne. 6)  then
	    write(6,13)run_name
	    go to 2
	  end if
	end if
	inquire(file=run_name,exist=exist)
	if (.not. exist)  then
	  write(6,14)run_name
	  run_name=' '
	  go to 1
	end if

	open (unit=1,file=run_name,status='old'
     &   ,form='unformatted',readonly)

	read (1)run_id,nevents,duration,stdur,mode,source,date,mjd,
     &   frjd,ra,dec,ut,st,azimuth,elevation,skyq,comms

*-- Hopefully the "date" read in will be >1,000,000 for years after Y2K
*-- so all date related logic should still work (i said hopefully) SJF 990906
        if(date.lt.940800) then 
	  call vaxdbl(duration)
        end if

 
*-- Now check to see whether pedestals already exist for this file
*-- Open the relevant pedestal file.
c-- Simplified MP920505.

	
        yy=int((float(date)+0.5)/10000.)
        write(yyc,'(i2.2)')mod(yy,100)
        if(date-yy*10000.lt.600)then
          hy='a'
        else
          hy='b'
        end if

*-- JB
	telsrc = run_name(1:2)
	if(telsrc.eq.'ge') then
*-- MAC  Change to use pthcal and pthntub instead of eminfdir.
           if (pthcal.eq.'default') pthcal=eminfdir
           if (pthntub.eq.'default') pthntub=eminfdir
           
*-- JQ
*-- JB 961209 changed back to teleid=2 for 11m, 1 for 10m.
	   teleid=2
	   call num_chan(teleid,date,n_pmts,m_pmts)
*MAC	   db_file=eminfdir(1:lnblnk(eminfdir))//'hrc'//yyc//hy//
*MAC     +                   '.cpeds'//char(0)
	   db_file=pthcal(1:lnblnk(pthcal))//'hrc'//yyc//hy//
     +                   '.cpeds'//char(0)
	         ! char(0) because strings must be terminated by '\0' in C
*-- JQ          ftubelist=eminfdir//'hrc.ntubelist'
*MAC           ftubelist=eminfdir(1:lnblnk(eminfdir))//'hrc'//yyc//hy//
*MAC     +                   '.ntubelist'//char(0)       
           ftubelist=pthntub(1:lnblnk(pthntub))//'hrc'//yyc//hy//
     +                   '.ntubelist'//char(0)       
	   db_status=chk_peds(run_name,date,db_file)
	   if (db_status.EQ.1) then
	      goto 101
	   else if (db_status.EQ.0) then
	      write(6,*)' CP ** Pedestal values already exist for ',
     +        run_name,'- halting execution.'
	      call exit(0)
C	   else if (db_status.EQ.2) then
C	      write(6,*)'  CP ** Error: Database not found or locked'
C	      write(6,*)'  CP ** Halting execution.'
C	      call exit(0)
	   endif

*-- Below commented out by JQ 
*          inquire(file=eminfdir//'hrc'//yyc//hy//'.cpeds',exist=exist)
*          if(exist)then 
*	    open(2,file=eminfdir//'hrc'//yyc//hy//'.cpeds',status='unknown')
*
* 941015 JB changed from 3:7 to 3:6
*
*	    pdid=run_name(3:6)
*	    call getnpeds(2,pdid,avg,rms,*101)
*-- Pedestal values for the file have been found .. halt execution
*	    close(2)
*	    write(6,*)'  CP ** Pedestal values already exist for ',run_name
*	    write(6,*)'  CP ** Halting execution.'
*            call exit(0)
*	  end if
	else if (telsrc.eq.'gt') then
*-- MAC Change to use pthcal and pthntub instead of infdir
           if (pthcal.eq.'default') pthcal=infdir
           if (pthntub.eq.'default') pthntub=infdir
	   teleid=1
	   call num_chan(teleid,date,n_pmts,m_pmts)
*debug	   print*,'npmts=',n_pmts
*debug	   print*,'date=',date
	   print*,' '
*MAC	   db_file=infdir(1:lnblnk(infdir))//'hrc'//yyc//hy//
*MAC     +                   '.cpeds'//char(0)
	   db_file=pthcal(1:lnblnk(pthcal))//'hrc'//yyc//hy//
     +                   '.cpeds'//char(0)
* 970506 JQ changed ftubelist to be half yearly.
*MAC           ftubelist=infdir(1:lnblnk(infdir))//'hrc'//yyc//hy//
*MAC     +                   '.ntubelist'//char(0)       
           ftubelist=pthntub(1:lnblnk(pthntub))//'hrc'//yyc//hy//
     +                   '.ntubelist'//char(0)       
	   db_status=chk_peds(run_name,date,db_file)
	   if (db_status.EQ.1) then
	      goto 101
	   else if (db_status.EQ.0) then
	      write(6,*)' CP ** Pedestal values already exist for ',
     +        run_name,' - halting execution.'
	      write(6,*)
	      call exit(0)
C	   else if (db_status.EQ.2) then
C	      write(6,*)'  CP ** Error: Database not found or locked'
C	      write(6,*)'  CP ** Halting execution.'
C	      call exit(0)
	   endif
*          inquire(file=infdir//'hrc'//yyc//hy//'.cpeds',exist=exist)
*          if(exist)then
*	    open(2,file=infdir//'hrc'//yyc//hy//'.cpeds',status='unknown')
*
* 941015 JB changed from 3:7 to 3:6
*
*	    pdid=run_name(3:6)
*	    call getnpeds(2,pdid,avg,rms,*101)
*-- Pedestal values for the file have been found .. halt execution
*	    close(2)
*	    write(6,*)'  CP ** Pedestal values already exist for ',run_name
*	    write(6,*)'  CP ** Halting execution'
*            call exit(0)
*	  end if
	else			! Telescope id = 10m
	   write(6,*),'  CP ** Assuming old 10m data from prefix: ',telsrc
*-- MAC Use pthcal and pthntub instead of oldinfdir
           if (pthcal.eq.'default') pthcal=oldinfdir
           if (pthntub.eq.'default') pthntub=oldinfdir
	   teleid=1
	   call num_chan(teleid,date,n_pmts,m_pmts)
	   print*,'npmts=',n_pmts
*MAC	   db_file=oldinfdir(1:lnblnk(oldinfdir))//'hrc'//yyc//hy//
*MAC     +                        '.cpeds'//char(0)
*MAC           ftubelist=oldinfdir(1:lnblnk(oldinfdir))//'hrc'//yyc//hy//   
*MAC     +                   '.ntubelist'//char(0)       
	   db_file=pthcal(1:lnblnk(pthcal))//'hrc'//yyc//hy//
     +                        '.cpeds'//char(0)
           ftubelist=pthntub(1:lnblnk(pthntub))//'hrc'//yyc//hy//   
     +                   '.ntubelist'//char(0)       
*--JQ           ftubelist=oldinfdir//'hrc.ntubelist'
	   db_status=chk_peds(run_name,date,db_file)
	   if (db_status.EQ.1) then
	      goto 101
	   else if (db_status.EQ.0) then
	      write(6,*)' CP ** Pedestal values already exist for ',
     +        run_name,' - halting execution'
	      call exit(0)
C	   else if (db_status.EQ.2) then
C	      write(6,*)'  CP ** Error: Database not found or locked'
C	      write(6,*)'  CP ** Halting execution'
C	      call exit(0)
	   endif
*          inquire(file=oldinfdir//'hrc'//yyc//hy//'.cpeds',exist=exist)
*          if(exist)then 
*	    open(2,file=oldinfdir//'hrc'//yyc//hy//'.cpeds',status='unknown')
*
* 941015 JB changed from 3:7 to 3:6
*
*	    pdid=run_name(3:6)
*	    call getnpeds(2,pdid,avg,rms,*101)
*-- Pedestal values for the file have been found .. halt execution
*	    close(2)
*	    write(6,*)'  CP ** Pedestal values already exist for ',run_name
*	    write(6,*)'  CP ** Halting execution.'
*            call exit(0)
	end if
*	end if
      

*-- Pedestal values for the file do not exist so continue execution
101     continue

*-- JQ	close(2)

	do i=1,n_pmts
	  sum(i)=0
	  sqrsum(i)=0
	end do
	evnt_cnt=0
	evnt_chr=0

	do k=1,nevents
          if(date.lt.940800) then
            read(1,end=777)code,time1,(pmts(chi),chi=1,n_pmts)
          else
            read(1,end=777)code,time1,time2,time3,(pmts(chi),
     +                          chi=1,n_pmts)
          endif
          if ((code.eq.1).or.(code.eq.2)) then
c	     type *,'event:',k,'pmts:',pmts
	     evnt_cnt=evnt_cnt+1
*
* MAC - Modified to only check inner 109 PMTs 
*
             do i=1,m_pmts
                if (pmts(i).ge.pcutoff) then
		   type *,'  CP ** Pedestal fluctuated above cutoff: ',pcutoff
		   type *,'  CP ** this event not counted in ped average.'
                   evnt_chr=evnt_chr+1
                   go to 776
                end if
             end do
	     do i=1,n_pmts
	        sum(i)=sum(i)+pmts(i)
	        sqrsum(i)=sqrsum(i)+pmts(i)*pmts(i)
             end do
          end if
776       continue
	enddo
777     continue
	if ((k-1).ne.nevents) then
	   type *,' '
	   type *,'  CP ** Event number mismatch ',k,nevents
	   type *,' '
	end if

	count=float(evnt_cnt-evnt_chr)
	if (count.lt.1.2) then
	   type *,'ERROR: Only ',evnt_cnt-evnt_chr,' events in ped avg'
	   type *,'Not able to calculate peds.  Exiting...'
	   call exit(1)
	end if

	do i=1,m_pmts
	  avg(i)=float(sum(i))/count
	  rms(i)=float(sqrsum(i))/count
	  rms(i)=rms(i)-avg(i)*avg(i)
	  rms(i)=sqrt(rms(i))
	  wrms(i)=rms(i)
	end do

*  the following code finds the 1/4, median, and 3/4 values for the 
*  pedvars, but only for the tubes which have peds above the minimum
*  cutoff (i.e., those with HV on) and stores the values in av_val

*
*  WARNING: Here rms() is replaced by a sorted version of rms().
*  wrms() is the array that should be written to the ped file
*
        imin=0
	do i=1,m_pmts	! Bubble sort: for each PMT (1-109)
          rsig=-1.0
*
* l=0 not allowed -- fix this later
*
          l=0
*MAC          do j=i,m_pmts	! From tube m_pmt-i+1 to tube 1
*MAC            m=m_pmts+1-j
          do m=m_pmts+1-i,1,-1	! From tube m_pmt-i+1 to tube 1
            if (rms(m) .gt. rsig)  then
                rsig=rms(m)	! Find the maximum
                l=m		! and save the index of rms()
            end if
          end do
          m=m_pmts+1-i		! set index m to i^th tube from end
*
* Modified by JB 941017 -- added condition if(l.ne.0)
*
          if(l.ne.0) then
            rms(l)=rms(m)	! If any tube > 0, swap the two
                                ! rms values so that max val floats
                                ! to rms(m), and rms(l) gets old rms(m)
                                ! in the next iteration, ignore rms(m)
          end if
          rms(m)=rsig
          if ((rms(m).lt.flowthresh).and.first) then
             imin=m
             first=.false.
          endif
        end do
*
* Actually av_val() seems to contain three values from the sorted 
* pedestal variances, and not the average for each bank.  I really
* don't know why VC does this.
*
        istep=nint(float(m_pmts-imin)/4.)
        do i=1,3
          av_val(i)=rms(istep*i+1+imin)
        end do


*-- JQ
*-- C database routines take care of lock file
***
***-- Before proceeding further check to see whether anyone else is 
***-- running "cpeds" by checking to see whether a lock file
***-- exists in the /tmp directory. If so then sleep for 20 seconds
***-- and then retry. Process aborts after 6 retries.
***
***5       inquire(file='/tmp/cpeds.lock',exist=exist)
***
***-- If no lock file then create one that is write accessible to all.
***        if(.not.exist) then
***           open(3,file='/tmp/cpeds.lock',status='unknown')
***           call system('chmod a+rw /tmp/cpeds.lock')
***
***        else

***-- Lock file present .. sleep for 20 seconds and then retry
***           write(6,*)
***     &      'cpeds locked by another user .. will retry in 20 s'
***           call system('sleep 20')
***           n_fail=n_fail+1
***           if(n_fail.gt.5)then
***             write(6,*)'  CP ** WARNING: cpeds waited for 2min'
***             write(6,*)'        process still locked, continuing'
***             write(6,*)'        anyway - ped values may be lost'
****            call exit(1)
***             goto 10
***           endif
***           goto 5
***        endif    !??????
***10	continue
***
***
*-- Open the output file in append mode.
*-- Get pedestals
***	if(telsrc.eq.'ge') then
***	  inquire(file=eminfdir//'hrc'//yyc//hy//'.cpeds',exist=exist)
***	  open(2,file=eminfdir//'hrc'//yyc//hy//'.cpeds'
***     &     ,access='append',status='unknown')
***	  syscmd = 'chmod 777 '//eminfdir//'hrc'//yyc//hy//'.cpeds'
***          fname=eminfdir//'hrc'//yyc//hy//'.cpeds'
***          ftubelist=eminfdir//'hrc.ntubelist'
***        else if(telsrc.eq.'gt') then
***	  inquire(file=infdir//'hrc'//yyc//hy//'.cpeds',exist=exist)
***	  open(2,file=infdir//'hrc'//yyc//hy//'.cpeds'
***     &      ,access='append',status='unknown')
***	  syscmd = 'chmod 777 '//infdir//'hrc'//yyc//hy//'.cpeds'
***          fname=infdir//'hrc'//yyc//hy//'.cpeds'
***          ftubelist=infdir//'hrc.ntubelist'
***	else 
***	  inquire(file=oldinfdir//'hrc'//yyc//hy//'.cpeds',exist=exist)
***   	  open(2,file=oldinfdir//'hrc'//yyc//hy//'.cpeds'
***     &     ,access='append',status='unknown')
***          syscmd = 'chmod 777 '//oldinfdir//'hrc'//yyc//hy//'.cpeds'
***          fname=oldinfdir//'hrc'//yyc//hy//'.cpeds'
***          ftubelist=oldinfdir//'hrc.ntubelist'
***	endif
***	if(.not. exist)	then	! If this file was just created
***	  call system(syscmd) 
***	end if


*-- Check for abnormal pedestals to avoid software crashes
	do i=1,n_pmts
	  if(avg(i).ge.1000.0) then
	    avg(i)=999.999
	    write(*,50) i
	  end if
	  if(rms(i).ge.1000.0) then
	    rms(i)=999.999
	    write(*,51) i
	  end if
	  if(wrms(i).ge.1000.0) then
	    wrms(i)=999.999
	    write(*,51) i
	  end if
	  if(wrms(i).lt.0.001) then
	    wrms(i)=0.001
	    write(*,51) i
	  end if
	enddo
	
	ped_events=evnt_cnt-evnt_chr

	db_status=write_peds(run_name,date,db_file,avg,wrms,ped_events,
     +	                     av_val,mode,n_pmts)

***	write (2,18)run_name,evnt_cnt-evnt_chr,av_val,mode
***	write (2,19)(avg(i),i=1,n_pmts)
***	write (2,19)(wrms(i),i=1,n_pmts)
***	close (unit=2)
	write(6,20)fname,run_name,date,evnt_cnt-evnt_chr,evnt_chr,av_val
*       write(6,40)comms

***c-- Now delete the lock file
***        call system('rm -f /tmp/cpeds.lock')


*
* JB 941028 -- Turn off all tubes whose pedestal variances differ from
* the median by more than a fixed fraction
* 
        lowerlimit = av_val(2)*flowthresh
        upperlimit = av_val(2)*fupthresh
        write(*,'('' Average variance: '',f7.3)')
     &   av_val(2)
        write(*,'('' Lower limit     : '',f7.3)')lowerlimit
        write(*,'('' Upper limit     : '',f7.3)')upperlimit
        nnt = 0
        do i=1,m_pmts
           if((wrms(i).gt.upperlimit).or.
     >          (wrms(i).lt.lowerlimit)) then
              if(nnt.lt.MAXOFF) then
                 nnt = nnt + 1
                 ntubesoff(nnt) = i
                 write(*,3) i,wrms(i)
              else
                 write(*,3) i,wrms(i)
                 write(*,*) 'CP ** ',MAXOFF,' tube limit was exceeded'
              end if
           end if
        end do
	write(*,*) ' '
        if(nnt.ge.1) then  
*	   call ntoff(ftubelist,run_name,nnt,ntubesoff)
	   db_status=chk_toff(run_name,nnt,ftubelist)
	   if (db_status.eq.1) then
	      db_status=write_toff(run_name,nnt,ntubesoff,ftubelist)
	   end if
        end if
	call exit(0)

*Formats
 3      format(' Variance for PMT ',i4,':',f7.3,' exceeds threshold')
 11     format('Enter run name: ',$)
 12     format(q,a6)
 13     format(' Run name must be 6 characters: ',a)
 14     format(' File does not exist: ',a)
 18     format(1x,a6,i7,3f7.3,4x,a3)
 19     format(10f7.3)
* 20    format(1x,a6,': ',i6,x,i7,x,3f7.3,4x,a3)
 20     format('  CP **       File: ',a30,' Run: ',a6,' Date: ',i6,/,
     +       ' Events processed: ',i6,23x,'  Above threshold: ',i6,/,
     +       '         Avg peds: ',3f7.3)
 40     format('         Comments:',/,5(80a,/),a4)
 50     format(' WARNING: Peds for PMT ',i3,' bad.  Substituting',
     >       1x,'999.999')
 51     format(' WARNING: Ped var. for PMT ',i3,' bad.  ',
     >       'Substituting 999.99')
 810    format('  CP ** Configuration file ',a,' not found.')
 820    format('  CP ** Section "par" in ',a,' not found.')
 825    format('  CP ** Section "cut" in ',a,' not found.')
 830    format(  '  A list of standard config files may be found in ',
     &       /,'  /home/sc/config.  Copy one to ',a,' in',
     &       /,'  your current working directory.',
     &       /,'Exiting...')
 840    format('  CP ** Unable to read enough lines from ',a)
 850    format('  CP ** Unable to read ',a)
 860    format('  CP ** Using ',a30,' for database path ')
 861    format('  CP ** Using ',a30,' for tubelist path')

	end
