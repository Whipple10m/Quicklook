***************************************************************************
*
      program fz2red
*
*     941017
*     JHB, ADK, JQ, MAC
*
*     10m version.
*     Program to convert GDF format data to old reduced data format.
*     This program has also grown to become a first pass data analysis
*     program, fixing up a number of time dependent problems in the
*     system, and writing a more compact, faster intermediate data 
*     format which is used in the remainder of the data analysis.
*
*     Note that all of the data integrity checks should be done (and
*     errors corrected) in this first pass analysis.  Thereafter, 
*     diagnostics such as the full ccd records, redundant times,
*     full tracking records, etc. are no longer needed and a slightly
*     expanded version of the old data format is more than adequate.
*
*     WARNING:  Due to binary number representation problems, the 
*     reduced files which are produced by fz2red2 are not transportable
*     across platforms (e.g. DEC to SUN).  You must rerun a local
*     version of fz2red on the zebra file.  This isn't much more painful
*     than recompiling GDF on the other system.
*
*     Based on fz2red.f by ADK.
*     This data format became the default format 940800
*
*     Usage:            fz2red2 [Filename] [date] [option]
*     Explicit Inputs:  filename   data filename, e.g. gt000123.fz	
*     			date       UT date of this data run, e.g. 940321
*                       option     f - interactive mode to fix header
*
*     Revision History:
*
*     Modified by JQ (941220) to write out the cumulative livetime with
*     each event instead of the actual event time. Also the nominal
*     oscillator value is used in all time calculations since the 
*     previous version of this program gets confused when there are 
*     missing second markers.
*
*     Modified by JB 950612 to calculate useful diagnostics and to
*     create a header file for each run containing important information
*     like the mode, source name, run number, weather conditions, etc.
*
*     Modified by JB 950625 to calculate first sidereal time, and to
*     write livetime-duation to header for automatic duration matching
*     in gparamdat
*
*     Modified by JQ 961205 to handle variable numbers of channels.
*
*     Modified by MAC 970904 (at least started) to use new GDF v.076 
*     which uses Fortran90.  First step: replace "include gdf.inc" calls
*     with "use gdf" and change gdf_ev10.trigger (etc.) to gdf_ev10%trigger
*     and likewise for all of the other fields.  I had to copy gdf.mod to
*     this directory for "use gdf" statement to work.  Added gdf$option
*     call to tell GDF the type of machine where this is running.
*
*     Modified by MAC 970916 to keep track of livetime from livetime
*     scaler as well and use that if the gate_open and gate_close values are
*     zero.  Changed livetime calculation from gate_open and
*     gate_close to be able to use the gdf_ev10 variables mark_gps, 
*     gate_open and gate_close which are now logical*4 variables instead of
*     integer*4.
*
*     Modified by MAC 970921 to identify pedestal events by the 
*     gdf_ev10%trigger record.  They are no longer taken in the frame
*     or "environmental" records but are stored in the same field as real
*     events.
*
*     Modified by MAC 971124 to offset elapsed time by the time of the 
*     first event because the elapsed time oscillator starts about 1 minute
*     before the run starts.
*
*     981007 MAC: Significant revision of this program to add more diagnostic
*     checks, printouts and, most importantly, a histogram file 
*     (gtnnnn_d.hbook)which contains plots of some distributions that can be
*     used for diagnostic checks of the system.  The changes implemented are
*     outlined below.  
*     (1) gtnnnn_d.hbook contains plots of tracking records, comparisons of 
*     the GPS clock and the elapsed time oscillator, trigger bits for the
*     events, HV, current, and singles rate records.
*     (2) By comparing the GPS and oscillator times, I can derive an accurate
*     estimate of the elapsed oscillator frequency.  This estimate is printed
*     out in each run.  The true oscillator frequency appears to be around
*     9.999992 MHz, not the 10 MHz we had assumed.  The elapsed oscillator 
*     times are corrected for this error and the GPS time conversion in 
*     TRUETIME use this frequency now so their times are more accurate.
*     (3) File gtnnnn.errors is printed out with all errors that were thought
*     to need immediate attention.  Too many ones to list here.
*     (4) GRS clock status bits are checked now, as are the individual bits,
*     to make sure none are stuck.  Fixed a bug in the subroutine TRUETIME
*     which converts GRS clock bits to times.  I was not using the sub-second
*     data, so the gpsutc variable was only good to 1 second.  Did not affect
*     MJD.
*     (5) HV status bits are checked now and tubes which are turned off
*     are flagged in gtnnnn.inf.  Not actually tested yet, because the 
*     HV records have not been written to the data as far as I can tell.
*     Current monitor values are also checked.
*     (6) Singles rates are checked and highest and lowest 5 PMTs are 
*     printed out in gtnnnn.inf.
*
*     981102 MAC: Put in a fix to allow N2 runs to be processed from the
*     DACQ.  For some reason, some of the N2 files do not have an event=1
*     and because of a fix to an older problem (data from previous runs being 
*     in subsequent data files), fz2red will not analyze data until it sees
*     event=1.  Put in a date switch for this and it prints out an error if
*     it does not see event=1.
*
*     991227 MAC: Change elapsed scaler frequency from 9.999992 MHz to 
*     9.99999227 MHz.
*
*     000718 PM : Minor modifications to formats and spacings
***************************************************************************
*
* GDF definitions
*
*970904 MAC
      use gdf
*
      implicit none
* PAW stuff
      integer PAW_NWORDS
      parameter (PAW_NWORDS = 2*1024*256) ! 2 MBytes
      integer h(PAW_NWORDS)
      common/ PAWC /
     +     h
*
* Parameters
*
      real*8 two16,two32,utst,convd
      parameter (two16 = 65536.0d0)
      parameter (two32 = 4294967296.0d0)
      parameter (utst=366.2422d0/365.2422d0)
      parameter (convd = 57.29577951)
      
      integer*4 max_chan
      parameter (max_chan=1000)           ! JQ 961205

      integer*4 scunit
      parameter (scunit = 13)
      integer*4 hunit
      parameter (hunit = 10)
      integer*4 eunit
      parameter (eunit = 11)

      integer MAXMINS
      parameter (MAXMINS = 200)
      
*
* GDF definitions
*

*old      include '/usr/local/gdf/v066/gdf.inc'
*
* Header contents
*
      character crunid*4
      integer*4 inevents
      integer*4 nevts
      real*8 dduration
      integer*4 istdur
      character cmode*3
      integer*4 mode
      character csource*20
      integer*4 idate
      real*8 dmjd
      real*8 fmjd
      real*8 dfrjd
      real*4 ra_rad, dec_rad
      real*4 rra
      real*4 rdec
      real*4 in_ra,in_dec
      real*4 ra_offset,ra_off2
      real*4 rra_tele
      integer*4 iut
      integer*4 ist
      real*4 razimuth
      real*4 relevation
      character cskyq*6
      character ccomms*404
      character longcomms*8000
      integer*4 comlength,namelength
      integer*4 igpsbeg(7)
      logical runid_exists
      logical disgpsmark
      data disgpsmark / .false. /
      logical missing_secmark		
      data missing_secmark / .false. /
*
* Date stuff
*
      character*6 cdate
      integer*4 defdate,dumdate
      data dumdate / 940100 /
*
* Event contents
*
      integer*4 icode
      real*8 dtime
      integer*2 sshort(max_chan)        ! JQ 961205      
      integer*4 nadc
      integer*4 npmt
      integer*4 teleid                  ! telescope id
      data teleid / 1 /                 ! 1 for 10m, 2 for 11m
      integer*4 chi                     ! channel i 
*
* Tracking variables
*
      real*8 sum_el		! average elevation
      real*8 ave_el
      real*8 trdev, sumtrdev, sumtrdev2, ntrack
      real*8 meantrdev, sigtrdev
      data ntrack / 0.0 /
      data sum_el / 0.0 /
      data sumtrdev / 0.0 /
      data sumtrdev2 / 0.0 /
      real*8	hra,hh,mra,mm,ss
      real*8	sidhh,sidmm,sidss
      real*8 ddec,dd,mdec
*
* CCD variables
*
      integer*4 nstar
      integer*2 xstar(3),ystar(3),istar(3)
      logical	firstccd
      data firstccd / .true. /
*
* HV variables
*
      logical firsthv
      data firsthv / .true. /
*
* Event rate variables
*
      real*8 sumerate, meanerate, sum2erate, sigerate, nmins
      data sumerate / 0.0 /
      data sum2erate / 0.0 /
      data nmins / 0.0 /
*
* Timing stuff
*
*
* NOTE: Oscillator frequency should be 20MHz (as assumed by gpr10)
* but is actually 10MHz
*
      real*8 nomoscfreq		! Nominal oscillator frequency
      data nomoscfreq / 9.99987462d6 /
      real*8 caloscfreq		! Calibrated osc. frequency
      data caloscfreq / 9.99987462d6 /
      integer*4 jmod
      real*8 oscfreq(MAXMINS+1)
      real*8 dmintime(MAXMINS+1)
      integer*4 evtpmin(MAXMINS+1)
      integer*4 evtlastmin
      real*8 dminute,dlastminute,firstev,lastev,tlastmin
      real*8 dsecond,dlastsecond,tlastsec
      real*8 timtick,lasttick,elapsed_time
      real*8 time_sc,last_sc
      integer*4 elsec_old,elns_old,elsec_off,elns_off,tmwarns
      integer*4 elsec_tmp,elns_tmp
      integer*4 first_sec,first_ns
      integer*4 minute,lastmin,minmarks,minmarks2,nroll,nrollmin
      integer*4 second,lastsec,secmarks,secmarks2,nrollsec
      integer*4 abssec, first_abssec	! absolute second from the 
					! GPS time, not from the GPS
					! second mark
      integer*2  gps(3)   ! bcd encoded GPS time
      integer*4  grs(3)   ! bcd encoded GPS time from True Time clock
      real*8     mjd      ! UTC time [mjd]
      real*8     gpsutc	  ! UTC time from GPS clock [sec]
      real*8     lgpsutc  ! GPS time last event
*     real*8     lgate_close    ! Gate close time at last event
      real*8	 lastgpsminute	! GPS time at last min. mark [sec]
      real*8	 gpsdev	  ! Deviation in GPS time between minute marks
      integer	 year	  ! Gregorian year
      character  timestr*80	! returned time string [hh:mm:ss.sss]
      character	 tstr*10	! temp. time string [hh:mm:ss.s]
      character  fsidstr*20,futstr*20,fnsidstr*20,lastsidstr*20
      character  firsttstr*80
      real*8	 ttime		! temp. time value
      integer	 itime		! temp. integer time value
      real*8	 timefix	! time diff. between first two evts.
      real*8 timedev, sumtdev, sumtdev2, meantdev, sigtdev
      real*8 maxtdev
      real*8 ntdev
      data sumtdev / 0.0 /
      data sumtdev2 / 0.0 /
      data ntdev / 0.0 /
      data maxtdev / 0.0 /
      real	theta
      real      utdec
      real      alt                     ! Telescope elevation returned
                                        ! by derot() subroutine
      real      az                      ! Telescope azimuth returned by
                                        ! derot subroutine
      real      stime                   ! sidereal time of first event
      real 	nstime			! nom. sidereal start time from
					! VHEGRO
      real	laststime		! sidereal time of last event
      real	sid_duration		! sidereal duration of run (i.e.
					! diff between first and last
      logical	first2evt		! First event in pass #2
      logical	first2min		! First min mrk event in pass #2
      data first2evt / 0.0 /
      data first2min / 0.0 /
      logical	fosccal			! Flag that the osc. is cal'd
      data 	fosccal / .false. /
      logical 	justcald 
      data	justcald / .false. /
*
* Variables used in diagnostic tests.
*
      logical    lbit0,lbit1,lbit2  !Set true if trigger bits are set.
      logical    newpmtoff
      integer*4  bt,nbt,itrk,ieltm,ndtm,idtm
      integer*4  ipeds,event_old
      integer*4  nsings,n_hv
      integer*4  nel_tics
      integer*4  ntic_chks,itics
      integer*4  ndelmjd        !Num. of times delmjd.gt.0.1 usec.
      integer*4  pmtsoff(max_chan)  !Which tubes turned off (based on HV)
      integer*4  npmtsoff       !Num. of tubes turned off.
      
      real*4     s,sx,sy,sxx,sxy,dtslope !Var. for estimating
                                !the 10MHz oscillator frequency
      real*4     grs_bits(3,24) !ave. of bits set for the GRS clock
      real*4     trk_devs(40)   !tracking deviations in 0.01 deg steps
      real*4     el_v_tm(61)    !Ave. elev. v. elapsed time of run
      real*4     az_v_tm(61)    !Ave. az. v. elapsed time of run
      real*4     ra_v_tm(61)    !Ave. RA v. elapsed time of run
      real*4     dec_v_tm(61)   !Ave. Dec. v. elapsed time of run
      real*4     ntrks_tm(61)   !Num. of tracking records in each bin
      real*4     deltime(80)    !Dist. of diffs in elapsed time b/w GPS
                                !and the elapsed time scaler in msec.
      real*4     dtm_v_tm(61)   !Ave. time diff. v. elapsed time of run
      real*4     ntm_v_tm(61)   !Num. evts. v. elapsed time of run
      real*4     avedeltime     !Average time diff. in msec.
      real*4     cal_osc_freq   !Calibrated oscillator freq. - not used
                                !for anything right now except diagnostics.
      real*4     ntrig_bits(10) !Num. of times trig bits for 10m events are
                                !set. 1=peds, 2=pst, 3=trig, 4=peds+pst,
                                !5=peds+trig, 6=pst+trig, 7=peds+pst+trig,
                                !8=pst but not trig, 9=trig but not pst,
                                !10=nothing
      real*4     ave_hv(max_chan) !Ave. measured HV for PMT
      real*4     ave_ai(max_chan) !Ave. Anode current.
      real*4     dev_hv(max_chan) !Std. dev. measured HV for PMT
      real*4     dev_ai(max_chan) !Std. dev Anode current.
      real*4     osc_dev(201)   !Distribution of osc. freq-10 MHz
      real*4     ave_sing_rate(max_chan) !Ave. singles rate for each ch.
      real*4     dev_sing_rate(max_chan) !Std dev singles rate for each ch.
      real*4     max_sing_rate(5,2)        !5 highest singles rates (rate,PMT)
      real*4     min_sing_rate(5,2)        !5 lowest singles rates (rate,PMT)
      real*4     sum_osc
      real*8     med_osc        !Median calc. osc. freq.
      real*8     ave_osc,dev_osc   !Ave and std. dev. of calc. osc. freq.
      real*8     dtics,nave_tics,navesq_tics
      real*8     nel_sec
      real*8     firstutc,frutc_old
      real*8     delmjd         !usec Diff. b/w decoded GPS time and DACQ.
      real*8     utc_old,utc_min_old
      real*8     elosc_freq        !Best guess of true freq. of elapsed scaler
      real*8     defosc_freq       !Value used by DACQ.
      real*8     nel_sec_old
      data nel_sec_old / 0.d0 /
*old      data elosc_freq / 9.999992d6 /
      data elosc_freq / 9.99999227d6 /
      data defosc_freq / 1.d7 /
      integer*4  nel_tics_old
      data nel_tics_old / -1 /
*
* Variables needed for livetime (JQ 941220)
*
      real*8 gopen,gclose,lastgopen,lastgclose,livetime,rlivetime
      integer*4 nrollopen,nrollclose
* Added a few more (MAC 970910)
      real*8 livetm_sc                     !Livetime from livetime scaler
      real*8 live_time
      integer*4 livesec_old,livens_old     !Previous event's live sec and ns
      integer*4 livesec_off,livens_off     !If there is a reset of the scaler
      integer*4 livesec_tmp,livens_tmp
      integer*4 lvwarns
      integer*4 etwarns
      integer*4 ltwarns
      data etwarns / 0 /
      data ltwarns / 0 /

*
* Kludge variables for converting new logical*4 declarations of gdf_ev10
* variables mark_gps,mark_open,mark_close,gate_open,gate_close
*
      logical*4 lmark_gps,lmark_open,lmark_close,lgate_open,lgate_close
      logical*4 trig_bits
      integer*4 itrig_bits
      integer*4 imark_gps,imark_open,imark_close,igate_open,igate_close
      equivalence (lmark_gps,imark_gps),(lmark_close,imark_close),
     >     (lmark_open,imark_open),(lgate_open,igate_open),
     >     (lgate_close,igate_close),(trig_bits,itrig_bits)

*
* Program variables
*
      integer u          ! I/O unit variable
      real*8 dabs	 ! Fortran library function
      integer*4 lnblnk   ! Fortran library function
      integer*4 i,j,ierr,ierr2,last,irun,iextra,k,l
      integer*2 zero(max_chan)                      ! JQ  961205
      character*64 filename
      character pathin*64, pathout*64
      logical exist,trackdone
      logical event1
      real*8 firstgpsutc
      real*8 oldfirstgps
      logical nominmarks
      integer*4	ncode8,ncode12,ncode6
      integer idbg
      data idbg / 0 /
      data zero / max_chan*0 /
      character options*20
      integer   noptions
      character optchar
      logical uselvgate  !True if use livetime from gate_open/gate_close
      logical usetmtick  !True if using elapsed time from gate_close
      logical   ffixhdr	! Flags that fz2red is being run in fix-header
                        ! mode
      logical first_lv

      data ffixhdr / .false. /
      data uselvgate /.true./
      data usetmtick /.true./
      data first_lv /.true./

      write(*,*)
     & '***********************  FZ2RED  **************************'

*
* Get filename and date from command line  Date added 940425 ADK
*
      call getarg(1,filename)
*
* 941017 JB -- changed from intarg to getarg due to bug with intarg
*
      call getarg(2,cdate)
      read(cdate,'(i6)') defdate
*     write(*,*)'defdate: ',defdate
*
* Prompt for file if not present
*
      if (filename .eq. ' ') then
         write(6,130)
         read(5,'(a)')filename
      end if
*
* Get options string
*
      call getarg(3,options)
*
* Parse "options" command line parameter
*
      noptions = lnblnk(options)
      do i=1,noptions
        optchar = options(i:i)
        if((optchar.eq.'f').or.(optchar.eq.'F')) then
          ffixhdr = .true.
        end if
      end do
*
* If no date argument, substitute a dummy value  ADK 940425
*
      if (defdate.eq.0) defdate = dumdate

*
* This one line for Y2K compliance. What's all the fuss about ? SJF 990831
*
      if (defdate.lt.800000) defdate = defdate + 1000000

*
* Check file's existence
*
*      filename = pathin(1:lnblnk(pathin))//filename
      inquire(file=filename,exist=exist)
      if (.not.exist) then
         write(6,140)filename(1:lnblnk(filename))
         call exit(1)
      end if
*
* Initialize header data
* Take runid from filename, in case there is no header information.
*                                                   ADK  940413
*
        
      write(*,*)' Filename: ',filename
      j = lnblnk(filename)
      crunid = filename(j-6:j-3)
      inevents = 0
      dduration = 0.
      istdur = 0
      cmode = '   '
      csource = '???       '
*
* Modified by JB 940426 to handle the case that no run
* header information is detected, in order that the date
* is still specified by the command line date
*   
      idate = defdate

* After this date, we no longer have the good old livetime oscillator.
* The elapsed oscillator frequency seems to have changed between 1997
* and 1998 (not surprising I suppose).  These are the current best guesses
* for both of those.
      if (idate.gt.970900) then
         uselvgate=.false.
         usetmtick=.false.
         disgpsmark=.true.
         write(*,158) elosc_freq
      endif


* Following lines added by JQ 961205:
      call num_chan(teleid,idate,nadc,npmt)

      dmjd = 0.
      dfrjd = 0.
      rra = 0.
      rdec = 0.
      iut = 0
      ist = 0
      razimuth = 0.
      relevation = 0.
      cskyq = '      '
      ccomms = '      '
      longcomms = '     '
      do i = 1,7
         igpsbeg(i) = 0
      end do
*
* Initialise livetime variables (JQ 941220)
*
      gopen=0.0d6
      gclose=0.0d6
      lastgopen=0.0d6
      lastgclose=0.0d6
      livetime=0.0d6
      nrollopen=0
      nrollclose=0
*
* Initialize Zebra
*
      call gdf$init('Z',ierr)
*
* All GDF routines return ierr=0 if successful.  Check this condition.
*
      call chkerr(ierr,'Zebra initialization')
* MAC 970905 - Tell GDF the system type so that it knows how to read
* the data properly.
      call gdf$option('AXP',.TRUE.,ierr)
      call chkerr(ierr,'System type initialization')
*
* Open data file and scratch file
*
      call gdf$open(10,filename,'R',ierr)
      call chkerr(ierr,'Opening input file')
      open(scunit,name='temp',form='unformatted',status='unknown')
*     open(scunit,form='unformatted',status='scratch')
*
* ** Now build up the header information.  Since we are forcing GDF data
* ** into the old format, some of this data is unavailable.  The other
* ** pieces may be in the run status information section, the tracking
* ** system information, the 10m frame information, or in the event 
* ** data itself. 
*
* ** The data structure (gdf_run) for run status information contains
* ** the VAX date and time (although this seems not to be working), the
* ** run number and type, the sky quality, the observers, and comments.
* ** In run 230 (file gd000230.fz), I have noted that there are three
* ** occurrences of gdf_run, and two are blank.  For our purposes here,
* ** we will take the last one which is not blank.
*
* ** The data structure (gdf_track(2)) for 10m tracking information
* ** contains source RA and DEC (in both J2000 and contemporary
* ** coordinates), the telescope azimuth and elevation, and the local
* ** sidereal time.  Take the first of these (if any) which is non-zero,
* ** since 10m azimuth and elevation information are correlated with
* ** the first sidereal minute.  The event data has been seen to have
* ** left-over events at the beginning, and such problems could plague
* ** our apprach on which gdf_track record to use, but we adopt it
* ** for now anyway.  Update this information until after event number
* ** 1 is reached.
*
* ** The data structure (gdf_fr10) for the 10m frame contains
* ** the time of the last UT second marker (from the 10m GPS clock)
* ** but Joachim comments that this might not be correct.  The 
* ** GPS second mark is also stored in gdf_ev10.
* ** The ADC values of 2 injected pedestal events are also stored in
* ** a 10m frame (gdf_fr10)
*
* ** The data structure (gdf_ev10) for the 10m event contains the event
* ** ADC values. It also contanins the GPS clock time.
* ** It also contains the latched oscillator value (gate_close), the 
* ** oscillator value of the last GPS second marker from the 10m, the 
* ** event number, and the live time counter. I have seen some events 
* ** (at the beginning of gd000230.fz, for example) with high counter 
* ** values before the counter is reset to 1. Skip all events until the
* ** first event is seen.
*
* ** Now implement this mess!
*
*
*
* Take one pass through data file, writing data to a scratch file.
* Second pass will be through scratch file. On the first pass, count
* the number of events and compile the header information.  Also
* compile a list of effective oscillator frequencies
* minute-by-minute
*
      do i = 1,MAXMINS
         oscfreq(i) = nomoscfreq   !Nominal oscillator value
      end do
      minmarks = 0
      secmarks = 0
      nroll = 0
      nrollmin = 0
      nrollsec = 0
      iextra = 0
      event1 = .false.
      trackdone = .false.
      runid_exists = .false.
*
* Fix the header by prompting for all header information if the
* f-option was specified
*
      if(ffixhdr) then
          write(*,*) '  FZ ** Please provide the missing information:'
          write(*,'('' RA now (hhmmss.s)    :'',$)')
          read(5,*)in_ra
          hh = aint(in_ra/10000.0)
          mm = aint((in_ra - hh*10000.0)/100.0)
          ss = (in_ra - hh*10000.0 - mm*100.0)
          rra = (hh+mm/60.0+ss/3600.0)*0.2617993878
          write(6,'(''  -> RA (rad.):'',f10.6)')rra
          write(*,'('' DEC now (+/-ddmmss.s):'',$)')
          read (5,*)in_dec
          dd = aint(in_dec/10000.0)
          mm = aint((in_dec - dd*10000.0)/100.0)
          ss = (in_dec - dd*10000.0 - mm*100.0)
          rdec = (dd+mm/60.0+ss/3600.0)/convd  
          write(6,'(''  -> DEC (rad.):'',f10.6)')rdec
          ra_rad = rra
          dec_rad = rdec
*
* Convert back using same code for rad to hhmmss etc as a check
*
          hra = rra*convd*24.0/360.0
          hh = aint(hra)
          mra = 60.0*(hra-hh)
          mm = aint(mra)
          ss = 60.0*(mra-mm)
          rra = hh*10000.+mm*100.+ss

          ddec = rdec*convd
          dd = aint(ddec)
          mdec = 60.0*(ddec-dd)
          mm = aint(mdec)
          ss = 60.0*(mdec-mm)
          rdec = dd*10000.+mm*100.+ss

          write(*,'('' Source name          :'',$)')
          read(5,'(a12)')csource(1:12)
          write(*,'('' Mode (0,1,2)         :'',$)')
          read(5,*)mode
          write(6,*) '  FZ ** Tracking information is now complete.'
*         trackdone = .true.
      end if
*
* First pass, main event loop.  Extract comments, run header info
* and calibrate the oscillator time.
*
      do while (ierr.eq.0)
         call gdf$read(10,' ',ierr)
         if ((gdf_run%new).and.(gdf_run%valid)) then
             if (gdf_run%run .ne. 0) then
               runid_exists = .true.
               idate = gdf_run%idate
*              year = gdf_run%year
*               write(*,*)'idate: ',idate   
               if (idate .eq. 0) idate = defdate
               year = 1900 + idate/10000
               irun = mod(gdf_run%run,10000)
* JQ 980804               write(crunid,'(i4.4)')irun
               if (gdf_run%sky_quality .eq. 1) then
                  cskyq = 'A'
               else if (gdf_run%sky_quality .eq. 2) then
                  cskyq = 'B'
               else if (gdf_run%sky_quality .eq. 3) then
                  cskyq = 'C'
               end if
*MAC970905               istdur = gdf_run%sid_length
               istdur = nint(gdf_run%sid_length(1))
               comlength = lnblnk(gdf_run%comment)
               longcomms = gdf_run%comment(1:comlength)
               ccomms = gdf_run%comment(1:400)
             end if ! gdf_run%run.ne.0
         end if ! gdf_run%new
*
* Get a sample of ccd data and display it
* 
         if (event1 .and. 
     &     (gdf_ccd(1)%new) .and. (gdf_ccd(1)%valid)) then
           if (firstccd) then
             do i=1,3
               nstar = gdf_ccd(1)%nstar
               if(i.lt.nstar) then
                 xstar(i) = gdf_ccd(1)%star(1,i)
                 ystar(i) = gdf_ccd(1)%star(2,i)
* MAC 971014
                 istar(i) = gdf_ccd(1)%intensity(i)
*MAC                 istar(i) = gdf_ccd(1)%star(3,i)
               endif
             end do
             firstccd = .false.
           endif
         endif  !ccd new
*
* Obtain HV records and store the values.  There should only be one
* record unless there is a change in HV status (e.g., tube turned off
* mid-run).  If there is a change, affected tubes will be flagged and
* printed out.
*
* DO NOT UNCOMMENT UNTIL THIS IS TESTED.
*m         type*,'1'
*m         type*,' '
         if (gdf_hv(1)%new.and.gdf_hv(1)%valid) then
            n_hv=n_hv+1
*debug            type*,'n_hv=',n_hv
*debug            do j=1,33
*debug               write(*,199) (nint(gdf_hv(1)%v_actual(k)), 
*debug     >              k=10*(j-1)+1,10*j)
*debug            enddo
            do j=1,npmt
               ave_hv(j)=ave_hv(j)+gdf_hv(1)%v_actual(j)
               dev_hv(j)=dev_hv(j)+(gdf_hv(1)%v_actual(j))**2
               ave_ai(j)=ave_ai(j)+gdf_hv(1)%i_anode(j)
               dev_ai(j)=dev_ai(j)+(gdf_hv(1)%i_anode(j))**2
* Loop through the status bits for the channel.
* bit 0
               if (iand(gdf_hv(1)%status(j),1).eq.1) then
                  if (abs(gdf_hv(1)%v_actual(j)-
     >                 gdf_hv(1)%v_set(j)).gt.20.) 
     >                 then
                     write(*,190) j,gdf_hv(1)%v_actual(j),
     >                    gdf_hv(1)%v_set(j)
                     write(eunit,190) j,gdf_hv(1)%v_actual(j),
     >                    gdf_hv(1)%v_set(j)
                  endif
               else
                  newpmtoff=.true.
                  do k=1,npmtsoff
                     if (pmtsoff(k).eq.j) newpmtoff=.false.
                  enddo
                  if (newpmtoff) then
                     npmtsoff=npmtsoff+1
                     pmtsoff(npmtsoff)=j
                  endif
                  if (abs(gdf_hv(1)%v_actual(j)
     >                 -gdf_hv(1)%v_set(j)).lt.20.) then
                     write(*,194) j,gdf_hv(1)%v_actual(j)
                     write(eunit,194) j,gdf_hv(1)%v_actual(j)
                  endif
               endif
*test* bit 1
*test               if (iand(ishft(gdf_hv(1)%status(j),-1),1).eq.1) then
*test                  write(*,192) j,
*test     >                 'HV ramping up, gain may be changing'
*test                  write(eunit,192) j,
*test     >                 'HV ramping up, gain may be changing'
*test               endif
*test* bit 2
*test               if (iand(ishft(gdf_hv(1)%status(j),-2),1).eq.1) then
*test                  write(*,192) j,
*test     >                 'HV ramping down, gain may be changing'
*test                  write(eunit,192) j,
*test     >                 'HV ramping down, gain may be changing'
*test               endif
*test* bit 3
*test               if (iand(ishft(gdf_hv(1)%status(j),-3),1).eq.1) then
*test                  if (iand(gdf_hv(1)%status(j),1).ne.1) then
*test                     write(*,198) j,'Hardware status is OFF',
*test     >                    ' but the software status is ON'
*test                     write(eunit,198) j,'Hardware status is OFF',
*test     >                    ' but the software status is ON'
*test                  endif
*test               else
*test                  if (iand(gdf_hv(1)%status(j),1).eq.1) then
*test                     write(*,198) j,'Hardware status is ON',
*test     >                    ' but the software status is OFF'
*test                     write(eunit,198) j,'Hardware status is ON',
*test     >                    ' but the software status is OFF'
*test                  endif
*test               endif
*test* bit 4
*test               if (iand(ishft(gdf_hv(1)%status(j),-4),1).eq.1) then
*test                  write(*,192) j,
*test     >                 'Violation of supply limit'
*test                  write(eunit,192) j,
*test     >                 'Violation of supply limit'
*test               endif
*test* bit 5
*test               if (iand(ishft(gdf_hv(1)%status(j),-5),1).eq.1) then
*test                  write(*,192) j,
*test     >                 'Violation of current limit'
*test                  write(eunit,192) j,
*test     >                 'Violation of current limit'
*test               endif
*test* bit 6
*test               if (iand(ishft(gdf_hv(1)%status(j),-6),1).eq.1) then
*test                  write(*,192) j,'Voltage error'
*test                  write(eunit,192) j,'Voltage error'
*test               endif
*test* bit 7
*test               if (iand(ishft(gdf_hv(1)%status(j),-7),1).eq.1) then
*test                  write(*,192) j,'Violation of voltage limit'
*test                  write(eunit,192) j,'Violation of voltage limit'
*test               endif
*test* bit 9
*test               if (iand(ishft(gdf_hv(1)%status(j),-9),1).ne.1) then
*test                  write(*,192) j,'Crate HV is off'
*test                  write(eunit,192) j,'Crate HV is off'
*test               endif
*test* bit 10
*test               if (iand(ishft(gdf_hv(1)%status(j),-10),1).ne.1) then
*test                  write(*,192) j,'Eeprom status is bad'
*test                  write(eunit,192) j,'Eeprom status is bad'
*test               endif
*test* bit 11
*test               if (iand(ishft(gdf_hv(1)%status(j),-11),1).ne.1) then
*test                  write(*,192) j,'Battery status is bad'
*test                  write(eunit,192) j,'Battery status is bad'
*test               endif
*test* bit 13
*test               if (iand(ishft(gdf_hv(1)%status(j),-13),1).ne.1) then
*test                  write(*,192) j,'24V status is bad'
*test                  write(eunit,192) j,'24V status is bad'
*test               endif
*test* bit 14
*test               if (iand(ishft(gdf_hv(1)%status(j),-14),1).eq.1) then
*test                  write(*,192) j,'Crate panic condition is active'
*test                  write(eunit,192) j,'Crate panic condition is active'
*test               endif
            enddo               !End of loop through PMTs
         endif                  !New gdf HV record.
*m         type*,'2'
*m         type*,' '
*
* Get first track event to obtain the initial ra, dec, az and el
* Also, accumulate tracking deviations for subsequent tracking 
* records.  Only do this, however, if the (true) first event has
* been detected.  This check is necessary due to early problems with
* the writer(?).
*
         if (event1 .and. 
     &        (gdf_track(1)%new) .and. (gdf_track(1)%valid)) then
           if (.not.trackdone)  then
             rra = gdf_track(1)%rasc_today
             rdec = gdf_track(1)%decl_today
					! Convert from rad to min RA
             ra_offset = gdf_track(1)%rasc_offset * 229.1831181
             rra_tele = gdf_track(1)%rasc_tele
             ra_off2 = (rra_tele-rra) * 229.1831181
             ra_rad = rra
             dec_rad = rdec
*
* Transform coordinates from radians to hhmmss.s, ddmmss.s
* for compatibility with old data format.
*
             hra = rra*convd*24.0/360.0
             hh = aint(hra)
             mra = 60.0*(hra-hh)
             mm = aint(mra)
             ss = 60.0*(mra-mm)
             rra = hh*10000.+mm*100.+ss

             ddec = rdec*convd
             dd = aint(ddec)
             mdec = 60.0*(ddec-dd)
             mm = aint(mdec)
             ss = 60.0*(mdec-mm)
             rdec = dd*10000.+mm*100.+ss

             razimuth = gdf_track(1)%azimuth
             relevation = gdf_track(1)%elevation
	     namelength = lnblnk(gdf_track(1)%source)
	     if (namelength.gt.20) then
               namelength = 20
             end if
             csource = gdf_track(1)%source(1:namelength)
             mode = gdf_track(1)%mode
*            if (event1)  then
               trackdone = .true.
*            end if
           end if ! .not. trackdone
           trdev = gdf_track(1)%deviation * convd
           sumtrdev = sumtrdev + trdev
	   sumtrdev2 = sumtrdev2 + (trdev*trdev)
           sum_el = sum_el + gdf_track(1)%elevation
           ntrack = ntrack + 1.0
* Additional diagnostics for tracking info.
           itrk = max(1,min(40,int(trdev/0.005)+1))
           trk_devs(itrk)=trk_devs(itrk)+1.
           ieltm=max(1,min(61,
     >          int(sngl(gdf_track(1)%utc-gdf_run%utc_start)*
     >          1440./0.5)+2))
           ntrks_tm(ieltm)=ntrks_tm(ieltm)+1.
           el_v_tm(ieltm)=el_v_tm(ieltm)+sngl(gdf_track(1)%elevation
     >          *convd)
           az_v_tm(ieltm)=az_v_tm(ieltm)+sngl(gdf_track(1)%azimuth
     >          *convd)
           ra_v_tm(ieltm)=ra_v_tm(ieltm)
     >          +sngl(gdf_track(1)%rasc_tele*convd)
           dec_v_tm(ieltm)=dec_v_tm(ieltm)
     >          +sngl(gdf_track(1)%decl_tele*convd)
         end if ! event1.and.gdf_track(1)%new...
*m         type*,'3'
*m         type*,' '
*
* In the first pass through the events we are mainly concerned with
* counting the events and with keeping track of the GPS second
* markers to calculate the oscillator frequency for each second.
* We also keep track of the oscillator values for each second.
* Here we have one variable (dlastsecond) for the last second marker,
* and an array (dsectime) containing the oscillator values of each
* second marker.  With some care, the dlastsecond variable could be
* eliminated.  In the interests of simplicity and getting the program
* running quickly, it has been retained.
*
         if ((gdf_ev10%new) .and. (gdf_ev10%valid)) then
*
* MAC 970911 here I set the logical*4 variables (some strange attempt at
* the equivalent of C's unsigned integers) in the GDF gdf_ev10
* record to some intermediate logical*4 variables which are equivalenced
* to integer*4 variables that we can actually use in our analysis.
*
            lmark_gps=gdf_ev10%mark_gps
            lmark_open=gdf_ev10%mark_open
            lmark_close=gdf_ev10%mark_close
            lgate_open=gdf_ev10%gate_open
            lgate_close=gdf_ev10%gate_close
            trig_bits=gdf_ev10%trigger
*
* Keep track of the number of time the various combinations of trigger bits 
* appear.
*
            if (defdate.gt.970901) then
               lbit0=.false.
               lbit1=.false.
               lbit2=.false.
               if (iand(ishft(itrig_bits,-gdf_trig_ped),1).eq.1) 
     >              lbit0=.true.
               if (iand(ishft(itrig_bits,-gdf_trig_pst),1).eq.1) 
     >              lbit1=.true.
               if (iand(ishft(itrig_bits,-gdf_trig_mul),1).eq.1) 
     >              lbit2=.true.
               if (lbit0) ntrig_bits(1)=ntrig_bits(1)+1.
               if (lbit1) ntrig_bits(2)=ntrig_bits(2)+1.
               if (lbit2) ntrig_bits(3)=ntrig_bits(3)+1.
               if (lbit0.and.lbit1) 
     >              ntrig_bits(4)=ntrig_bits(4)+1.
               if (lbit0.and.lbit2) 
     >              ntrig_bits(5)=ntrig_bits(5)+1.
               if (lbit1.and.lbit2) 
     >              ntrig_bits(6)=ntrig_bits(6)+1.
               if (lbit0.and.lbit1.and.lbit2) 
     >              ntrig_bits(7)=ntrig_bits(7)+1.
               if (lbit1.and..not.lbit2) 
     >              ntrig_bits(8)=ntrig_bits(8)+1.
               if (lbit2.and..not.lbit1) 
     >              ntrig_bits(9)=ntrig_bits(9)+1.
               if (.not.lbit0.and..not.lbit1.and..not.lbit2)
     >              ntrig_bits(10)=ntrig_bits(10)+1.
            endif
*
* event1 is initially false.  Skip events until event number is 1 or
* the first event has been detected.  This
* eliminates left over events at the beginning of the run.  After the first
* valid event is encountered set event1 = .true. and start processing all
* events.
*
            if (abs(defdate-940900).lt.30000) then
               if ((.not. event1) .and. ((gdf_ev10%event).ne.1)) 
     >              go to 998
            else
               if ((.not. event1) .and. ((gdf_ev10%event).lt.1)) 
     >              go to 998
            endif
            if ((.not.event1).and.(gdf_ev10%event.gt.1)) then
               write(*,141) gdf_ev10%event
               write(eunit,141) gdf_ev10%event
            endif
            event1 = .true.
* MAC - as of 970901 the pedestals are also stored as regular events,
* but identified by the field gdf_ev10%trigger.  So, for dates after
* 970901 add a check of this field before proceeding with converting
* this event.
*debug            write(*,'(''Trigger bits values = '',I4.4)') 
*debug     >           gdf_ev10%trigger
            if ((defdate.gt.970901).and.
     >           (iand(ishft(itrig_bits,-gdf_trig_ped),1)
     >           .eq.1)) goto 997
*

            inevents = inevents + 1
*
* MAC 970916 - Wisconsin clock time is now the standard
*
*old            if ((gdf_ev10%grs_clock(1).gt.0).or.
*old     >           (gdf_ev10%grs_clock(2).gt.0)) then
            if (defdate.gt.970901) then

               grs(1) = gdf_ev10%grs_clock(1)
               grs(2) = gdf_ev10%grs_clock(2)
               grs(3) = gdf_ev10%grs_clock(3)
               call truetime(grs,year,mjd,gpsutc,timestr)
* MAC test
               nbt=nbt+1
               do k=1,3
                  do j=1,24
                     bt=iand(ishft(grs(k),1-j),1)
                     grs_bits(k,j)=grs_bits(k,j)+float(bt)
                  enddo
               enddo
            else
               gps(1) = gdf_ev10%gps_clock(1)
               gps(2) = gdf_ev10%gps_clock(2)
               gps(3) = gdf_ev10%gps_clock(3)
               call gpstime(gps,year,mjd,gpsutc,timestr)
            endif
**
*	    if (inevents. lt. 15) then
*               write(*,'('' '',i3,'' gpsutc :'',f10.3,'' osctime :'',
*    &          f10.3)')inevents,gpsutc,gdf_ev10%gate_close/nomoscfreq
*           endif
**
*
* Initialize time stuff for the first event
*
            if (inevents .eq. 1) then
* UTC time of first event 
               firstutc=gdf_ev10%utc
               lastsec = imark_gps
*MAC970911               lastsec = gdf_ev10%mark_gps
               lasttick = igate_close
*MAC970911               lasttick = gdf_ev10%gate_close
               firstev = dfloat(igate_close)
*MAC970911               firstev = dfloat(gdf_ev10%gate_close)
               oldfirstgps = gpsutc	! Old way of calculate first
					! evt gps - should agree 
               if (firstev .lt. 0.) firstev = firstev + two32
 				! Compensate for event counter rolling 
				! over?
               first_abssec = gdf_ev10%utc

*MAC 970911 
               
               lastgopen=dfloat(igate_open)      
               lastgclose=dfloat(igate_close)    
*old               lastgopen=dfloat(gdf_ev10%gate_open)      ! JQ 941220
*old               lastgclose=dfloat(gdf_ev10%gate_close)    ! JQ 941220
*              lgate_close = gdf_ev10%gate_close
*
* Mod by JB 950628 - correct for possibility of lastgopen<0,
* lastgclose<0.  I believe this possibility was overlooked, although
* JQ might have had his reasons.
*
	       if(lastgopen.lt.0.)lastgopen = lastgopen + two32
               if(lastgclose.lt.0.)lastgclose = lastgclose + two32
*
* Mod by MAC 971124 for non-zero starting point for the livetime.
               first_sec=gdf_ev10%elapsed_sec
               first_ns=gdf_ev10%elapsed_ns
* Check that first livetime is not very large - if it is, put in offset.
               if ((gdf_ev10%live_sec
     >              +nint(dfloat(gdf_ev10%live_ns)/1.d9))
     >              .gt.1) then
* Kludge for 2^26 step jumps in oscillators during this period.
                  if ((defdate.gt.980429).and.
     >                 (defdate.lt.980523)) then
                     livesec_tmp=0
                     livens_tmp=0
                     do while ((gdf_ev10%live_sec+livesec_tmp
     >                    +nint(dfloat(
     >                    gdf_ev10%live_ns+livens_tmp)/1.d9))
     >                    .gt.1)
                        livesec_tmp=livesec_tmp-6
                        livens_tmp=livens_tmp-710886400
                        if (livens_tmp.lt.0) then
                           livesec_tmp=livesec_tmp-1
                           livens_tmp=livens_tmp+1000000000
                        endif
                     enddo
* Not a 2^26 screw up, it's something else, just set to zero.
                     if (abs(gdf_ev10%live_sec+livesec_tmp+
     >                    nint(dfloat(gdf_ev10%live_ns+
     >                    livens_tmp)/1.d9)).gt.1) then
                        livesec_tmp=0
                        livens_tmp=0
                        livesec_off=-gdf_ev10%live_sec
                        livens_off=-gdf_ev10%live_ns
                        type*,' FZ ******* WARNING *******'
                        type*, 'Live time >1s for 1st event 
     >                     - setting to zero'
*                       type*,' Setting first event livetime to zero'
                        type*,' **************************'
                     endif
                  else
                     type*,' FZ ******* WARNING *******'
                     type*,' Live time >1s for 1st event 
     >                     - setting to zero'
*                    type*,'Setting first event live time to zero'
                     type*,' **************************'
                     livesec_off=-gdf_ev10%live_sec
                     livens_off=-gdf_ev10%live_ns
                  endif  !End of check that date is not gt 980429
               endif !End of check that live time is not large for 1st event
            end if
***
* Record of live time from the live time scaler (MAC 970910)
* livesec_off and livens_off are non-zero only if there has been a reset
* or rollover of the live time scaler.
            livetm_sc=dfloat(gdf_ev10%live_sec+livesec_off)
     >           +dfloat(gdf_ev10%live_ns+livens_off)/1.d9
*debug            if (mod(inevents,20).eq.0) then
*debug               type*,'live_time,live_sec,live_ns',
*debug     >              livetm_sc,gdf_ev10%live_sec,gdf_ev10%live_ns
*debug            endif
            if (((gdf_ev10%live_sec+gdf_ev10%live_ns.lt.0).or.
     >          (livetm_sc.gt.3.d4)).and.first_lv) then
               type*,'  FZ ******* WARNING *******'
               type*,'  Live time value is screwy - watch your analysis'
               type*,'  Live time=',livetm_sc
               type*,'  Live sec,ns',gdf_ev10%live_sec,
     >              gdf_ev10%live_ns
               type*,'  FZ ***********************'
               first_lv=.false.
            endif
* Livetime stuff below (JQ 941220)
*
            gopen=dfloat(igate_open)
*MAC970911            gopen=dfloat(gdf_ev10%gate_open)
            if (gopen.lt.0) gopen=gopen+two32
            gopen=gopen+nrollopen*two32
            if (gopen.lt.lastgopen) then
               gopen=gopen+two32
               nrollopen=nrollopen+1
            endif
          
            gclose=dfloat(igate_close)
*MAC970911            gclose=dfloat(gdf_ev10%gate_close)
            if (gclose.lt.0) gclose=gclose+two32
            gclose=gclose+nrollclose*two32
            if (gclose.lt.lastgclose) then
               gclose=gclose+two32
               nrollclose=nrollclose+1
            endif


*MAC 970910 Add check that live time scaler values have not rolled 
* over.
            if(inevents.gt.1) then
               livetime=livetime+gclose-gopen 
*
* Kludge for screwy livetime oscillator values (they get shifted 
* by multiples of 2^26 oscillations.)
*
               livesec_tmp=0
               livens_tmp=0
               
               if ((defdate.gt.980429).and.(defdate.lt.980523)) then
* The live time went backward.
                  if (livetm_sc.lt.(dfloat(livesec_old)+
     >                 dfloat(livens_old)/1.d9)) then
* Try to correct it by adding multiples of 2^26 cycles/10^7 Hz until the 
* time is greater than the live time of the last event.
                     do while ((livetm_sc+dfloat(livesec_tmp)+
     >                    dfloat(livens_tmp)/1.d9)
     >                    .lt.(dfloat(livesec_old)+
     >                    dfloat(livens_old)/1.d9))
                        livesec_tmp=livesec_tmp+6
                        livens_tmp=livens_tmp+710886400
                        if (livens_tmp.gt.1000000000) then
                           livesec_tmp=livesec_tmp+1
                           livens_tmp=livens_tmp-1000000000
                        endif
                     enddo
* If the live time difference after correction is >0.5 sec, assume this
* was a real rollover and add to the offset values.
                     if ((livetm_sc-
     >                    dfloat(livesec_old-livesec_tmp)-
     >                    dfloat(livens_old-livens_tmp)/1.d9)
     >                    .gt.1.0d0) then
                        livesec_tmp=0             !Reset these if not used
                        livens_tmp=0
                        livesec_off=livesec_old
                        livens_off=livens_old
* Correction for funny case where a rollover occurred but there is 
* significant elapsed time despite this being the first event.
                        if (abs(gdf_ev10%live_sec).gt.0) then
                           if ((dfloat(gdf_ev10%live_sec)+
     >                          dfloat(gdf_ev10%live_ns)/1.d9)
     >                          .gt.0.d0) then
                              do while ((dfloat(gdf_ev10%live_sec+
     >                          livesec_tmp)+dfloat(
     >                          gdf_ev10%live_ns+livens_tmp)/1.d9)
     >                          .gt.6.d0)
                                 livesec_tmp=livesec_tmp-6
                                 livens_tmp=elns_tmp-710886400
                                 if (livens_tmp.lt.0) then
                                    livesec_tmp=livesec_tmp-1
                                    livens_tmp=livens_tmp+1000000000
                                 endif
                              enddo
                           else
                              do while ((dfloat(gdf_ev10%live_sec+
     >                          livesec_tmp)+dfloat(
     >                          gdf_ev10%live_ns+livens_tmp)/1.d9)
     >                          .lt.0.d0)
                                 livesec_tmp=livesec_tmp+6
                                 livens_tmp=livens_tmp+710886400
                                 if (livens_tmp.gt.0) then
                                    livesec_tmp=livesec_tmp+1
                                    livens_tmp=livens_tmp-1000000000
                                 endif
                              enddo
                           endif
                        endif
                        livetm_sc=dfloat(gdf_ev10%live_sec+
     >                       livesec_old+livesec_tmp)
     >                       +dfloat(gdf_ev10%live_ns+livens_old
     >                       +livens_tmp)/1.d9
                     else
                        livetm_sc=livetm_sc+dfloat(livesec_tmp)
     >                       +dfloat(livens_tmp)/1.d9
                     endif
                  else
* The livetime jumped by too much.
                     if ((livetm_sc-(dfloat(livesec_old)
     >                    +dfloat(livens_old)/1.d9)).gt.6.d0) then
                        do while ((livetm_sc-
     >                       dfloat(livesec_old-livesec_tmp)-
     >                       dfloat(livens_old-livens_tmp)/1.d9)
     >                       .gt.6.7108864d0)
                           livesec_tmp=livesec_tmp-6
                           livens_tmp=livens_tmp-710886400
                           if (livens_tmp.lt.0) then
                              livesec_tmp=livesec_tmp-1
                              livens_tmp=livens_tmp+1000000000
                           endif
                        enddo
* Weird
                        if (abs(livetm_sc-
     >                       dfloat(livesec_old-livesec_tmp)-
     >                       dfloat(livens_old-livens_tmp)/1.d9)
     >                       .gt.1.0d0) then
                           ltwarns=ltwarns+1
                           if (ltwarns<=10) then
                             type*,' FZ ******* WARNING *******'
                             type*,' Livetime jumped forward (not by 
     >                       multiple of 2^26) - probably inaccurate.'
                             type*,' **************************'
*                            type*,' FZ ******* WARNING *******'
*                            type*,' Live time jumped forward but not'
*                            type*,' by multiple of 2^26 - livetime is'
*                            type*,' probably inaccurate'
*                            type*,' **************************'
                           endif
                        endif
                        livetm_sc=livetm_sc+dfloat(livesec_tmp)
     >                       +dfloat(livens_tmp)/1.d9
                     endif
                  endif   !End if for check if live time went backward
               else
                  if (livetm_sc.lt.(dfloat(livesec_old)
     >                 +dfloat(livens_old)/1.d9)) then
                     if (abs(livetm_sc-dfloat(livesec_old)
     >                    -dfloat(livens_old)/1.d9).lt.10.) then
                        type*,'  FZ ******** WARNING **************'
                        type*,'  Live time scaler went backwards!'
                        type*,'  Live time may be inaccurate!'
                        type*,'  FZ *******************************'
                     else
* If there is more than 1 sec on the rolled over live time, there
* is some additional time gap that has not been accounted for, so 
* as a sort of fix, don't include the seconds in the live time calculation.
                        type*,'  FZ ********** WARNING *************'
                        type*,'  Live time scaler value rolled over.'
                        type*,'  Was the CAMAC crate reset?'
                        if (gdf_ev10%live_sec.gt.0) then
                           type*,'gdf_ev10 live_sec=',
     >                          gdf_ev10%live_sec
                           type*,'Some additional time gap'
                           type*,'Subtracting from the livetime ',
     >                          gdf_ev10%live_sec,' s'
                        endif
                        livesec_off=livesec_old
                        livens_off=livens_old
                        livetm_sc=dfloat(livesec_old)
     >                       +dfloat(gdf_ev10%live_ns+livens_old)/1.d9
                     endif
                     type*,'  ***********************************'
                  else !End if for check that live time has not decreased
* Weird
                     if ((livetm_sc-dfloat(livesec_old)-
     >                    dfloat(livens_old)/1.d9).gt.1.0d0) then
                        type*,' FZ *********** WARNING **************'
                        type*,' Livetime jumped forward >1s - may be ',
     >                          'inaccurate.'
                        type*,' *************************************'

*                        type*,'***** WARNING *****'
*                        type*,'Livetime jumped forward by more than'
*                        type*,'1 sec - livetime may be inaccurate.'
*                        type*,'*******************'
                     endif
                  endif
               endif         !End of check that UTdate is not gt 980429

               if ((gdf_ev10%live_sec.eq.livesec_old).and.
     >              (gdf_ev10%live_ns.eq.livens_old).and.
     >              (gdf_ev10%live_sec+gdf_ev10%live_ns).gt.0) then
                  lvwarns=lvwarns+1
                  if (lvwarns.le.5) then
                    type*,'  FZ *************WARNING******************'
                    type*,'  Livetime not increasing from livetime'
                    type*,'  scaler.  Livetime might not be accurate.'
                    type*,'  *****************************************'
                  endif
                  if (lvwarns.eq.5) then
                    type*,' **********************************',
     >                    '*********************'
                    type*,'  Too many warnings of this type to', 
     >                    ' continue displaying!'
                    type*,' **********************************',
     >                    '*********************'         
                  endif
               endif
* MAC 970912
               if ((livetime.eq.0).or.(.not.uselvgate)) then
                  if (uselvgate) then
                     type*,'  FZ ************  WARNING  **************'
                     type*,'  Livetime from gate_open/gate_close is 0.'
                     type*,'  Using the livetime scaler.'
                     type*,'  ****************************************'
                     uselvgate=.false.
                  endif
                  live_time=livetm_sc
               else
                  live_time=livetime
               endif

            endif  !End if for check that event number is gt 1
* The off correction is for rollovers, the tmp correction is for the
* 2^26 jumps starting 980429.
            livesec_old=gdf_ev10%live_sec+livesec_off+livesec_tmp
            livens_old=gdf_ev10%live_ns+livens_off+livens_tmp
               
* (inevents>1) to prevent negative livetime for 1st event.
*
*           livetime = dfloat(gdf_ev10%live_sec)*(20.0d6/nomoscfreq)
*
***
            if(inevents .eq. 2) then
*
* Patch the first gps time by using the second gps time and the
* time difference from the oscillator.  I think that this is rigorous,
* despite being a bit of a mess.
*
*              timefix = (gdf_ev10%gate_close - lgate_close)/nomoscfreq
               timefix = (gclose-lastgclose)/nomoscfreq
               firstgpsutc = gpsutc - timefix	! fix first gps time
               fmjd = mjd - timefix/86400.0	! fix first mjd
	       if(abs(firstgpsutc-oldfirstgps).gt.0.100) then
                 write(*,*)' FZ ** First gps time is incorrect,'
                 write(*,'(''     ** correcting from: '',
     &    f10.3,'' sec to: '',f10.3,'' sec'')')oldfirstgps,firstgpsutc
               end if
               firsttstr = timestr
               utdec = firstgpsutc/3600.0
               lgpsutc = firstgpsutc

            else if(inevents.gt.2) then
*
* Accumulate deviations between deltat for GPS clock and for oscillator
* after the first valid calibration of the oscilator has been obtained.
* This serves as a diagnostic only.
* 
*
              if(fosccal.and.(mod(inevents,50).eq.0)) then
*               timedev=(gdf_ev10%gate_close-lgate_close)/caloscfreq -
*    &           (gpsutc-lgpsutc)
                timedev=(gclose-lastgclose)/caloscfreq -
     &           (gpsutc-lgpsutc)
****
*
* Code for debugging timing problems, comment out if not wanted
*
*               if(idbg.lt.15) then
*                 idbg = idbg + 1
*                 write(*,'('' '',i4,'': dgps [sec]:'',f8.4,
*    &             '' dosc [sec]:'',f8.4, 
*    &             '' tdev:'',f8.4)')
*    &             inevents,gpsutc-lgpsutc,
*    &             (gclose - lastgclose)/caloscfreq,
*    &             timedev 
*               endif
*              
*
****
                sumtdev = sumtdev + timedev
                sumtdev2 = sumtdev2 + timedev*timedev
                ntdev = ntdev + 1.0
*               lgate_close = gdf_ev10%gate_close
* MAC 970911 Just a check
*                type*,'livetime from gates',livetime/nomoscfreq
*                type*,'livetime from scaler',livetm_sc
*                type*,'difference=',livetime/nomoscfreq-livetm_sc
*                type*,'*****'
              end if ! osccal and every 50th event
            end if ! nevents.gt.2

            lgpsutc = gpsutc 
            lastgopen=gopen
            lastgclose=gclose           
*
* Do minute handling
*
            abssec = gdf_ev10%utc
            second = imark_gps
*MAC970911            second = gdf_ev10%mark_gps

            if ((second .ne. lastsec).and.(.not.disgpsmark)
     >           .and.uselvgate) then
               secmarks = secmarks + 1
               dsecond = dfloat(second)
               if (dsecond .lt. 0.) dsecond = dsecond + two32
               if (secmarks .ne. 1) then
                  dsecond = dsecond + nrollsec*two32
                  if (dsecond .lt. dlastsecond) then
                     dsecond = dsecond + two32
                     nrollsec = nrollsec + 1	! No. times couter
						! rolled over?
                  end if
               end if 
               lastsec = second
               dlastsecond = dsecond

	       if(((secmarks/60)*60).eq.secmarks) then ! If new minute
                  minmarks = minmarks + 1
                  dminute = dsecond

                  if(minmarks.lt.MAXMINS) then
                    dmintime(minmarks) = dminute
                  else
	            write(*,*),'  FZ ** WARNING: Run exceeded 200 minutes'
*		    write(*,*),' oscillator frequency calibration'
*	            write(*,*),' is based on the last second of the'
*		    write(*,*),' 200th minute. '
                  end if
                  icode = 6	! Each new minute marker is taken to
				! be a code 6 event.
*old                  if ((gdf_ev10%grs_clock(1).gt.0).or.
*old     >                 (gdf_ev10%grs_clock(2).gt.0)) then
                  if (defdate.gt.970901) then
                     grs(1) = gdf_ev10%grs_clock(1)
                     grs(2) = gdf_ev10%grs_clock(2)
                     grs(3) = gdf_ev10%grs_clock(3)
                     call truetime(grs,year,mjd,gpsutc,timestr)
                  else
                     gps(1) = gdf_ev10%gps_clock(1)
                     gps(2) = gdf_ev10%gps_clock(2)
                     gps(3) = gdf_ev10%gps_clock(3)
                     call gpstime(gps,year,mjd,gpsutc,timestr)
                  endif
                  iextra = iextra + 1
                  ncode6 = ncode6 + 1
                  if(minmarks.gt.1) then
                    gpsdev = dabs(gpsutc - lastgpsminute - 60.0d0)
                    if(gpsdev.ge.0.5d0) then	! This means at least one
                      missing_secmark = .true.	! second mark is missing,
		      write(*,*)'  FZ ** WARNING: Missing second mark!!!'
                                                ! and osc calib is impossible
**
* The following lines were added by JB 950628 to calculate the
* oscilator frequency when possible
*
* Use the last oscfreq (Note oscfreq(1)=nomoscfreq) if a minute mark
* is dropped.
*
                      if(minmarks.lt.MAXMINS) then
                        oscfreq(minmarks) = oscfreq(minmarks-1)   
                      end if
                    else
                      if(minmarks.lt.MAXMINS) then
                       oscfreq(minmarks) = (dminute-dlastminute)/60.0d0
                       if(.not.fosccal) then
                         fosccal = .true.
                         justcald = .true.
                         caloscfreq = oscfreq(minmarks)
                       endif
                      endif
                    end if		
*
**
                    if(minmarks.lt.MAXMINS) then
                      evtpmin(minmarks) = ncode8 - evtlastmin
                    endif
                  end if !minmarks.gt.1
                  if(justcald) then
                     justcald = .false.
                     write(6,170) minmarks,gpsutc,evtpmin(minmarks),
     >                    oscfreq(minmarks)/1.0d06
                  else
                     if(minmarks.eq.1) then
                        write(6,171) minmarks,gpsutc,
     >                       evtpmin(minmarks),oscfreq(minmarks)/1.0d06
                     else
                        write(6,172) minmarks,gpsutc,evtpmin(minmarks),
     >                       oscfreq(minmarks)/1.0d06
                     end if
                  end if
*
* Calculate variance of event rate
*
                  if(minmarks.gt.1) then
                    sumerate = sumerate + evtpmin(minmarks)
                    sum2erate = sum2erate + 
     &               (evtpmin(minmarks)*evtpmin(minmarks))
                    nmins = nmins + 1.0
                  endif
*
                  lastgpsminute = gpsutc
                  evtlastmin = ncode8
*                 write(scunit)icode,dminute,zero
*
* ??? WARNING: Check this next line of code later
* 
* MAC 970912 
                  write(scunit)icode,elapsed_time,gpsutc,
     +                 max(0.d0,live_time), 
     +                 (zero(chi),chi=1,nadc)
*old                  write(scunit)icode,timtick,gpsutc,livetime, 
*old     +                 (zero(chi),chi=1,nadc)
                  dlastminute = dminute
               end if ! new minute mark
*              dlastminute = dminute
            end if ! new second mark
*
* Do event handling
*
					! Event time is given by the
					! osc. time when the event
					! gate is closed.
            timtick = dfloat(igate_close)
*MAC970911            timtick = dfloat(gdf_ev10%gate_close)
            if (timtick .lt. 0) timtick = timtick + two32
					! as usual compensate for the
					! scaler rolling over...
            timtick = timtick + nroll*two32
            if (timtick .lt. lasttick) then
               nroll = nroll + 1
               timtick = timtick + two32
            end if
            lasttick = timtick

*MAC 970912 - elapsed time is only filled after 970901
            time_sc=(dfloat(gdf_ev10%elapsed_sec+elsec_off-first_sec)
     >           +dfloat(gdf_ev10%elapsed_ns+elns_off-first_ns)/1.e9)
     >           *defosc_freq/elosc_freq
*debug            if (mod(inevents,20).eq.0) then
*debug               type*,'elapsed time, sec, ns',time_sc,
*debug     >              gdf_ev10%elapsed_sec,gdf_ev10%elapsed_ns
*debug            endif
            if (inevents.gt.1) then
*
* Kludge for screwy livetime oscillator values (they get shifted 
* by multiples of 2^26 oscillations.)
*
               elsec_tmp=0
               elns_tmp=0
               if ((defdate.gt.980429).and.
     >              (defdate.lt.980523)) then
* The live time went backward.
                  if (time_sc.lt.((dfloat(elsec_old)+
     >                 dfloat(elns_old)/1.d9)
     >                 *defosc_freq/elosc_freq)) then
* Try to correct it by adding multiples of 2^26 cycles/10^7 Hz until the 
* time is greater than the live time of the last event.
                     do while ((time_sc+(dfloat(elsec_tmp)+
     >                    dfloat(elns_tmp)/1.d9)
     >                    *defosc_freq/elosc_freq)
     >                    .lt.((dfloat(elsec_old)+
     >                    dfloat(elns_old)/1.d9)
     >                    *defosc_freq/elosc_freq))
                        elsec_tmp=elsec_tmp+6
                        elns_tmp=elns_tmp+710886400
                        if (elns_tmp.gt.1000000000) then
                           elsec_tmp=elsec_tmp+1
                           elns_tmp=elns_tmp-1000000000
                        endif
                     enddo
* If the live time difference after correction is 1 sec, assume this
* was a real rollover and add to the offset values.
                     if ((time_sc-
     >                    (dfloat(elsec_old-elsec_tmp)+
     >                    dfloat(elns_old-elns_tmp)/1.d9)
     >                    *defosc_freq/elosc_freq)
     >                    .gt.1.0d0) then
                        elsec_tmp=0             !Reset these if not used
                        elns_tmp=0
                        elsec_off=elsec_old+first_sec
                        elns_off=elns_old+first_ns
* Correction for funny case where a rollover occurred but there is 
* significant elapsed time despite this being the first event.
                        if (abs(gdf_ev10%elapsed_sec).gt.0) then
                           if ((dfloat(gdf_ev10%elapsed_sec)+
     >                          dfloat(gdf_ev10%elapsed_ns)/1.d9)
     >                          .gt.0.d0) then
                              do while ((dfloat(gdf_ev10%elapsed_sec+
     >                          elsec_tmp)+dfloat(
     >                          gdf_ev10%elapsed_ns+elns_tmp)/1.d9)
     >                          .gt.6.d0)
                                 elsec_tmp=elsec_tmp-6
                                 elns_tmp=elns_tmp-710886400
                                 if (elns_tmp.lt.0) then
                                    elsec_tmp=elsec_tmp-1
                                    elns_tmp=elns_tmp+1000000000
                                 endif
                              enddo
                           else
                              do while ((dfloat(gdf_ev10%elapsed_sec+
     >                          elsec_tmp)+dfloat(
     >                          gdf_ev10%elapsed_ns+elns_tmp)/1.d9)
     >                          .lt.0.d0)
                                 elsec_tmp=elsec_tmp+6
                                 elns_tmp=elns_tmp+710886400
                                 if (elns_tmp.gt.0) then
                                    elsec_tmp=elsec_tmp+1
                                    elns_tmp=elns_tmp-1000000000
                                 endif
                              enddo
                           endif
                        endif
                        time_sc=(dfloat(gdf_ev10%elapsed_sec+
     >                       elsec_old+elsec_tmp)
     >                       +dfloat(gdf_ev10%elapsed_ns+
     >                       elns_old+elns_tmp)/1.d9)
     >                       *defosc_freq/elosc_freq
                     else
                        time_sc=time_sc+(dfloat(elsec_tmp)
     >                       +dfloat(elns_tmp)/1.d9)
     >                       *defosc_freq/elosc_freq
                     endif
                  else
* The elapsed time jumped by too much.
                     if ((time_sc-((dfloat(elsec_old)
     >                    +dfloat(elns_old)/1.d9))
     >                    *defosc_freq/elosc_freq).gt.6.d0) then
                        do while ((time_sc-
     >                       (dfloat(elsec_old-elsec_tmp)+
     >                       dfloat(elns_old-elns_tmp)/1.d9)
     >                       *defosc_freq/elosc_freq)
     >                       .gt.6.7108864d0)
                           elsec_tmp=elsec_tmp-6
                           elns_tmp=elns_tmp-710886400
                           if (elns_tmp.lt.0) then
                              elsec_tmp=elsec_tmp-1
                              elns_tmp=elns_tmp+1000000000
                           endif
                        enddo
* Weird
                        if (abs(time_sc-
     >                       (dfloat(elsec_old-elsec_tmp)+
     >                       dfloat(elns_old-elns_tmp)/1.d9)
     >                       *defosc_freq/elosc_freq)
     >                       .gt.1.0d0) then
                           etwarns=etwarns+1
                           if (etwarns.le.10) then
                             type*,' FZ ******* WARNING *******'
                             type*,' Elapsed time jumped forward ', 
     >                         '(not by multiple of 2^26) - probably', 
     >                         ' inaccurate.'
*                            type*,'***** WARNING *****'
*                            type*,'Elapsed time jumped forward but not'
*                            type*,'by multiple of 2^26 - elaps time is'
*                            type*,'probably inaccurate'
*                            type*,'*******************'
                           endif
                           if (etwarns.eq.10) then
                              type*,' ********************************',
     >                             '************************'
                              type*,'  Too many warnings of this type ', 
     >                             'to continue displaying!'
                              type*,' ********************************',
     >                             '************************'
                           endif
                        endif
   
                        time_sc=time_sc+(dfloat(elsec_tmp)
     >                       +dfloat(elns_tmp)/1.d9)
     >                       *defosc_freq/elosc_freq
                     endif
                  endif   !End if for check if elapsed time went backward
               else
* This should never happen, but try to fix it just in case
                  if (time_sc.lt.(dfloat(elsec_old)
     >                 +dfloat(elns_old)/1.e9)
     >                 *defosc_freq/elosc_freq) then
                     if (abs(time_sc-(dfloat(elsec_old)
     >                    +dfloat(elns_old)/1.e9)
     >                    *defosc_freq/elosc_freq).lt.10.) then
                        type*,' FZ ************  WARNING  ************'
                        type*,' Elapsed time scaler went backwards!'
                        type*,' Run duration may be off !'
                        type*,' FZ ***********************************'
                     else
                        type*,' FZ ************  WARNING  ************'
                        type*,' Elapsed time scaler rolled over !'
                        type*,' Run duration may be off !'
                        type*,' FZ ***********************************'
                        elsec_off=elsec_old+first_sec
                        elns_off=elns_old+first_ns
*                        type*,'time_sc=',time_sc
                        time_sc=(dfloat(gdf_ev10%elapsed_sec+elsec_old)
     >                       +dfloat(gdf_ev10%elapsed_ns+elns_old)/
     >                       1.d9)*defosc_freq/elosc_freq
*                       type*,'first_sec,first_ns',first_sec,first_ns
*                       type*,'elsec_old,elns_old',elsec_old,elns_old
*                       type*,'elsec,elns',gdf_ev10%elapsed_sec,
*     >                    gdf_ev10%elapsed_ns
*                       type*,'time_sc=',time_sc
                     endif
                  else
* Weird
                     if (abs(time_sc-(dfloat(elsec_old)+
     >                    dfloat(elns_old)/1.d9)
     >                    *defosc_freq/elosc_freq).gt.1.0d0) then
                        type*,' FZ ***********  WARNING  ***********'
                        type*,' Elapsed time jumped forward >1s ',
     +                  '- may be inaccurate.'
*                        type*,' ***********************************'
                     endif
                  endif  
               endif !End of kludge check for 980429
               if ((gdf_ev10%elapsed_sec.eq.elsec_old)
     >              .and.(gdf_ev10%elapsed_ns.eq.elns_old)
     >              .and.(gdf_ev10%elapsed_ns+gdf_ev10%elapsed_sec
     >              .gt.0)) then
                  tmwarns=tmwarns+1
                  if (tmwarns.le.5) then
                     type*,' FZ **********  WARNING  ***********'
                     type*,' Elapsed time scaler not increasing.'
                     type*,' ***********************************'
                  endif
                  if (tmwarns.eq.5) then
                    type*,' **********************************',
     >                    '*********************'
                    type*,'  Too many warnings of this type to', 
     >                    ' continue displaying!'
                    type*,' **********************************',
     >                    '*********************'         
                  endif
               endif  !End if for check that elapsed time scaler is increasing
            endif !End if for check that inevents.gt.1
            elsec_old=gdf_ev10%elapsed_sec+elsec_off-first_sec
     >           +elsec_tmp
            elns_old=gdf_ev10%elapsed_ns+elns_off-first_ns
     >           +elns_tmp
 

* If timtick is zero, gate_close is 0 or malfunctioning, so
* in that case use elapsed time from the elapsed time scaler.

            if ((timtick.eq.0.).or.(.not.usetmtick)) then
               if (usetmtick) then
                  type*,' FZ ********  WARNING  ********'
                  type*,' Elapsed time from gate close is 0 ',
     >                  '- using elapsed time scaler'
                  usetmtick=.false.
               endif
               elapsed_time=time_sc
            else
               elapsed_time=timtick
            endif

*
* Compare elapsed time scaler and the GRS times as a diagnostic.  Also,
* compare our conversion of GRS bits to that of the data - it's not going
* to be exactly the same because we convert the GPS sub-second bits with
* a frequency slightly different than 10 MHz.
*
*            type*,'gdf_ev10%utc',gdf_ev10%utc
*            type*,'Our GRS time',mjd
            delmjd=(mjd-gdf_ev10%utc)*8.64d10
            if (dabs(delmjd).gt.1.d1) then
               ndelmjd=ndelmjd+1
*debug               type*,'mjd',mjd
*debug               type*,'grs',gdf_ev10%utc
*debug               type*,'delmjd (usec)',delmjd
            endif
*            type*,'elapsed time',time_sc*1.d3
*            type*,'elapsed UTC time',(gdf_ev10%utc-
*     >           firstutc)*8.64d7
            idtm=max(1,min(80,
     >           int(((gdf_ev10%utc-firstutc)*8.64d7-time_sc*1.d3)
     >           /1.d-1)+40))
*            type*,'idtm',idtm
            deltime(idtm)=deltime(idtm)+1.
            sx=sx+sngl((gdf_ev10%utc-firstutc)*1.44d3)
            sxx=sxx+(sngl((gdf_ev10%utc-firstutc)*1.44d3))**2
            sy=sy+sngl((gdf_ev10%utc-firstutc)*8.64d7-time_sc*1.d3)
            sxy=sxy+sngl((gdf_ev10%utc-firstutc)*1.44d3*
     >           ((gdf_ev10%utc-firstutc)*8.64d7-time_sc*1.d3))
            ieltm=max(1,min(60,
     >          int(sngl(gdf_ev10%utc-firstutc)*1440./0.5)+1))
            ntm_v_tm(ieltm)=ntm_v_tm(ieltm)+1.
            dtm_v_tm(ieltm)=dtm_v_tm(ieltm)+
     >           sngl((gdf_ev10%utc-firstutc)*8.64d7-time_sc*1.d3)
            avedeltime=avedeltime
     >           +sngl((gdf_ev10%utc-firstutc)*8.64d7-time_sc*1.d3)
            ndtm=ndtm+1

            icode = 8			! A code 8 is a normal event
            do i = 1,nadc
               sshort(i) = gdf_ev10%adc(i)
            end do
            call swaptubes(idate,sshort)! Correct for the entire history
					! of channel swapping errors
            ncode8 = ncode8 + 1

*            write(scunit)icode,timtick,sshort
*            write(scunit)icode,livetime,sshort   !JQ 941202
* MAC 970912
            write(scunit)icode,elapsed_time,gpsutc, 
     +                 max(0.d0,live_time), 
     +           (sshort(chi),chi=1,nadc)
*old            write(scunit)icode,timtick,gpsutc,livetime, !JB 950408
*old     +           (sshort(chi),chi=1,nadc) !JQ 961205

         end if  !End if for check that this is a ev10 record
 997     continue           !Go here if the ev10 record is a pedestal.
*m         type*,'4'
*m         type*,' '
*
* Example of how to 
* get minute marks from timtick rather than the gps second mark
* for the purpose of calculating the event rate.
*
* Think about this later.
*
*       mmark = int((timtick-firstev)/(60.0*nomoscfreq))

*
* Get injected pedestals from frame header
* 970920 MAC After 970901 pedestals are no longer taken as 
* frame records, but as real events with the gdf_ev10%trigger set
* value set to 0.  So, put in date check for this possibility.
* There is also no longer a second pedestal taken so it seemed
* easier just to put the date switch in.
*
         if (defdate.gt.970901) then
            if (event1.and.(gdf_ev10%new).and.(gdf_ev10%valid).and.
     >           (iand(ishft(itrig_bits,-gdf_trig_ped),1).eq.1)) then
               ipeds=ipeds+1
* Here keep track of the real oscillator frequency using the elapsed time 
* oscillator oscillations between GPS second marks (when pedestals are taken).
               nel_sec=(dfloat(gdf_ev10%elapsed_sec)+
     >              dfloat(gdf_ev10%elapsed_ns)/1.d9)
     >              *defosc_freq/elosc_freq
*m         type*,'4.1'
*m         type*,' '
               if (nel_sec_old.ne.0.d0) then
* Only use events with reasonable elapsed second values.
                  if ((nel_sec.gt.nel_sec_old).and.
     >                 (nel_sec.gt.0.d0).and.
     >                 (nel_sec_old.gt.0.d0)) then
                     dtics=(nel_sec-nel_sec_old)*elosc_freq
c*debug                     type*,'dtics el',dtics/
c*debug     >                    ((gdf_ev10%utc-utc_old)*8.64d4)
* Do not use event if last pedestal event was missed (although you could).
*m                     type*,'gdf_ev10%event',gdf_ev10%event
*m                     type*,'gdf_ev10%utc',gdf_ev10%utc
*m                     type*,'utc_old',utc_old
                     if (abs(gdf_ev10%utc-utc_old)*8.64d4
     >                    .lt.1.1d0) then
                        ntic_chks=ntic_chks+1
                        nave_tics=nave_tics+dtics
                        navesq_tics=navesq_tics+dtics*dtics
*save                        nave_tics=nave_tics+dtics/
*save     >                    ((gdf_ev10%utc-utc_old)*8.64d4)
*save                        navesq_tics=navesq_tics+dtics*dtics/
*save     >                    (((gdf_ev10%utc-utc_old)*8.64d4)**2)
                        itics=max(1,min(201,nint(dtics-
     >                       elosc_freq+100.5d0)))
                        osc_dev(itics)=osc_dev(itics)+1.
*m         type*,'4.2'
*m         type*,' '
                     endif
c*debug                  else
c*debug                     type*,'Elapsed time went backwards!'
c*debug                     type*,'nel_sec',nel_sec
c*debug                     type*,'nel_sec_old',nel_sec_old
*
* Calculate mean and variance of event rate
*
                     if (mod(ipeds,60).eq.0) then
                        minmarks=minmarks+1
                        if (minmarks.gt.1) then
                           sumerate = sumerate + 
     >                          float(gdf_ev10%event-event_old)
     >                          /(sngl(gdf_ev10%utc-utc_min_old)
     >                          *1440.)
                           sum2erate = sum2erate + 
     >                          (float(gdf_ev10%event-event_old)
     >                          /(sngl(gdf_ev10%utc-utc_min_old)
     >                          *1440.))**2
*m         type*,'4.3'
*m         type*,' '
                           nmins=nmins+1.0
                        endif
                        utc_min_old=gdf_ev10%utc
                        event_old=gdf_ev10%event
                     endif  !End of check that this is multiple of 60 ped evts 
                  endif  !End of check that elapsed seconds is new.
               endif !End of check that this is not the first elapsed sec.
*m         type*,'4.4'
*m         type*,' '
               utc_old=gdf_ev10%utc
               nel_sec_old=(dfloat(gdf_ev10%elapsed_sec)+
     >              dfloat(gdf_ev10%elapsed_ns)/1.d9)
     >              *defosc_freq/elosc_freq
*m         type*,'4.5'
*m         type*,' '
               icode = 1                ! Code 1 is a pedestal event
               iextra = iextra + 1
               do i = 1,nadc
                  sshort(i) = gdf_ev10%adc(i)
               end do
               call swaptubes(idate,sshort)
               write(scunit)icode,elapsed_time,gpsutc,
     +                 max(0.d0,live_time), 
     +              (sshort(chi),chi=1,nadc)
               ncode12=ncode12+1
            endif  !End of check for pedestal events.
*m         type*,'5'
*m         type*,' '
* Here the singles rate scalers should still be packed into the frame
* records, so this flag should be set for those records.
* Keep track of the average and std dev. of the singles rates.
            if ((event1).and.
     >           (gdf_fr10%new).and.(gdf_fr10%valid)) then
               if (frutc_old.gt.0.d0) then
                  nsings=nsings+1
                  if (gdf_ev10%utc-frutc_old.lt.1.1d0) then
                     do j=1,npmt
*save                     ave_sing_rate(j)=ave_sing_rate(j)+
*save     >                    gdf_fr10%scals(j)
*save     >                    /(sngl(gdf_ev10%utc-frutc_old)*8.64e4)
*save                     dev_sing_rate(j)=dev_sing_rate(j)+
*save     >                    (gdf_fr10%scals(j)
*save     >                    /(sngl(gdf_ev10%utc-frutc_old)*8.64e4))**2
                        ave_sing_rate(j)=ave_sing_rate(j)+
     >                       gdf_fr10%scals(j)
                        dev_sing_rate(j)=dev_sing_rate(j)+
     >                       (gdf_fr10%scals(j))**2
                     enddo
                  endif
               endif
               frutc_old=gdf_ev10%utc
            endif !End of check for frame records.
*m         type*,'6'
*m         type*,' '
         else ! Date switch
            if ((event1).and.
     >           (gdf_fr10%new).and.(gdf_fr10%valid)) then
               icode = 1                ! Codes 1 and 2 are long and
					! short peds.  Now code-1 means
					! the ADC vals for the first
					! random trigger and code-e mean
					! ADC vals for the second rnd
					! trigger
               iextra = iextra + 1
               do i = 1,nadc
                  sshort(i) = gdf_fr10%ped_adc1(i)
               end do
               call swaptubes(idate,sshort)
*	    type *,'code ',icode
*	    type *,sshort
*            write(scunit)icode,timtick,sshort
*            write(scunit)icode,livetime,sshort     !JQ 941202
               write(scunit)icode,elapsed_time,gpsutc,
     +                 max(0.d0,live_time), 
     +              (sshort(chi),chi=1,nadc)
*old            write(scunit)icode,timtick,gpsutc,livetime,     !JB 950408
*old     +             (sshort(chi),chi=1,nadc)                 !JQ 961205
               icode = 2
               iextra = iextra + 1
               do i = 1,nadc
                  sshort(i) = gdf_fr10%ped_adc2(i)
               end do
               call swaptubes(idate,sshort)
*	    type *,'code ',icode
*            write(scunit)icode,timtick,sshort
*            write(scunit)icode,livetime,sshort     !JQ 941202
               write(scunit)icode,elapsed_time,gpsutc,
     +              max(0.d0,live_time), 
     +              (sshort(chi),chi=1,nadc)
*old            write(scunit)icode,timtick,gpsutc,livetime,     !JB 950408
*old     +             (sshort(chi),chi=1,nadc)                 !JQ 961205
               ncode12 = ncode12 + 2
            end if  !End if for check that this is a frame event.
         endif  !End if for date switch tag.

998      continue
      end do			! End of the first pass, now go back
				! and add the oscillator calibration
				! info, etc. for the second pass.
      if(.not. runid_exists) then
        write(*,*) ' FZ ** SEVERE: Run header was not detected.'
        write(*,*) '    ** Run number is missing from the header file.'
        if(ffixhdr) then
          write(*,*) ' FZ ** Please provide the missing information:'
          write(6,'(''Run id (xxxx)      :'',$)')
          read(5,*)irun
          write(crunid,'(i4.4)')irun
*         print *,'10m RUN: ',crunid
          if (idate .eq. 0) idate = defdate
          year = 1900 + idate/10000
          write(6,'(''Sky quality (A,B,C):'',$)')
          read(5,*)cskyq
          write(6,*) ' FZ ** Run information is now complete'
          runid_exists = .true.
        else
          write(*,*) '    ** RE-RUN FZ2RED USING THE FIX-HEADER OPTION:'
          write(6,'(''       Type: fz2red2 '',a11,'' '',a6,'' f'')')
     &     filename(1:11),cdate(1:6)
          write(*,*) '       at the command line.'
        end if
      end if
*
* EOF reached
*
      call gdf$close(10,ierr)
*      call gdf$print(gdf_run,ierr2)

*******************************************************************************
* MAC 980821 Here we sum up the diagnostic information and write some of it out
* to an hbook file (eventually) and some to an errors file.
*******************************************************************************
      open(eunit,file='gt'//crunid//'.errors',status='unknown')
* Initialize histograms
      call HLIMIT(PAW_NWORDS)
* Book the histograms
*Tracking histograms
      call hbook1(1,'Tracking deviations',40,0.,0.2,0.)
      call hbook1(2,'Ave Elev v. time',61,-0.5,30.,0.)
      call hbook1(3,'Ave Az v. time',61,-0.5,30.,0.)
      call hbook1(4,'Ave RA v. time',61,-0.5,30.,0.)
      call hbook1(5,'Ave Dec v. time',61,-0.5,30.,0.)
*Time histograms
      call hbook1(11,'Elapsed time diffs.',80,-4.,4.,0.)
      call hbook1(12,'Time diff v. time',60,0.,30.,0.)
      call hbook1(13,'Osc freq. - 1e7',201,-100.5,100.5,0.)
*Trigger histograms
      call hbook1(21,'Trigger bits',10,-0.5,9.5,0.)
*HV histograms
      call hbook1(31,'Ave. HV values',npmt,0.5,float(npmt)+0.5,0.)
      call hbook1(32,'Ave. Anode current values',
     >     npmt,0.5,float(npmt)+0.5,0.)
      call hbarx(31)
      call hbarx(32)
*Singles rate histograms
      call hbook1(41,'Ave. singles rates',
     >     npmt,0.5,float(npmt)+0.5,0.)
      call hbarx(41)

*
* Average values of the GRS clock bits: sub-seconds and seconds should be
* close to 0.5 but if they are stuck, will have aves of 0 or 1.
*
      if (nbt.gt.0) then
* Compute the average bit values
         do k=1,2
            do j=1,24
               grs_bits(k,j)=grs_bits(k,j)/float(nbt)
            enddo
         enddo
* Stuck clock bits
         do j=1,24
            if ((grs_bits(1,j).lt.0.05)
     >           .or.(grs_bits(1,j).gt.0.95)) then
               write(*,150) j,1,grs_bits(1,j)
               write(eunit,150) j,1,grs_bits(1,j)
            endif
         enddo
         do j=1,7
            if ((grs_bits(2,j).lt.0.05)
     >           .or.(grs_bits(2,j).gt.0.95)) then
               write(*,150) j,2,grs_bits(2,j)
               write(eunit,150) j,2,grs_bits(2,j)
            endif
         enddo
         if (grs_bits(2,8).gt.0.01) then
            write(*,151) 8,2,grs_bits(2,8)
            write(eunit,151) 8,2,grs_bits(2,8)
         endif
* Status bits
         if (abs(grs_bits(3,17)-float(nbt)).gt.0.5) then 
            write(*,152) 17,3,nbt-nint(grs_bits(3,17)),
     >           'sec mark may be bad'
            write(eunit,152) 17,3,nbt-nint(grs_bits(3,17)),
     >           'sec mark may be bad'
         endif
         if (abs(grs_bits(3,18)-float(nbt)).gt.0.5) then 
            write(*,152) 18,3,nbt-nint(grs_bits(3,18)),
     >           'check external osc. stability'
            write(eunit,152) 18,3,nbt-nint(grs_bits(3,18)),
     >           'check external osc. stability'
         endif
         if (abs(grs_bits(3,19)-float(nbt)).gt.0.5) then 
            write(*,152) 19,3,nbt-nint(grs_bits(3,19)),
     >           'RS232 connection may be bad'
            write(eunit,152) 19,3,nbt-nint(grs_bits(3,19)),
     >           'RS232 connection may be bad'
         endif
         if (abs(grs_bits(3,20)-float(nbt)).gt.0.5) then 
            write(*,152) 20,3,nbt-nint(grs_bits(3,20)),
     >           'RS232 time does not match internal time'
            write(eunit,152) 20,3,nbt-nint(grs_bits(3,20)),
     >           'RS232 time does not match internal time'
         endif
         if ((grs_bits(3,21).gt.0.01).or.
     >        (grs_bits(3,22).gt.0.01).or.
     >        (grs_bits(3,23).gt.0.01).or.
     >        (grs_bits(3,24).gt.0.01)) then
            write(*,153) nint(grs_bits(3,21))
            write(eunit,153) nint(grs_bits(3,21))
         endif
      endif !End of check that there were some GRS times recorded.
*
* Tracking records
*
      if (ntrack.gt.0.0) then
         do j=1,40
            call hpak(1,trk_devs)
         enddo
         do j=1,61
            if (ntrks_tm(j).gt.0.) then
               el_v_tm(j)=el_v_tm(j)/ntrks_tm(j)
               az_v_tm(j)=az_v_tm(j)/ntrks_tm(j)
               ra_v_tm(j)=ra_v_tm(j)/ntrks_tm(j)
               dec_v_tm(j)=dec_v_tm(j)/ntrks_tm(j)
            endif
         enddo
         call hpak(2,el_v_tm)
         call hpak(3,az_v_tm)
         call hpak(4,ra_v_tm)
         call hpak(5,dec_v_tm)
      
         meantrdev = sumtrdev/ntrack
         sigtrdev = sumtrdev2/ntrack - meantrdev*meantrdev
         if(sigtrdev .lt. 0.00) then 
            sigtrdev = 0.00
         end if
         sigtrdev = sqrt(sigtrdev)
         ave_el = (sum_el/ntrack)*convd
         if ((meantrdev.gt.0.08).or.(sigtrdev.gt.0.08)) then
            write(eunit,160) min(99.999,meantrdev),
     >           min(99.999,sigtrdev)
         endif
      else
         write(*,*)' FZ ** SEVERE: No tracking records were detected.'
         write(*,*)'    ** Image derotation and elevation parameter'
         write(*,*)'    ** calculation will fail for this file.'
         write(eunit,161) 
         if(ffixhdr) then
         else
            write(*,*) '    ** RE-RUN FZ2RED WITH FIX-HEADER OPTION:'
            write(6,'(''     ** At the command line, type:'')')
            write(6,'(''     **     fz2red2 '',a11,'' '',a6,
     >           '' f'')') filename(1:11),cdate(1:6)
         end if
      endif       !End of check that there were tracking records

*
* Timing checks
*
* ndelmjd is > 0 if there are any events where the decoding of the GRS
* bits done here results in a time more than 1 usec different than 
* DACQs decoded GRS time.
      if (ndelmjd.gt.0) then
         write(*,154) ndelmjd
         write(eunit,154) ndelmjd
      endif
* Histogram of differences in elapsed time between GRS and elapsed time
* oscillator.
      call hpak(11,deltime)
      if (ndtm.gt.0) then
         avedeltime=avedeltime/float(ndtm)
         if ((avedeltime.gt.1.5).and.(defdate.gt.970901)) then
            write(eunit,157) avedeltime
         endif
      endif
* Estimate the true 10 MHz oscillator frequency based on time differences
* between GRS clock time and elapsed time scaler as a function of time.
* If 10 MHz is not the true frequency, plotting this time difference 
* versus time will yield roughly a line.  The dividing 10 MHz by the one
* plus the slope of this line will give the true frequency.  The 10 MHz
* frequency should be good to 1 part in 10^6 so we don't expect any big
* differences.
      s=float(inevents)
      dtslope=(s*sxy-sy*sx)/(s*sxx-sx*sx)
* Convert slope from msec/min to msec/msec.
      dtslope=dtslope/6.e4
      cal_osc_freq=sngl(elosc_freq)/(1.+dtslope)
      if ((abs(cal_osc_freq-sngl(elosc_freq)).gt.100.).and.
     >     (defdate.gt.970901))then
         write(eunit,156) cal_osc_freq
      endif
* Histogram the mean difference in elapsed time estimates as a function of 
* elapsed time.
      do j=1,60
         if (ntm_v_tm(j).gt.0.) then
            dtm_v_tm(j)=dtm_v_tm(j)/ntm_v_tm(j)
         endif
      enddo
      call hpak(12,dtm_v_tm)
* Average and std. dev. of 10 MHz oscillator frequency based on number of
* elapsed oscillator ticks between GPS second marks.
      if (ntic_chks.gt.0) then
         ave_osc=nave_tics/dfloat(ntic_chks)
         dev_osc=sqrt((navesq_tics-
     >        nave_tics*nave_tics/dfloat(ntic_chks))/dfloat(ntic_chks))
         if (abs(ave_osc-elosc_freq).gt.100.d0) then
            write(eunit,156) ave_osc
         endif
         j=0
         sum_osc=0.
         do while (sum_osc.lt.0.5*float(ntic_chks)
     >        .and.(j.lt.201))
            j=j+1
            sum_osc=sum_osc+osc_dev(j)
         enddo
         med_osc=elosc_freq-101.+float(j)
      endif
      call hpak(13,osc_dev)

      if (ltwarns>10) then 
         write(*,186)ltwarns
      endif
      if (etwarns>10) then
         write(*,185)etwarns
      endif
*
* Trigger bits checks.
*
      call hpak(21,ntrig_bits)
      if (ntrig_bits(4).gt.0.) then
         write(*,180) 'Ped and PST',nint(ntrig_bits(4))
         write(eunit,180) 'Ped and PST',
     >        nint(ntrig_bits(4))
      endif
      if (ntrig_bits(5).gt.0.) then
         write(*,180) 'Ped and Mult. Trig',
     >        nint(ntrig_bits(5))
         write(eunit,180) 'Ped and Mult. Trig',
     >        nint(ntrig_bits(5))
      endif
      if (ntrig_bits(7).gt.0.) then
         write(*,180) 'Ped, PST, and Mult. Trig',
     >        nint(ntrig_bits(7))
         write(eunit,180) 'Ped, PST, and Mult. Trig',
     >        nint(ntrig_bits(7))
      endif
      if (ntrig_bits(8).gt.0.) then
         write(*,180) 'PST but not Mult Trig',
     >        nint(ntrig_bits(8))
         write(eunit,180) 'PST but not Mult Trig',
     >        nint(ntrig_bits(8))
      endif
*
* HV records
*
      if (n_hv.gt.0) then
         do j=1,npmt
            ave_hv(j)=ave_hv(j)/float(n_hv)
            dev_hv(j)=sqrt(max(1.e-10,
     >           dev_hv(j)/float(n_hv)-ave_hv(j)*ave_hv(j)))
            ave_ai(j)=ave_ai(j)/float(n_hv)
            dev_ai(j)=sqrt(max(1.e-10,
     >           dev_ai(j)/float(n_hv)-ave_ai(j)*ave_ai(j)))
         enddo
      endif
      call hpak(31,ave_hv)
      call hpake(31,dev_hv)
      call hpak(32,ave_ai)
      call hpake(32,dev_ai)

*
* Singles rate records
*
      if (nsings.gt.0) then
         do j=1,5
            max_sing_rate(j,2)=-1.
            min_sing_rate(j,2)=9.9999e20
         enddo
         do j=1,npmt
            ave_sing_rate(j)=min(99999.,ave_sing_rate(j)/float(nsings))
            dev_sing_rate(j)=min(99999.,
     >           (sqrt(dev_sing_rate(j)/float(nsings)
     >           -ave_sing_rate(j)*ave_sing_rate(j))))
            if (ave_sing_rate(j).gt.max_sing_rate(5,2)) then
               k=5
               do while ((ave_sing_rate(j).gt.max_sing_rate(k,2))
     >              .and.(k.ge.1))
                  k=k-1
               enddo
               k=k+1
               do l=5,k+1,-1
                  max_sing_rate(l,1)=max_sing_rate(l-1,1)
                  max_sing_rate(l,2)=max_sing_rate(l-1,2)
               enddo
               max_sing_rate(k,1)=float(j)
               max_sing_rate(k,2)=ave_sing_rate(j)
            endif
            if (ave_sing_rate(j).lt.min_sing_rate(5,2)) then
               k=5
               do while ((ave_sing_rate(j).lt.min_sing_rate(k,2))
     >              .and.(k.ge.1))
                  k=k-1
               enddo
               k=k+1
               do l=5,k+1,-1
                  min_sing_rate(l,1)=min_sing_rate(l-1,1)
                  min_sing_rate(l,2)=min_sing_rate(l-1,2)
               enddo
               min_sing_rate(k,1)=float(j)
               min_sing_rate(k,2)=ave_sing_rate(j)
            endif
         enddo
         call hpak(41,ave_sing_rate)
         call hpake(41,dev_sing_rate)
      else
         write(eunit,200)
      endif !End of check that there are singles rates records.

* Close errors file
      close(eunit)
      call hrput(0,'gt'//crunid//'_d.hbook','T')
*******************************************************************************
* End of diagnostic checks.
*******************************************************************************

* MAC 970916 Modify to use the livetime scaler value if the gate_open
* and gate_close fields are missing. Convert from the cycles to seconds
* for the livetime if gate_open and gate_close are used.
*
      if (uselvgate) then
         live_time=live_time/caloscfreq
      endif

      if(nmins.gt.0) then
        meanerate = sumerate/nmins
        sigerate = sum2erate/nmins - meanerate*meanerate
        if(sigerate .lt. 0.00) then
          sigerate = 0.00
        end if
        sigerate = sqrt(sigerate)/sqrt(meanerate)
      else
        write(*,*)
        write(*,*)' FZ ** Rate and rms spread can not be calculated.'
      endif 

      if(ntdev.gt.0) then
        meantdev = sumtdev/ntdev
        sigtdev = sumtdev2/ntdev - meantdev*meantdev
        if(sigtdev .lt. 0.00) then
          sigtdev = 0.00
        end if
        sigtdev = sqrt(sigtdev)
      endif

      dmjd = gdf_run%utc_start	! JB 950408
*
* Calculate the nominal sidereal start time from the acquisition
* computer (VHEGRO) utc start time
*
      call mjd2sid(dmjd,nstime)	! calculate sidereal start time (nominal)
      call timefmt(nstime,fnsidstr,ist)
*
* Calculate the sidereal time of the last event and convert to a string
*
      call mjd2sid(mjd,laststime)
      call timefmt(laststime,lastsidstr,ist)
*
* Calculate the sidereal start time and convert to a string, and an 
* integer representation of the first sidereal minute.
*
      call mjd2sid(fmjd,stime)
      call timefmt(stime,fsidstr,ist)
      sid_duration = (laststime-stime)*60.0d0	! Sid duration [min]
*
* Convert the ut start time to a string.
*
      call timefmt(utdec,futstr,iut)
*
* Calculate az and alt that the telescope should be at for the first
* event time.  This is for comparison with the values given by the 
* first tracking record, which will be somewhat delayed.
*
      if((ntrack.gt.0.0).or.ffixhdr) then 
        call derot(ra_rad,dec_rad,dmjd,utdec,
     &   theta,alt,az,ttime)
      endif
*     write(6,'('' Run duration (elapsed UT from oscillator):'',f8.3)')
*    +  dduration

      open(hunit,File='gt'//crunid//'.inf',
     &  status='unknown')

      do u=6,hunit,hunit-6
       write(u,*)' '
       write(u,*)'Run Information: '
       write(u,*)'---------------  '
       write(u,*)'10m RUN                  : ',crunid
       write(u,'('' UT Date                  : '',i8)')19000000+idate
       write(u,'('' Source                   : '',a12)')csource(1:12)
       write(u,'('' Mode                     : '',i2)')mode
       write(u,*),' '
       write(u,*),'Tracking: '
       write(u,*),'-------- '
       write(u,'('' RA [hhmmss.s]            : '',f9.1)')rra
       write(u,'('' DEC [ddmmss.s]           : '',f9.1)')rdec
       write(u,'('' Ave. EL [deg]            : '',f8.3)')
     &   ave_el
       write(u,'('' Mean tracking error [deg]: '',f8.3)')meantrdev
       write(u,'('' RMS tracking error [deg] : '',f8.3)')sigtrdev
       if((meantrdev.gt.0.08).or.(sigtrdev.gt.0.08)) then
         write(u,*) ' FZ ** WARNING: Tracking errors are large enough'
         write(u,*) '       alpha plot might be significantly degraded'
       end if
       write(u,'('' AZ (track comp.) [deg]   : '',f8.3)')razimuth*convd
       write(u,'('' AZ (first evt calc) [deg]: '',f8.3)')az*convd
       write(u,'('' EL (track comp.) [deg]   : '',f8.3)')
     &  relevation*convd
       write(u,'('' EL (first evt calc) [deg]: '',f8.3)')alt*convd
       write(u,*)' '
       write(u,*)'Weather: '
       write(u,*)'-------  '
       write(u,*)'Sky quality              :     ',cskyq      
       write(u,'('' Mean event rate [min^-1] : '',f8.2)')meanerate
       write(u,'('' RMS event rate [sigma]   : '',f8.2)')sigerate 
       do i=1,3
        write(u,'('' CCD star at ('',i6,'','',i6'')  Amplitude:'',
     &   i6)')xstar(i),ystar(i),istar(i) 
       end do

* HV printout

       write(u,*),' '
       write(u,*),'High Voltage: '
       write(u,*),'------------ '

       if (n_hv.eq.0) then
          write(u,196) 
       else
          write(u,*) 'PMTs turned off according to HV records'
          write(u,*) '---------------------------------------'
          write(u,197) (pmtsoff(j),j=1,npmtsoff)
       endif
      
*
* Singles rates printout
*
       write(u,*),' '
       write(u,*),'Singles Rates: '
       write(u,*),'------------- '
       if (nsings.gt.0) then
          write(u,201)
          write(u,202)
          do j=1,5
             write(u,203) nint(max_sing_rate(j,1)),max_sing_rate(j,2),
     >            nint(min_sing_rate(j,1)),min_sing_rate(j,2)
          enddo
       else
          write(u,*) ' FZ ** No singles rate records in this run.'
       endif
*
       write(u,*),' '
       write(u,*),'Event and Run Timing: '
       write(u,*),'-------------------- '
*
* Actually this is the time string from the SECOND event not the first,
* but this is really only here
*
       write(u,'('' GPS [yr day hh:mm:ss.ss] : '',a20)')firsttstr(1:20)
       write(u,'('' Sidereal start time (GPS): '',a20)')fsidstr(1:20)
       write(u,'('' Sid. start time (VHEGRO) : '',a20)')fnsidstr(1:20)
       if(abs(stime - nstime).gt.5.0) then
        write(u,*) ' FZ ** WARNING: Mismatch in expected and actual '
        write(u,*) '    ** start time - VAXstation system clock may '
        write(u,*) '    ** be wrong (correct by typing sync at the '
        write(u,*) '    ** control prompt). ' 
       endif
       write(u,'('' Sid. time last evt (GPS) : '',a20)')lastsidstr(1:20)
       write(u,'('' Sid duration (GPS) [min] : '',f6.2)')sid_duration
       write(u,'('' Nom. sid run length [min]: '',f6.2)')
     &  gdf_run%sid_length
       write(u,'('' Nom. sid cycle time [min]: '',f6.2)')
     &  gdf_run%sid_cycle+gdf_run%sid_length
       write(u,'('' Track comp RA offset[min]: '',f6.2)')
     &  ra_offset
*      write(*,*)' ra_off2: ',ra_off2
*      write(u,'('' Offset for this run [min]: '',f6.2)')
*    &  ra_off2
       if(abs(ra_offset-gdf_run%sid_cycle-gdf_run%sid_length(1))
     &    .gt. 0.05) then
        if(idate.lt.950629) then
         write(u,*) ' FZ ** RELAX: Tracking computer wrote incorrect '
         write(u,*) '    ** offset before 950629 '
         write(u,*)
     &     '    ** Setting the RA offset to nominal value: 30min'
         ra_offset = 30.0
        else
         write(u,*) ' FZ ** WARNING: A sidereal time mismatch in on'
         write(u,*) '    ** and off runs is likely, since the cycle'
         write(u,*) '    ** time on VHEGRO does not match the tracking'
         write(u,*) '    ** computer''s RA offset. '
        end if
       end if
       write(u,'('' UTC start MJD (GPS)      : '',f12.5)')
     &  fmjd
       write(u,'('' UTC start MJD (VHEGRO)   : '',f12.5)')
     &  gdf_run%utc_start
       write(u,'('' UTC end [mjd]            : '',f12.5)')
     &  gdf_run%utc_end
       write(u,'('' UTC last event [mjd]     : '',f12.5)')gdf_run%utc
       write(u,'('' First sidereal minute    : '',i4)')ist
       write(u,'('' First UT minute          : '',i4)')iut 
*      write(u,'('' VAX time                 : '',f11.4)')gdf_run%itime
       write(u,*) ' '
       if (defdate.lt.970901) then
          write(u,'('' Mean GPS-osc. dev [ms]   : '',f10.3)')
     >         meantdev*1000.0
          write(u,'('' RMS GPS-osc. dev [ms]    : '',f10.3)')
     >         sigtdev*1000.0 
          if(sigtdev.gt.0.0015) then
             write(u,*) 
     >         ' FZ ** WARNING: GPS clock problem.  Time series'
             write(u,*) 
     >         '    ** analysis using GPS time should not be done.'
          endif
       else
          write(u,'('' Mean GPS-osc. dev [ms]   : '',f10.3)') 
     >         avedeltime
          write(u,'('' Nominal elapsed osc. freq. [Hz]   : '',
     >         e14.8)') elosc_freq
          write(u,'('' Ave. est. 10 MHz osc. freq. [Hz]  : '',
     >         e14.8)')ave_osc
          write(u,'('' Std. dev. of est. freq. [Hz]      : '',
     >        e14.8)')dev_osc
          write(u,'('' Median est. 10 MHz osc. freq. [Hz]: '',
     >         e14.8)') med_osc
          if (avedeltime.gt.1.5) then
             write(u,157) avedeltime
          endif
          if (abs(ave_osc-sngl(elosc_freq)).gt.100.) then
             write(u,156) ave_osc
          endif
       endif
*      write(u,'('' GPS first evt [sec>midnt]: '',f11.4)')firstgpsutc
       i = 1

*
       nominmarks = .true.
*      dduration = (lasttick - firstev) / nomoscfreq   !JQ 941220
* MAC 970916 Modified for the case when the gate_open and gate_close
* is missing.
       if (usetmtick) then
          dduration = (lasttick - firstev) / caloscfreq !JB 950628
       else
          dduration = time_sc
       endif

       write(u,'('' Number of gps second marks               : '',i4)')
     &  secmarks
       write(u,'('' Run duration (elapsed UT from oscillator): '',
     &  f8.3)')dduration
       write(u,'('' Run duration (elapsed GPS time [s])      : '',
     &  f8.3)')(gpsutc-firstgpsutc)
*
** Convert livetime scaler into sidereal minutes.
* Forget about the livetime scaler for now.  It seems to be unreliable
* and provides nothing other than the calculated livetime (sum of gate_close
* minus gate_open) JB 950625
*
*     livetime = dfloat(gdf_ev10%live_sec)*(20.0d6/nomoscfreq)  !8.3
*     write(u,'('' Live time (UT seconds)                   : '',f12.3)')
*    +  livetime
* MAC 970916
*old       write(u,'('' Live time (UT seconds)                   : '',
*    + f8.3,'' ('', f6.2,''%)'')') livetime/nomoscfreq,
*old     + f8.3,'' ('', f6.2,''%)'')') livetime/caloscfreq,	
*    +  (livetime/(nomoscfreq*(gpsutc-firstgpsutc)))*100.0
*old     +  (livetime/(caloscfreq*(gpsutc-firstgpsutc)))*100.0
       if (gpsutc-firstgpsutc.gt.0.) then
          write(u,'('' Live time (UT seconds)                   : '',
     +         f8.3,'' ('', f6.2,''%)'')') 
     +         max(0.d0,min(3.d4,live_time)), 
     +         (max(0.d0,min(3.d4,live_time))
     +         /(gpsutc-firstgpsutc))*100.0
       else
          write(u,'('' Live time (UT seconds)                   : '',
     +         f8.3,'' ('', f6.2,''%)'')') 
     +         max(0.d0,min(3.d4,live_time)),0.
       endif
       istdur = nint(max(0.d0,min(3.d4,live_time))*utst/60.0d0)
*old       istdur = nint(livetime*utst/60.0d0)
*      write(u,'('' Live time (SID minutes)       : '',i6)')
*    +  istdur
       write(u,*),' '
       write(u,*),'Breakdown of Events:'
       write(u,*),'-------------------'
       write(u,'('' Number of events (code 8)                : '',i6)')
     +  ncode8
       write(u,'('' Number of injected ped events            : '',i6)')
     +  ncode12
       write(u,'('' Number of minute mark events             : '',i6)')
     +  ncode6
       write(u,*)
*      write(u,*)
      end do
      close(hunit)
*
* Write fixed format header file
*
      open(hunit,File='gt'//crunid//'.hdr',
     &  status='unknown')
      write(hunit,'(''run_id           '',a4)'),crunid
      write(hunit,'(''date           '',i8)'),19000000+idate
      write(hunit,'(''source         '',a12)'),csource(1:12)
      write(hunit,'(''mode           '',i2)'),mode
      write(hunit,'(''ra             '',f9.1)')rra
      write(hunit,'(''dec            '',f9.1)')rdec
      write(hunit,'(''ave_el         '',f8.3)')ave_el
      write(hunit,'(''mean_trk_err   '',f8.3)'),meantrdev
      write(hunit,'(''rms_trk_err    '',f8.3)'),sigtrdev
      write(hunit,'(''sky_grade       '',a2)'),cskyq
      write(hunit,'(''mean_rate      '',f8.2)')meanerate
      write(hunit,'(''sigma_rate     '',f8.2)')sigerate
      write(hunit,'(''gps_sid_start  '',f8.4)')stime
      write(hunit,'(''gps_sid_dur    '',f6.2)')sid_duration
      write(hunit,'(''nom_sid_dur    '',f6.2)'),gdf_run%sid_length
      write(hunit,'(''trk_ra_offset  '',f6.2)')
     &  ra_offset
      write(hunit,'(''gps_mjd_start  '',f12.5)') fmjd
*     write(hunit,'(''nom_mjd_start  '',f12.5)'),gdf_run%utc_start
*     write(hunit,'(''nom_utc_end    '',f12.5)'),gdf_run%utc_end
      write(hunit,'(''first_sid_min  '',i6)')ist
      write(hunit,'(''first_ut_min   '',i6)')iut
      write(hunit,'(''mean_time_dev  '',f10.3)'),meantdev*1000.0
      write(hunit,'(''sig_time_dev   '',f10.3)'),sigtdev*1000.0
*     write(hunit,'(''first_evt_utc  '',f11.4)'),firstgpsutc
      write(hunit,'(''n_sec_marks    '',i6)'),secmarks
      write(hunit,'(''osc_duration   '',f8.3)'),dduration
      write(hunit,'(''gps_duration   '',f8.3)'),gpsutc-firstgpsutc
*     write(hunit,'(''livetime       '',f8.3)'),livetime/nomoscfreq
*old      write(hunit,'(''livetime       '',f8.3)'),livetime/caloscfreq
      write(hunit,'(''livetime       '',f8.3)'),
     +         max(0.d0,min(3.d4,live_time))
      write(hunit,'(''n_events       '',i6)'),ncode8
      write(hunit,'(''n_ped_events   '',i6)'),ncode12
      close(hunit)
*     write(*,*)'Comments: '
*     write(*,*)'-------- '
*
*     do while((i.lt.comlength).and.(i.lt.1600))
*       print *,longcomms(i:i+78)
*       i = i + 80
*     end do
*
*
*
*  Open output file and write header out.
*
*      open(20,file=pathout(1:lnblnk(pathout))//'gt'//crunid,
*     &  status='new',form='unformatted')
      open(20,file='gt'//crunid,
     &  status='unknown',form='unformatted')
*     &  status='new',form='unformatted')

      inevents = inevents + iextra
*
* Convert from IEEE floating point to VAX floating point for
* consistency with previous formats.
*
* 941017 JB -- abandon consistency with previous format. Switch on the date
* in the rest of the analysis code.  
*
*     call ieflt(razimuth)
*     call ieflt(relevation)
*     call ieflt(rra)
*     call ieflt(rdec)
*     call iedbl(dduration)
*     call iedbl(dmjd)
*     call iedbl(dfrjd)

* 
* Modified by JB 950615 - put livetime in place of elapsed osc.
* time in the header. 
*     write(20)crunid,inevents,dduration,istdur,cmode,csource,
*old      write(20)crunid,inevents,livetime/caloscfreq,istdur,cmode,
      write(20)crunid,inevents,
     +         max(0.d0,min(3.d4,live_time)),istdur,cmode,
     &         csource,idate,dmjd,dfrjd,rra,rdec,iut,ist,razimuth,
     &         relevation,cskyq,ccomms,igpsbeg
*
* Now make second pass through data.  Convert oscillator values to 
* event times and write them out.
*
      rewind(scunit)
      minmarks2 = 0
      nevts = 0
*
* ADK's old algorithm described below has been changed for the time
* being.
** Here we actually handle the data.  The toughest part is the time.
** All times are relative to the first event.  If we are before the
** first minute marker, then use the oscillator frequency as determined
** between markers 1 and 2; if we are in the middle of the run, use the
** time of the previous minute marker and the oscillator frequency for 
** the minute we're in; if we're after the last minute marker, then
** extrapolate from the time of the last minute using the last well-
** determined oscillator frequency.
*
900   continue
      nevts = nevts + 1
*     read(scunit,end=999)icode,timtick,sshort
      read(scunit,end=999)icode,timtick,gpsutc,rlivetime,   !JB 950408
     +             (sshort(chi),chi=1,nadc)                 !JQ 961205
      if (icode .eq. 6) then
         minmarks2 = minmarks2 + 1 
*        write(6,'('' Minute mark: '',i2,'' Osc. Freq.: '',f11.7,'' MHz'')')
*    &     minmarks2, oscfreq(minmarks2)/1.0d6

         if (minmarks2 .gt. 1) then
            dtime = tlastmin + 60.0d0
         else
            dtime = (timtick - firstev)/oscfreq(1)
         end if
         tlastmin = dtime
*
* Mod. by JB 950626 - If the first event is a minute marker, then 
* its gps time might be wrong, so replace it with the corrected value
* firstgpsutc
*
         if(first2min .and. (nevts .eq. 1)) then	
           first2min = .false.
           gpsutc = firstgpsutc
         endif
      else  !Non-time codes
*
* For non-time codes, convert the oscillator time and the livetime
* using the nominal oscillator frequency
* and record the gpstime in seconds.
*
             
*          dtime=timtick/nomoscfreq ! JQ 941220
*          rlivetime = rlivetime/nomoscfreq ! JB 950408

* MAC 970916
         if (usetmtick) then
            dtime=timtick/caloscfreq ! JB 950628
         else
            dtime=timtick
         endif
         if (uselvgate) then
            rlivetime = rlivetime/caloscfreq ! JB 950628
         endif
*
* Mod. by JB 950626 - If this is the first event, then 
* its gps time might be wrong, so replace it with the corrected value
* firstgpsutc
*
           if(first2evt .and. (icode .eq. 8)) then
             first2evt = .false.
             gpsutc = firstgpsutc
           endif
      end if
*      write(20)icode,dtime,sshort
      write(20)icode,dtime,gpsutc,rlivetime,     ! JB 950504
     +             (sshort(chi),chi=1,nadc)      ! JQ 961205
      go to 900
*
*
999   continue
      close(scunit)
      close(20)
      call gdf$exit(ierr)
      call exit(0)	! exit with status ok
*
* Formats
*
100   format(
     &     /,'  Paths for input and output files may be entered in a',
     &     /,'  file named "paths", with the input path on the first',
     &     /,'  line and the output path on the second line.  This',
     &     /,'  file was either not found or an error occurred',
     &     /,'  while reading it.',
     &     /)
110   format(' Enter path for input file: ',$)
120   format(' Enter path for output file: ',$)
130   format(' Enter input file name: ',$)
140   format(' File ',a,' does not exist.  Exiting...')

 141  format(' FZ ** First event has gdf_ev10%event = ',i7)

 150  format(' FZ ** GRS bit ',i2,' (bit range = 1-24) of set ',i1,
     >     ' appears stuck. Ave bit value = ',f4.2)
 151  format(' FZ ** GRS bit ',i2,' (bit range=1-24) of set ',i1,
     >     ' is non-zero (it should never be) = ',f4.2)
 152  format(' FZ ** GRS bit ',i2,' of set ',i1,
     >     ' is not 1 ',i9,' times: ',a)
 153  format(' FZ ** GRS bits 21-24 of set 3 are non-zero ',
     >     i7,' times.  Call Purdue in the morning.')
 154  format(' FZ ** Quicklook GRS - DACQ GRS >10 usec ',
     >     i7,' times.  Check bit conversions and MJD.')
 155  format(' FZ ** Cannot estimate 10 MHz oscillator true freq. ',
     >     'due to elapsed osc. problems.')
 156  format(' FZ ** Estd. osc. freq. is ',e14.8,
     >     ' Hz - very different from expected.')
 157  format(' FZ ** Ave. time difference is ',e13.7,
     >     ' - check GRS clock and elapsed time osc.')
 158  format('  FZ ** Using elapsed oscillator frequency of ',e14.8,
     >     ' Hz')

 160  format(' FZ ** Tracking deviations are large: mean dev = ',
     >     f6.3,' Sigm dev = ',f6.3,'.  Check the telescope tracking.')
 161  format(' FZ ** No tracking records were detected.  ',
     >     'Is there a problem there?')

 170  format(' Minute:',i3,' GPS [sec]:',f10.3,' Rate:',i5,
     >     ' Freq [MHz]:',f11.8,'->calib.')
 171  format(' Minute:',i3,' GPS [sec]:',f10.3,' Rate:',i5,
     >     ' Freq [MHz]:',f11.8,'<-nominal')
 172  format(' Minute:',i3,' GPS [sec]:',f10.3,' Rate:',i5,
     >     ' Freq [MHz]:',f11.8)

 180  format('  FZ ** ',a,' bits set ',i6,' times. This should',
     >     ' never happen.')

 185  format('  Elapsed time warning occurred ',i6,' times.')
 186  format('  Live time warning occurred    ',i6,' times.')

 190  format(' FZ ** PMT ',i3,' HV reads ',f7.1,' but should be ',
     >     f7.1,'.  Check this PMT.')
 191  format(' FZ ** Another HV record ',f4.2,' minutes into run.',
     >     ' Something changed during the run.')
 192  format(' FZ ** Status bit error for PMT ',i3,': ',a)
 193  format(' FZ ** PMT ',i3,' is disabled but HV = ',f7.1,
     >     ' which is too low to do anything.')
 194  format(' FZ ** PMT ',i3,' is disabled but HV = ',f7.1)
 196  format(' FZ ** WARNING: No HV records in this run!')
 197  format(1X,40(I3,1X))
 198  format(' FZ ** Status bit error for PMT ',i3,': ',a,a)
 199  format(1x,10(i5,1x))

 200  format(' FZ ** No singles rate scaler records this run.')
 201  format(1X,'Highest Rates',5X,' Lowest Rates')
 202  format(1X,'PMT     Rate',6X,'PMT     Rate')
 203  format(1X,i3,1X,f9.2,5X,i3,1X,f9.2)

 861  format(a,e14.7)
 862  format(a,i4)
 863  format(a,i12)
 864  format(a,f15.6)

      end
*
**************************************************************************


**************************************************************************
*
      subroutine chkerr(ierr,msg)
*
* ** This subroutine checks whether the argument ierr is 0, and if it
* ** isn't stops program execution and prints the message msg.
*
***************************************************************************
* MAC 970911
      use gdf
*
      implicit none

      integer*4 ierr,lnblnk
      character msg*45
*
***************************************************************************
*
      if (ierr.ne.0) then
         write(6,'(a,i2,a)')'Exiting...error ',ierr,' condition in ',
     &                                     msg(1:lnblnk(msg))
         call gdf$exit(ierr)
         call exit(1)
      end if

      return
      end
*
****************************************************************************

****************************************************************************
*
c
c--routine to return command line args as integers
c--user supplies:      npos = argument number (counting from 1)
c--routine returns:    inum = arg. as an integer
c--	mfc 22 oct 86
c
c	subroutine intarg(npos,inum)
c	character*100 cbuf
c	character*1 blank
c	dimension int(100)
c
c--invoke system subroutine to accept arg. in character format
c
c	call getarg(npos,cbuf)
c
c--determine extent of arg. by finding first blank
c
c	blank=' '
c	idex=index(cbuf,blank)
c	ilen=idex-1
c
c--convert from characters to integers by using internal read
c
c	inum=0
c	do 10 i=1,ilen
c	read(cbuf(i:i),101)int(i)
c101	format(i1)
c	inum=inum+int(i)*10**(ilen-i)
c10	continue
c
c	return
c	end
**************************************************************************
*
      SUBROUTINE GPSTIME(GPS,GYEAR,MJD,UTC,TIMESTR)

*
*     decode GPS time and return MJD in modified Julian days
*     Modified 941102 by JB.  Based on code by CA, MS, JR
*
      IMPLICIT NONE
*
*     modified Julian days for begining of each Gregorian year
*
      INTEGER      YEAR_MIN

      INTEGER      YEAR_MAX
      PARAMETER   (YEAR_MIN=1985)
      PARAMETER   (YEAR_MAX=1997)
      REAL*8       MJD_YEAR(YEAR_MIN:YEAR_MAX)

      DATA MJD_YEAR(1985) / 46066D0 /
      DATA MJD_YEAR(1986) / 46431D0 /
      DATA MJD_YEAR(1987) / 46796D0 /
      DATA MJD_YEAR(1988) / 47161D0 /
      DATA MJD_YEAR(1989) / 47527D0 /
      DATA MJD_YEAR(1990) / 47892D0 /
      DATA MJD_YEAR(1991) / 48257D0 /
      DATA MJD_YEAR(1992) / 48622D0 /
      DATA MJD_YEAR(1993) / 48988D0 /
      DATA MJD_YEAR(1994) / 49353D0 /
*
* Mod. by JB 950526 - Added dates from astronomical almanac.  Note
* that AA gives dates as the mjd of the 0th day of the year, where
* the above values are given for the 1st day of the year.
*    
      DATA MJD_YEAR(1995) / 49718D0 /
      DATA MJD_YEAR(1996) / 50083D0 /
      DATA MJD_YEAR(1997) / 50449D0 /      !JQ 970113, mjd in almanac + 1
   

C---- arguments
      INTEGER*2  GPS(3)   ! bcd encoded GPS time

      INTEGER    GYEAR    ! Gregorian year
      REAL*8     MJD      ! returned UTC time [mjd]
      INTEGER*4	 UTDATE	  ! returned date [yymmdd]
      REAL*8	 UTC	  ! returned UT time in seconds
      CHARACTER  TIMESTR*80	! returned time string [hh:mm:ss.sss]
      INTEGER*4  STATUS   ! returned status bits
      REAL*8	 DSEC

C---- misc
      INTEGER D1,D2,D3                          ! decimal digits
      INTEGER YEAR,DAY,HOUR,MIN,SEC,SEC03,SEC06 ! decoded time


c---- day of the year
      d1   = iand(ishft(gps(1), -6),15)
      d2   = iand(ishft(gps(1),-10),15)
      d3   = iand(ishft(gps(1),-14), 3)
      DAY  = d3*100 + d2*10 + d1

c---- hour of the day
      d1   = iand(      gps(1),    15)
      d2   = iand(ishft(gps(1),-4), 3)
      HOUR = d2*10 + d1

c---- minutes
      d1   = iand(ishft(gps(2), -9),15)
      d2   = iand(ishft(gps(2),-13), 7)
      MIN  =  d2*10 + d1

c---- seconds
      d1  = iand(ishft(gps(2),-2),15)
      d2  = iand(ishft(gps(2),-6), 7)
      SEC = d2*10 + d1

c---- msecs
      d1 = iand(ishft(gps(3), -6),15)
      d2 = iand(ishft(gps(3),-10),15)
      d3 = iand(ishft(gps(3),-14), 3)
      d3 = ior(d3,ishft(iand(gps(2),3),2))         ! this bcd digit spans words
      SEC03 = d3*100 + d2*10 + d1

c---- error code
      STATUS = iand(ishft(gps(3),-2),15)

c---- quarter msecs
      SEC06  = iand(gps(3),3) * 250

C---- check Gregorian year
      IF (GYEAR.LT.100) THEN
        YEAR = GYEAR + 1900
      ELSE
        YEAR = GYEAR
      END IF
*     IF (YEAR.LT.YEAR_MIN.OR.YEAR.GT.YEAR_MAX) THEN
*        WRITE(*,*) 'FZ ** Illegal Gregorian year: ',GYEAR
*        RETURN
*     ENDIF

      DSEC = float(SEC)+float(SEC03)/1.0d3+float(SEC06)/1.0d6
      UTC = float(HOUR)*3600.0d0+float(MIN)*60.0d0+
     & float(SEC)+float(SEC03)/1.0d3+float(SEC06)/1.0d6
      write(timestr,200)YEAR,DAY,HOUR,MIN,DSEC
C----
      IF (YEAR.GE.YEAR_MIN.AND.YEAR.LE.YEAR_MAX) THEN
        MJD =   MJD_YEAR(YEAR)
     .      + DFLOAT(DAY) - 1D0               ! -1D0 because 1. January
DAY=1
     .      + DFLOAT(HOUR ) /    24D0
     .      + DFLOAT(MIN  ) /  1440D0
     .      + DFLOAT(SEC  ) / 86400D0
     .      + DFLOAT(SEC03) / 86400D3
     .      + DFLOAT(SEC06) / 86400D6

      ENDIF
200   format(i4,x,i3,x,i2,':',i2,':',f9.6)

      RETURN
      END
**************************************************************************
*
      SUBROUTINE TRUETIME(GPS,GYEAR,MJD,UTC,TIMESTR)

*
*     decode GPS time from True Time clock and return MJD in modified 
*     Julian days
*     Modified 941102 by JB.  Based on code by CA, MS, JR
*     Modified 9709116 by MAC.  Converted from GPSTIME subroutine to
*     decode the Wisconsin clocks which have replaced the Michigan clocks.
*
      IMPLICIT NONE
*
*     modified Julian days for begining of each Gregorian year
*
      INTEGER      YEAR_MIN

      INTEGER      YEAR_MAX
      PARAMETER   (YEAR_MIN=1985)
      PARAMETER   (YEAR_MAX=2004)
      REAL*8       MJD_YEAR(YEAR_MIN:YEAR_MAX)

      DATA MJD_YEAR(1985) / 46066D0 /
      DATA MJD_YEAR(1986) / 46431D0 /
      DATA MJD_YEAR(1987) / 46796D0 /
      DATA MJD_YEAR(1988) / 47161D0 /
      DATA MJD_YEAR(1989) / 47527D0 /
      DATA MJD_YEAR(1990) / 47892D0 /
      DATA MJD_YEAR(1991) / 48257D0 /
      DATA MJD_YEAR(1992) / 48622D0 /
      DATA MJD_YEAR(1993) / 48988D0 /
      DATA MJD_YEAR(1994) / 49353D0 /
*
* Mod. by JB 950526 - Added dates from astronomical almanac.  Note
* that AA gives dates as the mjd of the 0th day of the year, where
* the above values are given for the 1st day of the year.
*    
      DATA MJD_YEAR(1995) / 49718D0 /
      DATA MJD_YEAR(1996) / 50083D0 /
      DATA MJD_YEAR(1997) / 50449D0 /      !JQ 970113, mjd in almanac + 1
      DATA MJD_YEAR(1998) / 50814D0 /
      DATA MJD_YEAR(1999) / 51179D0 /
      DATA MJD_YEAR(2000) / 51544D0 /
      DATA MJD_YEAR(2001) / 51910D0 /
      DATA MJD_YEAR(2002) / 52275D0 /
      DATA MJD_YEAR(2003) / 52640D0 /
      DATA MJD_YEAR(2004) / 53005D0 /
   

C---- arguments
      INTEGER*4  GPS(3)   ! bcd encoded GPS time

      INTEGER    GYEAR    ! Gregorian year
      REAL*8     MJD      ! returned UTC time [mjd]
      INTEGER*4	 UTDATE	  ! returned date [yymmdd]
      REAL*8	 UTC	  ! returned UT time in seconds
      CHARACTER  TIMESTR*80	! returned time string [hh:mm:ss.sss]
      INTEGER*4  STATUS   ! returned status bits
      REAL*8	 DSEC

C---- misc
      INTEGER D1,D2,D3                          ! decimal digits
      INTEGER YEAR,DAY,HOUR,MIN,SEC             ! decoded time

C* Data variables
      real*8 oscfreq      !In Hz
      data oscfreq /9.99999227d6/
*000103      data oscfreq /9.999992d6/
*old      data oscfreq /1.0d7/


c---- day of the year
      d1   = iand(gps(3),15)
      d2   = iand(ishft(gps(3),-4),15)
      d3   = iand(ishft(gps(3),-8), 3)
      DAY  = d3*100 + d2*10 + d1

c---- hour of the day
      d1   = iand(ishft(gps(2),-16),15)
      d2   = iand(ishft(gps(2),-20),15)
      HOUR = d2*10 + d1

c---- minutes
      d1   = iand(ishft(gps(2), -8),15)
      d2   = iand(ishft(gps(2),-12),15)
      MIN  =  d2*10 + d1

c---- seconds
      d1  = iand(gps(2),15)
      d2  = iand(ishft(gps(2),-4),15)
      SEC = d2*10 + d1

c---- fraction of seconds (down to 100 nsec)
      DSEC = dfloat(gps(1))/oscfreq

c---- error code
      STATUS = iand(ishft(gps(3),-16),15)

C---- check Gregorian year
      IF (GYEAR.LT.100) THEN
        YEAR = GYEAR + 1900
      ELSE
        YEAR = GYEAR
      END IF
*     IF (YEAR.LT.YEAR_MIN.OR.YEAR.GT.YEAR_MAX) THEN
*        WRITE(*,*) 'FZ ** Illegal Gregorian year: ',GYEAR
*        RETURN
*     ENDIF

      DSEC = dfloat(SEC)+DSEC
      UTC = dfloat(HOUR)*3600.0d0+dfloat(MIN)*60.0d0+DSEC
      write(timestr,200)YEAR,DAY,HOUR,MIN,DSEC
C----
      IF (YEAR.GE.YEAR_MIN.AND.YEAR.LE.YEAR_MAX) THEN
        MJD =   MJD_YEAR(YEAR)
     .      + DFLOAT(DAY) - 1D0       ! -1D0 because 1. January DAY=1
     .      + DFLOAT(HOUR ) /    24D0
     .      + DFLOAT(MIN  ) /  1440D0
     .      + DSEC          / 86400D0

      ENDIF
200   format(i4,x,i3,x,i2,':',i2,':',f10.7)

      RETURN
      END
*
* SWAPTUBES
* 941104 JB
*
* Fix errors in tube mapping throughout the years.
*
	subroutine swaptubes(utdate,sshort)
        
	implicit none

        integer*4 max_chan           ! JQ 961205
        parameter (max_chan=156)

        integer*2 sshort(max_chan)
        integer*2 temp
	integer*4 utdate

	if((utdate.gt.940800).and.(utdate.lt.941024)) then
	  temp = sshort(6)
          sshort(6) = sshort(9)
          sshort(9) = temp
        end if

	return

	end

*----------------------------------------------------------------------
* Subroutine TIMEFMT
* 950615
* JB
*
* Convert a time in fractional hours after midnight to a string in the
* format hh:mm:ss.s
* 
	subroutine timefmt(time,tstr,itime)
	implicit none
	real time
	real*8 dtime
	real*8 rh,rm,rs
	integer hh,mm
	character*20 tstr
        integer*4 itime
        
	dtime = time
	rh = aint(dtime)
        rm = aint((dtime-rh)*60.0d0)
        rs = (dtime-rh-rm/60.0d0)*3600.0d0
        itime = nint(rh*100.0+rm+1.0)
	hh = nint(rh+0.1) ! Avoid possible roundoff/type conversion	
	mm = nint(rm+0.1) ! problems by adding a bit before truncating
	write(tstr,'(i2,'':'',i2,'':''f4.1)')hh,mm,rs
  
        return
        end
*-----------------------------------------------------------------------
* Subroutine MJD2SID
* 950615
* JB
* 
* Convert fractional Julian day (mjd) to sidereal time in fractional
* hours after midnight
*

        subroutine mjd2sid(mjd,stime)

        implicit none

        real uttime
        real tt,tt0,stime
        real*8 mjd,dmjd

        dmjd = aint(mjd) 
        uttime = (mjd-dmjd)*24.0
        tt=(dmjd-51544.5)/36525.0      
        tt0=6.697374558d0+(2.400051336d3*tt)+(2.5862d-5*tt*tt)
        stime=mod(uttime*1.002737909+tt0-7.392333333,24.0)
        if(stime.lt.0.0)stime=stime+24.0
        return
        end

*-----------------------------------------------------------------------
*  Subroutine DEROT
*  950615
*  MP,JB
*
* Based on "position" and UT to ST conversion in Duffett-Smith
*
* rra		RA in radians
* rdec		DEC in radians
* mjd		
* uttime        decimal UT (hr. fraction after UT midnight)
* theta         angle by which FoV rotates
*               i.e. angle by which to de-rotate.
* alt           returned value of telescope alt (radians)
* azi           returned value of telescope azimuth (radians)
* stime         returned value of sidereal time
*
* Reference: R.M.Green "Spherical Astronomy" Cambridge U.P. 1985
*
*   See page 28 for alt-az to ra-dec conversion derived from one
*   application of the cosine law for the included angle A, one
*   application of the cosine law for the included angle H, and
*   the sine law: sinA/cosdelta = sinH/sinz to derive the arctan
*   relation needed to remove the sign ambiguity.
*
*   To derive the derotation angle, again refer to figure 2.3 and
*   again do the trick of applying the sine and cosine laws to
*   elim. the sign ambiguity: i.e.
*   sinphi = sindelta cosz + cosdelta sinz costheta
*   cosphi/sintheta = cosdelta/sinA
*   where theta is the angle PXZ between the decl. dir. and the alt.
*   direction
*

        subroutine derot(rra,rdec,mjd,uttime,theta,alt,azi,stime)

        implicit none

        real uttime,theta
        real rra,rdec,str,hr,alt,azi
        real tt,tt0,stime
        real arg, arga, argb
        real sk,ck
        real*8 mjd,dmjd
        real convd,pi,rlat,titor
        data convd/57.29577951/,pi/3.141592654/,rlat/.553065751136/,
     +  titor/.261799387/

        dmjd = aint(mjd)        ! Added by JB to truncate possible
                                ! fractional mjd
*
* convert decimal UT to decimal ST
*       tt=(mjd-51544.5)/36525.0
        tt=(dmjd-51544.5)/36525.0       ! JB 950411
        tt0=6.697374558d0+(2.400051336d3*tt)+(2.5862d-5*tt*tt)
        stime=mod(uttime*1.002737909+tt0-7.392333333,24.0)
        if(stime.lt.0.0)stime=stime+24.0
*       write(*,*)'Decimal ST',stime

*
* convert ST to radians
*
        str=stime*titor

*
* hr is the hour angle in radians
*
        hr=str-rra

*
* Now get the elevation
*
        arg = sin(rdec)*sin(rlat)+cos(rdec)*cos(rlat)*cos(hr)
        if(arg.ge.1.0) then
          arg = 1.0
        else if(arg.le.-1.0) then
          arg = -1.0
        endif
        alt=asin(arg)
*       alt=asin(sin(rdec)*sin(rlat)+cos(rdec)*cos(rlat)*cos(hr))
*
*       Compare this with CA's code (azelrot.f):
*       elev = PI/2.0 - acos(s_dec*s_lat + c_dec*c_lat*c_ang)
* Now get the azimuth
* Using atan2 will put azimuth in correct quadrant
*
        arga = cos(rdec)*cos(rlat)*sin(hr)*-1.0
        argb = (sin(rdec)-sin(rlat)*sin(alt))
        if(arga.ge.1.0) then
          arga = 1.0
        else if(arga.le.-1.0) then
          arga = -1.0
        endif
        if(argb.ge.1.0) then
          argb = 1.0
        else if(argb.le.-1.0) then
          argb = -1.0
        endif
        azi = atan2(arga, argb)
*
*       azi=atan2(-cos(rdec)*cos(rlat)*sin(hr),
*    +   (sin(rdec)-sin(rlat)*sin(alt)))
*
*       Compare this with CA's code (azelrot.f):
*       azimuth = atan2(-c_dec*s_ang, s_dec*c_lat-c_ang*c_dec*s_lat)+2PI
*       azimuth = dmod(azimuth,2PI)
*
* Alt is elevation, Azi is azimuth, both in radians.
*
        sk=cos(rlat)*sin(azi)*-1.0
        ck=cos(alt)*sin(rlat)-sin(alt)*cos(rlat)*cos(azi)
        if(sk.ge.1.0) then
          sk = 1.0
        else if(sk.le.-1.0) then
          sk = -1.0
        endif
        if(ck.ge.1.0) then
          ck = 1.0
        else if(ck.le.-1.0) then
          ck = -1.0
        endif
        theta=atan2(sk,ck)

        return
        end

**********************************************************************
        SUBROUTINE MODJULDAT(IC,MONTH,DAY,YEAR,DATE)
*
* Based on JULDAT, modified by JB to convert between MODIFIED
* Julian date and date, rather than Julian date.  Also modified
* to compile under DEC FORTRAN on the ALPHA.
*
* CONVERTS CALENDAR DATE TO MODIFIED JULIAN DATE  (WHEN IC .LE. 0)
* CONVERTS MODIFIED JULIAN DATE TO CALENDAR DATE  (WHEN IC .GT. 0)
*
      IMPLICIT NONE

      REAL*8 XMONTH,XDAY,XYEAR,DATE
      INTEGER DAY,YEAR,MONTH,MONTHD
      INTEGER I,IC,ND
*
* SEE FUNCTION YCOFF FOR COMMON DATA
*
      dimension   MONTHD (12)
      DATA MONTHD/31,28,31,30,31,30,31,31,30,31,30,31/

      XMONTH = MONTH
      XDAY = DAY
      XYEAR = YEAR
      IF( IC )1, 1, 6
*
* DATE TO MODIFIED JULIAN
*
1     ND = int(XDAY + 365.250D0 * XYEAR)
*1     ND = DAY + 365.250D0 * YEAR
      IF ((YEAR/4*4 - YEAR).NE.0)  GO TO 4
2     IF (MONTH.GT.2)  GO TO 4
3     ND = ND - 1
4     IF(MONTH.EQ.1)GO TO 40
      DO 5 I=2,MONTH
5     ND=ND+MONTHD (I-1)
40    DATE=1721044.50D0+dfloat(ND)
*
* Modified by JB 950620 - convert back fro JD to MJD since this
* is the expected output.  Note that there is an error in this
* conversion factor
*
      date = date -2400000.5
      RETURN
*
* MODIFIED JULIAN TO DATE
*
* Modified by JB 950620 - convert from MJD to JD since this is
* the expected input value in date
*
6     date = date+2400000.5
      ND = int(DATE - 1721044.50D0)
      YEAR = int(dfloat(ND) / 365.250D0)
*
* Modified by JB 950620 - IFIX causes link error
*
*     ND = ND - IFIX( 365.250D0 * YEAR )
      ND = ND - INT( 365.250D0 * dfloat(YEAR))
      IF ((YEAR/4*4 - YEAR).NE.0)  GO TO 10
7     IF( ND - 59 ) 9, 8, 10
8     MONTH = 2
      DAY = 29
      RETURN
9     ND = ND + 1
10    DO 11 I=1,12
      IF (ND.LE.MONTHD(I))  GO TO 12
11    ND=ND-MONTHD (I)
      I=12
12    MONTH = I
      DAY = ND
      RETURN
      END
*************************************************************************







































