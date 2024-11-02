*---------------------------------------------------------------------
*   GPARAMDAT_331.F
*   961104
*   JB, JQ, ADK, GV, MFC, MP, CA, etc.
*
*   NOTE: This is identical to gparamdat.f BUT it uses 331 PMTs in the
*     determination of max1, max2, max3 to match the full field trigger.
*
*   Use:
*     gparamdat_331 fname utdate n2id [options [hhmmss.s ddmmss.s] [xoff yoff]]
*
*       where options is a sequence of nonblank characters with the
*       following effect:
*       n  --  Produce a paw ntuple file of parameters
*       a  --  Produce an old ascii format parameter file
*       m  --  Perform a 2-dimensional mesh analysis
*       f  --  Perform a 2-dimensional V.Fomin like analysis (N/A)
*       v  --  Perform a 2-d analy. with supercuts + vc,rcl dist cut
*       u  --  Calculate parameters needed for a uv analysis
*       i  --  Run the program in interactive mode; higz event display
*       s  --  Use asymmetry parameter (costs some run time)
*       l  --  Calculate muon ring parameters and pe/dc (slow!)
*       x  --  Calculate position of most likely point of origin
*       p  --  Show percent complete
*       e  --  No image derotation
*       c  --  Enter object coordinates on command line following the
*              options string in the format ra: hhmmss.s, and 
*              dec ddmmss.s
*       o  --  Enter xoff,yoff offset for analysis
*
*   Example:
*       gparamdat_331 gt1234 940901 0232 nmuxs
*
*   If the nitrogen id is set to zero sky gains for the given UT date are
*   used instead of nitrogen gains.
* 
*   Modified by ADK 931129 dat.config, is read to determine paths of
*   calibration and data files and to determine picture and boundary
*   thresholds.
*  
*   Modified by JB 940421 to handle camera layout for 11m 
*   and to generate an ntuple for use with PAW
*
*   Modified by JB 941101 to include a muon parameter for picking
*   out muon rings, number of PMTs in the picture, and 2-fold and 3-fold
*   neighbor flags to the parameterized data ntuple.
*
*   Modified by JB 941208 to include Michael Punch's asymmetry parameter
*   for use in 2-d analysis.
*  
*   Modified by JB 950408 to read livetime, as well as absolute gps time
*   and oscillator time put in the data by fz2red2.
*
*   Modified by JB 950410 to include hadronicity parameter for UV
*   analysis
*
*   Modified by JB 950412 to perform image derotation and 2-d 
*   analysis
*
*   Modified by JB 960108 to allow command line specification of
*   RA and DEC for image derotation. 
*
*   Modified by JB 960109 to have 2-d cuts dependent on telescope
*   id (i.e. 10m or 11m)
*
*   Modified by JB 961206 all parameters are generated using derotated
*   camera coordinates.  Added command line argument to
*   specify an x and y offset for the analysis.  These are camera
*   coordinates aligned to the ra and dec axes.  Note that the .hdr
*   file now no longer contains the correct RA and DEC!
*
*   Modified by JB 961209 maximum number of tubes is now a dynamic
*   variable.  Databases include number of tubes.  No problems until
*   number of tubes exceeds 1000.
*
*   Modified by JQ 970506 to utilise new routines in C library which
*   allow any number of tubes to be switched off.
*
*   Modified by MAC 970922 to reduce the amount of time spent in 
*   the gimclean subroutine.  We only check neighbors of tubes listed
*   as being in the picture from the first loop through.
*
*   Modified by MAC 980130 to use 331 tubes in the calculation of 
*   max1, max2, and max3.
*
*------------------------------------------------------------------------

      program gparamdat

      implicit none

*
* Definitions
*
      integer*4 npmt,mpmt,ppmt
      common /num/ npmt,mpmt,ppmt
      integer*4 maxpmts
      integer*4 chi
      parameter(maxpmts=1000)
      integer*4 ndbg		! Line number at which you have trapped
                                ! an error - print detailed info after
                                ! this line.  This is a last resort if
                                ! dbx or dbg fail
      parameter(ndbg=13733)
      real*4	XSPACING,YSPACING  ! Spacing betw. pmts for int conv
      parameter(XSPACING=0.130)
      parameter(YSPACING=0.224)

      character filename*6,newname*50,gainlabel*2,runid*4,mode*3,
     &          source*20,skyq*6,comms*404,mdev*20,pads*1,trig*6,lfch*1
      character srcid*2
      character*8 runname
      character*11 infdir
      character*11 eminfdir
      character*9  oldinfdir
      data infdir/'/usr/dbgt/'/
      data eminfdir/'/usr/dbge/'/
      data oldinfdir /'/usr/db/'/

      integer	imode		! 0 - off, 1 - on
      integer   isky
      character*80 cdate,cn2id
      character*80 rachar,decchar
      character*80 xoffchar,yoffchar
      real*4 xoff,yoff
      integer*4	ievt
      integer*4 utdate,n2id,calibrators(10),nevents,stdur,date,ut,st,
     &          softcode,code,nunit,dunit,prunit,infunit,i,j,k,l,
     &          pmts(maxpmts),munit,tubesoff(maxpmts),ptubesoff(maxpmts)
      integer*4 nevtd100 	! number of events divided by 100
      integer*4 remevt
      integer*4 noff_on,noff_off,noff_pr
      integer*4 stubesoff(10)	! Short list of tubes off for consistency
				! with old .pdat format.
* 
* Note tubesoff must not be an adjustable array since it is used in a common
* block.  So replaced parameter array length with explicit definition.
* One more reason this should be rewritten in C.
*
      integer*4 size,eventsin,eventsout
      integer*4 vetoflag,outflag
      integer*2 short(maxpmts)
      integer*4 gpsbeg(7),iver,nver
*
* Restored by JB 940503
*
      character*4 pdid
      real*4 peds(maxpmts),pedvar(maxpmts),
     &       gains(maxpmts),gainvar(maxpmts),
     &       event(maxpmts),revent(maxpmts),
     &       ra,dec,azimuth,elevation,length,width,miss,dist,azwidth,
     &       sinalpha,frac2,frac3,hadron,arcness,asymmetry,phi,
     &       muon,xgamma,ygamma,azwidthp,missp,distp,asymp
      real*4 pictpmt(maxpmts)
      real*4 pictx(maxpmts),picty(maxpmts)
      integer pictnum(maxpmts)
      real*4 arc,gain,mugain,soal,muskew,arclen,xave,yave,rave
      integer npts
      integer nring
      data nring /0/
      real*4 xtrip(200),ytrip(200)
      real*4 smiss,sdist,sazwidth,ssinalpha,sxmean,symean
      real*4 sphi,sasymmetry
      integer*4 irevent(maxpmts),ievent(maxpmts) ! Integer events for viewing
      real*4 wxdeg(maxpmts),wydeg(maxpmts),wradius(maxpmts)
      real*4 grevent(maxpmts)
      common /wpict/ grevent 
      real*4 max1,max2,max3,tmp	! Maximum three signals after gimclean
      integer*4 loc1,loc2,loc3,nb2,nb3,npicture
      real*4 gx1,gx2,gx3,gy1,gy2,gy3
      real*4 maxsep,sep
      real*4 nbthresh
      real*4 rmax1,rmax2,rmax3	! Maximum three signals in raw data
      integer*4 rloc1,rloc2,rloc3
      real xmean, ymean

      real*8 duration,mjd,frjd,osctime,gpsutc,livetime,time,first7
      real*8 prduration,minduration
      real*4 ltime
      real*8 firsttime		! Time of first event - for use with
				! old data format
      real*8 oldtime		! Time of last event for calculating
				! deltat
      integer*4 bits(0:31),nbr(0:23,maxpmts)

      logical exist,long,start,gpsyes
      logical fnb_trig		! Set this to true to calculate the 
                                ! neighbor trigger condition.
*      data fnb_trig /.false./
      data fnb_trig /.true./
      logical rotcoords
      data rotcoords /.true./

      integer*4 yy
      character yyc*2, hy*1
      integer*4	dbg		! Debug pointer in case we have a
				! bad crash and the debugger gets lost
      logical valid

      character*80	hdron,hdroff
      real*8 sid_start_off,sid_start_on,diff_sid_start,trk_ra_offset
      character*80	dummystr,astr,bstr

*
* C library functions
*
      integer*4 get_coords
      integer*4 read_peds
      integer*4 read_n2gains
      integer*4 read_toff
      integer*4 status
      integer*4 ped_npmt, ped_mpmt, ped_eventcnt
      character*2	ped_mode 
      character*80	ped_file
      character*80	hrc_name
      character*80	n2_file
      character*80	ntoff_file
      integer*4		n2_pedid, n2_mpmt 
      character*4	n2_code
      integer*4 nbr_list(7,maxpmts) ! new camera has 7 neighbours in places
      common /neigh/ nbr_list
*
* TEMPORARY
*
      integer*4 tnbr_list(7*maxpmts)
********************************************************************
*
* Hard-coded cut values
*
      real*4 lengthlo,lengthhi
      real*4 widthlo,widthhi
      real*4 distlo,disthi
      real*4 sinalphahi
      real*4 sizelo
      real*4 triglo
      real*4 lwidhi,llenhi
      real*4 azwidhi
      real*4 ltimehi
      common / cuts /
     +  lengthlo,lengthhi,
     +  widthlo,widthhi,
     +  distlo,disthi,
     +  sinalphahi,
     +  sizelo,
     +  triglo,
     +  lwidhi,llenhi,
     +  azwidhi,
     +  ltimehi

      data lengthlo /0.16/, lengthhi /0.45/	! supercuts
      data widthlo /0.073/, widthhi /0.160/	! supercuts
      data distlo /0.51/, disthi /2.4/		! supercuts
      data sinalphahi /0.174/			! supercuts
      data sizelo /0.0/
      data triglo /45.0/
*     data lwidhi /0.16/		! approx. Vacanti et al. 1991
      data lwidhi /0.26/                ! Vacanti et al. 1991, zone 1
*     data llenhi /0.35/		! approx. Vacanti et al. 1991
      data llenhi /0.50/		! Fomin et al. 1994
      data azwidhi /0.17/		! approx. Vacanti et al. 1991
      data ltimehi /1620.0/		! Approx. 27min duration

      real*4 mumax1hi
      data mumax1hi /200.0/		! Maximum value of Max1 for a
                                        ! muon ring
      integer*4 munpiclo		! Minumum number of pixels for
*     data munpiclo /5/			! a good muon ring
      data munpiclo /12/                ! a good muon ring
*
* Variables from gcutdat.f needed to parse dat.config
*
      character pthdummy*80
      integer*4 codedo(12),maxdo,trigmult
      real*8 defmaxdur
      logical	wait,veto
      character waitstat*3,vetostat*3
*********************************************************************
*
* Counters for tabulating events passing cuts
*
      integer NALPBINS
      parameter (NALPBINS=18)
      real*4 DALPHA
      parameter (DALPHA=5.0)
      integer*4 alphahist(NALPBINS)
      integer*4 npraw
      integer*4 nptrig
      integer*4 npshape
      integer*4 nporient
      integer*4 npboth
      common /npass/ ievt,npraw,nptrig,npshape,nporient,npboth
      common /ahist/ alphahist
      data npraw /0/, nptrig /0/, npshape /0/, nporient /0/,
     &  npboth /0/
*********************************************************************
*
* Variables added by JB for 2-d analysis
*
      real	hh,mm,ss,dd,era,arcmm,arcss
      real	rra,rdec
      real	racmd,deccmd		! RA and DEC from command line
      real	radeg, decdeg		! RA and DEC in degrees
      real	utdec
      real	firstutdec, lastutdec
      real	alt			! Telescope elevation returned
                                        ! by derot() subroutine
      real	az			! Telescope azimuth returned by
                                        ! derot subroutine
      real	stime			! Decimal sidereal time returned
                                        ! by derot() subroutine
      character*12   hfname
      integer   NXGRID,NYGRID
      parameter (NXGRID=30, NYGRID=30)
*331  parameter (NXGRID=48, NYGRID=48)
*91      parameter (NXGRID=30, NYGRID=30)
*     parameter (NXGRID=21, NYGRID=21)
      integer	ALPHAID
      parameter (ALPHAID=51)
      integer   MESHID			! HBOOK histogram id for mesh
      parameter (MESHID=11)
      real*4    xmesh(NXGRID),ymesh(NYGRID)	! x,y coords of mesh
      real*4    nmesh(NXGRID,NYGRID)		! Number of counts on grid
      integer   tmesh(NXGRID,NYGRID)		! Temporary number of counts
      character chmesh(NXGRID,NYGRID)
      integer   ngrid
      common /meshv/ xmesh,ymesh,nmesh,tmesh,chmesh
      real*4    frac_count
      integer   ix,iy
      real      theta, ctheta, stheta, mstheta, mctheta
      real*4    xt,yt
      real*4    rwxdeg(maxpmts),rwydeg(maxpmts)	! Rotated camera coords
      common /sinfo/rra,rdec,mjd
      common /rwcoords/rwxdeg,rwydeg
      real*8	convd, pi, rlat
      data convd/57.29577951/,pi/3.141592654/,rlat/.553065751136/

*********************************************************************
*
* Variables added by JB for Ntuple 
*
*     character parmfmt*20		! Format of parameterized output
*					! files:
      character telsrc*2
      character options*20		! Command line options
					! n -- NTUPLE FILE
					! a -- ASCI FILE (old
					!      parameterize format)
                                        ! 2 -- 2-d analysis
                                        ! u -- UV analysis
                                        ! i -- Interactive with evt disp

      integer   noptions
      character optchar
      logical	fntuple
      logical   fasci
      data	fntuple /.FALSE./, fasci /.TRUE./

      integer N_TELE
      parameter (N_TELE = 2)
      integer N_PARAM
      parameter (N_PARAM = 38)
      integer N_TIME
      parameter (N_TIME = 4)

      integer*4 tele_id
      character*3 tele_name(N_TELE)
      common / tele /
     +  tele_id,
     +  tele_name

      logical      gr_param_valid(N_TELE)          ! data valid
      real gr_param_value(N_PARAM,N_TELE)          ! actual value
      character*8  gr_param_name (N_PARAM)         ! parameter names
      
      logical      gr_time_valid(N_TIME)
      real*8 gr_time_value(N_TIME,N_TELE) 
      character*8  gr_time_name(N_TIME)

      common / gr_param /
     +    gr_param_valid,
     +    gr_param_value,
     +    gr_param_name,
     +    gr_time_value,
     +    gr_time_name

      integer        nt_unit	! Ntuple file logical descripter 
      integer        nt_ierr
      character*8    nt_dir	! Ntuple directory
      character*10   nt_topdir
      character*80   nt_fname	! Ntuple file name

      common / nt /
     +     nt_unit,
     +     nt_ierr,
     +     nt_dir,
     +     nt_topdir,
     +     nt_fname

      integer	NT_ID
      parameter (NT_ID = 100)

      integer	PAW_NWORDS
      parameter (PAW_NWORDS = 2*1024*256)	! 2 MB
      integer	h(PAW_NWORDS)
      common/ PAWC /
     +     h

      integer
     +     TELE_10,
     +     TELE_11
*
      parameter
     +    (TELE_10 = 1,
     +     TELE_11 = 2)
*
      integer
     +     PI_AZWIDTH,
*    +     PI_CONCENT,
     +     PI_WIDTH,
     +     PI_LENGTH,
     +     PI_DISTANCE,
     +     PI_SINALPHA,
     +     PI_MISS,
*    +     PI_RPEAK,
     +     PI_XCENTR,
     +     PI_YCENTR,
     +     PI_PHI,
*    +     PI_AXINTER,
     +     PI_TOTSIG,
     +     PI_MAX1,
     +     PI_LOC1,
     +     PI_MAX2,
     +     PI_LOC2,
     +     PI_MAX3,
     +     PI_LOC3
*
      integer
     +     PI_RMAX1,
     +     PI_RLOC1,
     +     PI_X1,
     +     PI_Y1,
     +     PI_RMAX2,
     +     PI_RLOC2,
     +     PI_X2,
     +     PI_Y2,
     +     PI_HADRON,
     +     PI_NB2,
     +     PI_NB3,
     +     PI_NPICTURE,
     +     PI_ASYMM,
     +     PI_VETO,
     +     PI_EL,
     +     PI_FRAC2,
     +     PI_MAXSEP,
     +     PI_XGAM,
     +     PI_YGAM,
     +     PI_GAIN,
     +     PI_ARC,
     +     PI_MUON,
     +     PI_UTC,
     +     PI_DELTAT,
     +     PI_LIVETIME,
     +     PI_PHASE

      parameter
     +    (PI_AZWIDTH = 1,
*    +     PI_CONCENT = 2,
     +     PI_WIDTH = 2,
     +     PI_LENGTH = 3,
     +     PI_DISTANCE = 4,
     +     PI_SINALPHA = 5,
     +     PI_MISS = 6,
*    +     PI_RPEAK = 8,
     +     PI_XCENTR = 7,
     +     PI_YCENTR = 8,
     +     PI_PHI = 9,
*    +     PI_AXINTER = 12,
     +     PI_TOTSIG = 10,
     +     PI_MAX1 = 11,
     +     PI_LOC1 = 12,
     +     PI_MAX2 = 13,
     +     PI_LOC2 = 14,
     +     PI_MAX3 = 15,
     +     PI_LOC3 = 16)

      parameter
     +    (PI_RMAX1 = 17,
     +     PI_RLOC1 = 18,
     +     PI_X1 = 19,
     +     PI_Y1 = 20,
     +     PI_RMAX2 = 21,
     +     PI_RLOC2 = 22,
     +     PI_X2 = 23,
     +     PI_Y2 = 24,
     +     PI_HADRON = 25,
     +     PI_NB2 = 26,
     +     PI_NB3 = 27,
     +     PI_NPICTURE = 28,
     +     PI_ASYMM = 29,
     +     PI_VETO = 30,
     +     PI_EL = 31,
     +     PI_FRAC2 = 32,
     +     PI_MAXSEP = 33,
     +     PI_XGAM = 34,
     +     PI_YGAM = 35,
     +     PI_GAIN = 36,
     +     PI_ARC = 37,
     +     PI_MUON = 38,
     +     PI_UTC = 1,
     +     PI_DELTAT = 2,
     +     PI_LIVETIME = 3,
     +     PI_PHASE = 4)

      data gr_param_name(PI_AZWIDTH) / 'azwidth'  /
*     data gr_param_name(PI_CONCENT) / 'concent'  /
      data gr_param_name(PI_WIDTH) / 'width'    /
      data gr_param_name(PI_LENGTH) / 'length'   /
      data gr_param_name(PI_DISTANCE) / 'distance' /
      data gr_param_name(PI_SINALPHA) / 'sinalpha' /
      data gr_param_name(PI_MISS) / 'miss'     /
*     data gr_param_name(PI_RPEAK) / 'rpeak'    /
      data gr_param_name(PI_XCENTR) / 'xcentr'   /
      data gr_param_name(PI_YCENTR) / 'ycentr'   /
      data gr_param_name(PI_PHI) / 'phi'  /
*     data gr_param_name(PI_AXINTR) / 'axinter'  /
      data gr_param_name(PI_TOTSIG) / 'totsig'   /
      data gr_param_name(PI_MAX1) / 'max1' /
      data gr_param_name(PI_LOC1) / 'loc1' /
      data gr_param_name(PI_MAX2) / 'max2' /
      data gr_param_name(PI_LOC2) / 'loc2' /
      data gr_param_name(PI_MAX3) / 'max3' /
      data gr_param_name(PI_LOC3) / 'loc3' /
      data gr_param_name(PI_RMAX1) / 'rmax1' /
      data gr_param_name(PI_RLOC1) / 'rloc1' /
      data gr_param_name(PI_X1) / 'x1' /
      data gr_param_name(PI_Y1) / 'y1' /
      data gr_param_name(PI_RMAX2) / 'rmax2' /
      data gr_param_name(PI_RLOC2) / 'rloc2' /
      data gr_param_name(PI_X2) / 'x2' /
      data gr_param_name(PI_Y2) / 'y2' /
      data gr_param_name(PI_HADRON) / 'hadron' /
      data gr_param_name(PI_NB2) / 'nb2' /
      data gr_param_name(PI_NB3) / 'nb3' /
      data gr_param_name(PI_NPICTURE) / 'npicture' /
      data gr_param_name(PI_ASYMM) / 'asymm' /
      data gr_param_name(PI_VETO) / 'veto' /
      data gr_param_name(PI_EL) / 'el' /
      data gr_param_name(PI_FRAC2) / 'frac2' /
      data gr_param_name(PI_MAXSEP) / 'maxsep' /
      data gr_param_name(PI_XGAM) / 'xgam' /
      data gr_param_name(PI_YGAM) / 'ygam' /
      data gr_param_name(PI_GAIN) / 'gain' /
      data gr_param_name(PI_ARC) / 'arc' /
      data gr_param_name(PI_MUON) / 'muon' /

      data gr_time_name(PI_UTC) / 'utc' /
      data gr_time_name(PI_DELTAT) / 'deltat' /
      data gr_time_name(PI_LIVETIME) / 'livetime' /
      data gr_time_name(PI_PHASE) / 'phase' /

      data tele_name(TELE_10) / 't10' /
      data tele_name(TELE_11) / 't11' /

      data nt_dir    /   'GRNTPL' /
      data nt_topdir / '//GRNTPL' /
      data nt_unit /9/
*      data nt_unit   /84/

      logical HEXIST 
*************************************************************************
*
*  This section contains variables added for during creation of
*  paramdat, 9311 by ADK
*
      integer*4 NPATHS    !Number of paths read from dat.config
      parameter(NPATHS=4)

      character cfgfile*10
      parameter(cfgfile='dat.config')

      integer*4 lnblnk   !library function

      logical pad_save
      real*4 prpeds(maxpmts),prpedvar(maxpmts)
      real*4 diffped,sign
      real*4 sig_thresh,nbr_thresh

      character name1*6,name2*6,prid*4,padstat*3,prname*6
*
* Mod. by JB 940503 pthpr*80 -> pthpr*90
*
      character path(NPATHS)*80, pthdat*80, pthcal*80, pthntub*80,
     &          pthout*80,pthpr*90,line*80
      equivalence (path(1),pthdat), (path(2),pthcal),
     &            (path(3),pthntub), (path(4),pthout)
      logical pad,eof


      logical fview
      logical fmesh
      logical ffomin
      logical fsuperd
      logical f2d
      logical fxy
      logical fuv
      logical fasymm
      logical floose
      logical fmuon
      common /flags/ fview,fmesh,ffomin,fsuperd,f2d,fxy,fuv,fasymm,
     &               floose,fmuon
      logical fpct
      logical fderot
      logical fcoords
      logical foffset

      character cdummy*10
      logical pass_cuts
      logical pass_shape
      logical pass_fscuts
      logical pass_superd
      logical pass_mscuts
      logical pass_loose
      logical pass_mesh
      logical fpassed
      character passmsg*200
      logical nextpass
      logical nextloose
      logical nextevt
      logical nextmu
      logical showraw
      logical pairexist
      data nextpass / .false. /
      data nextloose / .false. /
      data nextevt / .true. /
      data nextmu / .true. /
      data showraw / .false. /
      data pairexist / .true. /
      character cevtnum*5
      integer kwtype
*
**************************************************************************
*
      common /wcoords/wxdeg,wydeg
      common /info/peds,pedvar,prpeds,prpedvar,gains,tubesoff,nbr
      common /bytes/bits
*
* Initializations
      data softcode/2/,nunit/1/,dunit/2/,munit/3/,prunit/22/
      data infunit/13/
      data mdev/'parameterize.msg'/,pads/' '/
      data bits/'00000001'X,'00000002'X,'00000004'X,'00000008'X,
     +          '00000010'X,'00000020'X,'00000040'X,'00000080'X,
     +          '00000100'X,'00000200'X,'00000400'X,'00000800'X,
     +          '00001000'X,'00002000'X,'00004000'X,'00008000'X,
     +          '00010000'X,'00020000'X,'00040000'X,'00080000'X,
     +          '00100000'X,'00200000'X,'00400000'X,'00800000'X,
     +          '01000000'X,'02000000'X,'04000000'X,'08000000'X,
     +          '10000000'X,'20000000'X,'40000000'X,'80000000'X/
      lfch=char(10)

*
* JB 940510
* Initialize PAW Memory
*
      call HLIMIT(PAW_NWORDS)

*
* Initialize variable trig to " " as it is needed to keep the first 10
* records 80 characters wide (UNIX ignores null characters .fixup JR 900703)
*
      trig=' '
      start=.true.
*
*  Before doing anything else, check for the configuration file.  Quit
*  if it's not there.
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
* Parameterize section of cfgfiledat.config found.
*
* Lines beginning with a space or a # ar comment lines.  Skip these,
* and read the paths for data files, calibration files, the 
* hrc.ntubelist file, and the output files.  The line after the
* directory paths contains the picture and boundary threshold
* values.  The next line determines whether software padding is
* on or off; if it is on, then one more path is read to determine
* the location of the file "pairs".
*
* Loop until the proper number of non-comment lines are read
*

      i = 1
      do while(i.le.NPATHS)
	 call getline(10,line,eof)
         path(i) = line
*	 write(*,866)'paths: ',i,path(i)
	 i = i + 1
      end do
      
*
* After the paths are read get the picture and boundary thresholds.
*
      call getline(10,line,eof)
      read(line,*)sig_thresh,nbr_thresh
*     read(line,*,err=720)sig_thresh,nbr_thresh
*
* Check if we want to do software padding.  If we do, then
* we need one last path to use to find the file "pairs".
*
* Modified by JB 941029 -- try to locate pairs file even if no padding.
* this is still needed to ensure that the same tubes are turned off in
* the on and off runs for on/off analysis.
*

      call getline(10,line,eof)
      if ((line(1:1).eq.'y') .or. (line(1:1).eq.'Y')) then
         pad = .true.
         padstat = ' ON'
      else
         pad = .false.
         padstat = 'OFF'
      end if

*     if (pad) then
      call getline(10,line,eof)
      pthpr = line
*     write(*,861)'Padding path ',pthpr
*
* Mod. by JB 940503 - added dummy path to elim problem with appending
* / beyond array bounds which occurs when compiling with new f77
*
* Mod. by JB 941029 - pairs path must be specified even if no padding
*
*     else
*       pthpr = '.'
*     end if
*

*     Finally read number of tubes we wish to use in when parameterizing
*     images

      call getline(10,line,eof)
      read(line,*)ppmt

* We have all we need from the configuration file...close it.
*
      close(10)
*---------------------------------------------------------------------
*
* Re-open the configuration file, this time get cut values from the
* cutdat section
*
      open(10,file=cfgfile,status='old',readonly)
*
* The cutdat section of cfgfile is denoted by a line with
* "cut" as the first three characters.  Find this.
*
      line = '   '
      do while (line(1:3) .ne. 'cut')
         call getline(10,line,eof)
         if (eof) then
           write(6,825)cfgfile
           write(6,830)cfgfile
           call exit(1)
         end if
      end do
*
* Skip over lines with paths of input and output files
*
      call getline(10,pthdummy,eof)
      call getline(10,pthdummy,eof)
*
* Get array of codes to include in analysis (not used)
*
      call getline(10,line,eof)
      read(line,*)codedo
      do i = 1,12
         if (codedo(i) .ne. 0) then
            maxdo = i
         end if
      end do
*
* Get the default maximum duration for which to analyze
*
      call getline(10,line,eof)
      read(line,*)defmaxdur
*
* Find out if all data before the first code 7 should be ignored
*
      call getline(10,line,eof)
      if ((line(1:1).eq.'y') .or. (line(1:1).eq.'Y')) then
         wait = .true.
         waitstat = ' ON'
      else
         wait = .false.
         waitstat = 'OFF'
      end if
*
* Find out if we should kill events triggering the veto.
*
      call getline(10,line,eof)
      if ((line(1:1).eq.'y') .or. (line(1:1).eq.'Y')) then
         veto = .true.
         vetostat = ' ON'
      else
         veto = .false.
         vetostat = 'OFF'
      end if
*
* Read width, length, dist, alpha, trigger, and size cut values
*
      call getline(10,line,eof)
      read(line,*)widthlo,widthhi

      call getline(10,line,eof)
      read(line,*)lengthlo,lengthhi

      call getline(10,line,eof)
      read(line,*)distlo,disthi

      call getline(10,line,eof)
      read(line,*)sinalphahi

      call getline(10,line,eof)
      read(line,*)trigmult,triglo	! Temp. disable trig multiplicity
      if(trigmult.ne.2) then
        write(*,*) ' PD ** Trigger multiplicity can not be changed.'
        write(*,*) ' PD ** Set to default multiplicity: 2'
      endif

      call getline(10,line,eof)
      read(line,*)sizelo
*
* Done with configuration file.  Close it, and stop execution if
* we have encountered an end of file along the way
*
      close(10)
*
*------------------------------------------------------------------
*
* We haven't been checking to see if we got new data when needed; do
* so now.  All data should have been read without generating an
* end of file.
*
      if (eof) then
        write(*,840)cfgfile
        write(*,830)cfgfile
        call exit(1)
      end if
*
* Getting here means we read enough lines
* Insert a trailing slash in path names
*
      do i = 1,NPATHS
	 j = lnblnk(path(i))
	 if (path(i)(j:j) .ne. '/') path(i)(j+1:j+1) = '/'
      end do
*     write(*,*)'pthpr: ',pthpr
      j = lnblnk(pthpr)
      if (pthpr(j:j) .ne. '/') pthpr(j+1:j+1) = '/'
*
* Done with the parameterize part of cfgfile...resume with
* the "traditional" part of parameterization.  In using path names,
* we use fortran library function LNBLNK which returns the index of
* the last non-blank character of its (string) argument.
*
*-- Get parameters from command line
*
* JB - intarg() does not work with DEC Fortran 3.1, the problem is
* currently unknown.  Use getarg instead
*
      call getarg(1,filename)
      call getarg(2,cdate)
      call getarg(3,cn2id)
      call getarg(4,options)
      read(cdate,*)utdate
      read(cn2id,*)n2id
      fmuon = .false.
      fview = .false.
      fmesh = .false.
      ffomin = .false.
      fsuperd = .false.
      f2d = .false.
      fxy = .false.
      fuv = .false.
      fntuple = .false.
      fasci = .false.
      fasymm = .false.
      fpct = .false.
      fderot = .true.
      fcoords = .false.
      foffset = .false.
*
* Parse "options" command line parameter
*
      noptions = lnblnk(options)
      do i=1,noptions
        optchar = options(i:i)
        if((optchar.eq.'n').or.(optchar.eq.'N')) then
          fntuple = .true.
        else if((optchar.eq.'a').or.(optchar.eq.'A')) then
          fasci = .true.
        else if((optchar.eq.'f').or.(optchar.eq.'F')) then
          ffomin = .true.
        else if((optchar.eq.'v').or.(optchar.eq.'V')) then
          fsuperd = .true.
        else if((optchar.eq.'m').or.(optchar.eq.'M')) then
          fmesh = .true.
        else if((optchar.eq.'u').or.(optchar.eq.'U')) then
          fuv = .true.
        else if((optchar.eq.'i').or.(optchar.eq.'I')) then
          fview = .true.
        else if((optchar.eq.'s').or.(optchar.eq.'S')) then
          fasymm = .true.
        else if((optchar.eq.'p').or.(optchar.eq.'P')) then
          fpct = .true.
        else if((optchar.eq.'l').or.(optchar.eq.'L')) then
          fmuon = .true.
        else if((optchar.eq.'x').or.(optchar.eq.'X')) then
          fxy = .true.
        else if((optchar.eq.'e').or.(optchar.eq.'E')) then
          fderot = .false.
        else if((optchar.eq.'c').or.(optchar.eq.'C')) then
          fcoords = .true.
        else if((optchar.eq.'o').or.(optchar.eq.'O')) then
          foffset = .true.
        end if
      end do
      if(ffomin.or.fmesh.or.fsuperd) then
        f2d = .true.
      else
        f2d = .false.
      endif
*
* Added 951031 to allow specification of ra and dec on command line
* if the tracking records were missing.
*
      xoff = 0.0 
      yoff = 0.0
      if(fcoords) then
        call getarg(5,rachar)
        read(rachar,*)racmd
        call getarg(6,decchar)
        read(decchar,*)deccmd
	if(foffset) then
	  call getarg(7,xoffchar) 
          read(xoffchar,*)xoff
	  call getarg(8,yoffchar)
          read(yoffchar,*)yoff
        endif
      else
        if(foffset) then
          call getarg(5,xoffchar)
          read(xoffchar,*)xoff
          call getarg(6,yoffchar)
	  read(yoffchar,*)yoff
        endif
      endif
*MAC 970918 What is this?
*      if(foffset) then
*      endif
       
*     if(cview.eq.'view') then
*       fview = .true.
*     endif
*     fmesh = .false.
*     if(cmesh.eq.'mesh') then
*       fmesh = .true.
*     endif
*
* Open gt????.inf file in append mode to record important information
* about this data file.  gt????.inf was probably previously created
* by fz2red2
*
      open(infunit,file=filename//'.inf',
     &  access='append',status='unknown')
      write(infunit,490)
      write(infunit,*)' '
*
* Initialize hplot for disp_event
*
      if(fview) then
        call IGINIT(0)  ! Initialize HIGZ
        call IGWKTY(kwtype)	! Get workstation type
*       call IGSSE(6,kwtype)	! Errors writen to unit 6, set ws type	
*       call HPLINT(1)  ! Inititalize HPLOT
        call HPLINT(kwtype)  ! Inititalize HPLOT
        call HPLCAP(0)  ! No metafile
      end if

      write(*,490)
*PM 000804: Fixed fmt 890 to show leading zeros in dates after 991231
      write(*,890)filename,utdate,n2id
      write(infunit,890)filename,utdate,n2id
      write(*,867)sig_thresh,nbr_thresh
      write(infunit,867)sig_thresh,nbr_thresh
      write(6,912)
      write(infunit,912)
      write(6,913)2,triglo,sizelo,distlo,disthi,lengthlo,lengthhi,
     &  widthlo,widthhi,sinalphahi
      write(infunit,913)2,triglo,sizelo,distlo,disthi,lengthlo,lengthhi,
     &  widthlo,widthhi,sinalphahi
      if(pad) then
        write(*,*)' PD ** Padding ON'
        write(*,*)' '
      else
        write(*,*)' PD ** Padding OFF'
        write(*,*)' '
      endif
*
* Get telescope ID from filename
*
      telsrc=filename(1:2)
*
* If none, ask for them
*
      if(filename.eq.' ')then	! TEMPORARY! should be null, but not
                                ! allowed
	 write(6,100)
	 read(5,'(a)')filename
	 write(6,110)
	 read(5,*)utdate
	 write(6,130)
	 read(5,*)n2id
      endif
      pdid=filename(3:6)        ! JB - was 3:7, but this doesn't make
                                ! sense

*
* This one line for Y2K compliancy. What's all the fuss about? SJF 990831
*
      if (utdate.lt.800000) utdate = utdate + 1000000

*
* Check file's existence
*
      inquire(file=pthdat(1:lnblnk(pthdat))//filename,exist=exist)
      if(.not.exist)then
	 write(6,90)pthdat(1:lnblnk(pthdat))//filename
	 stop
      endif

*
* Filename is suposed to be like cr1234
*
      newname=pthout(1:lnblnk(pthout))//filename//'.pdat'
*
* Ntuple filename is like em1234.rz
*
*     nt_fname = pthout(1:lnblnk(pthout))//filename//'.rz'
      nt_fname = filename//'.rz'
*      write(*,861)' PD ** Creating ntuple file: ',nt_fname
      write(6,400)filename

*-- Get grid of coordinates and neighboring tubes
      if(telsrc.eq.'ge') then	! Trick the program into setting
				! the camera type to that which was
				! on the 10m on 930101, i.e. the camera
				! which is now on the 11m.
 	write(*,*)'GP ** 11m data file'
        tele_id = TELE_11
       else
        tele_id = TELE_10
       endif
       call num_chan(tele_id,utdate,mpmt,npmt) 
       if((ppmt.le.0).or.(ppmt.gt.npmt)) then
          ppmt=npmt
       endif
       wxdeg(1) = 1.
       wxdeg(2) = 2.
       wxdeg(3) = 3.
       wxdeg(4) = 4.
       wydeg(1) = 1.
       wydeg(2) = 2.
       wydeg(3) = 3.
       wydeg(4) = 4.
       tnbr_list(1) = 1
       tnbr_list(2) = 2
       status = get_coords(tele_id,utdate,wxdeg,wydeg,wradius,nbr_list)
       if(status.gt.0) then
         write(*,*)' PD ** Failed to generate coordinates, aborting'
	 call exit(1)
       endif
* MAC 970922 Stop using this
*mac       call make_nbrmask(nbr)

*
* Load these coordinates into rotated camera coordinates which are
* actually used by ghillakpar
*
*      do ix=1,mpmt		
* MAC 970922 does not seem like we use this at all
*save       do ix=1,npmt		
*save         rwxdeg(ix) = wxdeg(ix) 
*save         rwydeg(ix) = wydeg(ix)
*save       end do 
*
* ** If we're doing padding, we need to know what file to pair with this
* ** one.  Get that information.
*
* Modified by JB 941029 - need to know what file to pair with this one
* even if we're not doing padding.  This is necessary to make sure same
* tubes are turned off in on and off runs
*
*     if (pad) then
*	 write(*,*)'Doing padding'
         inquire(file=pthpr(1:lnblnk(pthpr))//'hrc.pairlist',
     &                                               exist=exist)
         if (.not.exist) then
            write(6,850)pthpr(1:lnblnk(pthpr))//'hrc.pairlist'
            if (pad) then
              write(*,*) ' PD ** Aborting gparamdat.'
              call exit(1)	! exit and return a bad status
            endif
            imode = 1
         else
           open(1,file=pthpr(1:lnblnk(pthpr))//'hrc.pairlist',
     &      status='unknown')
3001       continue
           read(1,3002,end=3009)name1,name2
3002       format(a6,1x,a6)
           pairexist = .true.
           if (filename.eq.name1) then
              prid=name2(3:6)
              prname = name2(1:6)
              imode = 1
              go to 3008
           end if
           if (filename.eq.name2) then
              prid=name1(3:6)
              prname = name1(1:6)
              imode = 0
              go to 3008
           end if
           go to 3001
3009       continue
           pairexist = .false.
           write(*,*)
     &      ' PD ** WARNING! Cannot find corresponding pair file.'
           write(*,*)' PD **         Tubes turned off may not match.'
           write(infunit,*)
     &      ' PD ** WARNING! Cannot find corresponding pair file.'
           write(infunit,*)
     &      ' PD **          Tubes turned off may not match.'
           if (pad) then
             write(*,*) ' PD ** Aborting gparamdat.'
             write(infunit,*) ' PD ** Aborting gparamdat.'
             call exit(1)
           end if
3008       continue
           close(1)
        end if
*     end if
*
* Now determine which pedestal file should be accessed (divided into
* six months at a time, all of 88 being incorporated)   JR  901107
* Simplified MP920505
*
      yy=int((float(utdate)+0.5)/10000.)
      write(yyc,'(i2.2)')mod(yy,100)
      if (utdate-yy*10000.lt.600)then
	 hy='a'
      else
	 hy='b'
      end if
*
* Modified by JB 950602 
* if pthcal is set to "default" in dat.config then use the default
* values
* 
      if(pthcal(1:7).eq.'default') then
        write(*,*)' PD ** Using default paths for calibration database'
        write(*,*)' '
        if(telsrc.eq.'ge') then	
          pthcal = eminfdir
        else if(telsrc.eq.'gt') then
          pthcal = infdir
        else
          pthcal = oldinfdir
        endif
      endif
      if(pthntub(1:7).eq.'default') then
        if(telsrc.eq.'ge') then	
          pthntub = eminfdir
        else if(telsrc.eq.'gt') then
          pthntub = infdir
        else
          pthntub = oldinfdir
        endif
      endif
*
* Many calibration files; set j equal to last non-blank character 
* calibration file path to save typing
*
      j = lnblnk(pthcal)
*
* Get code pedestals
*
* Mod. by JB 940503, use (injected or code 1,2) cpeds not npeds
*

      inquire(file=pthcal(1:j)//'hrc'//yyc//hy//'.cpeds',
     &        exist=exist)
      if (.not.exist)then
	 write(6,425)filename
	 stop
      end if
*
* 961209 - use new C libraries
*
      ped_file = pthcal(1:j)//'hrc'//yyc//hy//'.cpeds' 

*     write(*,864)'  PD ** Read pedestals for run: ',pdid,
*    +   ' from file: ',pthcal(1:j)//'hrc'//yyc//hy//'.cpeds'
*     write(infunit,864)'  PD ** Read pedestals for run: ',pdid,
*    +   ' from file: ',pthcal(1:j)//'hrc'//yyc//hy//'.cpeds'

      write(*,864)'  PD ** Read pedestals for run: ',pdid,
     +   ' from file: ', ped_file(1:lnblnk(ped_file))
      write(infunit,864)'  PD ** Read pedestals for run: ',pdid,
     +   ' from file: ', ped_file(1:lnblnk(ped_file))

      ped_npmt = npmt
      ped_mpmt = mpmt
C      status = read_peds(pdid,utdate,ped_file,peds,pedvar,ped_eventcnt,
C     &                   ped_npmt,ped_mode,ped_mpmt)  JQ 961227
      status = read_peds(filename,utdate,ped_file,peds,pedvar,
     &                   ped_eventcnt,ped_npmt,ped_mode,ped_mpmt)

      if (status.gt.0) then
        write(*,*)' PD ** Could not read the pedestal file:',
     &   ped_file(1:lnblnk(ped_file))
        call exit(1)
      endif
*pm      write(*,*) ' PD ** Pedestals:'
*pm      write(6,895)(peds(i),i=1,npmt)
*pm      write(*,*) ' PD ** Pedestal variances:'
*pm      write(6,895)(pedvar(i),i=1,npmt)
*
      write(infunit,*)' PD ** Pedestals:'
      write(infunit,895)(peds(i),i=1,npmt)
      write(infunit,*)' PD ** Pedestal variances:'
      write(infunit,895)(pedvar(i),i=1,npmt)
*
* Get peds for pair if doing padding
*
*     write(*,*)'get peds for pair if padding'
      if (pad) then
         ped_file = pthcal(1:j)//'hrc'//yyc//hy//'.cpeds' 
         ped_npmt = npmt
         ped_mpmt = mpmt
C         status = read_peds(prid,utdate,ped_file,prpeds,prpedvar,
C     &           ped_eventcnt,ped_npmt,ped_mode,ped_mpmt) JQ 961227
         status = read_peds(prname,utdate,ped_file,prpeds,prpedvar,
     &           ped_eventcnt,ped_npmt,ped_mode,ped_mpmt)
         if (status.gt.0) then
           write(*,*)' PD ** Could not read the pedestal file:',
     &      ped_file(1:lnblnk(ped_file))
           call exit(1)
         endif
         write(*,861),
     +   '  PD ** Read pedestals from other file in padding pair: ',prid
*pm         write(*,*)' PD ** Pedestals:'
*pm         write(6,895)(prpeds(i),i=1,npmt)
*pm         write(*,*)' PD ** Pedestal variances:'
*pm         write(6,895)(prpedvar(i),i=1,npmt)
*
         write(infunit,861),
     +   '  PD ** Read pedestals from other file in padding pair: ',prid
         write(infunit,*)' PD ** Pedestals:'
         write(infunit,895)(prpeds(i),i=1,npmt)
         write(infunit,*)' PD ** Pedestal variances:'
         write(infunit,895)(prpedvar(i),i=1,npmt)
      end if
*
* Determine threshold for neighbor trigger.  This threshold is 
* a multiple (sig_thresh) of the average of the pedvars for all of
* the pmts
* 
* Modified 950615 by JB to make sure threshold is the same for on
* and off runs by using the maximum of the two pedvars (as in
* software padding)
*
*mac
      if (fnb_trig) then
         nbthresh = 0.0
         do i=1,ppmt
            nbthresh = nbthresh + amax1(pedvar(i),prpedvar(i))
         end do
         nbthresh = nbthresh*sig_thresh/float(ppmt)
         write(*,'(/,''  PD ** Neighbor trigger threshold [d.c.]: '',
     &        f7.3)')nbthresh
         write(infunit,
     >        '(/,''  PD ** Neighbor trigger threshold [d.c.]: '',
     &        f7.3)')nbthresh
      endif
*
* Get sky gains
*
15    if(n2id.eq.0)then
	 write(*,*)' PD ** Get sky gains since no nitrogen id is given'
	 write(infunit,*)' PD ** Get sky gains since no nitrogen id is given'
         open(1,file=pthcal(1:j)//'hrc.skygains',status='unknown')
*        open(1,file=pthcal(1:j)//'hrc.skygains',status='old')
	 call getskygains(1,utdate,calibrators,gains,gainvar,*20)
	 write(*,861)'Skygains filename: ',pthcal(1:j)//'hrc.skygains'
	 write(infunit,861)'Skygains filename: ',pthcal(1:j)//'hrc.skygains'

	 gainlabel='SG'
*
* Not found
*
      else
*
* Get nitrogen gains
*
         write(*,*)' '
	 write(*,*)' PD ** Getting N2 gains from file: ',
     &    pthcal(1:j)//'hrc'//yyc//hy//'.n2gains' 
         write(infunit,*)' '
	 write(infunit,*)' PD ** Getting N2 gains from file: ',
     &    pthcal(1:j)//'hrc'//yyc//hy//'.n2gains' 
	 inquire(file=pthcal(1:j)//'hrc'//yyc//hy//'.n2gains',
     &                                              exist=exist)
	 if (.not.exist)then
	    write(6,450)filename
	    write(infunit,450)filename
	    stop
	 end if
*
* Modified by JB 961208 - Don't create if it doesn't exist.
*
         n2_file = pthcal(1:lnblnk(pthcal))//'hrc'//yyc//hy//'.n2gains'
     &                     //char(0)  ! JQ 961227
         n2_mpmt = mpmt
         status = read_n2gains(n2id,utdate,n2_file,gains,gainvar,
     &                         n2_code,n2_pedid,n2_mpmt)

	 if (status.gt.0) then
	    write(*,*)' PD ** Cannot read nitrogen gains from file:',
     &        n2_file(1:lnblnk(n2_file))
	    call exit(1)
         endif

 	 write(*,861)'  Nitrogen gain filename: ',n2_file(1:lnblnk(n2_file))
         write(*,*)' '
*pm	 write(*,*)' PD ** Gains: '
*pm	 write(6,895)(gains(i),i=1,npmt)
*pm	 write(*,*)' PD ** Gain variances: '
*pm	 write(6,895)(gainvar(i),i=1,npmt)
*
	 write(infunit,*)' PD ** Gains: '
	 write(infunit,895)(gains(i),i=1,npmt)
	 write(infunit,*)' PD ** Gain variances: '
	 write(infunit,895)(gainvar(i),i=1,npmt)

	 gainlabel='N2'
	 calibrators(1)=n2id
      endif

*
* Open I/O files
*
* Given the pair filename, open the file and read the header to 
* insure that the duration of the two files match for cuts within
* gparamdat
*
      if(pairexist) then
          if(name1.eq.name2) then	! Save some time...
            prduration = 1.0d20	! Force minimum duration to be
				! given by the other file 
          else
*
* First see if a file of the format gt????.hdr exists.  This is the
* new ascii header file which is generated by fz2red2 and contains
* all kinds of useful data, in a more easily accessible format than
* the old binary header.  Since this program should still work on
* the old data format, use the binary header information for the
* most important thing - on/off duration matching, and use the ascii
* header for warning about sidereal time mismatches.  
*
            hdron = pthdat(1:lnblnk(pthdat))//filename//'.hdr'
            hdroff = pthdat(1:lnblnk(pthdat))//prname//'.hdr'

            inquire(file=hdron,exist=exist)
            if(exist) then
              inquire(file=hdroff,exist=exist)
              if(exist) then
                open(10,file=hdron,status='old',readonly)
                line='                                        '
                do while (line(1:13) .ne. 'gps_sid_start')
                  line='                                        '
                  call getline(10,line,eof)
                  if(eof) then
                    goto 2006 ! Not found, force break
                  endif
                end do
                read(line(15:30),*)sid_start_on
                close(10)

                open(10,file=hdroff,status='old',readonly)
                line='                                        '
                do while (line(1:13) .ne. 'gps_sid_start')
                  line='                                        '
                  call getline(10,line,eof)
                  if(eof) then
                    goto 2006 ! Not found, force break
                  endif
                end do
                read(line(15:30),*)sid_start_off
                do while (line(1:13) .ne. 'trk_ra_offset')
                  line='                                        '
                  call getline(10,line,eof)
                  if(eof) then
                    goto 2006 ! Not found, force break
                  endif
                end do
                read(line(15:30),*)trk_ra_offset
                close(10)
                 
                diff_sid_start = 60.0*abs(sid_start_off-sid_start_on)

                write(6,
     &           '(''  PD ** Sidereal start time diff. for runs: '',
     &           a6,'' '',a6,'': '',f8.3,
     &           '' min '')')filename,prname,diff_sid_start
                write(6,'(''  PD ** Tracking computer RA-offset :'',
     &            f8.3,'' min '')')trk_ra_offset
*
                write(infunit,
     &           '(''  PD ** Sidereal start time diff. for runs: '',
     &           a6,'' '',a6,'': '',f8.3,
     &           '' min '')')filename,prname,diff_sid_start
                write(infunit,
     &            '(''  PD ** Tracking computer RA-offset :'',
     &            f8.3,'' min '')')trk_ra_offset
*
                if(abs(diff_sid_start-trk_ra_offset).gt.0.1) then
                  write(6,
     &             '(''  PD ** WARNING: Sidereal time mismatch! '')')
                  write(infunit,
     &             '(''  PD ** WARNING: Sidereal time mismatch! '')')
                endif
              endif
            endif

2006        continue

*
* Modified by JB 961208 - Don't create if does not exist
*
*           open(prunit,file=pthdat(1:lnblnk(pthdat))//prname,
*    +        form='unformatted',status='old')
            open(prunit,file=pthdat(1:lnblnk(pthdat))//prname,
     +        form='unformatted',status='unknown')
*
* Read header of pair file, all but prduration overwriten later
*
            read(prunit,err=3010)
     +       runid,nevents,prduration,stdur,mode,source,date,mjd,
     +       frjd,ra,dec,ut,st,azimuth,elevation,skyq,comms,gpsbeg
            if(utdate.lt.940800) then
              call vaxdbl(prduration)
            endif
            close(prunit)
            go to 3011
3010        continue
            rewind(prunit)
            read(prunit)	! No gps, read without gpsbeg
     +       runid,nevents,prduration,stdur,mode,source,date,mjd,
     +       frjd,ra,dec,ut,st,azimuth,elevation,skyq,comms
            if(utdate.lt.940800) then
              call vaxdbl(prduration)
            endif
            close(prunit)
3011        continue
          endif 
      else
          prduration = 1.0d20	! Force minimum duration to be
				! given by the other file 
      endif
*
* Now open the data file to be parameterized
*
* Modified by JB 961208 - Don't create if doesn't exist
*
      open(dunit,file=pthdat(1:lnblnk(pthdat))//filename,
     &  form='unformatted',status='unknown')
*     open(dunit,file=pthdat(1:lnblnk(pthdat))//filename,
*    &  form='unformatted',status='old')
      if(fasci) then
        open(nunit,file=newname(1:lnblnk(newname)),status='unknown')
      end if 
*     write(*,*)'opened files'
*
*   Create and open Ntuple file
*
      if(fntuple) then 
        call gr_ntuple('N')
       end if 
*
* Read in the header record. Output the header records
* Try to read gps time of first event...if error, then don't
* read it     ADK  930525
*
      read(dunit,err=901)
     +   runid,nevents,duration,stdur,mode,source,date,mjd,
     +   frjd,ra,dec,ut,st,azimuth,elevation,skyq,comms,gpsbeg
      gpsyes=.true.
      go to 902

901   continue    !gps time not in header
      rewind(dunit)
      read(dunit)
     +   runid,nevents,duration,stdur,mode,source,date,mjd,
     +   frjd,ra,dec,ut,st,azimuth,elevation,skyq,comms
      gpsyes=.false.

902   continue    !Header read, one way or the other

* Convert floating and double precision from VAX to IEEE format
      if(utdate.lt.940800) then
        call vaxflt(ra)
        call vaxflt(dec)
        call vaxflt(azimuth)
        call vaxflt(elevation)
        call vaxdbl(duration)
        call vaxdbl(mjd)
        call vaxdbl(frjd)
      endif
*
* If the command line option 'c' has been chosen, override the ra
* and dec in the header with the values specified on the command line
*
      if(fcoords) then
        ra = racmd
        dec = deccmd
      endif
      
*
* Calculate the minimum duration for ON or OFF runs
*
      if(duration.le.prduration) then
         minduration = duration
      else
         minduration = prduration
      endif
*
* Now set the high livetime cut value equal to the minimum of
* the nominal value or minduration
*
      if(minduration.lt.ltimehi) ltimehi = minduration
      write(6,'(/,''  PD ** Duration of '',a6,11x,'': '',f9.3,
     & '' seconds'')')filename(1:6),duration 
      if(prduration.lt.10000.0) then
       write(6,'(''  PD ** Duration of pair file '',a6,'' : '',
     & f9.3,'' seconds'')')
     &  prname(1:6),prduration
      endif
      write(6,'(''  PD ** Cutting to '',18x,'': ''f9.3,
     & '' seconds '')')ltimehi
*
      write(infunit,'(/,''  PD ** Duration of '',a6,11x,'': '',
     &  f9.3,'' seconds'')')filename(1:6),duration 
      if(prduration.lt.10000.0) then
       write(infunit,'(''  PD ** Duration of pair file '',a6,'' : '',
     &  f9.3,'' seconds'')')prname(1:6),prduration
      endif
      write(infunit,'(''  PD ** Cutting to '',18x,'': '',f9.3,
     & '' seconds '')')ltimehi
*
* Convert back to numbers again
* 
      if(skyq(1:1).eq.'A') then
        isky = 1 
      else if (skyq(1:1).eq.'B') then
        isky = 2
      else if (skyq(1:1).eq.'C') then
        isky = 3
      endif
*----------------------------------------------------------------
* 
* Calculate some things for image derotation in the 2-d analysis
*
* First convert ra,dec from hhmmss.,ddmmss. to degrees:
*
       hh=aint(ra/10000.)
       mm=aint((ra-hh*10000)/100.)
       ss=ra-hh*10000-mm*100.
       radeg=(hh+mm/60.+ss/3600.)*15.0
 
       dd=aint(dec/10000.)
       arcmm=aint((dec-dd*10000.)/100.)
       arcss=dec-dd*10000.-arcmm*100.
       decdeg=dd+arcmm/60.0+arcss/3600.
*
* Now convert ra,dec to radians
*
        write(*,*)
        write(6,'(''  PD ** Source        :  '',a20)')source
        rra=radeg/convd
        write(6,'(''        RA (hhmmss.s) :  '',f8.1,''   RA (deg) : '',
     &   f8.4)')ra,radeg
        rdec=decdeg/convd
        write(6,'(''        DEC (ddmmss.s):  '',f8.1,''   DEC (deg): '',
     &   f8.4)')dec,decdeg
        write(6,'(''        MJD           : '',f9.3,/)')mjd+frjd
*
        write(infunit,*)
        write(infunit,'(''  PD ** Source        :  '',a20)')source
        rra=radeg/convd
        write(infunit,'(''        RA (hhmmss.s) :  '',f8.1,
     &    ''   RA (deg) : '',f8.4)')ra,radeg
        rdec=decdeg/convd
        write(infunit,'(''        DEC (ddmmss.s):  '',f8.1,
     &    ''   DEC (deg): '',f8.4)')dec,decdeg
        write(infunit,'(''        MJD           :   '',f9.3,/)')mjd+frjd
*
* Initialize arrays for 2-d "mesh" analysis
*
* Note:  If the range is from -1.5 to 1.5 in 30 0.1 deg steps, then
* the first xvalue is -1.45 the mean value of the two bin edges
* -1.5,-1.4.  Think about it for a minute, and you will see
* what I'm getting at.
*
        do ix=1,NXGRID
          xmesh(ix) = -0.05*dfloat(nxgrid) + 0.05 + dble(ix-1)*(0.100)
*         xmesh(ix) = -1.00 + dble(ix-1)*(0.100)
          ymesh(ix) = xmesh(ix)
        end do
*
* Initialize HBOOK stuff
*
        if(.not.fntuple) then	      ! added by JB 950523
          call HCDIR('//PAWC',' ')    
        endif

*	if (HEXIST(ALPHAID)) then     ! previous histogram?
*         call HDELET(ALPHAID)        ! delete it
*       end if
*	call HBOOK1(ALPHAID,source//runid,36,0.0,90.0,0.)
        if(f2d) then
	  if (HEXIST(MESHID)) then     ! previous histogram?
            call HDELET(MESHID)        ! delete it
          end if
          call HBOOK2(MESHID,source,NXGRID,-0.05*float(nxgrid),
     &         0.05*float(nxgrid),
     &     NYGRID,-0.05*float(nygrid),0.05*float(nygrid),0.)
*          call HBOOK2(MESHID,source,NXGRID,-1.5,1.5,
*     &     NYGRID,-1.5,1.5,0.)
*         write(*,*)'Initialized histograms'
*         call HBOOK2(MESHID,source,NXGRID,-1.0,1.0,
*    &     NYGRID,-1.0,1.0,0.)
        endif
        hfname = filename//'.h2d'
*        
*
*----------------------------------------------------------------

*
* Tubes off section re-written to use C libraries - JQ 970506
*
*
      runname = filename(1:2)//runid//char(0)
      
      j=lnblnk(pthntub)
      ntoff_file = pthntub(1:j)//'hrc'//yyc//hy//'.ntubelist'//char(0) 

      status = read_toff(runname,noff_on,tubesoff,ntoff_file)
      noff_pr=noff_on


*
* Modified by JB 941029.  Force same tubes off in on/off pair if pairs
* file exists.  This is the case even if padding is off.
*
*      close(4)
*     write(*,'(''  PD ** Tubes off for run '',a6,'':'')')srcid//runid 
      write(*,'(''  PD ** Tubes off for run '',a6,'':'')')filename
      write(*,'(''        '',15i4)')(tubesoff(i),i=1,noff_on)
*
      print *
*     write(infunit,'(''  PD ** Tubes off for run '',a6,'':'')')
*    +     srcid//runid 
      write(infunit,'(''  PD ** Tubes off for run '',a6,'':'')')
     +     filename
      write(infunit,'(''        '',15i4)')(tubesoff(i),i=1,noff_on)
      if (pairexist) then
         status = read_toff(prname,noff_off,ptubesoff,ntoff_file)
         i = 1
         do while (i.le.noff_off) ! For ea. tube of pair
            j = 1		! Determine whether it is a duplicate
            do while (j.le.noff_pr.and.(tubesoff(j).ne.ptubesoff(i)))
               j = j + 1
            enddo               ! Either found a matching tube and quit 
                                ! or no match found.
            if (j.gt.noff_pr) then
               tubesoff(j) = ptubesoff(i) ! no match so append tube.
               noff_pr=noff_pr+1
            endif
         i = i + 1
         enddo
      endif
      write(*,'(''  PD ** Tubes off for pair run '',a6,'':'')')prname 
      write(*,'(''        '',15i4)')(ptubesoff(i),i=1,noff_off)
      write(*,*)' '
      write(*,*)' PD ** Combined tubes off:'
      write(*,'(''        '',15i4)')(tubesoff(i),i=1,noff_pr)
*
      write(infunit,'(''  PD ** Tubes off for pair run '',
     &  a6,'':'')')prname 
      write(infunit,'(''        '',15i4)')(ptubesoff(i),i=1,noff_off)
      write(infunit,*)' '
      write(infunit,*)' PD ** Combined tubes off:'
      write(infunit,'(''        '',15i4)')(tubesoff(i),i=1,noff_pr)
      do i=1,10		! Generate short list of tubes off for .pdat header
			! Only include 10 tubes for consistency with old fmt
        stubesoff(i) = tubesoff(i)
      end do
*
* softcode identifies the type of formats to be used
* Also put the thresholds into this first record  ADK 931129
* softcode = 1 : first version
* softcode = 2 : second version @ FLWO. 
* Record 9 contains the list of tubes turned off.
*
      if(fasci) then
        write(nunit,190)softcode,sig_thresh,nbr_thresh,padstat,
     +    (pads,i=1,58)
        write(nunit,200)runid,source,mode,(pads,i=1,50)
        write(nunit,210)date,ut,st,(pads,i=1,63)
        write(nunit,220)idint(mjd),frjd,(pads,i=1,60)
        write(nunit,230)elevation,azimuth,(pads,i=1,67)
        write(nunit,240)nevents,duration,(pads,i=1,61)
        write(nunit,250)pdid,gainlabel,calibrators,(pads,i=1,22)
        write(nunit,260)ifix(ra),ifix(dec),trig,skyq,(pads,i=1,50)
        write(nunit,270)stubesoff,(pads,i=1,39)
* 
* One blank record for future info
* Use last record for GPS time of first ev
*
        if (gpsyes) then
	 write(nunit,281)gpsbeg,(pads,i=1,51)
        else
	 write(nunit,280)(pads,i=1,79)
        end if
      end if 				! If writing asci param file

*
* After 881101 (run id 603) use integer*2 to store the events
* Condition on the date takes care of reset run ids
*
      if (date.ge.881101) then
	 long=.false.
      else
	 long=.true.
      endif

*------------
*-- MAIN LOOP
*------------
*     write(*,*)'main loop'
      eventsin = 0
      eventsout = 0
      nevtd100 = nevents/100
*
* Loop over the events
*
      do ievt=1,nevents
         if(fpct) then
           remevt = mod(ievt,nevtd100)
           if(remevt.eq.0) then
             write(*,*)(ievt/nevtd100),'%'
           endif 
         endif
*
* Read events in correct format
*        write(*,*)ievt
         if(utdate.lt.940800) then	! JB 950408
	   if(long)then
              read(dunit,end=778)code,time,(pmts(chi),chi=1,npmt)
	   else
              read(dunit,end=778)code,time,(short(chi),chi=1,npmt)
	      do j=1,mpmt
	         pmts(j)=short(j)
	      end do
	   endif
           eventsin = eventsin + 1
*          
* Convert from VAX to IEEE double precision
*
	   call vaxdbl(time)
	   if(ievt.eq.1) then
             firsttime = time
           endif
           ltime = time - firsttime
           gpsutc = ltime + frjd*(24.0d0*60.0d0*60.0d0)
           utdec = gpsutc/3600.0
         else 
	   if(long)then
             read(dunit,end=778)code,osctime,gpsutc,livetime,
     &                          (pmts(chi),chi=1,npmt)
             utdec = gpsutc/3600.0 ! Convert to dec. hrs. for derot.
             time = osctime	! Prior to 950410 was livetime
             ltime = livetime	! Convert from double to float
	   else
             read(dunit,end=778)code,osctime,gpsutc,livetime,
     &                          (short(chi),chi=1,npmt)
             utdec = gpsutc/3600.0 ! Convert to dec. hrs. for derot.
             time = osctime	! Prior to 950410 was livetime
             ltime = livetime
	     do j=1,mpmt
	       pmts(j)=short(j)
	     end do
	   endif
*
* Modified by JB 951016 - We totally rely on livetime for duration
* matching, but for some dates the Michigan memory modules which record
* the gate_open and gate_close were not working.  For these dates,
* set the livetime to the elapsed real time.
*
           if((utdate.eq.951016).or.(utdate.eq.950308)) then
	     if(ievt.eq.1) then
               firsttime = time
             endif
             ltime = time - firsttime
           endif
           eventsin = eventsin + 1
         end if
*        if (ievt.gt.ndbg) write(*,*)'1' 
*
* Remember the time of the first code 7
*
	 if ((start).and.(code.eq.7)) then
	    start = .false.
	    first7 = time
	 end if
*
*  Special codes:  1 = Short pedestal
*                  2 = Long pedestal
*                  5 = atomic clock on the half minute
*                  6 = UT minute marker (from WWVB clock)
*                  7 = Sidereal minute marker
*  Parameterize all other codes
*
	 if((code.ne.1).and.(code.ne.2).and.(code.ne.5).and.
     &                          (code.ne.6).and.(code.ne.7)) then
*           if (ievt.gt.ndbg) write(*,*)'2' 
	    if(i.eq.1) then 
               oldtime = time
            endif
*
* Set vetoflag if the veto shield fired
* 950530 TCW says veto should be in 118,119,120
* 970922 MAC - no more muon veto with the new focus box
*
*           if (ievt.gt.ndbg) write(*,*)'3' 
	    if(utdate.lt.961201) then
	      if ((pmts(119).gt.50).or.(pmts(120).gt.50)) then
	        vetoflag = 1
	      else
	        vetoflag = 0
	      end if
            else if (utdate.lt.970901) then
	      if ((pmts(152).gt.50).or.(pmts(153).gt.50)) then
	        vetoflag = 1
	      else
	        vetoflag = 0
	      end if
	    endif
*
* MAC 970923 Rather large change in philosophy here.  Because of the
* large increase in PMTs and the desire for more speed in the analysis
* we are going to clean the image first to eliminate the tubes which
* are not in the picture, and subsequently only loop through those
* tubes which are in the picture.  Hopefully this will save a lot of
* unnecessary checks of PMT value, etc.
*
*
* Clean up the image; zero tubes without appreciable signal, subtract
* peds, apply gains.  Do software padding if desired.
*
*           if (ievt.gt.ndbg) write(*,*)'5' 
            npicture = 0
* MAC 970922 Do not pad if this is a tracking run, the random number
* stuff takes a long while.

            if (filename.eq.prname) then
               pad_save=pad
               pad=.false.
            else
               pad_save=pad
            endif
	    call gimclean(pmts,sig_thresh,nbr_thresh,pad,
     &       event,npicture,pictpmt,pictx,picty,pictnum)
            pad=pad_save
            if (fview) then
              do k=1,npmt
                ievent(k) = int(event(k))
                if(ievent(k).lt.0) ievent(k) = 0
              end do
            end if
*
* Image cleaned...all tubes not in picture have been zeroed
* Set a flag if the outer 18 tubes contain any picture
* Obsolete 961209
*
*           if (ievt.gt.ndbg) write(*,*)'6' 
*	    tmp = 0.
*	    do j = 92,109
*	       tmp = tmp + abs(event(j))
*	    end do
*	    if (tmp.gt.1.) then
*	       outflag = 1
*	    else
*	       outflag = 0
*	    end if
*
* Added by JB 940917 to calculate the parameters rmaxsig1, rloc1,...
* the three highest "raw" adc signals, ie. the difference between
* the raw adc values and the pedestal values. 
* Note that for a 2/91 trigger rloc2, rmaxsig2 are much more representative
* of the location and the signal on the triggering phototube than loc2,
* and maxsig2
*
* irevent()   Array of integer PMT values with ped. subtracted
* revent()    Array of float PMT values with ped. subtracted
* grevent()   Array of PMT values with ped. subtraction and gain correction
*
*           if (ievt.gt.ndbg) write(*,*)'4' 
            
            if (fmuon.or.fuv.or.fview) then
               do k=1,npmt
                  irevent(k)=pmts(k)-peds(k)
                  if(irevent(k).lt.0) irevent(k) = 0
                  revent(k) = float(irevent(k))
                  grevent(k) = revent(k)*gains(k) 
               enddo
            else
*mac               do l=1,npicture
*mac                  k=pictnum(l)
               do k=1,npmt
                  irevent(k)=pmts(k)-peds(k)
                  if(irevent(k).lt.0) irevent(k) = 0
                  revent(k) = float(irevent(k))
*                  grevent(k) = revent(k)*gains(k) 
               enddo
            endif
* Start here
	    call nmax(revent,ppmt,rmax1,rloc1,rmax2,rloc2,
     +       rmax3,rloc3)
*
* Zero Hillas parameters
*
	    length=0.0
	    width=0.0
	    miss=0.0
	    dist=0.0
	    azwidth=0.0
            sinalpha=0.0
	    frac2=0.0
	    frac3=0.0
	    size=0
	    max1=0.
	    max2=0.
	    max3=0.
	    loc1=0
	    loc2=0
	    loc3=0
	    hadron=0.
            muon=0.
            nb2=0
            nb3=0
	    asymmetry = 0.
            phi = 0.
            
*           if (ievt.gt.ndbg) write(*,*)'6.5' 
* MAC 970922 Only look if more than 2 pmts pass the picture or boundary
* threshold - eliminates a loop within the program
 	    if(fnb_trig) then
               call nb_trigger(revent,nbthresh,nb2,nb3)
            endif

*
* Calculate Hillas Parameters
*
*
* Run derot at this point to calculate the telescope elevation
* for every event, regardless of whether 2-d analysis is being
* done.   In this version of the program the camera coordinates
* are also transformed at this point.  This can do absolutely no
* harm, and should change nothing in the analysis.  It also allows
* an ra and dec offset on the command line. 
*
            call derot(utdec,theta,alt,az,stime) 
            if (.not.fderot) theta = 0.0
	    if (ievt.eq.1) then
	      write(6,'(''  PD ** xoffset,yoffset (deg) : '',f8.4,
     &         f8.4)')xoff,yoff
*             write(*,*)'elevation',alt*convd
              write(6,'(/,''  PD ** First utdec (dec.hr)  :'',f9.4,
     &         /,''        Sidereal time (hr)    :'',f9.4,
     &         /,''        Derotation angle (deg):'',f9.4,
     &         /,''        Telescope AZ (deg)    :'',f9.4,
     &         /,''        Telescope EL (deg)    :'',f9.4)')
     &         utdec,stime,theta*convd,az*convd,alt*convd
*
	      write(infunit,'(''  PD ** xoffset,yoffset (deg) : '',
     &         f8.4,f8.4)')xoff,yoff
              write(infunit,
     &         '(/,''  PD ** First utdec (dec.hr)  :'',f9.4,
     &         /,''        Sidereal time (hr)    :'',f9.4,
     &         /,''        Derotation angle (deg):'',f9.4,
     &         /,''        Telescope AZ (deg)    :'',f9.4,
     &         /,''        Telescope EL (deg)    :'',f9.4)')
     &         utdec,stime,theta*convd,az*convd,alt*convd
               firstutdec = utdec
            endif

            ctheta = cos(theta)
            stheta = sin(theta)
            mstheta = -1.0*stheta 
            mctheta = -1.0*ctheta
* MAC - to save CPU time, only loop through tubes in the picture for
* rotating the coordinates.
*mac            do ix=1,npmt
            do iy=1,npicture
               if ((pictnum(iy).gt.maxpmts).or.
     >              (pictnum(iy).le.0)) then
                  type*,'pictnum=',pictnum(iy)
               endif
               ix=pictnum(iy)
              xt = wxdeg(ix)
              yt = wydeg(ix)
              rwxdeg(ix) = xt*ctheta+yt*mstheta ! Derotate camera coords
              rwydeg(ix) = xt*stheta+yt*ctheta
            end do 

*           xt=0.0
*           yt=0.0
* 
	    call ghillakpar(event,code,npicture,pictnum,
     &           length,width,miss,
     &           dist,azwidth,frac2,frac3,size,loc1,loc2,loc3,
     &           max1,max2,max3,xmean,ymean,hadron,asymmetry,
     &           muon,phi,valid,xoff,yoff,theta,xgamma,ygamma,
     &           azwidthp,missp,distp,asymp)
*
* Calculate maximum separation between two pixels for UV or
* Visible filter analysis
*
	      if(loc1.ge.1.and.loc1.le.npmt) then
                gx1 = wxdeg(loc1)
                gy1 = wydeg(loc1)
	      else
                gx1 = 1000.0
                gy1 = 1000.0
	      end if
	      if(loc2.ge.1.and.loc2.le.npmt) then
                gx2 = wxdeg(loc2)
                gy2 = wydeg(loc2)
	      else
                gx2 = 1000.0
                gy2 = 1000.0
	      end if
	      if(loc3.ge.1.and.loc3.le.npmt) then
                gx3 = wxdeg(loc3)
                gy3 = wydeg(loc3)
	      else
                gx3 = 1000.0
                gy3 = 1000.0
	      end if
              maxsep = 0.0
              sep = sqrt((gx1-gx2)*(gx1-gx2)+(gy1-gy2)*(gy1-gy2))
              if(sep.gt.maxsep) maxsep = sep
              sep = sqrt((gx1-gx3)*(gx1-gx3)+(gy1-gy3)*(gy1-gy3))
              if(sep.gt.maxsep) maxsep = sep
              sep = sqrt((gx3-gx2)*(gx3-gx2)+(gy3-gy2)*(gy3-gy2))
              if(sep.gt.maxsep) maxsep = sep
*
* Load ntuple variables with the appropriate parameter values
* and write one entry in the ntuple file
*
*           if (ievt.gt.ndbg) write(*,*)'7.1' 
	    if(dist.gt.1.0d-6) then
	          sinalpha = (miss/dist)
	    else
	          sinalpha = 1.0
	    endif
*
* Calculate the muon arc parameter and the gain (pe/dc ratio)
*
            if(fmuon.or.(fuv.and.(muon.lt.0.1).and.
     &       (nb3.eq.1).and.(npicture.gt.munpiclo).and.
     &       (max1.lt.mumax1hi))) then
	     call findarc(event,rave,xave,yave,arc,gain,mugain,soal,
     &        muskew,arclen,npicture,pictpmt,pictx,picty,pictnum,xtrip,
     &        ytrip,npts)
	     nring = nring + 1
            endif
*           if (ievt.gt.ndbg) write(*,*)'7.2' 
*
* Modified by JB 950710 - Moved this code to more logical place,
* before the first call to ghillakpar - shouldn't make any difference
*
*           if (ievt.gt.ndbg) write(*,*)'7.3' 
*           call derot(utdec,theta,alt,az,stime) 
*           if (ievt.eq.1) then
*             write(6,'(/,''  PD ** First utdec (dec.hr)  :'',f9.4,
*    &         /,''        Sidereal time (hr)    :'',f9.4,
*    &         /,''        Derotation angle (deg):'',f9.4,
*    &         /,''        Telescope AZ (deg)    :'',f9.4,
*    &         /,''        Telescope EL (deg)    :'',f9.4)')
*    &         utdec,stime,theta*convd,az*convd,alt*convd
*           endif

*
* If the 2-d analysis is being done, store the values of the parameters
* which change with the changing grid coordinates.  Restore these
* values at the end of the 2-d analysis loop.
*
*           if (ievt.gt.ndbg) write(*,*)'8' 
            if(f2d) then
              smiss = miss
              sdist = dist
              sazwidth = azwidth
              ssinalpha = sinalpha
              sxmean = xmean
              symean = ymean
              sphi = phi
              sasymmetry = asymmetry
            endif
*           if (ievt.gt.ndbg) write(*,*)'9' 
            if(fntuple) then
              gr_param_value(PI_LENGTH,tele_id) = length
              gr_param_value(PI_WIDTH,tele_id) = width
              gr_param_value(PI_MISS,tele_id) = miss
              gr_param_value(PI_DISTANCE,tele_id) = dist
              gr_param_value(PI_AZWIDTH,tele_id) = azwidth
              gr_param_value(PI_TOTSIG,tele_id) = size
	      gr_param_value(PI_SINALPHA,tele_id) = sinalpha
	      gr_param_value(PI_XCENTR,tele_id) = xmean
	      gr_param_value(PI_YCENTR,tele_id) = ymean
              gr_param_value(PI_PHI,tele_id) = phi
              gr_param_value(PI_MAX1,tele_id) = max1
              gr_param_value(PI_LOC1,tele_id) = loc1
              gr_param_value(PI_MAX2,tele_id) = max2
              gr_param_value(PI_LOC2,tele_id) = loc2
              gr_param_value(PI_MAX3,tele_id) = max3
              gr_param_value(PI_LOC3,tele_id) = loc3
              gr_param_value(PI_RMAX1,tele_id) = rmax1
              gr_param_value(PI_RLOC1,tele_id) = rloc1
	      if(rloc1.ge.1.and.rloc1.le.npmt) then
                gr_param_value(PI_X1,tele_id) = wxdeg(rloc1)
                gr_param_value(PI_Y1,tele_id) = wydeg(rloc1)
	      else
                gr_param_value(PI_X1,tele_id) = 1000.0
                gr_param_value(PI_Y1,tele_id) = 1000.0
	      end if
              gr_param_value(PI_RMAX2,tele_id) = rmax2
              gr_param_value(PI_RLOC2,tele_id) = rloc2
	      if(rloc2.ge.1.and.rloc2.le.npmt) then
                gr_param_value(PI_X2,tele_id) = wxdeg(rloc2)
                gr_param_value(PI_Y2,tele_id) = wydeg(rloc2)
	      else
                gr_param_value(PI_X2,tele_id) = 1000.0
                gr_param_value(PI_Y2,tele_id) = 1000.0
	      end if
              gr_param_value(PI_HADRON,tele_id) = hadron
              gr_param_value(PI_NB2,tele_id) = nb2
              gr_param_value(PI_NB3,tele_id) = nb3
              gr_param_value(PI_NPICTURE,tele_id) = npicture
              gr_param_value(PI_ASYMM,tele_id) = asymmetry
              if (utdate.lt.961201) then
                gr_param_value(PI_VETO,tele_id) = pmts(118)+pmts(119)
     &            + pmts(120)
              else
                gr_param_value(PI_VETO,tele_id) = pmts(152)+pmts(153)
     &            + pmts(154)
              endif
              gr_param_value(PI_EL,tele_id) = alt*convd
              gr_param_value(PI_FRAC2,tele_id) = frac2
              gr_param_value(PI_MAXSEP,tele_id) = maxsep
              gr_param_value(PI_XGAM,tele_id) = xgamma
              gr_param_value(PI_YGAM,tele_id) = ygamma
              gr_param_value(PI_GAIN,tele_id) = gain
              gr_param_value(PI_ARC,tele_id) = arc
              gr_param_value(PI_MUON,tele_id) = muon
              gr_param_valid(tele_id) = valid 
*             gr_time_value(PI_UTC,tele_id) = time
              gr_time_value(PI_UTC,tele_id) = gpsutc
              gr_time_value(PI_PHASE,tele_id) = osctime
*             gr_time_value(PI_LIVETIME,tele_id) = livetime
              gr_time_value(PI_LIVETIME,tele_id) = ltime
	      gr_time_value(PI_DELTAT,tele_id) = (time-oldtime)
	      oldtime=time
              gr_time_valid(tele_id) = .TRUE.     ! TEMPORARY
            end if	! end if ntuple
*
* Do 2-d mesh analysis
*
*           if (ievt.gt.ndbg) write(*,*)'10' 
*
            if(f2d) then
               floose = pass_loose(length,width,dist,size,max2,ltime)
*              if (ievt.gt.ndbg) write(*,*)'11' 
               if(floose) then
	        ngrid = 0
*               write(*,*),'azwidth:',azwidth
*
* Modified 961206 - coordinate transformation has already been done
*
*               ctheta = cos(theta)		
*               stheta = sin(theta)
*               mstheta = -1.0*stheta 
*               mctheta = -1.0*ctheta
*               if (ievt.gt.ndbg) write(*,*)'12' 
*               do ix=1,npmt
*                 xt = wxdeg(ix)
*                 yt = wydeg(ix)
*                 rwxdeg(ix) = xt*ctheta+yt*mstheta ! Derotate camera coords
*                 rwydeg(ix) = xt*stheta+yt*ctheta
*               end do 
*               if (ievt.gt.ndbg) write(*,*)'13' 
	        do ix=1,NXGRID
                 do iy=1,NYGRID
                  xt = xmesh(ix) + xoff
                  yt = ymesh(iy) + yoff
                  if(fmesh) then   ! if azwidth-mesh, save some time
*
* Modified by JB - bug in oazwidth - temporarily use rehillakpar
*
	           call rehillakpar(event,code,npicture,pictnum,
     &                    length,width,miss,
     &              dist,azwidth,xmean,ymean,asymmetry,phi,valid,
     &              xt,yt,rotcoords)
                  else if(ffomin.or.fsuperd) then
	           call rehillakpar(event,code,npicture,pictnum,
     &                    length,width,miss,
     &              dist,azwidth,xmean,ymean,asymmetry,phi,valid,
     &              xt,yt,rotcoords)
                  end if
                  if(fmesh) then
		   if(pass_mesh(length,width,miss,dist,azwidth,
     &                          ltime)) then
                    tmesh(ix,iy) = 1
	            ngrid = ngrid+1
                    if(fasymm) then
                      if(asymmetry.gt.0.0) then
                        chmesh(ix,iy) = '+'
                        tmesh(ix,iy) = 2
                      else
                        chmesh(ix,iy) = '-'
                      endif
                    endif
                   else
                    tmesh(ix,iy) = 0
                    if(fasymm) then
                      chmesh(ix,iy) = '0'
                    endif
                   end if
                  else if(ffomin) then
		   if(pass_fscuts(length,width,miss,dist,asymmetry,
     &                            ltime)) then
                    tmesh(ix,iy) = 1
	            ngrid = ngrid+1
                    if(fasymm) then
                      if(asymmetry.gt.0.0) then
                        chmesh(ix,iy) = '+'
                        tmesh(ix,iy) = 2
                      else
                        chmesh(ix,iy) = '-'
                      endif
                    endif
                   else
                    tmesh(ix,iy) = 0
                    if(fasymm) then
                      chmesh(ix,iy) = '0'
                    endif
                   endif
                  else if(fsuperd) then
		   if(pass_superd(length,width,miss,dist,asymmetry,
     &                            ltime)) then
                    tmesh(ix,iy) = 1
	            ngrid = ngrid+1
                    if(fasymm) then
                      if(asymmetry.gt.0.0) then
                        chmesh(ix,iy) = '+'
                        tmesh(ix,iy) = 2
                      else
                        chmesh(ix,iy) = '-'
                      endif
                    endif
                   else
                    tmesh(ix,iy) = 0
                    if(fasymm) then
                      chmesh(ix,iy) = '0'
                    endif
                   endif
                  end if
                 end do
                end do
                if(ngrid.gt.0) then
                  frac_count = 1.0d0/dfloat(ngrid)
                  do ix=1,NXGRID
                   do iy=1,NYGRID
                    if(tmesh(ix,iy).gt.0) then
                     if(fasymm) then
                       if(tmesh(ix,iy).gt.1) then
                         nmesh(ix,iy) = nmesh(ix,iy) + 1.0
                         call HFILL(MESHID,xmesh(ix),ymesh(iy),1.0)
                       endif
                     else
                      nmesh(ix,iy) = nmesh(ix,iy) + 1.0
                      call HFILL(MESHID,xmesh(ix),ymesh(iy),1.0)
*                     nmesh(ix,iy) = nmesh(ix,iy) + frac_count
*                     call HFILL(MESHID,xmesh(ix),ymesh(iy),frac_count)
                     endif
                    endif
                   end do
                  end do
                endif
               endif
*
* Restore values of position dependent parameters 
*
               miss = smiss
               dist = sdist
               azwidth = sazwidth
               sinalpha = ssinalpha
               xmean = sxmean
               ymean = symean
               phi = sphi
               asymmetry = sasymmetry
*              if (ievt.gt.ndbg) write(*,*)'14' 
            endif ! end if f2d

*           if (ievt.gt.ndbg) write(*,*)'15' 
            passmsg = 'Passed cuts '	! Moved out of if(fview)
                                        ! 950529 by JB
            fpassed = pass_cuts(length,width,miss,dist,sinalpha,
     &         size,max2,ltime,passmsg)
            if(fview) then
*               if (ievt.gt.ndbg) write(*,*)'16' 
                if((floose.and.nextloose).or.
     &           (fpassed.and.nextpass).or.nextevt.or.
     &           (nextmu.and.(muon.lt.0.1).and.(nb3.eq.1).and.
     &           (npicture.gt.munpiclo).and.(max1.lt.mumax1hi))) then
*                 if (ievt.gt.ndbg) write(*,*)'17' 
903               if(showraw) then
                    call show109(6,irevent)
                    call disp_event(revent,xmean,ymean,length,
     &               width,phi,xgamma,ygamma,nextmu,
     &                xtrip,ytrip,npts,xave,yave,rave,
     &                   npicture,pictnum)
                  else
                    call show109(6,ievent)
                    call disp_event(event,xmean,ymean,length,
     &               width,phi,xgamma,ygamma,nextmu,
     &               xtrip,ytrip,npts,xave,yave,rave,npicture,
     &               pictnum)
                  end if ! show raw
                  write(*,*)'  '
*                 if (ievt.gt.ndbg) write(*,*)'18' 
                  write(6,900)size,length,width,dist,sinalpha
		  if(size.gt.0) then
                    arcness = 1000.0*length/float(size)
                  else
                    arcness = 0.0
                  end if ! size 
*                 if (ievt.gt.ndbg) write(*,*)'19' 
                  write(6,905)arcness,muon,nb2,nb3,npicture
                  write(6,906)azwidth,hadron,phi*convd,asymmetry
                  write(6,908)frac2,maxsep
                  if(fxy) then
                    write(6,907)azwidthp,missp,distp,asymp
                  endif
                  write(*,*)passmsg
* Added by MAC 000114
                  if (max1.gt.0.) then
                     write(6,914) 100.*length/max1
                  else
                     write(6,914) 0.
                  endif
                  write(6,915) 
                  do i=1,npicture
                     write(6,916) pictnum(i),nint(pictpmt(i)),
     >                    pictx(i),picty(i)
                  enddo
*		  if(fmuon) write(*,*)'.'
*                 if((fmuon.or.fuv).and.nextmu) then
*
* Modified by JB 961105
*
*                 if(fmuon) then
*                   write(*,*)' '
*                   write(*,*)'arc =',arc
*                   write(*,*)'gain =',gain
*                   write(*,*)'rave =',rave
*                   write(*,*)'xave,yave =',xave,yave
*                   write(*,*)'mugain =',mugain
*                   write(*,*)'soal = ',soal
*                   write(*,*)'muskew = ',muskew
*                   write(*,*)'arclen = ',arclen
*                 endif
*
* Output the results of the 2-d analysis
*
*                 if (ievt.gt.ndbg) write(*,*)'20' 
                  if(f2d) then
                    write(*,*)' '
                    write(6,'(''UT [hrs]:'',f10.5,'' Rot. angle [deg]'',
     &               f6.2)')utdec,theta*convd
                    write(*,*)' '
*                   if (ievt.gt.ndbg) write(*,*)'21' 
                    if(fasymm) then
                      do iy=NYGRID,1,-1
*                       write(6,'(30a2)')(chmesh(i,iy),i=1,NXGRID)
*                       write(6,'(21a2)')(chmesh(i,iy),i=1,NXGRID)
                      end do
                    else
                      do iy=NYGRID,1,-1
                        write(6,'(30i2)')(tmesh(i,iy),i=1,NXGRID)
*                       write(6,'(21i2)')(tmesh(i,iy),i=1,NXGRID)
                      end do
                    endif ! fasymm
                    write(*,*)' '
                  end if ! f2d
*                 if (ievt.gt.ndbg) write(*,*)'22' 
                  write(*,*)
     & '1-Next, 2-Pass, 3-Muon, 4-2d, 5-Raw/Cleaned, 6-Write, 7-Quit: '
                  read(5,'(a)')cdummy
                  if(cdummy.eq.'1') then
                    nextevt = .true.
                    nextpass = .false.
                    nextmu = .false.
                    nextloose = .false.
                  else if(cdummy.eq.'2') then
                    nextpass = .true.
                    nextevt = .false.
                    nextmu = .false.
                    nextloose = .false.
                  else if(cdummy.eq.'3') then
                    nextpass = .false.
                    nextevt = .false.
                    nextmu = .true.
                    nextloose = .false.
                  else if(cdummy.eq.'4') then
                    nextpass = .false.
                    nextevt = .false.
                    nextmu = .false.
                    nextloose = .true.
                  else if(cdummy.eq.'5') then
                    if(showraw) then
                      showraw = .false.
                      write(*,*)' Toggled to CLEANED image display'
                      write(*,*)' '
                    else
                      showraw = .true.
                      write(*,*)' Toggled to RAW image display'
                    endif
                    nextevt = .true.
                    nextpass = .false.
                    nextmu = .false.
                    goto 903
                  else if(cdummy.eq.'6') then
                    nextevt = .true.
                    nextpass = .false.
                    nextmu = .false.
                    write(cevtnum,'(i5.5)')i
                    open(8,file='evt'//cevtnum//'.dat',status='unknown')
                    write(8,'(f7.2)')(event(j),j=1,npmt)
                    close(8)
                  else if(cdummy.eq.'7') then
                    call IGEND	
                    call exit(0)
                  endif ! cdummy
                endif ! end if show an event
            endif ! end if fview
*
* Now write one row of the ntuple 
*
*           if (ievt.gt.ndbg) write(*,*)'23' 
            if(fntuple) then
              if(valid) then
                call gr_ntuple('W') 
*	        write(*,*)'wrote ntuple entry'
              endif
            end if 			! end if(fntuple)
	 else   !Event was code 1,2,5,6,7; set parameters to 0
*           if (ievt.gt.ndbg) write(*,*)'24' 
	    length=0.0
	    width=0.0
	    miss=0.0
	    dist=0.0
	    azwidth=0.0
            sinalpha=0.0
	    frac2=0.0
	    frac3=0.0
	    size=0
	    max1=0.
	    max2=0.
	    max3=0.
	    loc1=0
	    loc2=0
	    loc3=0
	    vetoflag=0
	    outflag=0
	 end if
*
* Write new record to asci output file
*
         if(valid) then
           if(fasci) then
*            if (ievt.gt.ndbg) write(*,*)'25' 
  	     write(nunit,300,err=58)
     &        code,time,length,width,miss,dist,azwidth,frac2,
     &        frac3,size,loc1,loc2,loc3,max1,max2,max3,vetoflag,outflag
*	   write(*,*)'wrote new .pdat entry'
           end if
           eventsout = eventsout + 1
         end if
	 go to 40
58       type *,' PD ** Bad write to output file...event dumped.'
40       continue
      end do
778   continue
*----------------
*-- END MAIN LOOP
*----------------
       write(6,'(/,''  PD ** Last utdec (dec.hr)   :'',f9.4,
     & /,''        Sidereal time (hr)    :'',f9.4,
     & /,''        Derotation angle (deg):'',f9.4,
     & /,''        Telescope AZ (deg)    :'',f9.4,
     & /,''        Telescope EL (deg)    :'',f9.4)')
     & utdec,stime,theta*convd,az*convd,alt*convd
*
       write(infunit,'(/,''  PD ** Last utdec (dec.hr)   :'',f9.4,
     & /,''        Sidereal time (hr)    :'',f9.4,
     & /,''        Derotation angle (deg):'',f9.4,
     & /,''        Telescope AZ (deg)    :'',f9.4,
     & /,''        Telescope EL (deg)    :'',f9.4)')
     & utdec,stime,theta*convd,az*convd,alt*convd
      lastutdec = utdec
      close(dunit)
      if(fasci) then
        close(nunit)
      end if
*
*   Print info and then close Ntuple file
*
      if(fntuple) then
*       write(*,*)'closing ntuple'
*       call gr_ntuple('P')
        call gr_ntuple('C')
      end if
*
* Save HBOOK file with alpha and 2-d histograms
*
*     if(f2d) then
*       write(*,*)'hrout histogram'
        call HRPUT(MESHID,hfname,'T')
*       call HRPUT(ALPHAID,hfname,'T')
*     endif
*
      write(6,'(/,''  PD ** Sid. dur. of '',a6,''   : '',
     &  f8.3,'' minutes'')')filename,
     +   (time-first7)*366.2422/(60.0*365.2422)
      write(6,'(''        Livetime (UT seconds) : '',f8.3)'),ltime
*
      write(infunit,'(/,''  PD ** Sid. dur. of '',a6,''   : '',
     &  f8.3)')filename,
     &    (time-first7)*366.2422/(60.0*365.2422)
      write(infunit,
     &   '(''        Livetime (UT seconds) : '',f8.3)'),ltime
      
      write(6,'(''        Number of Events      :'',i9)'),nevents
      write(6,'(''        Number of Muon Rings  :'',i9)'),nring
      write(infunit,'(''        Number of Events      :'',i9)'),
     & nevents
      write(infunit,'(''        Number of Muon Rings  :'',i9)'),
     & nring
*
* Save pedvars, difference in pedvars, nitrogen gains and derotated camera
* coordinates for use in quicklook distributions generated by paw
*
      open(10,
     & file=pthcal(1:lnblnk(pthcal))//'pgvec/'//filename//'.pgvec',
     & status='unknown',err=89)
      call derot((lastutdec+firstutdec)/2.0,theta,alt,az,stime)
      if(.not.fderot) theta = 0.0
      ctheta = cos(theta)		
      stheta = sin(theta)
      mstheta = -1.0*stheta 
      mctheta = -1.0*ctheta
      do i=1,npmt
        xt = wxdeg(i)
        yt = wydeg(i)
*       rwxdeg(i) = xt*mctheta+yt*stheta ! Derotate camera coords
*       rwydeg(i) = xt*mstheta+yt*mctheta
        rwxdeg(i) = xt*ctheta+yt*mstheta ! Derotate camera coords
        rwydeg(i) = xt*stheta+yt*ctheta
      end do 
      do i=1,npmt
        diffped = pedvar(i)*pedvar(i)-prpedvar(i)*prpedvar(i)
        if(diffped.gt.0.0) then
          sign = 1.0
          diffped = sqrt(diffped)
        else
          sign = -1.0
          diffped = sign*sqrt(sign*diffped)
        endif
        write(10,'(f11.3,f11.3,f11.3,f11.3,f11.3,f11.3)')rwxdeg(i),
     &   rwydeg(i),peds(i),pedvar(i),diffped,gains(i)
      end do 
      close(10)
89    continue
*
* Update number of events in header block.
* NOTE: you write 79 characters on each record in the parameter
* file, but it must be opened with recl=80. One extra character
* (line feed) is inserted by Fortran.
*
      if(fasci) then
        open(nunit,file=newname,access='direct',status='unknown',
     +    form='formatted',recl=80)
*
* It is important to write the line feed at the end of each record.
*
        write(nunit,241,rec=6)eventsout,duration,(pads,i=1,61),lfch
        close(nunit)
      end if
*
* Message
*      write(6,490)
      if(fasci) then
        write(6,501)'  PD ** Created parameterized ascii file : ',
     &   newname
      endif
      if(fntuple) then
        write(6,501)'  PD ** Created parameterized ntuple file: ',
     &   nt_fname
      endif
*     write(6,490)
*     call newheader(6,runid,eventsout,duration,stdur,mode,source,
*    &   date,mjd,frjd,ra,dec,ut,st,azimuth,elevation,skyq,comms)
      write(6,290)iver,nver
*
* Data base file is opened in append mode.  Old entries are not 
* automatically replaced but most be manually removed.
*
* The columns of the database correspond to the following quantities
* (1) source name
* (2) file name and run number
* (3) mode (2-TRK,1-ON,0-OFF)
* (4) UT date
* (5) average elevation
* (6) sky quality (1-A,2-B,3-C)
* (7) run duration (after livetime cut)
* (8) number of raw events (after livetime cut)
* (9) number of events after software trigger
* (10) number of events after trigger+shape
* (11) number of events after trigger+orientation
* (12) number of events after trigger+shape+orientation
*
      open(munit,file=pthcal(1:lnblnk(pthcal))//'hrc.dbase',
     &  access='append',status='unknown')
      if(duration.gt.ltimehi) duration = ltimehi
      write(munit,910)source(1:11),filename,imode,utdate,alt*convd,
     &  duration,isky,npraw,nptrig,npshape,nporient,npboth
      close(munit)

*
* Alpha data base writen to in append mode.
* 
* The columns of the database correspond to the following quantities
* (1) source name
* (2) file name and run number
* (3) mode (2-TRK,1-ON,0-OFF)
* (4) duration
* (5) sky quality
* (6,7...24) alpha(1)-alpha(18)
*
      open(munit,file=pthcal(1:lnblnk(pthcal))//'hrc.alpha',
     &  access='append',status='unknown')
      write(munit,911)source(1:11),filename,imode,
     &  duration,isky,(alphahist(i),i=1,18)
      close(munit)
*
* Auxiliary file is opened in append mode. It contains the header
* block for each parameters file.
*
      open(munit,file=pthout(1:lnblnk(pthout))//mdev,access='append',
     &  status='unknown')
      write(munit,490)
      write(munit,500)newname
      write(munit,490)
      if (eventsin.ne.nevents) then
         write(6,'('' PD ** Event # mismatch: '',i6,x,i6)')
     &                                           eventsin,nevents
	 write(munit,'(''**Event # mismatch: '',i6,x,i6)')
     &                                           eventsin,nevents
      end if
      write(munit,190)softcode,sig_thresh,nbr_thresh,padstat,
     &                                             (pads,i=1,58)
      write(munit,200)runid,source,mode,(pads,i=1,50)
      write(munit,210)date,ut,st,(pads,i=1,63)
      write(munit,220)idint(mjd),frjd,(pads,i=1,60)
      write(munit,230)elevation,azimuth,(pads,i=1,67)
      write(munit,240)nevents,duration,(pads,i=1,61)
      write(munit,250)pdid,gainlabel,calibrators,(pads,i=1,22)
      write(munit,260)ifix(ra),ifix(dec),trig,skyq,(pads,i=1,50)
      write(munit,270)stubesoff,(pads,i=1,39)
      if (gpsyes) then
	 write(munit,281)gpsbeg,(pads,i=1,51)
      else
	 write(munit,280)(pads,i=1,79)
      end if
      write(munit,290)iver,nver

      close(infunit)
      call IGEND	! Close HIGZ
      call exit(0)	! Exit with no error
      stop
*
* Handle ERROR conditions
*
20       write(*,150)utdate
	 stop
*
* Format Statements 
*
90    format(' PD ** File does not exist: ',a)
100   format(' Enter filename (zzxxxx): ',$)
110   format(' Enter UT date (yymmdd): ',$)
130   format(' Enter nitrogen gains id: ',$)
140   format(' PD ** Cannot find pedestals for given id:  ',i6,1x,i4)
150   format(' PD ** Cannot find skygains for given date: ',i6)
160   format(' PD ** Cannot find n2 gains for given date: ',i6,1x,i4)
165   format(' PD ** Tubes turned off: ',15i4)
190   format(i2,2x,f5.2,2x,f5.2,2x,a3,58a)
200   format(a4,1x,a20,1x,a3,50a)
210   format(i6,1x,i4,1x,i4,63a)
220   format(i5,1x,f13.11,60a)
230   format(f5.2,1x,f6.2,67a)
240   format(i5,1x,f12.6,61a)
241   format(i5,1x,f12.6,62a)
250   format(a4,1x,a2,10(1x,i4),22a)
260   format(i7,1x,i7,1x,a6,1x,a6,50a)
270   format(10i4,39a)
280   format(79a)
281   format(7(xi3),51a)
290   format('  PD ** TUBE COORDINATE VERSION ',i1,' used',
     &     /,'        NEIGHBOR VERSION ',i1,' used',/)
300   format(i2,xf12.6,7(xf5.3),xi6,3(xi3),3(xf5.0),xi1,xi1)
400   format('  PD ** Processing file ',a)
425   format('  PD ** Unable to open Pedestal file for ',a)
450   format('  PD ** Unable to open Nitrogen gain for ',a)
490   format(80('-'))
500   format(' File: ',a)
501   format(/,a,15a)
810   format('  PD ** Configuration file ',a,' not found.')
820   format('  PD ** Section "par" in ',a,' not found.')
825   format('  PD ** Section "cut" in ',a,' not found.')
830   format(  '  A list of standard config files may be found in ',
     &       /,'  /usr/sc/config.  Copy one to ',a,' in',
     &       /,'  your current working directory.',
     &       /,'Exiting...')
840   format(' PD ** Unable to read enough lines from ',a)
850   format(' PD ** Unable to read ',a)
861   format(a,12a)
862   format(a,f12.6)
863   format(a,i6)
864   format(a,a4,a,a20)
866   format(a,i2,a20)
867   format('  PD ** Picture threshold :',f5.2,'  Boundary threshold :'
     & ,f5.2,/)
890   format(' ***************************  GPARAMDAT  ***************
     &************',/,' File: ',a12,2x,'Date: ',i6.6,2x,'N2 ID: ',i6,/)
895   format(10f7.3)
900   format(' SIZE:',i5,1x,'LENGTH:',f7.2,1x,'WIDTH:',f7.2,1x,
     & 'DIST:',f7.2,1x,'SIN ALPHA:',f7.2)
905   format(' LENGTH/TOTSIG:',f7.4,' MUONICITY:',f7.4,' NB2:',i2,
     & ' NB3:',i2,' NPICT:',i4)
906   format(' AZWIDTH:',f7.4,' HADRON:',f7.4,' PHI:',f7.2,
     & ' ASYMM:',f7.4)
907   format(' AZWIDTHP:',f7.4,' MISSP:',f7.4,' DISTP:',f7.2,
     & ' ASYMMP:',f7.4)
908   format(' FRAC2:',f7.4,' MAXSEP:',f11.4)
910   format(a11,1x,a6,1x,i1,1x,i6,1x,f5.1,1x,f7.1,1x,i1,1x,
     &  i6,1x,i6,1x,i6,1x,i6,1x,i6)
911   format(a11,1x,a6,1x,i1,1x,f7.1,1x,i1,1x,18(1x,i4))
912   format('  PD ** Standard cut values: ')
913   format('    ',i2,'/91 > ',f6.1,
     &     /,'     size > ',f6.1,
     &     /,'    ',f5.3,' < dist   < ',f5.3,
     &     /,'    ',f5.3,' < length < ',f5.3,
     &     /,'    ',f5.3,' < width  < ',f5.3,
     &     /,'    sin(alpha) < ',f5.3,/)
 914  format(' length/max1=',f8.5)
 915  format(' PMT',1x,' P.H.',1x,'  Xpos',1x,'  Ypos')
 916  format(1x,i3,1x,i5,1x,f6.3,1x,f6.3)
      end


*-----------------------------------------------------------------------
*       radec2azel
*
*	Converts ra and dec (in radians) to az and el (in radians)
*       given the sidereal time (in degrees)
*
*     subroutine radec2azel(ra,dec,sid,az,el)
*
*      implicit none
*      real*4 ra,dec,sid,az,el,ha
*      real*4 sindec,cosdec
*      real*4 sinha,cosha
*      real*4 cosphi,sinphi	! phi is the observers latitude
*      data cosphi / ??? /, sinphi / ??? /
*
*      ha = sid - ra 
*      cosha = cos(ha)
*      sinha = sin(ha)
*      cosdec = cos(dec)
*      sindec = sin(dec)
*
*      arg = cosdec*cosha
*      az = atan((cosdec*sinha)/(arg*sinphi-sindec*cosphi))
*      el = asin(sindec*sinphi+arg*cosphi)
*
*      return
*      end
*	

*-----------------------------------------------------------------------
*       Function PASS_MESH
*       JB
*       950410
*
*	Returns true if the image parameters pass the width/azwidth cuts
*       (as used in the U.Mich analysis) plus a width/length dependent
*       in the spirit of RCL and VC.
*
      logical function pass_mesh(length,width,miss,dist,azwidth,
     +  livetime)

      implicit none

      real*4 length,width,miss,dist,azwidth
      real*4 livetime
      real*4    disp
      real*4	DISPC
*
* 11m
*
*     parameter	(DISPC=0.3)	! Appropriate for 11m only!
*
* 10m
*
*     parameter	(DISPC=0.2125)
*

      real*4 lengthlo,lengthhi
      real*4 widthlo,widthhi
      real*4 distlo,disthi
      real*4 sinalphahi
      real*4 sizelo
      real*4 triglo
      real*4 lwidhi,llenhi
      real*4 azwidhi
      real*4 ltimehi
      common / cuts /
     +  lengthlo,lengthhi,
     +  widthlo,widthhi,
     +  distlo,disthi,
     +  sinalphahi,
     +  sizelo,
     +  triglo,
     +  lwidhi,llenhi,
     +  azwidhi,
     +  ltimehi
*
* Modified by JB 960109 - pass telescope id needed for cuts which
* is needed since cuts depend on which telescope
*
      integer N_TELE
      parameter (N_TELE = 2)
      integer*4 tele_id
      character*3 tele_name(N_TELE)
      common / tele /
     +  tele_id,
     +  tele_name
      integer
     +     TELE_10,
     +     TELE_11
      parameter
     +    (TELE_10 = 1,
     +     TELE_11 = 2)
*
      if(tele_id.eq.TELE_11) then
        DISPC = 0.3
      else
        DISPC = 0.2125
      endif
*
      if(length.gt.0.0) then
        if(tele_id.eq.TELE_11) then
          disp = abs( dist - 1.63 + 1.63*(width/length) )
        else
*         disp = abs( dist - 1.87 + 1.87*(width/length) )
          disp = abs( dist - 1.7 + 1.7*(width/length) )
*         disp = abs( dist - 2.0 + 2.0*(width/length) )
        endif
      else
        disp = 1.0e6
      endif

      if(azwidth.gt.0.) then
        if((width/azwidth.gt.0.925d0).and.
     +    (disp.lt.DISPC).and.
     +    (livetime.lt.ltimehi)) then
          pass_mesh = .true.
        else
          pass_mesh = .false.
        endif
      else
        pass_mesh = .false.
      endif

*     if(azwidth.lt.azwidhi) then
*       pass_mesh = .true.
*     else
*       pass_mesh = .false.
*     endif

      return
      end

*---------------------------------------------------------------------
*       PASS_LOOSE
*	JB
*       950410
*
*	Returns true if the parameters satisfy some loose shape and
*       trigger cuts.  These cuts are only used to filter data 
*       for the purpose of speeding up the 2-d analysis.
*
      logical function pass_loose(length,width,dist,size,max2,livetime)

      implicit none

      real*4 length,width,dist
      integer*4 size
      real*4 max2
      real*4 livetime

      real*4 lengthlo,lengthhi
      real*4 widthlo,widthhi
      real*4 distlo,disthi
      real*4 sinalphahi
      real*4 sizelo
      real*4 triglo
      real*4 lwidhi,llenhi
      real*4 azwidhi
      real*4 ltimehi
      common / cuts /
     +  lengthlo,lengthhi,
     +  widthlo,widthhi,
     +  distlo,disthi,
     +  sinalphahi,
     +  sizelo,
     +  triglo,
     +  lwidhi,llenhi,
     +  azwidhi,
     +  ltimehi
      logical fview
      logical fmesh
      logical ffomin
      logical fsuperd
      logical f2d
      logical fxy
      logical fuv
      logical fasymm
      logical floose
      logical fmuon
      common /flags/ fview,fmesh,ffomin,fsuperd,f2d,fxy,fuv,fasymm,
     &               floose,fmuon


      if(ffomin) then
        if((length.gt.lengthlo).and.(width.gt.widthlo).and.
     &    (length.lt.llenhi).and.(dist.lt.disthi).and.
     &    (width.lt.lwidhi).and.
     &    (livetime.lt.ltimehi).and.
     &    (max2.gt.triglo).and.(size.gt.sizelo)) then
          pass_loose = .true.
        else
          pass_loose = .false.
        endif
      else if(fsuperd.or.fmesh) then
        if(((length.gt.lengthlo).and.(length.lt.lengthhi)).and.
     &    ((width.gt.widthlo).and.(width.lt.widthhi)).and.
     &    (dist.lt.disthi).and.
     &    (max2.gt.triglo).and.
     &    (livetime.lt.ltimehi).and.
     &    (size.gt.sizelo)) then
          pass_loose = .true.
        else
          pass_loose = .false.
        endif
      endif
      return
      end
*
*---------------------------------------------------------------------
*       PASS_SHAPE
*       JB
*       941201
*
*	Returns true if the image parameters pass the shape cuts
*
      logical function pass_shape(length,width,
     &  size,max2,livetime)

      implicit none

      real*4 lengthlo,lengthhi
      real*4 widthlo,widthhi
      real*4 distlo,disthi
      real*4 sinalphahi
      real*4 sizelo
      real*4 triglo
      real*4 lwidhi,llenhi
      real*4 azwidhi
      real*4 ltimehi
      common / cuts /
     +  lengthlo,lengthhi,
     +  widthlo,widthhi,
     +  distlo,disthi,
     +  sinalphahi,
     +  sizelo,
     +  triglo,
     +  lwidhi,llenhi,
     +  azwidhi,
     +  ltimehi

      real*4 length,width
      integer*4 size
      real*4 max2
      real*4 livetime

      if(((length.gt.lengthlo).and.(length.lt.lengthhi)).and.
     &  ((width.gt.widthlo).and.(width.lt.widthhi)).and.
     &  (max2.gt.triglo).and.(size.gt.sizelo).and.
     &  (livetime.lt.ltimehi)) then
        pass_shape = .true.
      else
        pass_shape = .false.
      endif

      return
      end

*-----------------------------------------------------------------------
*       Function PASS_CUTS
*       JB
*       941201
*
*       Returns true if parameter pass supercuts + trigger cut
*
      logical function pass_cuts(length,width,miss,dist,
     &  sinalpha,size,max2,livetime,passmsg)

      implicit none

      integer NALPBINS
      parameter (NALPBINS=18)
      real*4 DALPHA
      parameter (DALPHA=5.0)
      real*8	convd
      data convd/57.29577951/
      real*4 length,width,miss,dist,sinalpha,livetime
      integer*4 size
      real*4 max2
      real*4 lengthlo,lengthhi
      real*4 widthlo,widthhi
      real*4 distlo,disthi
      real*4 sinalphahi
      real*4 sizelo
      real*4 triglo
      real*4 lwidhi,llenhi
      real*4 azwidhi
      real*4 ltimehi
      common / cuts /
     +  lengthlo,lengthhi,
     +  widthlo,widthhi,
     +  distlo,disthi,
     +  sinalphahi,
     +  sizelo,
     +  triglo,
     +  lwidhi,llenhi,
     +  azwidhi,
     +  ltimehi
      logical fview
      logical fmesh
      logical ffomin
      logical fsuperd
      logical f2d
      logical fxy
      logical fuv
      logical fasymm
      logical floose
      logical fmuon
      common /flags/ fview,fmesh,ffomin,fsuperd,f2d,fxy,fuv,fasymm,
     &               floose,fmuon
      logical ftrigger
      logical fdist
      logical fshape
      logical forient
      logical fboth
      integer*4 ievt
      integer*4 npraw
      integer*4 nptrig
      integer*4 npshape
      integer*4 nporient
      integer*4 npboth
      common /npass/ ievt,npraw,nptrig,npshape,nporient,npboth
      integer*4 alphahist(NALPBINS)
      common /ahist/ alphahist
      real*4 ralpha
      integer ialpha
      character*200 passmsg
      character*7 fail1,fail2,fail3,fail4,fail5,fail6
      integer*4 lnblnk
      real*8 dsin

      ftrigger = .false.
      fdist = .false.
      fshape = .false.
      forient = .false.
      fboth = .false.
      if(livetime.lt.ltimehi) then
       npraw = npraw+1
       if((max2.gt.triglo).and.(size.gt.sizelo)) then
        ftrigger = .true.
        nptrig = nptrig+1
        if((dist.gt.distlo).and.(dist.lt.disthi)) then
          fdist = .true.
          if(((length.gt.lengthlo).and.(length.lt.lengthhi)).and.
     &       ((width.gt.widthlo).and.(width.lt.widthhi))) then
            fshape = .true.
            npshape = npshape+1
*
* Update alpha histogram for all events passing shape cuts
*
            dsin = dble(sinalpha)
            if(dsin.lt.-1.0d0) then
              dsin = -1.0d0
            else
              if(dsin.gt.1.0d0) dsin = 1.0d0
            endif
            ralpha = asin(dsin)*convd	! Calculate alpha in deg.
*           call HFILL(ALPHAID,ralpha,0.0,1.0)
            ialpha = int(ralpha/DALPHA)+1
            if((ialpha.gt.0).and.(ialpha.le.NALPBINS)) then
              alphahist(ialpha) = alphahist(ialpha) + 1
            endif
          end if
          if(sinalpha.lt.sinalphahi) then
            forient = .true.
            nporient = nporient+1
          end if
          if(fshape.and.forient) then
            fboth = .true.
            npboth = npboth+1
          end if ! shape
         end if ! dist
       end if ! trigger
      end if ! livetime

      if(.not. fboth) then
        if(fview) then
          if(.not. ((length.gt.lengthlo).and.(length.lt.lengthhi))) then
           fail1 = 'length '
*          failstr = failstr(1:lnblnk(failstr))//'length '
          else
           fail1 = ' '
          end if
          if(.not. ((width.gt.widthlo).and.(width.lt.widthhi))) then
           fail2 = 'width '
          else
           fail2 = ' '
          end if
          if(.not. ((dist.gt.distlo).and.(dist.lt.disthi))) then
           fail3 = 'dist '
          else
           fail3 = ' '
          end if
          if(.not. (sinalpha.lt.sinalphahi)) then
           fail4 = 'alpha '
          else
           fail4 = ' '
          end if
          if(.not. (size.gt.sizelo)) then
           fail5 = 'size '
          else
           fail5 = ' '
          end if
          if(.not. (max2.gt.triglo)) then
           fail6 = 'trig '
          else
           fail6 = ' '
          end if
          passmsg = 'Failed cuts: '
     &     //fail1(1:lnblnk(fail1))//' '//fail2(1:lnblnk(fail2))
     &     //' '//fail3(1:lnblnk(fail3))//' '//fail4(1:lnblnk(fail4))
     &     //' '//fail5(1:lnblnk(fail5))
     &     //' '//fail6(1:lnblnk(fail6))
         end if
         pass_cuts = .false.
      else
         pass_cuts = .true.
      end if
      return
      end

*-----------------------------------------------------------------------
*	NB_TRIGGER
*
*	Set flags nb2, nb3 to indicate that the signals in at least
*       two or three adjacent tubes are above threshold
*
      subroutine nb_trigger(event,nbthresh,nb2,nb3)

      integer*4 npmt,mpmt,ppmt
      common /num/ npmt,mpmt,ppmt
      integer*4 maxpmts
      parameter (maxpmts=1000)
      real*4 event(maxpmts)
      real*4 nbthresh
      integer*4 nb2,nb3
      integer evtdisc(maxpmts)
      integer disc_sum
      integer*4 nbr_list(7,maxpmts)
      common /neigh/ nbr_list
      integer i,j
      real*4 store
*
*  nbr_list modified 941116 by JB - PMT #110 now flags the absence of
*  a neighbor tube.  PMT #0 was previously used to flag no neighbor,
*  but this is an illegal index of an array and caused an
*  error.  Of course this is only an error in some versions of fortran.
*  (Yet another in a growing list of reasons why our code should be
*  written in C)
*
*  970922 MAC modified to eliminate zeroing of tube 110 (we have many more
*  than that now.  And eliminate the first loop through the PMTs.
*

C       print*,'Entering loop'
C       do i=1,6
C          do j=1,mpmt
C             if (nbr_list(i,j).GT.0.0) print *,'Nbr > 0'
C          enddo
C       enddo
C       print*,'Exiting loop'
C       print*,'Exiting loop'


      ndisc = 0
      nb2 = 0
      nb3 = 0
*      store = event(110)
*      event(110) = 0.0
      do i=1,npmt
        if (event(i).gt.nbthresh) then
           disc_sum = 7
           do j=1,7
              if (event(nbr_list(j,i)).gt.nbthresh) 
     >             disc_sum = disc_sum + 1
           end do
           if(disc_sum.ge.8) then
              nb2 = 1
              if(disc_sum.ge.9) then
                 nb3 = 1
                 return         ! no more to be learned so return
              end if ! if nb3
           end if               ! if nb2
        endif
      enddo
*mac      do i=1,npmt
*mac        if (event(i).gt.nbthresh) then
*mac          evtdisc(i) = 1 
*mac          ndisc = ndisc + 1
*mac        else
*mac          evtdisc(i) = 0
*mac        end if
*mac      end do
*mac      if (ndisc.ge.2) then
*mac        do i=1,npmt
*mac          if(evtdisc(i).eq.1) then
*mac            disc_sum = 6
*mac            do j=1,6
*mac              disc_sum = disc_sum + evtdisc(nbr_list(j,i))
*mac            end do
*mac            if(disc_sum.ge.7) then
*mac              nb2 = 1
*mac              if(disc_sum.ge.8) then
*mac                nb3 = 1
*mac                return	! no more to be learned so return
*mac              end if ! if nb3
*mac            end if ! if nb2
*mac          end if ! if evtdisc
*mac        end do ! do for each pmt 
*mac      end if ! if ndisc>=2
*      event(110) = store
      return
      end
*-----------------------------------------------------------------------
*     IMCLEAN
*
* Clean up image
*
* (1) zero tubes which are below thresholds
* 
* I/O parameters
*     pmts : ADC counts from HRC
*     sig_thresh : Number of sigma for picture
*     nbr_thresh : Number of sigma for boundary
*     pad : True to do padding
*     event : Cleaned, flat-fielded HRC event (Whipple labelling)
*
* Other variables of interest:
*     wxdeg,wydeg : grid of coordinates (Whipple labelling)
*
*
*     941115 Modified by JB to return npicture, the number of tubes
*     in the picture or boundary.
*
*     950629 Modified by JB to return an array of picture PMTs
*
*     950710 Modified by JB to zero slightly negative tubes
*
*     970922 Modified by MAC to only loop through tubes in picture
*     to determine the tubes passing the neighbor criteria.  Also,
*     reduce the number of loops through the pmts.
*

      subroutine gimclean(pmts,sig_thresh,nbr_thresh,pad,event,
     &  npicture,pictpmt,pictx,picty,pictnum)

      implicit none

      integer*4 npmt,mpmt,ppmt
      common /num/ npmt,mpmt,ppmt
      integer*4 maxpmts
      parameter(maxpmts=1000)

* Input/output variables
      real*4 event(maxpmts)
      real*4 pictpmt(maxpmts)
      real*4 pictx(maxpmts),picty(maxpmts)
      real*4 sig_thresh,nbr_thresh
      integer pictnum(maxpmts),npicture
      integer*4 pmts(maxpmts)
      logical pad

* Variables used in this program
      real*4 gasdev    !Gaussian deviation function
      real*4 secnds    !System time function
      character syscmd*80  !System call
      integer*4 i,j,k
      real*4 size
      integer inpict(maxpmts)
      integer*4 bits(0:31),nbr(0:23,maxpmts)
*mac      integer*4 picture(0:23),boundary(0:23)
*Replace bit masks with these to make things clearer and fit with my scheme.
      logical picture(maxpmts)
      integer*4 l,m
      real*4 peds(maxpmts),pedvar(maxpmts),gains(maxpmts)
      real*4 prpeds(maxpmts),prpedvar(maxpmts)
      integer*4 tubesoff(maxpmts)
      real*4 maxpedvar(maxpmts),diffpedvar(maxpmts),x
      integer*4 nboundary,subpict
      logical inboundary
*
* Stuff for random number generator
*
      integer*4 iseed,idum
*mac      data iseed / 12345678 /
      data idum / -123456 /

      common /info/peds,pedvar,prpeds,prpedvar,gains,tubesoff,nbr
      common /bytes/bits
    
*
* PMT coordinates and list of neighboring tubes for each pmt.
*
      real*4 wxdeg(maxpmts),wydeg(maxpmts)
      integer*4 nbr_list(7,maxpmts)
      common /wcoords/wxdeg,wydeg
      common /neigh/ nbr_list
*
* Data variables
*
      logical first
      data first/.true./
*
      npicture = 0 
*
*  Do padding as the first thing and subtract off the noise pedestals.
*  Also, reset some fields used later
*
      if (pad) then
* Set the seed for the random number generator as the seconds since 
* midnight.  For some reason, gasdev does not like seeds bigger than
* 10000, so remove that highest value from the seed.
         if (first) then
            iseed = nint(secnds(0.0))
            iseed = iseed -int(float(iseed)/10000.)*10000
            first=.false.
            type*,' '
            type*,' PD ** Padding seed =',iseed
            type*,' '
         endif
         do i = 1,npmt
*
            picture(i)=.false.
*
            maxpedvar(i)=amax1(pedvar(i),prpedvar(i))
*            diffpedvar(i)=sqrt(amax1(0.,
*     &           (maxpedvar(i)**2 - pedvar(i)**2 )))
* MAC 970922 - why loop through twice for this
*mac         end do

*mac         do k = 1,npmt
               x = gasdev(idum)
*save            x = gasdev(iseed)
*            event(i) = float(pmts(i))+diffpedvar(i)*x
            event(i) = float(pmts(i))+
     >           x*sqrt(amax1(0.,
     &           (maxpedvar(i)**2 - pedvar(i)**2 )))
     >           -peds(i)
         end do
      else
         do k = 1,npmt
*
            picture(k)=.false.
*
            event(k) = float(pmts(k))-peds(k)
         end do
      end if


*
* Topological threshold application and gains subtraction
* Initialize boundary and picture elements to zero
*
*mac      do k=0,23
*mac         picture(k)=0
*mac         boundary(k)=0
*mac      enddo
*
* Subtract noise pedestals
*
*mac      do k=1,npmt
*mac         event(k)=event(k)-peds(k)
*mac      enddo
*
* Turn off tubes from tubelist, this loop was moved up to here in
*  case of zero or near-zero tube variances, which would cause tube
*  to be in "picture" when it shouldn't, suggested MFC. MP911021
* Initialization of k included..left out after loop moved SF930329
* Moved here; was before calc of event.  Final fix, I hope ADK 930331
* Order of statements inside DO WHILE reversed  ADK 930331
*
      k=1
      do while (k.le.maxpmts.and.tubesoff(k).ne.0)
         event(tubesoff(k))=0.0
         k=k+1
      enddo

*
* must use maxpedvar in picture/boundary decisions if
* padding is switched on.  Fixed, MFC 940820
*
* Decide if tube should be in picture (=> over some # of sigma)
*
      if(pad)then
         do k=1,ppmt
            if(event(k).gt.sig_thresh*maxpedvar(k))then
*mac               l=mod(k-1,32)
*mac               m=(k-1)/32
*mac               picture(m)=picture(m).or.bits(l)
               npicture = npicture + 1
               inpict(npicture)=k
               picture(k)=.true.
            endif
         enddo
      else
         do k=1,ppmt
            if(event(k).gt.sig_thresh*pedvar(k))then
*mac               l=mod(k-1,32)
*mac               m=(k-1)/32
*mac               picture(m)=picture(m).or.bits(l)
               npicture = npicture + 1
               inpict(npicture)=k
               picture(k)=.true.
            endif
         enddo
      endif
*
* Decide if tube should be in boundary
*  (=> at least one neighbour in picture
*         and tube over some # of sigma)
*
      nboundary=0
      if(pad)then
*mac        do k=1,npmt
*mac           if(event(k).gt.nbr_thresh*maxpedvar(k))then
*mac              j=0
*mac              inboundary = .false.
*mac              do while (j.lt.24)
*mac                 if((nbr(j,k).and.picture(j)).ne.0)then
*mac                    l=mod(k-1,32)
*mac                    m=(k-1)/32
*mac                    boundary(m)=boundary(m).or.bits(l)
*                   inboundary = .true.
*mac                    j=23
*mac                 endif
*mac                 j=j+1
*mac              end do
*             if(inboundary) npicture=npicture+1
*mac           end if
*mac        end do
         do k=1,npicture
            do j=1,7
               l=nbr_list(j,inpict(k))
               if ((event(l).gt.nbr_thresh*maxpedvar(l)).and.
     >              (.not.picture(l)).and.
     >              (l.le.ppmt)) then
                  picture(l)=.true.
                  nboundary=nboundary+1
                  inpict(npicture+nboundary)=l
               endif
            enddo
         enddo
      else
*mac        do k=1,npmt
*mac           if(event(k).gt.nbr_thresh*pedvar(k))then
*mac              j=0
*mac              inboundary=.false.
*mac              do while (j.lt.24)
*mac                 if((nbr(j,k).and.picture(j)).ne.0)then
*mac                    l=mod(k-1,32)
*mac                    m=(k-1)/32
*mac                    boundary(m)=boundary(m).or.bits(l)
*                   inboundary = .true.
*mac                    j=23
*mac                 endif
*mac                 j=j+1
*mac              end do
*             if(inboundary) npicture=npicture+1
*mac           end if
*mac        end do
         do k=1,npicture
            do j=1,7
               l=nbr_list(j,inpict(k))
               if ((event(l).gt.nbr_thresh*pedvar(l)).and.
     >              (.not.picture(l)).and.
     >              (l.le.ppmt)) then
                  picture(l)=.true.
                  nboundary=nboundary+1
                  inpict(npicture+nboundary)=l
                  if (l.eq.0) then
                     type*,'pmt=',inpict(k)
                     type*,'inbr=',j
                     type*,'pnbr=',l
                  endif
               endif
            enddo
         enddo
      end if
* Picture comprises all picture and boundary tubes
*mac      do k=0,23
*mac         picture(k)=picture(k).or.boundary(k)
*mac      end do

* npicture is all tubes in picture and boundary
      npicture=npicture+nboundary
* apply gains  - only if the tube is in the picture.
*mac      do k=1,npmt
*mac         if (picture(k)) then
*save            event(k)=event(k)*gains(k)
*mac            event(k)=max(0.0,event(k)*gains(k))
*mac            if (event(k).gt.0.) then
*mac               npicture = npicture+1
*mac               pictpmt(npicture) = event(k)
*mac               pictx(npicture) = wxdeg(k)
*mac               picty(npicture) = wydeg(k)
*mac               pictnum(npicture) = k
*mac            endif
*mac         else
*mac            event(k)=0.
*mac         endif
*mac      end do
      subpict=0
      do k=1,npicture
         l=inpict(k)
         event(l)=max(0.0,event(l)*gains(l))
         if (event(l).gt.0.) then
            pictpmt(k-subpict) = event(l)
            pictx(k-subpict) = wxdeg(l)
            picty(k-subpict) = wydeg(l)
            pictnum(k-subpict) = l
         else
            subpict=subpict+1
         endif
      enddo
      npicture=npicture-subpict
* Turn off tubes below thresholds.
*mac      do k=1,npmt
*
* Added by JB 950710 - zero negative tubes which have somehow made it
* into the picture.
*
*save         if(event(k).lt.0.0) event(k)=0.0
*
*mac         l=mod(k-1,32)
*mac         m=(k-1)/32
*mac         if ((bits(l).and.picture(m)).eq.0) then
*save         if (.not.picture(k)) then
*save           event(k)=0.0
*save         else 
*save           npicture = npicture+1
*save           pictpmt(npicture) = event(k)
*save           pictx(npicture) = wxdeg(k)
*save           picty(npicture) = wydeg(k)
*save           pictnum(npicture) = k
*save         endif
*save      end do

*
* NEGATIVE SIZE
* Events are flagged if their size after calibration is negative.
* This is necessary to deal with pedestals' oscillation and could be
* removed at a later stage.
*
      size=0.0
      call nsum(event,npmt,size)
      if(size.le.0.0)size=-1.

      return
      end

************************************************************************
*     Subroutine REHILLAKPAR                        
*     JB
*     950412
*
*     Re-parameterize at a new origin for 2-d analysis, but only
*     redo necessary calculations to save time.
*
*
	subroutine rehillakpar(event,icode,npicture,pictnum,
     &     length,width,miss,
     &     dist,azwidth,xmean,ymean,asymmetry,phi,valid,x0,y0,
     &     rotcoords)
 
	implicit none
* Parameters
        real*8	PI
        parameter(PI=3.141592654)
	integer*4 maxpmts
	parameter(maxpmts=1000)
        real*4	HADR2		! Hadronicity radius^2
        parameter(HADR2=0.25)
* Input/Output variables
        integer npicture,pictnum(maxpmts)
*
        logical valid
        integer*4 npmt,mpmt,ppmt
        common /num/ npmt,mpmt,ppmt
        integer*4 i,j
	integer*4 icode
	real*4 event(maxpmts),rwxdeg(maxpmts),rwydeg(maxpmts)
        real*4 wxdeg(maxpmts),wydeg(maxpmts)
	real*4 length,width,miss,dist,azwidth
        real*4 muon,hadron,asymmetry,phi,x0,y0
	real*4 dx,dy,dmean2
	real*4 sumsig,sumxsig,sumx2sig,sumysig,sumy2sig,sumxysig
	real*4 xmean,x2mean,ymean,y2mean,xymean
	real*4 sdevx2,sdevy2,sdevxy
        real*4 x3mean,y3mean,x2ymean,xy2mean
        real*4 sumx3sig,sumy3sig,sumx2ysig,sumxy2sig
        real*4 sdevx3,sdevy3,sdevx2y,sdevxy2
        real*4 hyp,cosphi,sinphi,cos2phi,sin2phi
	real*4 d,z,u,v
        real*4 xi,yi,si
	real*4 wxbyev,wybyev,xmean2,ymean2,meanxy
        real*4 xmean3,ymean3
        real*4 wx2byev,wy2byev
        real*4 arg,argl
        real*4 arga,argb	! Result of intermed. calculation to
				! speed computation
        
	real*4 sign
        logical	rotcoords
	common /rwcoords/rwxdeg,rwydeg
        common /wcoords/wxdeg,wydeg
        logical fview
        logical fmesh
        logical ffomin
        logical fsuperd
        logical f2d
        logical fxy
        logical fuv
        logical fasymm
        logical floose
        logical fmuon
        common /flags/ fview,fmesh,ffomin,fsuperd,f2d,fxy,fuv,fasymm,
     &                 floose,fmuon

	sumsig=0.0
	sumxsig=0.0
	sumysig=0.0
	sumx2sig=0.0
	sumy2sig=0.0
	sumxysig=0.0
	sumx3sig=0.0
        sumy3sig=0.0
        sumx2ysig=0.0
        sumxy2sig=0.0
  
        valid = .true.
*mac	do i=1,npmt
*mac	  if (event(i).gt.0.0) then
        do j=1,npicture
           i=pictnum(j)
            if(rotcoords) then
              xi = rwxdeg(i) - x0
              yi = rwydeg(i) - y0
            else
              xi = wxdeg(i) - x0
              yi = wydeg(i) - y0
            endif
            si = event(i)
	    wxbyev = xi*si
            wx2byev = xi*wxbyev
	    wybyev = yi*si
            wy2byev = yi*wybyev
	    sumsig = sumsig + si
	    sumxsig = sumxsig + wxbyev
	    sumx2sig = sumx2sig + wx2byev
	    sumysig = sumysig + wybyev
	    sumy2sig = sumy2sig + wy2byev
	    sumxysig = sumxysig + xi*wybyev
*---------------------------------------------------------------------
* The following code is needed to calculate the asymmetry param.
*
            if(fasymm) then
              sumx3sig = sumx3sig + xi*wx2byev		! asymm
              sumy3sig = sumy3sig + yi*wy2byev		! asymm
              sumx2ysig = sumx2ysig + xi*xi*wybyev	! asymm
              sumxy2sig = sumxy2sig + yi*yi*wxbyev	! asymm
            endif
*
*---------------------------------------------------------------------
*mac	  endif
	enddo

	if (sumsig.le.0.0) then
          valid = .false.
	  icode=99
	  return
	endif

	xmean = sumxsig / sumsig
	x2mean = sumx2sig / sumsig
	ymean = sumysig / sumsig
	y2mean = sumy2sig / sumsig
	xymean = sumxysig / sumsig
	xmean2 = xmean * xmean
	ymean2 = ymean * ymean
	meanxy = xmean * ymean
 
	sdevx2 = x2mean - xmean2
	sdevy2 = y2mean - ymean2
	sdevxy = xymean - meanxy

*--------------------------------------------------------------------
* The following code has been added to calculate the asymmetry param
*
        if(fasymm) then
	  x3mean = sumx3sig /sumsig	! asymm
          y3mean = sumy3sig / sumsig	! asymm
          x2ymean = sumx2ysig / sumsig	! asymm
          xy2mean = sumxy2sig / sumsig	! asymm
          xmean3 = xmean2 * xmean		! asymm
          ymean3 = ymean2 * ymean		! asymm
          sdevx3 = x3mean-3.0*x2mean*xmean+2.0*xmean3	! Not same as MP 
	  sdevy3 = y3mean-3.0*y2mean*ymean+2.0*ymean3	! " "
	  sdevx2y = x2ymean-2.0*xymean*xmean+2.0*xmean2*ymean
     &            -x2mean*ymean
	  sdevxy2 = xy2mean-2.0*xymean*ymean+2.0*xmean*ymean2
     &            -xmean*y2mean
        endif
*
* Mod. by JB 940511 to avoid sqrt of neg. number
*
 	arg = xmean2 + ymean2
        if(arg.gt.0.0) then
          dist = SQRT(arg)
        else
          dist = 0.0
*         valid = .false.
        endif

	d = max(1.e-15,sdevy2 - sdevx2)
        arg = max(1.e-30,d*d + 4.0*sdevxy*sdevxy)
	if(arg.gt.0.0) then
	  z = SQRT(arg)
	else
          z = 0.0
        endif

	if(z.gt.0.0)then
	  u=1+d/z
	  v=2-u
          arg = (u*xmean2+v*ymean2)/2.0 - meanxy*(2.0*sdevxy/z)
          if(arg.ge.0.0) then
 	    miss=SQRT(arg)
          else
            miss = 0.0
*           valid = .false.
          end if
	else
	  miss = dist
        endif
*--------------------------------------------------------------------
* Begin asymm calculation
*
        if(fasymm) then
          arga = (d+z)*ymean+2.0*sdevxy*xmean
          if(arga.ge.1.0) then
            arga = 1.0
          else if(arga.le.-1.0) then
            arga = -1.0
          endif
          argb = 2.0*sdevxy*ymean-(d-z)*xmean
          if(argb.ge.1.0) then
            argb = 1.0
          else if(argb.le.-1.0) then
            argb = -1.0
          endif
          if((arga.eq.0.0).and.(argb.eq.0.0)) then
            phi = PI/2.0
          else
            phi = atan2(arga, argb)
          endif

*
* Put in yet another patch to eliminate a floating point error -
* this time floating underflow.  What gives?  Is it fortran or
* OSF/1 which is so perverse?
*
          if(abs(phi).lt.1.0e-8) phi = 0.0

          cosphi = cos(phi)
          sinphi = sin(phi)
*         phi = atan2((d+z)*ymean+2*sdevxy*xmean ,
*    &          2*sdevxy*ymean-(d-z)*xmean )
*         cosphi = cos(phi)
*         sinphi = sin(phi)
          cos2phi = cosphi*cosphi
          sin2phi = sinphi*sinphi

*
* Put in patch for floating underflow error, until a better
* solution is derived.
*
          if(abs(sdevx2y).lt.1.0e-10) sdevx2y = 0.
          if(abs(sdevx3).lt.1.0e-10) sdevx3 = 0.
          if(abs(sdevxy2).lt.1.0e-10) sdevxy2 = 0.
          if(abs(sdevy3).lt.1.0e-10) sdevy3 = 0.
*
          arg = sdevx3*cosphi*cos2phi+3.0*sdevx2y*cos2phi*sinphi+
     &      3.0*sdevxy2*cosphi*sin2phi+sdevy3*sinphi*sin2phi
          
          argl = length*length*length
          if(argl.gt.0.0) then
            asymmetry = arg/argl
            if(asymmetry.ge.0.0) then
              sign = 1.0
            else
              sign = -1.0
            endif 
            asymmetry = abs(asymmetry)
            asymmetry = sign*(asymmetry**(0.33333333))
          else
            if (arg.lt.0.0) then
              asymmetry = -1.0e10	! minus infinity
            else 
	      asymmetry = 1.0e10	! plus infinity
            endif
          endif ! argl.gt.0.0
        endif ! asymm
*
* Note that d and z are now defined differently than the values used
* to calculate miss and asymmetry
*
	d = y2mean - x2mean
*
* Modified by JB 950710 to avoid sqrt of neg number
*	z = SQRT(d*d+4*xymean*xymean)
        arg = d*d+4.0*xymean*xymean
        if(arg.gt.0.0) then
	  z = SQRT(arg)
        else
          z = 0.0	
        endif
        arg = x2mean+y2mean-z
	if(arg.gt.0.0) then
	  azwidth = SQRT(arg/2.0)
        else
*	  write(*,*)'ERROR: can't evaluate azwidth = sqrt',arg,'/2.0'
          azwidth = 0.0		! TEMPORARY: Should flag invalid event?
        endif
 
	return
	end

************************************************************************
*   Subroutine ghillakpar
*   M. Punch, G. Vacanti, J. Buckley
*   900105-950410
*   
*   Routine to calculate the Hillas parameters for cleaned images.
*   This routine is the heart and soul of the Whipple data analysis,
*   please treat it with respect and even a little bit of awe
*   and reverence.
*
*   Revision History (incomplete) :
*
*   Use of Akerlof Azwidth (Akwidth?) implemented MP 910112
*
*   Simpler Miss formula suggested CA, revised MP, implemented 901101
*
*   Simpler Width and Length formulae suggested CA, implemented 901017
*
*   Yet another little bug fixed by MP 900418- Generalize the case of
*       the horizontal line.
*
*   Bug fixed by RCL 900307-The case of the horizontal line image was
*       never considered.
*
*   Bug fixed by GV 900215
*   This version takes events in WHIPPLE format and parameterises them.
*   Max is called here, its results are passed to calling procedure.
*   event code is passed to this subroutine. If something goes wrong,
*   this is changed to 99.
*
*   G. Vacanti introduces the common statement: coordinates are
*   computed only in the main program.
*
*   Modified by M. Punch to make it faster April, 1989
*   mjl 10 dec 87
*  
*   Modified by JB 940917 to also return xmean, ymean, and valid
*   
*   Modified by JB to fix divide by zeros, sqrt of negative numbers
*   and to calculate the following parameters:
*	asymm		Third moment or skew as in Punch Ph.D. thesis
*                       with some corrections to the formulae by JB
*       hadron          Hadronicity parameter as used by MCC, M.Urban
*                       et al. in UV analysis
*       muon            New parameter by MCC,JB to pick out large
*                       rings needed for the muon calibration
*
*   Modified by JB to calculate parameters about the origin
*   specified by the variables x0,y0.  This is done partially for
*   2-d analysis and in anticipation of the day when the encoder 
*   values are stored for each event. (x0=0,y0=0 for standard analysis)
*
*   Modified by JB 950706 - now pass ghillakpar the derotation angle 
*   theta, and calculate the most likely point of origin of the gamma
*   ray: xgam,ygam as well as the azwidth, miss and dist taking
*   this point as the new origin:
*   azwidthp,missp,distp
*
*   Modified by JB 961206 - derotation is always done before this
*   routine is invoked
*

	subroutine ghillakpar(event,icode,npicture,pictnum,
     &       length,width,miss,
     &       dist,azwidth,frac2,frac3,isize,loc1,loc2,loc3,
     &       maxsig1,maxsig2,maxsig3,xmean,ymean,hadron,asymmetry,
     &       muon,phi,valid,x0,y0,theta,xgam,ygam,azwidthp,missp,
     &       distp,asymp)
 
	implicit none
* Parameters
        real*4	PI
        parameter(PI=3.141592654)
        integer*4 ndbg		! Line number at which you have trapped
                                ! an error - print detailed info after
        parameter(ndbg=13733)
	integer*4 maxpmts
	parameter(maxpmts=1000)
        real*4	HADR2		! Hadronicity radius^2
        parameter(HADR2=0.25)
        real*4	mrtd
        parameter(mrtd=0.057296)
* Input/output variables
        integer npicture,pictnum(maxpmts)
*
        logical valid
        logical validp
        integer*4 npmt,mpmt,ppmt
        common /num/ npmt,mpmt,ppmt
        integer*4 i,j
	integer*4 loc1,loc2,loc3,isize,icode
        integer*4 icodep
	real*4 event(maxpmts),rwxdeg(maxpmts),rwydeg(maxpmts)
        real*4 wxdeg(maxpmts),wydeg(maxpmts)
	real*4 maxsig1,maxsig2,maxsig3
	real*4 length,width,miss,dist,azwidth,frac2,frac3
        real*4 muon,hadron,asymmetry,phi,x0,y0
        real*4 xgam1p,xgam2p,ygam1p,ygam2p,xgam1,xgam2,ygam1,ygam2
        real*4 azwidthp,missp,distp,dummyp,asymp
        real*4 phip,lengthp,widthp,xmeanp,ymeanp
        real*4 displ
        real*4 xgam,ygam
        real*4 theta,ctheta,stheta
	real*4 dx,dy,dmean2
	real*4 sumsig,sumxsig,sumx2sig,sumysig,sumy2sig,sumxysig
	real*4 xmean,x2mean,ymean,y2mean,xymean
	real*4 sdevx2,sdevy2,sdevxy
        real*4 x3mean,y3mean,x2ymean,xy2mean
        real*4 sumx3sig,sumy3sig,sumx2ysig,sumxy2sig
        real*4 sdevx3,sdevy3,sdevx2y,sdevxy2
        real*4 hyp,cosphi,sinphi,cos2phi,sin2phi
	real*4 d,z,u,v
        real*4 xi,yi,si
	real*4 wxbyev,wybyev,xmean2,ymean2,meanxy
        real*4 xmean3,ymean3
        real*4 wx2byev,wy2byev
        real*4 arg,argl 	! Result of intermed. calculation to
        real*4 arga,argb        ! speed computation
        
        real*4 axslope	! TEMPORARY
	real*4 sign

        real*4	maj_axis,min_axis,x,y,rad,xavg,yavg,out

	common /rwcoords/rwxdeg,rwydeg
*	common /wcoords/wxdeg,wydeg
        logical fview
        logical fmesh
        logical ffomin
        logical fsuperd
        logical f2d
        logical fxy
        logical fuv
        logical fasymm
        logical floose
        logical fmuon
        common /flags/ fview,fmesh,ffomin,fsuperd,f2d,fxy,fuv,fasymm,
     &                 floose,fmuon
        logical rotcoords
*
* Modified by JB 960109 - pass telescope id needed for cuts which
* is needed since cuts depend on which telescope
*
        integer N_TELE
        parameter (N_TELE = 2)
        integer*4 tele_id
        character*3 tele_name(N_TELE)
        common / tele /
     +    tele_id,
     +    tele_name
        integer
     +       TELE_10,
     +       TELE_11
        parameter
     +      (TELE_10 = 1,
     +       TELE_11 = 2)
*

	rotcoords = .true.
*	rotcoords = .false.

	sumsig=0.0
	sumxsig=0.0
	sumysig=0.0
	sumx2sig=0.0
	sumy2sig=0.0
	sumxysig=0.0
	sumx3sig=0.0
        sumy3sig=0.0
        sumx2ysig=0.0
        sumxy2sig=0.0
        muon=0.0
  
        valid = .true.
*mac	do i=1,npmt
*mac	  if (event(i).gt.0.0) then
        do j=1,npicture
           i=pictnum(j)
            xi = rwxdeg(i) - x0
            yi = rwydeg(i) - y0
*           xi = wxdeg(i) - x0
*           yi = wydeg(i) - y0
            si = event(i)
	    wxbyev = xi*si
            wx2byev = xi*wxbyev
	    wybyev = yi*si
            wy2byev = yi*wybyev
	    sumsig = sumsig + si
	    sumxsig = sumxsig + wxbyev
	    sumx2sig = sumx2sig + wx2byev
	    sumysig = sumysig + wybyev
	    sumy2sig = sumy2sig + wy2byev
	    sumxysig = sumxysig + xi*wybyev
*
* The following code is needed to calculate the asymmetry param.
*
            if(fasymm) then
              sumx3sig = sumx3sig + xi*wx2byev		
              sumy3sig = sumy3sig + yi*wy2byev	
              sumx2ysig = sumx2ysig + xi*xi*wybyev
              sumxy2sig = sumxy2sig + yi*yi*wxbyev
            endif
*
*mac	  endif
	enddo


	if (sumsig.le.0.0) then
          valid = .false.
	  icode=99
	  return
	endif

	xmean = sumxsig / sumsig
	x2mean = sumx2sig / sumsig
	ymean = sumysig / sumsig
	y2mean = sumy2sig / sumsig
	xymean = sumxysig / sumsig
	xmean2 = xmean * xmean
	ymean2 = ymean * ymean
	meanxy = xmean * ymean
 
	sdevx2 = x2mean - xmean2
	sdevy2 = y2mean - ymean2
	sdevxy = xymean - meanxy

*
* The following code has been added to calculate the asymmetry param
*
        if(fasymm) then
	  x3mean = sumx3sig /sumsig
          y3mean = sumy3sig / sumsig	
          x2ymean = sumx2ysig / sumsig	
          xy2mean = sumxy2sig / sumsig
          xmean3 = xmean2 * xmean	
          ymean3 = ymean2 * ymean	
          sdevx3 = x3mean-3.0*x2mean*xmean+2.0*xmean3	! Not same as MP! 
	  sdevy3 = y3mean-3.0*y2mean*ymean+2.0*ymean3	! " "
	  sdevx2y = x2ymean-2.0*xymean*xmean+2.0*xmean2*ymean
     &            -x2mean*ymean
	  sdevxy2 = xy2mean-2.0*xymean*ymean+2.0*xmean*ymean2
     &            -xmean*y2mean
         endif
*
*
* Mod. by JB 941114 Calculate muonicity
* Note: parameter name has been changed from "hadron" to "muon: 950410
* since the historical definition of the "hadronicity" parameter is
* different.
* As of 950706 only use muonicity to select out big rings when running
* in interactive mode, or when fmuon flag is set to calculate muon
* gain parameter.
*
        if(fview.or.fmuon.or.fuv) then
*mac	  do i=1,npmt
*mac	    if (event(i).gt.0.0) then
          do j=1,npicture
             i=pictnum(j)
              dx = rwxdeg(i)-xmean
              dy = rwydeg(i)-ymean
*             dx = wxdeg(i)-xmean
*             dy = wydeg(i)-ymean
              dmean2 = dx*dx + dy*dy
              if(dmean2.le.HADR2) then
                muon = muon + event(i)
              endif
*mac	    endif
	  enddo
          muon = muon/sumsig
        endif

 	arg = xmean2 + ymean2
        if(arg.gt.0.0) then
          dist = SQRT(arg)
        else
          dist = 0.0
*         valid = .false.
        endif

*
* Simpler formulae for length & width suggested CA 901019
*
	d = sdevy2 - sdevx2
	arg = d*d + 4.0*sdevxy*sdevxy
	if(arg.gt.0.0) then
	  z = SQRT(arg)
	else
          z = 0.0
        endif
*       length = SQRT((sdevx2+sdevy2+z)/2.0)
*       width = SQRT((sdevy2+sdevx2-z)/2.0)
*
* Mod. by JB 940511 to avoid sqrt of neg. number
*
        arg = (sdevx2+sdevy2+z)/2.0
        if(arg.gt.0.0) then
 	  length = SQRT(arg)
        else
*          write(*,*)'ERROR: cant evaluate length'
          length = 0.0
        endif
*
* Mod. by JB 940511 to avoid sqrt of neg. number
*
 	arga = (sdevy2+sdevx2-z)/2.0
 	if(arga.gt.0.0) then
 	  width = SQRT(arga)
 	else
*         write(*,*)'ERROR: cant evaluate width = sqrt',arga
          width = 0.0 		
*         valid = .false.
        endif

*
* Simpler formula for miss introduced CA, MP 901101
* Revised MP 910112
*
	if(z.gt.0.0)then
	  u=1+d/z
	  v=2-u
*
* Mod. by JB 940511 to avoid sqrt of neg. number
*
          arga = (u*xmean2+v*ymean2)/2.0 - meanxy*(2.0*sdevxy/z)
          if(arga.ge.0.0) then
 	    miss=SQRT(arga)
          else
* 	    write(*,*)'ERROR: cant evaluate miss = sqrt',arg
            miss = 0.0
*           valid = .false.
          end if
	else
	  miss = dist
        endif

*
* Old method of calculating phi
*
*       arga = 2.0*sdevxy*ymean-(d+z)*xmean
*       if(arga.ne.0.0) then
*         axslope = (2.0*sdevxy*xmean+(d+z)*ymean)/arga
*         hyp = sqrt(axslope*axslope+1.0)
*         if(hyp.gt.0.0) then
*           cosphi = 1.0/hyp
*           sinphi = cosphi*axslope
*         else
*           valid = .false.
*           return
*         endif
*       else
*         axslope = 1.0e10      ! infinity
*         cosphi = 0.0
*         sinphi = 1.0
*       endif
*
* New method of calculating phi based on hadronicity code
*
        arga = (d+z)*ymean+2.0*sdevxy*xmean
        if(arga.ge.1.0) then
          arga = 1.0
        else if(arga.le.-1.0) then
          arga = -1.0
        endif
        argb = 2.0*sdevxy*ymean-(d-z)*xmean
        if(argb.ge.1.0) then
          argb = 1.0
        else if(argb.le.-1.0) then
          argb = -1.0
        endif
        if((arga.eq.0.0).and.(argb.eq.0.0)) then
          phi = PI/2.0
        else
          phi = atan2(arga, argb)
        endif
        cosphi = cos(dble(phi))
        sinphi = sin(dble(phi))

*
* Calculate most probable x-y point of origin for gamma-ray
* shower.
*
        if(fxy) then
*
* First tranform the centroid into the "event frame"
*
*         xavg = xmean*cosphi + ymean*sinphi
*         yavg = ymean*cosphi - xmean*sinphi
          if(length.gt.0.0) then
            if(tele_id.eq.TELE_11) then
              displ = abs( 1.63 - 1.63*(width/length))
            else
*             displ = abs(1.87 - 1.87*(width/length)) 
              displ = abs(1.7 - 1.7*(width/length)) 
            endif
*
* Next, calculate the two possible points of origin
*
*           xgam1p = xavg + displ ! two possible pts along the major axis
*           xgam2p = xavg - displ
*           ygam1p = yavg
*           ygam2p = yavg
*
* Now transform these back from the event frame to the camera frame
*
*           xgam1 = xgam1p*cosphi - ygam1p*sinphi 
*           xgam2 = xgam2p*cosphi - ygam2p*sinphi
*           ygam1 = ygam1p*cosphi + xgam1*sinphi
*           ygam2 = ygam2p*cosphi + xgam2*sinphi
*
* Method #2
*
            xgam1p = displ*cosphi
            xgam2p = -displ*cosphi
            ygam1p = displ*sinphi
            ygam2p = -displ*sinphi
            xgam1 = xgam1p + xmean
            ygam1 = ygam1p + ymean
            xgam2 = xgam2p + xmean
            ygam2 = ygam2p + ymean
*
* For one of these two points, recalculate the asymmetry.  At the 
* other point, it must just have the same value but opposite sign.
* Also recalculate alpha, and the azwidth at one of the two points
*
            widthp = width	! Protect variables from being
            lengthp = length 	! changed by rehillakpar
            missp = miss		! but still pass original values
            distp = dist		! in case they are needed.
            xmeanp = xmean
            ymeanp = ymean 
            phip = phi
            icodep = icode
            azwidthp = azwidth
            asymp = asymmetry
            validp = valid
				! Recalculate hillas parameters and
				! DO rotate the coordinates

	    call rehillakpar(event,icodep,npicture,pictnum,
     &           lengthp,widthp,missp,
     &           distp,azwidthp,xmeanp,ymeanp,asymp,phip,validp,
     &           xgam1,ygam1,rotcoords)
*
*           ctheta = cos(theta)
*           stheta = sin(theta)
            if(asymp.gt.0.0) then
*
* Modified 961206 - derotation is already done, don't do it again!
*
*             xgam = xgam1*ctheta-ygam1*stheta
*             ygam = xgam1*stheta+ygam1*ctheta
              xgam = xgam1
              ygam = ygam1
            else
*             xgam = xgam2*ctheta-ygam2*stheta
*             ygam = xgam2*stheta+ygam2*ctheta
              xgam = xgam2
              ygam = ygam2
            endif
          else
            displ = 100.0 ! Set to ridiculous value to flag invalid
            xgam = 10000.0
            ygam = 10000.0
          endif
        endif ! fxy
*
*       Begin hadronicity calculation using formula of
*       Chantell, Urban, et al.
*
*       Note: phi calculation also needed for asymmetry, moved
*       out of this block of code.  sumsig has already been
*       calculated, commented out.

        hadron = 0.0
        if(fuv) then
          xavg = xmean*cosphi + ymean*sinphi
          yavg = ymean*cosphi - xmean*sinphi
*         sumsig = 0
          out = 0
          maj_axis=mrtd*2.6*(10.2-0.34*(dist/mrtd))
          min_axis=mrtd*2.6*(6.6-0.22*(dist/mrtd))
*mac          do i = 1,npmt
*
* Mod. by JB 950706 - only do the following calculation of x,y
* event(i).gt.0 to save time
*
*mac            if(event(i).gt.0.0)then
          do j=1,npicture
             i=pictnum(j)
              x = rwxdeg(i)*cosphi + rwydeg(i)*sinphi
              y = rwydeg(i)*cosphi - rwxdeg(i)*sinphi
              if((maj_axis.gt.0.0).and.(min_axis.gt.0.0)) then
                arg = ((xavg-x)/maj_axis)**2 +
     &                ((yavg-y)/min_axis)**2 
                if(arg.gt.0.0) then
                  rad = SQRT(arg) 
                else
                  rad = 0.0
                endif
              else
                rad = 1.0e6
              endif
              if(rad .gt. 1)then
                out = out + event(i)
              endif
*mac            endif
          enddo
          hadron = out / sumsig
        endif
*
* Begin asymm calculation
*
        asymmetry = 0.0
        if(fasymm) then
          cos2phi = cosphi*cosphi
          sin2phi = sinphi*sinphi
          arg = sdevx3*cosphi*cos2phi+3.0*sdevx2y*cos2phi*sinphi+
     &      3.0*sdevxy2*cosphi*sin2phi+sdevy3*sinphi*sin2phi
          argl = length*length*length
          if(argl.gt.0.0) then
            asymmetry = arg/argl
            if(asymmetry.ge.0.0) then
              sign = 1.0
            else
              sign = -1.0
            endif 
            asymmetry = abs(asymmetry)
            asymmetry = sign*(asymmetry**(0.33333333))
          else
            if (arg.lt.0.0) then
              asymmetry = -1.0e10	! minus infinity
            else 
	      asymmetry = 1.0e10	! plus infinity
            endif
          endif 
        endif
*
* Begin Akerlof azwidth cut, MP 910112
*
* Note that d and z are now defined differently than the values used
* to calculate miss and asymmetry
*
	d = y2mean - x2mean
        arga = d*d+4.0*xymean*xymean
        if(arga.lt.0.0) arga = 0.0
	z = sqrt(arga)
*
* Mod. by JB 940511 to avoid sqrt of neg. number
*
        arg = x2mean+y2mean-z
	if(arg.gt.0.0) then
	  azwidth = SQRT(arg/2.0)
        else
*	  write(*,*)'ERROR: can't evaluate azwidth = sqrt',arg,'/2.0'
          azwidth = 0.0		! TEMPORARY: Should flag invalid event?
        endif
 
	call nmax(event,ppmt,maxsig1,loc1,maxsig2,loc2,maxsig3,loc3)
 
        if(sumsig.gt.0.0) then
	  frac2 = (maxsig1+maxsig2) / sumsig
	  frac3 = (maxsig1+maxsig2+maxsig3) / sumsig
        else
          frac2 = 0.0
          frac3 = 0.0
        end if

	call nmax(event,ppmt,maxsig1,loc1,maxsig2,loc2,maxsig3,loc3)

	isize=nint(sumsig)
 
	return
	end

*---------------------------------------------------------------------
*   
*     Function pass_superd
*     950117
*     J.Buckley
*
*     2-d supercuts plus some ideas from rcl and vc on limiting the
*     dist parameter based on length and width
* 
      logical function pass_superd(length,width,miss,dist,asymmetry,
     +  livetime)

      implicit none

      real*4	length,width,miss,dist,sinalpha,asymmetry,livetime
      real*4	disp
      real*4	DISPC
*
* JB reoptimized VC's value:
*
* 11m
*
*     parameter	(DISPC=0.3)	! FOR 11m ONLY!
*
* 10m
*
*     parameter	(DISPC=0.2125)
*
      real*4 lengthlo,lengthhi
      real*4 widthlo,widthhi
      real*4 distlo,disthi
      real*4 sinalphahi
      real*4 sizelo
      real*4 triglo
      real*4 lwidhi,llenhi
      real*4 azwidhi
      real*4 ltimehi
      common / cuts /
     +  lengthlo,lengthhi,
     +  widthlo,widthhi,
     +  distlo,disthi,
     +  sinalphahi,
     +  sizelo,
     +  triglo,
     +  lwidhi,llenhi,
     +  azwidhi,
     +  ltimehi
      logical fview
      logical fmesh
      logical ffomin
      logical fsuperd
      logical f2d
      logical fxy
      logical fuv
      logical fasymm
      logical floose
      logical fmuon
      common /flags/ fview,fmesh,ffomin,fsuperd,f2d,fxy,fuv,fasymm,
     &               floose,fmuon
*
* Modified by JB 960109 - pass telescope id needed for cuts which
* is needed since cuts depend on which telescope
*
      integer N_TELE
      parameter (N_TELE = 2)
      integer*4 tele_id
      character*3 tele_name(N_TELE)
      common / tele /
     +  tele_id,
     +  tele_name
      integer
     +     TELE_10,
     +     TELE_11
      parameter
     +    (TELE_10 = 1,
     +     TELE_11 = 2)
*
      if(tele_id.eq.TELE_11) then
        DISPC = 0.3
      else
        DISPC = 0.2125
      endif
*

      if(length.gt.0.0) then
*
* JB Reoptimized original RCL cut values
*       disp = abs( dist - 2.0 + 2.0*(width/length) )
        if(tele_id.eq.TELE_11) then
          disp = abs( dist - 1.63 + 1.63*(width/length) )
        else 
*         disp = abs( dist - 1.87 + 1.87*(width/length) ) 
          disp = abs( dist - 1.7 + 1.7*(width/length) ) 
        endif
      else
        disp = 1.0e6
      endif
      if(dist.gt.1.0d-6) then
        sinalpha = (miss/dist)
      else
        sinalpha = 1.0
      endif

      if((disp.lt.DISPC).and.
     &  ((length.gt.lengthlo).and.(length.lt.lengthhi)).and.
     &  ((width.gt.widthlo).and.(width.lt.widthhi)).and.
     &  (livetime.lt.ltimehi).and.
     &  (sinalpha.lt.sinalphahi)) then
*    &  (.not.(fasymm.and.(asymmetry.lt.0.0)))) then
        pass_superd = .true.
      else
        pass_superd = .false.
      endif
      return         
      end 
*---------------------------------------------------------------------
*   
*     Function pass_fscuts
*     950117
*     J.Buckley
*
*     Multiparameter cuts based on those used in the false source
*     method by Fomin et al., 1994, Astroparticle Physics.
*     These cuts are a bit more logical than supercuts from a shower-
*     development point of view and include a dist dependent alpha cut
*     as well as a cut on length/width as a function of distance (an 
*     attempt to build in the expected dependence of the shower shape on
*     impact parameter).  Parameter cuts are combined with a hypersphere-
*     like method.  Michael Punch's asymmetry parameter is
*     also revived and included here in an attempt to eliminate the 
*     two-fold ambiguity in the source direction for a given shower.
*     Also assumed here is an initial supercut-like dist cut, since
*     truncated showers are probably of limited use.
*
      logical function pass_fscuts(length,width,miss,dist,asymmetry)

      implicit none

      real*4	length,width,miss,dist,sinalpha,asymmetry
      real*4	GC,DISPC,LENGTHC,MISSC,DISTC,SUMC  
      parameter	(GC=1.8,DISPC=0.4,LENGTHC=0.5,MISSC=0.23,
     &  DISTC=1.1,SUMC=0.8)
      real*4	disp,misscut,arg,arg1,arg2,sum

      if(width.gt.0.0) then
        disp = abs((length/width)/(max(1.e-6,(1.0+GC*dist)-1.0d0)))
      else
        disp = 1.0e6
      endif
      if(dist.gt.0.0) then
        arg = MISSC/dist
        misscut = MISSC/sqrt(1.0 + arg*arg)
      else
        misscut = 1.0e-6
      endif
      if(misscut.gt.0.) then
        arg1 = miss/misscut
      else
        arg1 = 1.0e6
      endif
      arg2 = disp/DISPC
      sum = arg1*arg1 + arg2*arg2
      if((length.lt.LENGTHC).and.(sum.lt.SUMC)) then
        pass_fscuts = .true.
      else
        pass_fscuts = .false.
      endif
      return         
      end 

*---------------------------------------------------------------------------

      subroutine make_nbrmask(nbr_mask)
*
*
      implicit none

      integer*4 npmt,mpmt,ppmt
      common /num/ npmt,mpmt,ppmt
      integer*4 maxpmts
      parameter(maxpmts=1000)

      integer*4 bits(0:31)
      integer*4 nbr_list(7,maxpmts)
      common /neigh/ nbr_list
      integer*4 nbr_mask(0:23,maxpmts)
      integer*4 i,j,index1,index2
      common /bytes/ bits
 
      do j = 1,npmt
         do i = 0,23
            nbr_mask(i,j) = 0
         end do
      end do

      do i = 1,npmt
         do j = 1,7
            if (nbr_list(j,i) .ne. 0) then
               index1 = (nbr_list(j,i) - 1) / 32
               index2 = mod(nbr_list(j,i)-1,32)
               nbr_mask(index1,i) =
     &            nbr_mask(index1,i) .or. bits(index2)
            end if
         end do
      end do

      return
      end

*---------------------------------------------------------------------------
*
*  GR_NTUPLE
*  JB
*  940511
*
*  Write a PAW ntuple file of image parameters for each event.
*  Based on granite.for by Joachim Rose
*
        subroutine gr_ntuple(chopt)

	character*(*) chopt

        character telsrc*2
        integer N_TELE
        parameter (N_TELE = 2)
        integer N_PARAM
        parameter (N_PARAM = 38)
        integer N_TIME
        parameter (N_TIME = 4)

        integer*4 tele_id
        character*3 tele_name(N_TELE)
        common / tele /
     +    tele_id,
     +    tele_name

        logical      gr_param_valid(N_TELE)          ! data valid
        real gr_param_value(N_PARAM,N_TELE)          ! actual value
        character*8  gr_param_name (N_PARAM)         ! parameter names
      
        logical      gr_time_valid(N_TIME)
        real*8 gr_time_value(N_TIME,N_TELE) 
        character*8  gr_time_name(N_TIME)

        common / gr_param /
     +    gr_param_valid,
     +    gr_param_value,
     +    gr_param_name,
     +    gr_time_value,
     +    gr_time_name

        integer        nt_unit	! Ntuple file logical descripter 
        integer	       nt_ierr
        character*8    nt_dir	! Ntuple directory
        character*10   nt_topdir
        character*80   nt_fname	! Ntuple file name

        common / nt /
     +     nt_unit,
     +     nt_ierr,
     +     nt_dir,
     +     nt_topdir,
     +     nt_fname

        integer	NT_ID
        parameter (NT_ID = 100)

        integer   PAW_NWORDS
        parameter (PAW_NWORDS = 2*1024*256)       ! 2 MB
        integer   h(PAW_NWORDS)
        common/ PAWC /
     +     h

        integer
     +     TELE_10,
     +     TELE_11
*
        parameter
     +    (TELE_10 = 1,
     +     TELE_11 = 2)
*
        logical HEXIST 
*---------------------------------------------------------------------------
*     open and initialize ntuple file
*---------------------------------------------------------------------------
	if (index(chopt,'N').ge.1) then
         call HROPEN(                ! open RZ file
     +    nt_unit,                   ! logical unit
     +    nt_dir,                    ! top dirctory name
     +    nt_fname,                  ! file name
     +    'N',1024,nt_ierr)          ! new file
        
         if (nt_ierr.ne.0) then
          write(*,*)'*** Failed to open n-tuple file.'
	  return
	 end if
	 if (HEXIST(NT_ID)) then     ! previous n-tuple?
          call HDELET(NT_ID)         ! delete it
         end if
         call HBNT(NT_ID,'GRAS',' ') ! create a new one

*
* Modified by JB 941020 -- Only include entries in ntuple for telescopes
* which are operating.  This is a self-descriptive data format after all!
*
*	 do i=1,N_TELE
         i = tele_id 
	  call HBNAME(               ! define parameter valid
     +      NT_ID,
     +      nt_dir,
     +      gr_param_valid(i),
     +      tele_name(i)//'ipv:L')

          do j=1,N_PARAM             ! for all parameters
            call HBNAME(             ! define parameter values
     +        NT_ID,
     +        nt_dir,
     +        gr_param_value(j,i),
     +        tele_name(i)//gr_param_name(j)//
     +        ':R')
          end do		     ! for all params.

          call HBNAME(               ! define time valid
     +      NT_ID,
     +      nt_dir,
     +      gr_time_valid(i),
     +      tele_name(i)//'tiv'//
     +      ':L')
          do j=1,N_TIME              ! for all time parameters
            call HBNAME(             ! define parameter values
     +        NT_ID,
     +        nt_dir,
     +        gr_time_value(j,i),
     +        tele_name(i)//gr_time_name(j)//
     +        ':R*8')
          end do		     ! for all time params.
*	 end do			     ! for all teles.

         call HCDIR('//PAWC',' ')    ! back to default
*---------------------------------------------------------------------------
*     write data
*---------------------------------------------------------------------------
      else if (index(chopt,'W').ge.1) then
        call HCDIR(nt_topdir,' ')
        call HFNT(NT_ID)
        call HCDIR('//PAWC',' ')
*---------------------------------------------------------------------------
*     close file
*---------------------------------------------------------------------------
      else if (index(chopt,'C').ge.1) then
        call HCDIR(nt_topdir,' ')        ! change directory
*       write(*,*)'hrout ntuple'
        call HROUT(NT_ID,ICYCLE,' ')     ! store n-tuple in Rz file
        call HREND(nt_dir)               ! close file
        close(nt_unit)                   ! close logical unit
        call HCDIR('//PAWC',' ')
*---------------------------------------------------------------------------
*     print info?
*---------------------------------------------------------------------------
      else if (index(chopt,'P').ge.1) then
	write(*,*)'ntuple summary:'
        call HCDIR(nt_topdir,' ')
        if (HEXIST(NT_ID)) then
          call HPRNT(NT_ID)
        else
*         write(*,*) 'No n-tuple defined.'
        end if
        call HCDIR('//PAWC',' ')	! Added by JB 950523
*---------------------------------------------------------------------------
*     print invalid option
*---------------------------------------------------------------------------
      else
        write(*,*)'Invalid option(s): '//CHOPT
      end if

      return
      end 

*---------------------------------------------------------------------------

      subroutine make_nbrmask2(nbr_mask)
*
* ** Make mask of neighboring tubes from look-up table.  This
* ** version is for a partial outer ring of 18 one-inch tubes.
* ** The code is identical to make_nbrmask; only the data statment
* ** setting nbr_list is different.
*                                         ADK 940323
*
      implicit none

      integer*4 npmt,mpmt,ppmt
      common /num/ npmt,mpmt,ppmt
      integer*4 maxpmts
      parameter(maxpmts=1000) 

      integer*4 bits(0:31)
      integer*4 nbr_list(7,maxpmts)
      common /neigh/ nbr_list
      integer*4 nbr_mask(0:23,maxpmts)
      integer*4 i,j,index1,index2

      common /bytes/ bits
*
***************************************************************************
*
      do j = 1,npmt
         do i = 0,3
            nbr_mask(i,j) = 0
         end do
      end do

      do i = 1,npmt
         do j = 1,6
            if (nbr_list(j,i) .ne. 0) then
               index1 = (nbr_list(j,i) - 1) / 32
               index2 = mod(nbr_list(j,i)-1,32)
               nbr_mask(index1,i) =
     &            nbr_mask(index1,i) .or. bits(index2)
            end if
         end do
      end do

      return
      end

*-----------------------------------------------------------------------
*-- GV @ FLWO		SHOW109		! NOV 17, 1988 !
* 
* Modified by JB 941020 -- changed from real to real*4 for compatibility
* with gparamdat.f
*
*-----------------------------------------------------------------------
*-- Displays REAL HRC events in an hexagonal pattern.
*-- A few changes from PWK version (rtpmt109).
*-- Numbering convention (zones 0, 1 and 2):
*--
*--            16    17    18
*--
*--          15     6     7    19
*--
*--      14     5     1     2     8                        
*--
*--         13     4     3     9                  
*--
*--            12    11    10                        
*-----------------------------------------------------------------------

	subroutine r4show109(iunit,event)

	integer*4 maxpmts
	parameter (maxpmts=1000)
	real*4 event(maxpmts)
        integer iunit

	write(iunit,601)
     +	event(104),event(105),event(106),event(107)
	write(iunit,602)
     +	event(82),event(83),event(84),event(85),event(86),event(87)
	write(iunit,603)
     +	event(103),event(81),event(54),event(55),event(56),event(57),
     +	event(58),event(88),event(108)
	write(iunit,604)
     +	event(80),event(53),event(32),event(33),event(34),event(35),
     +	event(59),event(89)
	write(iunit,605)
     +	event(102),event(79),event(52),event(31),event(16),event(17),
     +	event(18),event(36),event(60),event(90),event(109)
	write(iunit,606)
     +	event(78),event(51),event(30),event(15),event(6),event(7),
     +	event(19),event(37),event(61),event(91)
	write(iunit,607)
     +	event(101),event(77),event(50),event(29),event(14),event(5),
     +	event(1),event(2),event(8),event(20),event(38),event(62),
     +	event(92)
	write(iunit,606)
     +	event(76),event(49),event(28),event(13),event(4),event(3),
     +	event(9),event(21),event(39),event(63)
	write(iunit,605)
     +	event(100),event(75),event(48),event(27),event(12),event(11),
     +	event(10),event(22),event(40),event(64),event(93)
	write(iunit,604)
     +	event(74),event(47),event(26),event(25),event(24),event(23),
     +	event(41),event(65)
	write(iunit,603)
     +	event(99),event(73),event(46),event(45),event(44),event(43),
     +	event(42),event(66),event(94)
	write(iunit,602)
     +	event(72),event(71),event(70),event(69),event(68),event(67)
	write(iunit,608)
     +	event(98),event(97),event(96),event(95)
601	format(t19,4(f5.1,6x),/)
602	format(t21,6(f5.1,1x),/)
603	format(t13,9(f5.1,1x),/)
604	format(t16,8(f5.1,1x),/)
605	format(t7,11(f5.1,1x),/)
606	format(t10,10(f5.1,1x),/)
607	format(t1,13(f5.1,1x),/)
608	format(t19,4(f5.1,6x))

	return
	end
*-----------------------------------------------------------------------
*-- GV @ FLWO		SHOW109		! NOV 17, 1988 !
*-----------------------------------------------------------------------
*-- Displays HRC events in an hexagonal pattern.
*-- A few changes from PWK version (rtpmt109).
*-- Numbering convention (zones 0, 1 and 2):
*--
*--            16    17    18
*--
*--          15     6     7    19
*--
*--      14     5     1     2     8                        
*--
*--         13     4     3     9                  
*--
*--            12    11    10                        
*-----------------------------------------------------------------------

	subroutine  show109(iunit,event)
*
* Modified by JB 941024 -- changed to int*4
*
	integer*4 event(120)
	write(iunit,601)
     +	event(104),event(105),event(106),event(107)
	write(iunit,602)
     +	event(82),event(83),event(84),event(85),event(86),event(87)
	write(iunit,603)
     +	event(103),event(81),event(54),event(55),event(56),event(57),
     +	event(58),event(88),event(108)
	write(iunit,604)
     +	event(80),event(53),event(32),event(33),event(34),event(35),
     +	event(59),event(89)
	write(iunit,605)
     +	event(102),event(79),event(52),event(31),event(16),event(17),
     +	event(18),event(36),event(60),event(90),event(109)
	write(iunit,606)
     +	event(78),event(51),event(30),event(15),event(6),event(7),
     +	event(19),event(37),event(61),event(91)
	write(iunit,607)
     +	event(101),event(77),event(50),event(29),event(14),event(5),
     +	event(1),event(2),event(8),event(20),event(38),event(62),
     +	event(92)
	write(iunit,606)
     +	event(76),event(49),event(28),event(13),event(4),event(3),
     +	event(9),event(21),event(39),event(63)
	write(iunit,605)
     +	event(100),event(75),event(48),event(27),event(12),event(11),
     +	event(10),event(22),event(40),event(64),event(93)
	write(iunit,604)
     +	event(74),event(47),event(26),event(25),event(24),event(23),
     +	event(41),event(65)
	write(iunit,603)
     +	event(99),event(73),event(46),event(45),event(44),event(43),
     +	event(42),event(66),event(94)
	write(iunit,602)
     +	event(72),event(71),event(70),event(69),event(68),event(67)
	write(iunit,608)
     +	event(98),event(97),event(96),event(95)
601	format(t19,4(i4,8x),/)
602	format(t21,6(i4,2x),/)
603	format(t13,9(i4,2x),/)
604	format(t16,8(i4,2x),/)
605	format(t7,11(i4,2x),/)
606	format(t10,10(i4,2x),/)
607	format(t1,13(i4,2x),/)
608	format(t19,4(i4,8x))

	return
	end

*-----------------------------------------------------------------------
*
* DISP_EVENT.F
* JQ 941116
*
* Displays event in x-window using hbook, higz
*

      subroutine disp_event(event,xmean,ymean,length,width,phi,
     &  xgam,ygam,nextmu,xtrip,ytrip,ntrips,xave,yave,rave,
     &  npicture,pictnum)

      implicit none

      integer NPAWC
      parameter (NPAWC=2*1024*256)
      integer*4 npmt,mpmt,ppmt
      common /num/ npmt,mpmt,ppmt
      integer*4 maxpmts
      parameter(maxpmts=1000)
      integer   NXGRID,NYGRID			! JB 950414
      parameter (NXGRID=30, NYGRID=30)		! JB 950414
*     parameter (NXGRID=21, NYGRID=21)		! JB 950414

      integer i,max_adc,min_adc,NT,IBN,j
      integer iver,event_int(maxpmts)
      integer npicture,pictnum(maxpmts)

      real radius1,radius2,event(maxpmts),evmax,sevent(maxpmts)
      real*4 xmean,ymean,length,width,phi
      real*4 xave,yave,rave
      real*4 xgam,ygam
      logical nextmu,inpict
      real XNDC,YNDC,XWC,YWC
      real*4    rwxdeg(maxpmts),rwydeg(maxpmts)	! JB 950414
      real*4	wxdeg(maxpmts),wydeg(maxpmts)
      common /rwcoords/rwxdeg,rwydeg		! JB 950414
      common /wcoords/wxdeg,wydeg
      character event_ch(maxpmts)*3

      integer H(NPAWC)
      common/PAWC/H

      integer	ix,iy,npts
      real*4	xpts((NXGRID+NYGRID)*2+3)
      real*4    ypts((NXGRID+NYGRID)*2+3)
      integer	npm
      real*4	xpm(1),ypm(1)
      data	npm /1/

      real*4    xmesh(NXGRID),ymesh(NYGRID)	! JB 950414	
      real*4    nmesh(NXGRID,NYGRID)		! JB 950414
      integer   tmesh(NXGRID,NYGRID)		! JB 950414
      character chmesh(NXGRID,NYGRID)		! JB 950414
      real*4	xtrip(200),ytrip(200)
      integer 	ntrips
      common /meshv/ xmesh,ymesh,nmesh,tmesh,chmesh	! JB 950414
      logical fview
      logical fmesh
      logical ffomin
      logical fsuperd
      logical f2d
      logical fxy
      logical fuv
      logical fasymm
      logical floose
      logical fmuon
      common /flags/ fview,fmesh,ffomin,fsuperd,f2d,fxy,fuv,fasymm,
     &               floose,fmuon

*     call getcoords(941116,pmtx,pmty,iver)

      call ISFACI(6)
      call ISFAIS(1)
      call ISTXAL(2,3)
      call ISCHH(.06)

      radius1=.059
      radius2=.115
      evmax=0.

      do i=1,npmt
         if (event(i).gt.evmax)evmax=event(i)
         event_int(i)=event(i)
         if (event_int(i).gt.999) event_int(i)=999
         if (event_int(i).lt.-99) event_int(i)=-99
         call IZITOC(event_int(i),event_ch(i))
      end do

      do i=1,379
         sevent(i)=event(i)*(radius1-.01)/max(1.e-20,evmax)
      end do
      do i=380,490
         sevent(i)=event(i)*(radius2-.01)/max(1.e-20,evmax)
      end do

      call HPLFRA(-1.7,1.7,-1.7,1.7,' ')
C      call HPLFRA(-2.4,2.4,-2.4,2.4,' ')
      if(f2d.and.floose.and.(.not. nextmu)) then
       do i=1,379
         inpict=.false.
         do j=1,npicture
            if (i.eq.pictnum(j)) inpict=.true.
         enddo
         if (inpict) then
            call ISPLCI(3)	! Set line color to light gray
            call IGARC(rwxdeg(i),rwydeg(i),radius1,radius1,0.,360.)
            call ISPLCI(1)
         else
            call IGARC(rwxdeg(i),rwydeg(i),radius1,radius1,0.,360.)
         endif
         call IGARC(rwxdeg(i),rwydeg(i),0,sevent(i),0.,360.)
         call ITX(rwxdeg(i),rwydeg(i),event_ch(i))
       end do
       do i=380,npmt
         call IGARC(rwxdeg(i),rwydeg(i),radius2,radius2,0.,360.)
         call IGARC(rwxdeg(i),rwydeg(i),0,sevent(i),0.,360.)
         call ITX(rwxdeg(i),rwydeg(i),event_ch(i))
       end do
      else
       do i=1,379
         inpict=.false.
         do j=1,npicture
            if (i.eq.pictnum(j)) inpict=.true.
         enddo
         if (inpict) then
            call ISPLCI(3)	! Set line color to light gray
            call IGARC(wxdeg(i),wydeg(i),radius1,radius1,0.,360.)
            call ISPLCI(1)
         else
            call IGARC(wxdeg(i),wydeg(i),radius1,radius1,0.,360.)
         endif
         call IGARC(wxdeg(i),wydeg(i),0,sevent(i),0.,360.)
         call ITX(wxdeg(i),wydeg(i),event_ch(i))
       end do
       do i=380,npmt
         call IGARC(wxdeg(i),wydeg(i),radius2,radius2,0.,360.)
         call IGARC(wxdeg(i),wydeg(i),0,sevent(i),0.,360.)
         call ITX(wxdeg(i),wydeg(i),event_ch(i))
       end do
      endif

      if((nextmu.or.fmuon).and.(.not.(floose.and.f2d))) then
        call ISMKSC(0.7)  ! Set marker scale factor
        call ISMK(24)
        do ix=1,ntrips    
          xpm(1) = xtrip(ix)
          ypm(1) = ytrip(ix)
          npm=1
          call IPM(npm,xpm,ypm)
        end do

        call ISMKSC(1.5)
        call ISMK(28)
        xpm(1) = xave
        ypm(1) = yave
        npm=1
        call IPM(npm,xpm,ypm)
        call IGARC(xave,yave,rave,rave,0.,360.)
      endif

      if(f2d.and.floose.and.(.not. nextmu)) then
        call ISPLCI(7)	! Set line color to light gray
        npts = 1
        do ix=1,NXGRID
          xpts(npts) = xmesh(ix)
          xpts(npts+1) = xmesh(ix)
          ypts(npts) = ymesh(1)
          ypts(npts+1) = ymesh(NYGRID)
          npts = npts + 2
        end do
        do iy=1,NYGRID
          ypts(npts) = ymesh(iy)
          ypts(npts+1) = ymesh(iy)
          xpts(npts) = xmesh(1)
          xpts(npts+1) = xmesh(NXGRID)
          npts = npts + 2
        end do
        npts = npts - 1
        call IML(npts,xpts,ypts)
        call ISPLCI(1)	! Set line color back to black
     
*       call ISPMCI(2)    ! Set the polymarker color index 
                       
        if(fxy) then
          xpm(1) = xgam     ! Put a mark at the most likely shower
          ypm(1) = ygam     ! origin position.
          npm=1
          call ISMKSC(3.5)  ! Set marker scale factor
          call ISMK(28)
          call IPM(npm,xpm,ypm)
        endif

        call ISMKSC(1.7)  ! Set marker scale factor
        do ix=1,NXGRID
          do iy=1,NYGRID
            if(tmesh(ix,iy).eq.2) then
              call ISMK(28)	! Set polymark to big plus
	      xpm(1) = xmesh(ix)
              ypm(1) = ymesh(iy)
              npm = 1
              call IPM(npm,xpm,ypm)
            else if(tmesh(ix,iy).eq.1) then
              call ISMK(24)	! Set polymark to big plus
	      xpm(1) = xmesh(ix)
              ypm(1) = ymesh(iy)
              npm = 1
              call IPM(npm,xpm,ypm)
            endif
          end do
        end do
*       call ISPMCI(1)    ! Set the polymarker color index 
        call ISMKSC(1.)  ! Set marker scale factor
      endif
      call IUWK(0,1)

      return
      end


*-----------------------------------------------------------------------
*  Subroutine DEROT   
*  950604
*  MP,JB
*
* Based on "position" and UT to ST conversion in Duffett-Smith
*
* uttime	decimal UT (hr. fraction after UT midnight)
* theta 	angle by which FoV rotates
*               i.e. angle by which to de-rotate.
* alt		returned value of telescope alt (radians)
* azi		returned value of telescope azimuth (radians)
* stime		returned value of sidereal time
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

	subroutine derot(uttime,theta,alt,azi,stime)

        implicit none

        real uttime,theta
	real rra,rdec,str,hr,alt,azi
	real tt,tt0,stime
        real arg, arga, argb
        real sk,ck
	real*8 mjd,dmjd
	common /sinfo/rra,rdec,mjd
        real convd,pi,rlat,titor
	data convd/57.29577951/,pi/3.141592654/,rlat/.553065751136/,
     +  titor/.261799387/

        dmjd = aint(mjd)	! Added by JB to truncate possible 
                                ! fractional mjd
*
* convert decimal UT to decimal ST
*	tt=(mjd-51544.5)/36525.0
	tt=(dmjd-51544.5)/36525.0 	! JB 950411
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
*       write(*,*)'Hr angle (rad)',hr
*       write(*,*)'Position (rad)',rra,rdec

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
*	azi=atan2(-cos(rdec)*cos(rlat)*sin(hr),
*    +	 (sin(rdec)-sin(rlat)*sin(alt)))
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
*-----------------------------------------------------------------------
*  Subroutine FINDARC   
*  950604
*  JB
*
*  Try to find muon arcs using the procedure of Altherr and Seixas
*  (NIM, A317 (1992) 335-338) but no neural net!
*
*  The algorithm goes as follows:
*  For each unique triplet of PMTs (N*(N-1)*(N-2)/6 in total)
*  calculate the center and radius of the uniquely defined circle.
*  OR the bits for PMTs which are elements of this triplet into 
*  the appropriate element (determined by the ring center) of a 
*  2-d array of bitmasks, and increment the corresponding element
*  of the array ntriplets(x,y). Also increment the cumulative ring
*  radius and cumulative x and y-center position arrays rad(x,y)
*  xsum(x,y), ysum(x,y).  Divide xsum, ysum, radsum by ntriplets in
*  the end.  Then do a boxcar average over
*  neighboring elements in the array ntriplets, and x,y,rad arrays.
*  Find the maximum of the boxcar averaged array ntripbox(x,y)
*  and calculate the values xave,yave,rave from these.  Calculate
*  a bitmask (for which tubes are in the arc) by ORing together the
*  bitmasks of the peak element and its neighbors.  From the
*  resulting bitmask calculate the signal sum for these pmts
*  and divide by the number of pmts and the radius squared.  This
*  gives a value proportional to the pe/dc ratio for this arc.  
*  A comparison of this value to the value derived for big rings
*  (from the hadronicity) can be combined with limits on the ring
*  radius to derive the liklihood that this event is in fact a muon
*  arc.
*
*
	subroutine findarc(event,rave,xave,yave,arc,gain,mugain,soal,
     &   muskew,arclen,npict,pictpmt,pictx,picty,pictnum,xtrip,
     &   ytrip,npts)

	implicit none

	real*4 DX,DY
        parameter(DX=0.25)
        parameter(DY=0.25)
        integer NX
        integer NY
        parameter(NX=24)
        parameter(NY=24)
        integer NXD2
        integer NYD2
        parameter(NXD2=12)
        parameter(NYD2=12)
        real*4 RSQMAX
        parameter(RSQMAX=4.0)
        real*4 RSATSQ
        parameter(RSATSQ=1.31312)
	real*4 RSAT
        parameter(RSAT=1.1459)
	real*4 PIXSPACE
        parameter(PIXSPACE=0.25)
	real*4 DRAD
	parameter(DRAD=0.26)
        real*4 PI
        parameter(PI=3.141592654)
        real*4 TWOPI
        parameter(TWOPI=6.283)
        real*4 CONVD
        parameter(CONVD=57.29577951)
        real*4 SQRT12
        parameter(SQRT12=3.464)
        integer*4 npmt,mpmt,ppmt
        common /num/ npmt,mpmt,ppmt
	integer*4 maxpmts
	parameter(maxpmts=1000)
	real*4 event(maxpmts)
        real*4 arc
        real*4 denom,numx,numy
        real*4 pictpmt(maxpmts)
        real*4 pictx(maxpmts),picty(maxpmts)
        real*4 pictxsq(maxpmts),pictysq(maxpmts)
        integer pictnum(maxpmts)
        integer*4 npict
        real*4 x0,y0,x1,x2,x3,y1,y2,y3
        real*4 x1sq,x2sq,x3sq,y1sq,y2sq,y3sq
	integer i,j,k,l,m,i1,i2,i3,n,ix,iy
        integer l1,l2,l3,m1,m2,m3
        integer ntrip(NX,NY)
        integer nbox
        integer ixmax,iymax,nmax	
        real*4 dnmax
        real*4 xsum(NX,NY)
        real*4 ysum(NX,NY)
        real*4 rsum(NX,NY)
        real*4 xave,yave,rave,rsq
        real*4 phsum,gain
        integer nph
        integer*4 bits(0:31),pmtmask(NX,NY,0:23)
        integer*4 totpmtmask(0:23)
        common /bytes/bits
        logical fview
        logical fmesh
        logical ffomin
        logical fsuperd
        logical f2d
        logical fxy
        logical fuv
        logical fasymm
        logical floose
        logical fmuon
        common /flags/ fview,fmesh,ffomin,fsuperd,f2d,fxy,fuv,fasymm,
     &               floose,fmuon
        real*4	xtrip(200),ytrip(200)	! Save points for display only
        integer npts
        real*4 grevent(maxpmts)
        common /wpict/ grevent 
        real*4 wxdeg(maxpmts),wydeg(maxpmts)
        common /wcoords/ wxdeg,wydeg
	real*4 r1,r2,r1sq,r2sq,rpsq 
	real*4 dxa,dya,stot,sarc,mugain,soal,soal2
	real*4 arclen,wsum,wxphi,wphi,wphisq,wt,phi,philo,phihi
	real*4 aanulus
        real*4 apix	! Area of a hexagonal pixel
	real*4 arcfrac1,arcfrac2
	real*4 xcs,ycs,dphi
	real*4 muskew
        integer	npix,narc
	integer fint

	apix = 0.866 * PIXSPACE*PIXSPACE
	fint = 0	! Set equal to zero to eliminate integration
			! over annulus.
  
        if(npict.lt.3) then
          return
        endif
        npts = 1
	do i=1,NX
          do j=1,NY
            ntrip(i,j) = 0
            xsum(i,j) = 0.0
            ysum(i,j) = 0.0
            rsum(i,j) = 0.0
          end do
        end do  
	do i=1,npict	! could save some time here by doing this once
			! in main program
          pictxsq(i) = pictx(i)*pictx(i)
          pictysq(i) = picty(i)*picty(i)
        end do
        n=0

*
* Loop through the N*(N-1)*(N-2)/(3*2*1) triplets, where N=npict
*
        do i1=1,npict
          l1=mod(pictnum(i1)-1,32)	! Calculate the bit index
          m1=(pictnum(i1)-1)/32		! and word index of this pmt

          do i2=i1+1,npict
            l2=mod(pictnum(i2)-1,32)
            m2=(pictnum(i2)-1)/32

            do i3=i2+1,npict
              l3=mod(pictnum(i3)-1,32)
              m3=(pictnum(i3)-1)/32
              n=n+1

              x1 = pictx(i1)
              x2 = pictx(i2)
              x3 = pictx(i3)
              y1 = picty(i1)
              y2 = picty(i2)
              y3 = picty(i3)

              x1sq = pictxsq(i1)
              x2sq = pictxsq(i2)
              x3sq = pictxsq(i3)
              y1sq = pictysq(i1)
              y2sq = pictysq(i2)
              y3sq = pictysq(i3)

*
* Only calculate the little bit you need to get to the
* next condition - there's a good chance you won't make it through
* the whole mess.
*
              denom = 2.0*((x2-x3)*(y1-y2)-(x1-x2)*(y2-y3))
              if(abs(denom).gt.1.0e-10) then
                numx = ((x2sq-x3sq+y2sq-y3sq)*(y1-y2)
     &                 -(x1sq-x2sq+y1sq-y2sq)*(y2-y3))
                x0 = numx/denom
                ix = NXD2 + nint(x0/DX)
                if((ix.gt.0).and.(ix.le.NX)) then
                  numy = ((x1sq-x2sq+y1sq-y2sq)*(x2-x3)
     &                   -(x2sq-x3sq+y2sq-y3sq)*(x1-x2))
                  y0 = numy/denom
                  iy = NYD2 + nint(y0/DY)
                  if((iy.gt.0).and.(iy.le.NY)) then
                    rsq = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)
                    if((rsq.gt.0.0).and.(rsq.lt.RSQMAX)) then
                      ntrip(ix,iy) = ntrip(ix,iy) + 1
*
* Accumulate averages at the
* appropriate bin
*
                      rsum(ix,iy) = rsum(ix,iy) + sqrt(rsq)
                      if(fview) then
                        xsum(ix,iy) = xsum(ix,iy) + x0
                        ysum(ix,iy) = ysum(ix,iy) + y0
                        if(npts.lt.200) then
                          xtrip(npts) = x0
                          ytrip(npts) = y0
                          npts = npts + 1 
                        endif
                      endif
*
* Set the bits in pmtmask corresponding to the pmts indexed by
* i1,i2,i3
*
                      pmtmask(ix,iy,m1)=pmtmask(ix,iy,m1).or.bits(l1)
                      pmtmask(ix,iy,m2)=pmtmask(ix,iy,m2).or.bits(l2)
                      pmtmask(ix,iy,m3)=pmtmask(ix,iy,m3).or.bits(l3)
                    end if ! rsq
                  end if ! iy
                end if ! ix
              endif ! denom
            end do ! i3
          end do ! i2
	end do ! i1
        nmax = 0
        do ix=1,NX-1
          do iy=1,NY-1
            nbox = ntrip(ix,iy)+ntrip(ix+1,iy)+ntrip(ix+1,iy+1)+
     &             ntrip(ix,iy+1)
            if(nbox.gt.nmax) then
              nmax = nbox
              ixmax = ix
              iymax = iy
            end if
          end do
        end do
        if(nmax.gt.0) then
          dnmax = float(nmax)
          if(fview) then
            xave = (xsum(ixmax,iymax)+xsum(ixmax+1,iymax)+
     &             xsum(ixmax+1,iymax+1)+xsum(ixmax,iymax+1))/dnmax
            yave = (ysum(ixmax,iymax)+ysum(ixmax+1,iymax)+
     &             ysum(ixmax+1,iymax+1)+ysum(ixmax,iymax+1))/dnmax
          endif
          rave = (rsum(ixmax,iymax)+rsum(ixmax+1,iymax)+
     &           rsum(ixmax+1,iymax+1)+rsum(ixmax,iymax+1))/dnmax
          do i=0,3
            totpmtmask(i)=pmtmask(ixmax,iymax,i).or.
     &                    pmtmask(ixmax+1,iymax,i).or.
     &                    pmtmask(ixmax+1,iymax+1,i).or.
     &                    pmtmask(ixmax,iymax+1,i)
          end do
          phsum = 0.0
          nph = 0
          do i=1,npict
            k=pictnum(i)
            l=mod(k-1,32)	
            m=(k-1)/32
*
* If tube has contributed to the arc, include its pulseheight in the sum
*
            if((bits(l).and.totpmtmask(m)).gt.0) then
              phsum = phsum + pictpmt(i)
              nph = nph + 1
            endif
          end do
          if(npict.gt.3) then
            arc = 6.0*dnmax/float(npict*(npict-1)*(npict-2))
          else
            arc = 0.0
          endif
          if((rave.gt.0.0).and.(nph.gt.0)) then
*
* Modified 961009 by JB - Acutually while the light per unit
* pathlength is proportional to theta^2 or r^2, the pathlength
* increases as theta^-1 or r^-1, and the light per pixel is proportional
* to the light per azimuth angle divided by r, so actually the gain should
* does not depend on radius to first order.
*
*           gain = RSATSQ*phsum/(rave*rave*nph)
            gain = phsum/float(nph)
*
* Added 961104 by JB - Calculate amount of light in an annulus about
* the ring center, approx. correct for fraction of area covered by
* pixels and for the linear dependence of the total amount of light
* on radius.
*
	    if (fint.eq.1) then
              r2 = rave + DRAD
              r1 = rave - DRAD
              if (r1.lt.0.0) then
                r1 = 0.0
              endif
              r1sq = r1*r1
              r2sq = r2*r2
	      aanulus = PI*(r2sq-r1sq)
*
* First calculate arclength...
*
              arclen = 0.0
              wsum = 0.0
	      xcs = 0.0
              ycs = 0.0
              do i=1,npict	! For all pixels contrib to muon
                k=pictnum(i)	! arc, calculate the centroid of
                l=mod(k-1,32)	! the light distribution rel. to the
                m=(k-1)/32	! ring center position
                if((bits(l).and.totpmtmask(m)).gt.0) then
                  dxa = pictx(i) - xave
                  dya = picty(i) - yave
                  wt = pictpmt(i)
                  wsum = wsum + wt
                  xcs = xcs + dxa*wt
                  ycs = ycs + dya*wt
                end if 
              end do
              if (wsum.gt.0.0) then
                xcs = xcs/wsum	! Calculate angle of vector from
                ycs = ycs/wsum	! ring center to centroid of
	        wphi = PI+atan2(ycs,xcs)	! light distribution in ring
                muskew = sqrt(xcs*xcs+ycs*ycs)/rave
                wphisq = 0.0
                do i=1,npict	! For all pixels contrib to muon
                  k=pictnum(i)	! arc, calculate the rms arclen
                  l=mod(k-1,32) ! w.r.t. the centroid angle.
                  m=(k-1)/32	
                  if((bits(l).and.totpmtmask(m)).gt.0) then
                    dxa = pictx(i) - xave
                    dya = picty(i) - yave
                    phi = PI+atan2(dya,dxa)
                    wt = pictpmt(i)
                    dphi = amod(abs(phi-wphi),TWOPI)
                    dphi = min(dphi,TWOPI-dphi)
*                   write(*,*)'wt,phi,wphi,phi-wphi,dphi: ',
*    &               int(wt),phi*CONVD,wphi*CONVD,(phi-wphi)*CONVD,
*    &               dphi*CONVD
                    wphisq = wphisq + (wt*dphi*dphi)
                  end if 
                end do	! for all pixels in the arc
                wphisq = wphisq/wsum
	        wphisq = sqrt(wphisq)
                arclen = wphisq
                philo = wphi - arclen
                phihi = wphi + arclen
              else
	        wphi = 0.0
                philo = 0.0
                phihi = 0.0
                arclen = 1.0e20	! CHECK THIS
              endif
*
* Now calculate light in arc and correct for the fraction
* lost over the edge of the camera
*
              npix = 0
	      narc = 0
              stot = 0.0
              sarc = 0.0
              do i=1,npmt
                dxa = (wxdeg(i) - xave)
                dya = (wydeg(i) - yave)
                rpsq = dxa*dxa + dya*dya
                if ((rpsq.gt.r1sq).and.(rpsq.lt.r2sq)) then
                  phi = PI+atan2(dya,dxa)
                  dphi = amod(abs(phi-wphi),TWOPI)
                  dphi = min(dphi,TWOPI-dphi)
                  if (abs(dphi).lt.arclen) then
                    narc = narc + 1
                    sarc = sarc + grevent(i)
                  endif
                  npix = npix + 1
                  stot = stot + grevent(i)
                endif
              end do
              if (npix.gt.0) then
                arcfrac1 = (apix * float(npix))/max(1.e-20,aanulus)
                mugain = (1.0/36.0)*stot * (RSAT/rave) / 
     >               max(1.e-20,arcfrac1)
                if(arclen.gt.1.0e-20) then
*
* Convert from rms arclen to equivalent interval or sqrt(12)*rms
* also convert from radians to deg.
*
	          arclen = arclen * CONVD * SQRT12
                  soal2 = (stot/arclen)*(RSAT/rave) 
*
* arcfrac2 is the fraction of the area in the arc from philo to phihi
* which falls on the camera.  Note that since dN/dphi scales like
* the radius, a correction for ring radius is also required.
*
                  arcfrac2 = (apix * float(narc)) / 
     &                 (max(1.e-20,aanulus)*(phihi-philo)/TWOPI)
                  soal = 
     &             10.0*(sarc/arclen)*(RSAT/rave)/
     >                 max(1.e-20,arcfrac2)
                endif
              else
                mugain = 0.0
              endif ! if (npix)
*             write(*,*)'wphi          : ',wphi*CONVD
*             write(*,*)'philo,phihi : ',philo*CONVD,
*    &                 phihi*CONVD
*             write(*,*)'arclen (deg): ',arclen
*             write(*,*)'narc        : ',narc
*             write(*,*)'arcfrac2    : ',arcfrac2
*             write(*,*)'soal        : ',soal
*             write(*,*)'npix        " ',npix
*             write(*,*)'arcfrac1    : ',arcfrac1
*             write(*,*)'mugain      : ',mugain
*             write(*,*)'muskew      : ',muskew
            endif ! if(fint)
          else
            gain = 0.0
            mugain = 0.0
          endif
*
* As always, set values to some overflow which flag an ivalid
* event
*
        else	
          xave = 1000.0
          yave = 1000.0
          rave = 1000.0
          gain = 1.0e10
          mugain = 1.0e10
          arc = 0.0
        endif

        return
        end 
