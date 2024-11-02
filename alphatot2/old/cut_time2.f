      REAL FUNCTION cut_time(mode)
*********************************************************
*                                                       *
* This file was generated by HUWFUN.                    *
*                                                       *
*********************************************************
*
*     Ntuple Id:      100  
*     Ntuple Title:   GRAS
*     Creation:       25/05/95 10.09.41
*
*********************************************************
* MC 960209
* Changed cuts so that upper bounds can be used on size, 
* max1, and max2.  Currently these are applied as or
* cuts.  That is, the event is kept if 
* ((t10totsig.lt.size(2)).or.(t10max1.lt.max1(2)).or.
* (t10max2.lt.max2(2)))  
*
      LOGICAL         CHAIN
      CHARACTER*128   CFILE
      INTEGER IDNEVT,NCHEVT,ICHEVT
      REAL VIDN1,VIDN2,VIDN3,VIDN(10)
*
      COMMON /PAWIDN/ IDNEVT,VIDN1,VIDN2,VIDN3,VIDN
      COMMON /PAWCHN/ CHAIN, NCHEVT, ICHEVT
      COMMON /PAWCHC/ CFILE
*
*--   Ntuple Variable Declarations
*
      REAL t10azwidth,t10width,t10length,t10distance,t10sinalpha,t10miss
     + ,t10xcentr,t10ycentr,t10phi,t10totsig,t10max1,t10loc1,t10max2
     + ,t10loc2,t10max3,t10loc3,t10rmax1,t10rloc1,t10x1,t10y1,t10rmax2
     + ,t10rloc2,t10x2,t10y2,t10hadron,t10nb2,t10nb3,t10npicture
     + ,t10asymm,t10veto,t10el,t10frac2,t10maxsep,t10xgam,t10ygam
     + ,t10gain,t10arc,t10muon
      DOUBLE PRECISION t10utc,t10deltat,t10livetime,t10phase
      LOGICAL t10ipv,t10tiv
*
      COMMON /PAWCR4/ t10ipv,t10azwidth,t10width,t10length,t10distance
     + ,t10sinalpha,t10miss,t10xcentr,t10ycentr,t10phi,t10totsig,t10max1
     + ,t10loc1,t10max2,t10loc2,t10max3,t10loc3,t10rmax1,t10rloc1,t10x1
     + ,t10y1,t10rmax2,t10rloc2,t10x2,t10y2,t10hadron,t10nb2,t10nb3
     + ,t10npicture,t10asymm,t10veto,t10el,t10frac2,t10maxsep,t10xgam
     + ,t10ygam,t10gain,t10arc,t10muon,t10tiv
      COMMON /PAWCR8/ t10utc,t10deltat,t10livetime,t10phase
*

*
*--   Enter user code here
* 
* This selection function actually cuts the data, and generates 
* statistics for each stage of cutting.  The resulting statistics
* have the following meaning:
*    
* stats(?,1) - Number of raw events after duration cut alone
* stats(?,2) - Number of events after duration and NB3 cuts
* stats(?,3) - Number of events after Trig, i.e. MAX12+Size+NB3+Dur+Length/Size
* stats(?,4) - Number of events after Trig+Dist+Length+Width
* stats(?,4) - Number of events after Trig+Dist+Length+Width+Asymm+OFF-Alpha
* stats(?,5) - Number of events after Trig+Dist+Asymm+ON-Alpha
* stats(?,7) - Number of events after Trig+Dist+length+Width+Asymm+ON-Alpha
* stats(?,9) - Number of upo's within radius cut rcut
*
      logical cuts(13),lastevt
      vector stats,durs,cv
      vector rcutsq
      vector lwfcor
      real delx, dely, delr2
*     vector ra,dec
      real mode,alphadeg
      integer hist_idi,i
      character*26 mess(2)
      character*12 meshname 
      character*10 name
      
      if idnevt.eq.1 then
         mess(1)='2nd Pass - Cutting ON run '
         mess(2)='2nd Pass - Cutting OFF run'
         i=1
         modei=mode
         name=cfile(1:10)
         write(6,100)mess(modei),i,name,durs(i)
         meshname = cfile(1:6)//'.h2d'
         call HRGET(11,meshname,'A') 
         write(6,101)meshname
      else if (ichevt.eq.1) then
         i=i+1
         nchel=.true.
         write(modei,200)name,(stats(modei,j),j=1,7)
         name=cfile(1:10)
         write(6,100)mess(modei),i,name,durs(i)
         meshname = cfile(1:6)//'.h2d'
         call HRGET(11,meshname,'A') 
         write(6,101)meshname
      else if (idnevt.eq.nchevt) then
         lastevt=.true.
      endif
      

      cuts(1)=(t10distance.gt.cv(1)).and.(t10distance.lt.cv(2))
      cuts(2)=(t10width.gt.cv(3)).and.(t10width.lt.cv(4))
      cuts(3)=(t10length.gt.cv(5)).and.(t10length.lt.cv(6))
      cuts(4)=(t10totsig.gt.cv(7)).and.(t10max1.gt.cv(9)).and.
     +     (t10max2.gt.cv(11))
      cuts(5)=(t10sinalpha.lt.cv(14))
      cuts(6)=((t10totsig.lt.cv(8)).or.(t10max1.lt.cv(10)).or.
     +     (t10max2.lt.cv(12)))
      cuts(7)=(t10sinalpha.gt.cv(15)).and.(t10sinalpha.lt.cv(16))
      cuts(8)=((t10length/t10totsig).gt.cv(17)).and.
     +        ((t10length/t10totsig).lt.cv(18))
      cuts(9)=t10asymm.gt.cv(19)
      cuts(10)=t10nb3.gt.cv(13)
      cuts(11)=t10livetime.lt.durs(i)
      cuts(12)=t10sinalpha.le.1.0
      cuts(13)=cuts(8).and.cuts(10).and.cuts(4).and.cuts(6).and.cuts(11)
*
* Accumulate 2-d histograms if Trig+Length+Width Cuts are satisfied
*
      if cuts(13).and.cuts(2).and.cuts(3) then
         call HFILL(16+modei,t10xgam,t10ygam,1.)
         delx = t10xgam - t10xcentr
         dely = t10ygam - t10ycentr
         delx = delx*lwfcor(1)
         dely = dely*lwfcor(1)
         delr2 = delx*delx + dely*dely
         if (delr2.lt.rcutsq) then
           stats(modei,9) = stats(modei,9) + 1.0 
         endif
      endif
*
* If Cuts: Trig(Length/Size+NB3+Size+MAX12+Duration) + Dist + Width
* + Length + Sinalpha<1 
* Then accumulate the shape statistics, if additionally the asymmetry
* cut is satisfied, then accumulate the alpha histogram
* 
      if cuts(13).and.cuts(1).and.cuts(2).and.cuts(3)
     +                                     .and.cuts(12) then
*        alphadeg=asin(0.999999*t10sinalpha)*(45.0/atan(1.0))
         stats(modei,4)=stats(modei,4)+1	! Trig+Dist+Shape
         alphadeg=asin(t10sinalpha)*(45.0/atan(1.0))
         if cuts(9) then
           call HFILL(modei,alphadeg,0.,1.)
         endif
*
* JB 950817
* Temporarily dispense with 20-90, only show 20-65
*
*        if (alphadeg.gt.20.) stats(modei,5)=stats(modei,5)+1 !shape20-90
*
* Now add in Alpha and Asymm cuts to Trig+Dist+Shape to get the
* total cuts for the Alpha-ON and Alpha-OFF regions
*
         if cuts(5).and.cuts(9) stats(modei,7)=stats(modei,7)+1
         if cuts(7).and.cuts(9) then
           stats(modei,8)=stats(modei,8)+1	! JB 950720
           stats(modei,5)=stats(modei,8)
         endif
      endif
*
* Stuff 2-d histograms
*
*     if cuts(13).and.cuts(2).and.cuts(3) then	! JB 950720
*         call HFILL(20+modei,t10xgam,t10ygam,1.)
*     endif

      if cuts(11) then				! Raw - duration only
         stats(modei,1)=stats(modei,1)+1
      endif
      
      if (cuts(10).and.cuts(11)) then		! Raw+NB3
        stats(modei,2)=stats(modei,2)+1
      endif

      if cuts(13) then				! MAX1,2+Size+NB3+dur+L/S
         stats(modei,3)=stats(modei,3)+1
      endif

      if cuts(13).and.cuts(1).and.cuts(5).and.cuts(9) then
         stats(modei,6)=stats(modei,6)+1	 ! Trig+Dist+Alpha+Asymm
      endif
 
      if lastevt then
        write(modei,200)cfile(1:10),(stats(modei,j),j=1,7)
      endif

 100  format(a26,1x,i4,1x,a10,1x,'Cutting to',f11.3,' seconds.')
 101  format('Opening 2-d hist file: ',1x,a12)
 200  format(a10,7(1x,i9))
* 
      END





