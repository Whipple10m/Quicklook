* 
*     ULCOMIS.F
*     JB
*     951017
*
*     Modified 951017 by JB to calculate "exact" Li and Ma upper limit
*   
*     PAW COMIS program to calculate upper limits according to method of
*     Helene et al. solving for ULCOUNTS in the equation: 
*
*     Probability of detecting DIFF excess counts =
*     (1-CONFIDENCE) = 
*       I[(ULCOUNTS-DIFF)/SIGMA] / I[-DIFF/SIGMA]
* 
*     and calculating significance according to the Li and Ma equation.
*
*     Implicit Inputs:
*
*        vec	stats(2,8)	Vec containing num evts. on and off
*                               after cuts (e.g. 15<alpha<65, alpha<10)
*        vec    roffon(1)	Vector containing the ratio of counts
*                               in the off alpha interval (15,65) to
*                               the on alpha interval (0,10)
*        vec    efarea(1)	Vector containing effective area of
*                               in units of 10^8 cm^2
*     Implicit Outputs:
*
*        vec    sigoo(1)	Li and Ma signif (in sigma) for ON/OFF
*        vec    uloo(1)		Helene 3sigma upper lim. for ON/OFF
*        vec    rateoo(1)	Gamma ray rate in gammas/min
*        vec    fluxoo(1)       Flux in gammas cm^-2 sec^-1 for ON/OFF
*	 vec    ulfloo(1)       Flux upper limit for ON/OFF
*
*        vec    sigtr(1)	Li and Ma signif (in sigma) for Tracking
*        vec    ultr(1)		Helene 3sigma upper lim. for Trk analysis
*        vec    ratetr(1)	Gamma ray rate in gammas/min for Trk
*        vec    fluxtr(1)       Flux in gammas cm^-2 sec^-1 for Trk 
*	 vec    ulfltr(1)       Flux upper limit for Tracking analysis
*
      subroutine ulcomis()
*
      implicit none
*
      vector	stats,roffon
      vector    sigoo,uloo,rateoo,fluxoo,ulfloo
      vector    sigtr,ultr,ratetr,fluxtr,ulfltr
      vector    efarea
      vector    tdurcut
*
      real xarg
      real RTBIS
      real diff, on, off, ratio, alpha, sigma, time
      real sigma1,sigma2,signif1,signif2
      real arg1,arg2
      real ulcounts,ulflux
      real inttarget, x1, x2
      real intnormal
      real CONFIDENCE
      data CONFIDENCE / 0.999 /
      real AEFF
      real XACC                 ! Accurace to iterate to
      data XACC / 0.01 /        ! Iterate to 0.01 sigma
      logical success
      logical fapprox

      fapprox = .true.
      AEFF=efarea(1)*1.0e-04 	! Effective Area (cm^2)/ 1x10^12
*
* Do ON/OFF analysis if there were OFF runs
*
      if(stats(2,7).gt.0.1) then
        on = stats(1,7)		! Do ON/OFF analysis first
        off = stats(2,7)
*       ratio = roffon(1)
        ratio = 1.0
        time = tdurcut(1)/60.0
*       write(*,*)'on,off,ratio,time:',on,off,ratio,time
        alpha = 1.0/ratio
        diff = on-alpha*off
        sigma1 = sqrt(max(1.e-20,alpha*(on+off)))
        if (sigma1.gt.0.0001) then
          signif1 = diff/sigma1
        else
          signif1 = diff/0.0001
        endif
*       write(*,*)'alpha,diff,sigma1,signif1:',alpha,diff,
*    &   sigma1,signif1
        if (on+off.gt.0.0001) then
           arg1 = on*log(max(1.e-30,((1.+alpha)/alpha)*(on/(on+off))))
           arg2 = off*log(max(1.e-30,(1.+alpha)*(off/(on+off))))
        else
           arg1 = on*log(max(1.e-30,((1.+alpha)/alpha)*(on/0.0001)))
           arg2 = off*log(max(1.e-30,(1.+alpha)*(off/0.0001)))
        endif
        signif2 = sqrt(2.)*sqrt(max(1.e-30,abs(arg1+arg2)))
*       if((arg1+arg2).gt.0.) then
*         signif2 = sqrt(2.)*sqrt(arg1+arg2)
*         fapprox = .false.
*       else
*         write(*,*) 'Can not do exact Li and Ma significance'
*         fapprox = .true.
*       endif
        if(signif2.gt.0.) then
          if(diff.lt.0.0) then
            signif2 = -1.0*signif2
          endif
          sigma2 = (on-alpha*off)/signif2
          fapprox = .false.
        else
          write(*,*) 'Can not do exact Li and Ma significance'
          fapprox = .true.
        endif
*
* Use approximate Li and Ma formula
*
        if(fapprox) then
          sigma = sigma1
        else
          sigma = sigma2
        endif
        if (abs(sigma).gt.0.00001) then
           inttarget = (1.0-CONFIDENCE)*(intnormal(-1.0*diff/sigma))
        else
           inttarget = (1.0-CONFIDENCE)*(intnormal(-1.0*diff/0.00001))
        endif
        x1 = 0.0
        x2 = 3.0	! Check that root is bracketed by [x1,x2]
*       write(*,*)'sigma,inttarget,x1,x2',sigma,inttarget,x1,x2
        call ZBRACR(intnormal,x1,x2,inttarget,success)
*       write(*,*)'success,x1,x2',success,x1,x2
        if(success) then ! Now solve for root intnormal(xarg)=inttarget
*         write(*,*)'success'
          xarg = RTBIS(intnormal,x1,x2,inttarget,XACC)
*         write(*,*)'x1,x2,inttarget,xarg',x1,
*    &     x2,inttarget,xarg
          ulcounts = sigma*xarg + diff
          if (time.lt.0.00001) then
             time=0.00001
          endif
          ulflux = (ulcounts/(AEFF*time*60.0))
*         write(*,*)'ulcounts,ulflux',ulcounts,ulflux
          if(fapprox) then
            sigoo(1) = signif1
          else
            sigoo(1) = signif2
          endif
          uloo(1) = ulcounts/time
          ulfloo(1) = ulflux
          rateoo(1) = diff/time
          fluxoo(1) = (diff/(AEFF*time*60.0))
        else
          write(*,*)' ERROR: Could not find root within the interval:',
     &      '[',x1,',',x2,']'
          write(*,*)' Edit the source code for this program to change'
          write(*,*)' the values of the variables: x1,x2'
        endif
      endif
*
* Do alpha analysis no matter what
*
      on = stats(1,7)
      off = stats(1,8)
      ratio = roffon(1)
      time = tdurcut(1)/60.0
*     write(*,*)'on,off,ratio,time:',on,off,ratio,time
      alpha = 1.0/ratio
      diff = on-alpha*off
      sigma1 = sqrt(max(1.e-30,alpha*(on+off)))
      if (sigma1.gt.0.00001) then
         signif1 = diff/sigma1
      else
         signif1 = diff/0.00001
      endif
*     write(*,*)'alpha,diff,sigma1,signif1',alpha,diff,sigma1,signif1
      if (on+off.gt.0.0001) then
         arg1 = on*log(max(1.e-30,((1.+alpha)/alpha)*(on/(on+off))))
         arg2 = off*log(max(1.e-30,(1.+alpha)*(off/(on+off))))
      else
         arg1 = on*log(max(1.e-30,((1.+alpha)/alpha)*(on/0.0001)))
         arg2 = off*log(max(1.e-30,(1.+alpha)*(off/0.0001)))
      endif
      signif2 = sqrt(2.)*sqrt(max(1.e-30,abs(arg1+arg2)))
      if(signif2.gt.0.) then
        if(diff.lt.0.0) then
          signif2 = -1.0*signif2
        endif
        sigma2 = (on-alpha*off)/signif2
        fapprox = .false.
      else
        write(*,*) 'Can not do exact Li and Ma significance'
        fapprox = .true.
      endif
*
*     if((arg1+arg2).gt.0.) then
*         signif2 = sqrt(2.)*sqrt(arg1+arg2)
*         fapprox = .false.
*     else
*	  write(*,*) 'Can not do exact Li and Ma significance'
*         fapprox = .true.
*     endif
*       if(signif2.gt.0.) then
*         sigma2 = (on-alpha*off)/signif2
*         fapprox = .false.
*     else
*         write(*,*) 'Can not do exact Li and Ma significance'
*         fapprox = .true.
*     endif
*
* Determine whether to use approximate Li and Ma formula
*
      if(fapprox) then
        sigma = sigma1
      else
        sigma = sigma2
      endif
      if (abs(sigma).gt.0.00001) then
         inttarget = (1.0-CONFIDENCE)*(intnormal(-1.0*diff/sigma))
      else
         inttarget = (1.0-CONFIDENCE)*(intnormal(-1.0*diff/0.00001))
      endif
      x1 = 0.0
      x2 = 3.0	! Check that root is bracketed by [x1,x2]
      call ZBRACR(intnormal,x1,x2,inttarget,success)
      if(success) then ! Now solve for root intnormal(xarg)=inttarget
        if (time.lt.0.00001) then
           time=0.00001
        endif
        xarg = RTBIS(intnormal,x1,x2,inttarget,XACC)
        ulcounts = sigma*xarg + diff
        ulflux = (ulcounts/(AEFF*time*60.0))
        if(fapprox) then
          sigtr(1) = signif1
        else
          sigtr(1) = signif2
        endif
        ultr(1) = ulcounts/time
        ulfltr(1) = ulflux
        ratetr(1) = diff/time
        fluxtr(1) = (diff/(AEFF*time*60.0))
      else
        write(*,*)' ERROR: Could not find root within the interval:',
     &    '[',x1,',',x2,']'
        write(*,*)' Edit the source code for this program to change'
        write(*,*)' the values of the variables: x1,x2'
      endif
      return
      end
*
*     Function intnormal
*     JB,ADK
*     950501
*
*     Use the complement of the error function to calculate the
*     integral of the normal distribution
*
      real function intnormal(x)
      real	x,erfc
*
*     write(*,*)'intnormal x',x
      intnormal = (erfc(x/sqrt(2.)) / 2.)
      return
      end
*
* Numerical Recipes Functions 
* Included here to make the program self-contained.
*
* WARNING: Some functions have been modified slightly to fix 
* problems arrising from new DEC fortran compiler, among other
* things.  
*
      REAL FUNCTION ERFC(X)
      REAL X,GAMMP,GAMMQ
      IF(X.LT.0.)THEN
        ERFC=1.+GAMMP(.5,X**2)
*       write(*,*)'gammp ',gammp(0.5,X**2)
      ELSE
        ERFC=GAMMQ(.5,X**2)
*       write(*,*)'gammq ',gammq(0.5,X**2)
      ENDIF
      RETURN
      END
*
      REAL FUNCTION GAMMP(A,X)
      REAL A,X
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END
*
      REAL FUNCTION GAMMQ(A,X)
      REAL A,X,GAMSER,GLN
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      RETURN
      END
*
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      INTEGER ITMAX
      REAL EPS,GAMMLN
      PARAMETER (ITMAX=100,EPS=3.E-7)
      REAL GAMSER,A,X,GLN,SUM,DEL,AP
      INTEGER N
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
*
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      INTEGER ITMAX
      REAL EPS
      PARAMETER (ITMAX=100,EPS=3.E-7)
      REAL GAMMCF,A,X,GLN,GAMMLN
      REAL GOLD,A0,A1,B0,B1,FAC,AN,ANA,ANF,G
      INTEGER N
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
*
      FUNCTION GAMMLN(XX)
      DOUBLE PRECISION COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
*
*     Subroutine zbracr
*     JB
*
*     Modified version of Numerical recipes subroutine to find a 
*     root of an arbitrary function.  Bracket root of f(x) = FVAL
*     rather than f(x) = 0
*
      SUBROUTINE ZBRACR(FUNC,X1,X2,FVAL,SUCCES)
      REAL FACTOR
      INTEGER NTRY
      PARAMETER (FACTOR=1.6,NTRY=50)
      REAL F1,F2,X1,X2,FVAL
      INTEGER J
      LOGICAL SUCCES
c--*******************************************************************
c-- Eliminating source of warning message during compilation JR 900108
	EXTERNAL FUNC
c--*******************************************************************
      IF(X1.EQ.X2)PAUSE 'You have to guess an initial range'
      F1=FUNC(X1)-FVAL
      F2=FUNC(X2)-FVAL
      SUCCES=.TRUE.
      DO 11 J=1,NTRY
        IF(F1*F2.LT.0.)RETURN
        IF(ABS(F1).LT.ABS(F2))THEN
          X1=X1+FACTOR*(X1-X2)
          F1=FUNC(X1)-FVAL
        ELSE
          X2=X2+FACTOR*(X2-X1)
          F2=FUNC(X2)-FVAL
        ENDIF
11    CONTINUE
      SUCCES=.FALSE.
      RETURN
      END
*
*
*     Subroutine rtbis
*     JB
*
*     Modified version of Numerical recipes subroutine to find a 
*     root of an arbitrary function.  Here solve for f(x) = FVAL
*     rather than f(x) = 0
*
*     Find root of FUNC(X)-FVAL=0 in the interval [X1,X2] 
*     return when accuracy on X is XACC
*
      FUNCTION RTBIS(FUNC,X1,X2,FVAL,XACC)
      REAL X1,X2,FVAL,XACC
      REAL DX,FMID,XMID
      INTEGER J,JMAX
      PARAMETER (JMAX=40)
c--*******************************************************************
c-- Eliminating source of warning message during compilation JR 900108
        EXTERNAL FUNC
c--*******************************************************************
      FMID=FUNC(X2)-FVAL
      F=FUNC(X1)-FVAL
      IF(F*FMID.GE.0.) PAUSE 'Root must be bracketed for bisection.'
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        FMID=FUNC(XMID)-FVAL
        IF(FMID.LE.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
      PAUSE 'too many bisections'
      END

