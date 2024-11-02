      subroutine calc_sigs()
      
      vector stats,sigmas
      double precision dpstats(2,7)
      double precision arg
      integer i

      do i=1,7
         dpstats(1,i)=stats(1,i)
         dpstats(2,i)=stats(2,i)
         arg = dpstats(1,i)+dpstats(2,i)
         if (arg.gt.0.) then
           arg = SQRT(arg)
           if (arg.gt.0.) then
             sigmas(i)=(dpstats(1,i)-dpstats(2,i))/arg
           else
             sigmas(i)=0.0
           endif
         else
           sigmas(i)=0.0
         endif
      enddo

      return

      end
      
