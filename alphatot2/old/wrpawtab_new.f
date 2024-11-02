      subroutine wrpawtab()
      vector cv,stats,sigmas,tduron,tduroff
      vector ondurs,offdurs,durs,tdurcut
      real alpha,excess(7)
      integer on(7),off(7),i,inton(7),intoff(7)
      integer thison(7),thisoff(7),k
      double precision dpon,dpoff
      character*10 onfile,offfile
      integer statontmp,statofftmp

C     Write Header/Cuts
      write(20,*)'PAWtab.out'
      write(20,*)' '
      write(20,*)'Cuts:'
      write(20,100)cv(1),cv(2)  !distance
      write(20,110)cv(3),cv(4)  !width
      write(20,120)cv(5),cv(6)  !length
      write(20,130)cv(7)        !size
      write(20,140)cv(8),cv(9)  !trigger
      write(20,150)cv(10)       !nbr3
      alpha=asin(cv(11))*(45.0/atan(1.0))
      write(20,160)alpha        !alpha
      write(20,*)' '
      write(20,*)' '
      write(20,300)                 ! lots of '-'
      write(20,310)                 ! raw, trig ....
      write(20,300)                 ! lots of '-'
C
C
 10   read(1,200,end=50) onfile,(on(i),i=1,7)
      read(2,200,end=50) offfile,(off(i),i=1,7)
      k=k+1
      do i=1,7
         thison(i)=on(i)-inton(i)
         thisoff(i)=off(i)-intoff(i)
         dpon=thison(i)
         dpoff=thisoff(i)
         excess(i)=(dpon-dpoff)/sqrt(dpon+dpoff)
         inton(i)=on(i)
         intoff(i)=off(i)
      enddo
      write(20,210)onfile,(thison(i),i=1,7)
      write(20,210)offfile,(thisoff(i),i=1,7)
      write(20,220)(excess(i),i=1,7)
      Write(20,230)ondurs(k),offdurs(k),durs(k)
      goto 10

 50   continue

      write(20,300)                 ! lots of '-'
      write(20,210)'Total ON',(stats(1,i),i=1,7)
      write(20,210)'Total OFF',(stats(2,i),i=1,7)
      write(20,220)(sigmas(i),i=1,7)
      Write(20,240)tduron(1),tduroff(1),tdurcut(1)      

     

 100  format(1x,f5.3,' < distance < ',f5.3)
 110  format(1x,f5.3,' < width < ',f5.3)
 120  format(1x,f5.3,' < length < ',f5.3)
 130  format(1x,'Size > ',f6.0)
 140  format(1x,'1/91 > ',f4.0,',  2/91 > ',f4.0)
 150  format(1x,'Nbr3 > ',f3.1)
 160  format(1x,'Alpha < ',f6.2,' degrees.')
 300  format(79('-'))
 310  format(14x,'Raw',7x,'Nbr3',5x,'Trigger',5x,'Shape',3x,'Sh2065',
     +       3x,'Orient',3x,'Gammas')

 200  format(a10,7(1x,i9))
 210  format(a10,3(i9,2x),3(i7,2x),i6)
 220  format('Excess',2x,3(5x,f6.2),1x,4(3x,f6.2))
 230  format('ON duration:',f10.3,4x,'OFF duration:',f10.3,4x,
     +'Both cut to:',f10.3)
 240  format('Total ON duration:',f10.3,2x,'Total OFF duration:',f10.3,
     + 2x,'Cut to:',f10.3)
 250  format(10x,3(i9,2x),3(i7,2x),i6)


      return
      end
