      subroutine ntoff(filename,runname,nnt,ntubesoff) 
*
* ** Update hrc.ntubelist, the standard list of tubes to be turned off in
* ** analysis.
*
* 951017 JB - Modified to wait for the lock file to dissappear, but
* ** eventually ignore it, since it shouldn't take more than 2min
* ** for another program to access the lock file.
*
      implicit none

      integer*4 MAXOFF
      parameter (MAXOFF=15)   !Maximum number of tubes to be turned off
                              ! Note that present 80 character record length
                              ! and a79 line format only can support
                              ! 17 tubes turned off.
      character*80 filename    ! Path and filename for ntoff list
                              ! e.g. /home/dbgt/hrc.ntubelist
      character*6 runname
      integer*4 nnt,ntubesoff(MAXOFF)	! New tubes to be turned off
     

      integer*4 irec,lnblnk
      integer*4 i,nt,inew,tubesoff(MAXOFF)

      logical exist

      character line*80,last*80,path*40,sel*1,run*6,lfch*1
      integer n_fail

*****************************************************************************
*
* ** Define the line feed character.  We need to put it at the end of
* ** record or life is not good.
*
      lfch = char(10)
*
* ** First check for a lock file to see whether someone else is running
* ** the program "ntoff".  Problems will occur when 2 or more processes
* ** access the same file in write mode.  Do this even though this new
* ** version includes a path.  This is easier than figuring out which
* ** path another user may be using.
*
      n_fail = 0
 5    inquire(file='/tmp/ntoff.lock',exist=exist)
*
* ** If no lock then create one and make it write accessible to all; if
* ** the lock file is present then halt execution
*
      if (.not.exist) then
         open(2,file='/tmp/ntoff.lock',status='unknown')
         call system('chmod a+rw /tmp/ntoff.lock')
      else
         write(6,*)'NT ** hrc.ntubelist locked by another user...'
         write(6,*)'      try again in 20sec'
         call system('sleep 20')
         n_fail=n_fail+1
         if(n_fail.gt.5)then
           write(6,*)' NT ** WARNING: ntoff waited its turn for 2min '
           write(6,*)'       but will now continue despite the lock'
           goto 10
         endif
         goto 5
      end if
*
10    continue
      write(6,'('' NT ** Filename: '',a)')
     &    filename(1:lnblnk(filename))
      
      inquire(file=filename(1:lnblnk(filename)),exist=exist)
      if (.not.exist) then
         write(6,'('' NT ** Creating file: '',a)')
     &    filename(1:lnblnk(filename))
      end if
*
* ** The file exists or should be created.  Open it in direct accesss mode.
*
      open(10,file=filename(1:lnblnk(filename)),form='formatted',
     &   access='direct',status='unknown',recl=80)
*
* ** If this is a newly created file, make it accessible by
* ** everyone
*
      if(.not.exist) then
        call system('chmod 777 '//filename(1:lnblnk(filename)))
      endif
*
* ** Initialize variable last, which is the last entry examined, to
* ** be empty.  Fortran probably does this for us, but what the hey.
*
      line = ' '
*
* ** Top of loop.  Get the name of the data run of interest.
* ** If the name doesn't at least have the right number of characters,
* ** ask for another one.
*
      last = line
      if (lnblnk(runname) .ne. 6) then
        write(*,'(''NT ** '',a6,''is not a valid run name. Exiting.'')')
     &   runname
        goto 999
      end if
*
* ** Now read the file a line at a time, looking for the run of 
* ** interest.  If we find it, then display it with user options.  If
* ** it is not found, then make a new entry.
*
* ** Top of file search loop
*
      irec = 0
40    continue
         irec = irec + 1
         read(10,rec=irec,fmt=120,err=99)line
*        read(10,120,end=99,rec=irec)line
         read(line,125)run,nt,(tubesoff(i),i=1,nt)
*        write(*,121)line
*
* ** Read another line if we this isn't the right one
*
         if (run .ne. runname) go to 40
      continue
*
* ** End of search loop
* ** At this point we've found the entry for the file of interest.  Display
*
50    continue
      write(*,*)' NT ** Old tubelist entry was found:'
      write(*,120)line
      write(*,*)' NT ** Exiting'
      goto 999

99    continue
      sel = 'n'    		!Make a new entry
      line(1:6) = runname	!Put filename into line
      do i = 7,79
         line(i:i) = ' '	!Pad the rest of the line with blanks
      end do
      if (nnt .gt. MAXOFF) then
         write(*,150)nnt
         nnt = MAXOFF
      end if
      write(line(7:10),170)nnt
      do i = 1,nnt
         write(line(7+i*4:10+i*4),170)ntubesoff(i)
      end do
      line(80:80) = lfch
      write(*,*)'NT ** Added new entry to hrc.tubelist:'
      write(6,119)line
      write(10,121,rec=irec)line

999   continue
      close(10)
      call system('rm -f /tmp/ntoff.lock')
      return
*
* ** Formats
*
119   format(x,a80)
120   format(a79)
121   format(a80)
125   format(a6,1x,i3,15(1x,i3))	! Note -- this must match MAXOFF
150   format('NT ** Sorry...max number of tubes is ',i2)
170   format(xi3)

      end
