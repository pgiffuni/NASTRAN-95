SUBROUTINE delscr
     
!     THIS SUBROUTINE IS CALLED AT THE BEGINNING OF EACH FUNCTIONAL
!     MODULE TO PHYSICALLY DELETE ALL SCRATCH FILES USED BY A
!     PREVIOUS FUNCTIONAL MODULE
 
 INTEGER*2       iunit
 
 COMMON /dsunit/ iunit(220)
 COMMON /xfist / ifist(2)
 COMMON /xpfist/ ipfist
 INCLUDE 'DSIOF.COM'
 
 DATA mask / 32767 /
 
 nfiles = ifist(2) - ipfist
 IF (nfiles == 0) RETURN
 istr = ipfist + 1
 iend = ifist(2)
 DO  i = istr, iend
   ifile = ifist(2*i+1)
   IF (ifile < 301 .OR. ifile > 320) CYCLE
   ifilex = 0
   CALL geturn (ifile)
   IF (ifilex == mask) CYCLE
   CALL dbmmgr ( 7 )
 END DO
 RETURN
END SUBROUTINE delscr
