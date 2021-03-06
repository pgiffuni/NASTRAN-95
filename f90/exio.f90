SUBROUTINE exio
     
!     THE MAIN PURPOSE OF THIS MODULE IS TO COPY DATA BETWEEN THE
!     RESIDENT SOF AND AN EXTERNAL TAPE OR DISK FILE.  AS AN EXTRA
!     ADDED ATTRACTION, IT WILL ALSO APPEND AN EXTERNAL SOF (CREATED BY
!     SOME OTHER NASTRAN RUNS) TO THE RESIDENT SOF AND COMPRESS THE
!     RESIDENT SOF.
 
!     OPTIONS ARE -
 
!     (1) DUMP (RESTORE) THE ENTIRE SOF TO (FROM) AN EXTERNAL FILE.
!         INTERNAL FORM ONLY.  THIS IS THE MOST EFFICIENT MEANS TO SAVE
!         OR RECOVER A BACKUP COPY OF THE SOF, EXCEPT FOR SYSTEM UTILITY
!         PROGRAMS.
 
!     (2) COPY SELECTED ITEMS BETWEEN THE SOF AND AN EXTERNAL FILE.
 
!     (3) CHECK THE EXTERNAL FILE AND PRINT OUT A LIST OF ALL SUBSTRUC-
!         TURES AND ITEMS ON IT ALONG WITH THE DATE AND TIME EACH WAS
!         CREATED.
 
!     (4) APPEND AN EXTERNAL SOF TO THE RESIDENT SOF.
 
!     (5) COMPRESS THE RESIDENT SOF. (PLACE ALL ITEMS IN CONTIGUOUS
!         BLOCKS ON THE SOF AND ELIMINATE ALL EMBEDDED FREE BLOCKS)
 
!     FEBRUARY 1974
 
 INTEGER :: FORMAT,exte,BLANK,head1,head2,dry,device,uname,  &
     pos,bcds(2,10),inbcds(2,5)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / dry,xmach,device(2),uname(2),FORMAT(2),mode(2),  &
     pos(2),datype(2),names(10),pdate,ptime
 COMMON /system/ sysbuf,nout
 COMMON /output/ head1(96),head2(96)
 EQUIVALENCE     (inte,bcds(1,7)),(exte,bcds(1,8)), (device(1),inbcds(1,1))
 DATA    BLANK / 4H       /
 DATA    bcds  / 4HSOFI   ,4HN     , 4HSOFO   ,4HUT    ,  &
     4HREST   ,4HORE   , 4HCHEC   ,4HK     ,  &
     4HCOMP   ,4HRESS  , 4HAPPE   ,4HND    ,  &
     4HINTE   ,4HRNAL  , 4HEXTE   ,4HRNAL  ,  &
     4HREWI   ,4HND    , 4HNORE   ,4HWIND  /
 
 DO  i = 1,96
   head2(i) = BLANK
 END DO
 loop30:  DO  i = 1,5
   DO  j = 1,10
     IF (inbcds(1,i) /= bcds(1,j)) CYCLE
     inbcds(2,i) = bcds(2,j)
     CYCLE loop30
   END DO
 END DO loop30
 
 DO  i = 1,2
   head2(i   ) = mode(i)
   head2(i+ 3) = FORMAT(i)
   head2(i+ 6) = device(i)
   head2(i+ 9) = uname(i)
   head2(i+12) = pos(i)
 END DO
 
!     INTERNAL FORMAT - GINO I/O IS USED FOR DATA WHICH WILL BE READ OR
!                       WAS WRITTEN ON THE SAME HARDWARE.
 
 IF (FORMAT(1) == inte) CALL exio1
 
!     EXTERNAL FORMAT - FORTRAN I/O IS USED FOR DATA WHICH WILL BE READ
!                       OR WAS WRITTEN ON A DIFFERENT MACHINE.
 
 IF (FORMAT(1) == exte) CALL exio2
 
!     CHECK VALIDITY OF FORMAT TO ASCERTAIN WHETHER EITHER EXIO1 OR
!     EXIO2 WAS CALLED.
 
 IF (FORMAT(1) == inte .OR. FORMAT(1) == exte) RETURN
 WRITE  (nout,50) uwm,FORMAT
 50 FORMAT (a25,' 6333, ',2A4,' IS AN INVALID FORMAT PARAMETER FOR ',  &
     'MODULE EXIO.')
 dry = -2
 RETURN
END SUBROUTINE exio
