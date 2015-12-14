SUBROUTINE exio2
     
!     EXIO2 COPIES SUBSTRUCTURE ITEMS BETWEEN THE SOF AND AN EXTERNAL
!     TAPE USING FORTRAN FORMATTED IO.  THE TAPE COULD HAVE BEEN CREATED
!     OR COULD BE READ ON A DIFFERENT BRAND OF COMPUTER.
 
 
 LOGICAL :: univac
 INTEGER :: dry      ,xblk     ,uname    ,pos      ,  &
     UNIT     ,fort     ,num(32)  ,sofin    , sofout   ,rewi     ,eqf
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm
 COMMON /BLANK / dry      ,xblk     ,device(2),uname(2) ,  &
     formt(2) ,mode(2)  ,pos(2)   ,datype(2),  &
     names(10),UNIT     ,univac   ,lbuf     , iadd
 COMMON /system/ sysbuf   ,nout     ,x1(36)   ,nbpc     , nbpw     ,ncpw
 DATA    fort  , sofin    ,sofout   ,rewi     ,eqf      /  &
     4HFORT, 4HSOFI   ,4HSOFO   ,4HREWI   ,4HEOF    /
 DATA    num   / 2H1   , 2H2      ,2H3      ,2H4      ,2H5      ,  &
     2H6   , 2H7      ,2H8      ,2H9      ,2H10     ,  &
     2H11  , 2H12     ,2H13     ,2H14     ,2H15     ,  &
     2H16  , 2H17     ,2H18     ,2H19     ,2H20     ,  &
     2H21  , 2H22     ,2H23     ,2H24     ,2H25     ,  &
     2H26  , 2H27     ,2H28     ,2H29     ,2H30     , 2H31  , 2H32     /
 
!     INITIALIZE
 
 nogo = 0
 
!     DECODE FORTRAN UNIT
 
 IF (uname(1) /= fort) GO TO 20
 DO  i = 1,32
   UNIT = i
   IF (uname(2) == num(UNIT)) GO TO 30
 END DO
 20 nogo = 1
 CALL page2 (-2)
 WRITE (nout,6356) uwm,uname
 
!     DECODE MODE OF OPERATION
 
 30 iomode = 0
 IF (mode(1) == sofout) iomode = 1
 IF (mode(1) == sofin ) iomode = 2
 IF (iomode  >      0) GO TO 40
 nogo = 1
 CALL page2 (-2)
 WRITE (nout,6338) uwm,mode
 
!     IF ERRORS THEN QUIT
 
 40 IF (nogo == 0) GO TO 50
 dry = -2
 GO TO 300
 
!     SET POSITION AND UNIVAC FLAGS
 
 50 univac = .true.
 IF (xblk <= 0) xblk = 3960
 xblk = xblk - MOD(xblk,132)
 lbuf = xblk/ncpw
 IF (MOD(xblk,ncpw) /= 0) lbuf = lbuf + 1
 iadd = 2
 IF (pos(1) == rewi) iadd = 1
 IF (pos(1) ==  eqf) iadd = 3
 
!     BRANCH ON MODE OF OPERATION
 
 SELECT CASE ( iomode )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 200
 END SELECT
 
!     SOFOUT
 
 100 CALL exo2
 GO TO 300
 
!     SOFIN
 
 200 CALL exi2
 
!     NORMAL MODULE COMPLETION
 
 300 RETURN
 
!     MESSAGE TEXT
 
 6338 FORMAT (a25,' 6338, ',2A4,' IS AN INVALID MODE PARAMETER FOR ',  &
     'MODULE EXIO')
 6356 FORMAT (a25,' 6356, ',2A4,' IS AN INVALID UNIT FOR MODULE EXIO,',  &
     ' EXTERNAL FORMAT')
END SUBROUTINE exio2
