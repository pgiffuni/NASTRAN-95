SUBROUTINE inptt3
     
!     THIS ROUTINE READS MATRIX DATA FROM AN INPUT TAPE, WRITTEN IN
!     ROCKWELL INTERNATIONAL COMPANY'S CUSTOMARY FORMAT, INTO NASTRAN
!     GINO MATRIX BLOCK.
!     (THE RI DATA IS IN A COMPACT FORTRAN-FORMATTED CODED FORM, DOUBLE
!     PRECISION, WHCIH APPEARS TO HAVE QUITE WIDESPREAD ACCEPTANCE IN
!     THE AEROSPACE FIELD, AND PARTICULARY IN MARSHALL SPACE FLIEGHT
!     CENTER (MSFC) AREA)
 
!     WRITTEN ORIGINALLY BY MEL MARTENS, ROCKWELL INTERNATIONAL, SPACE
!     DIVISION (213) 922-2316, AND MODIFIED UP TO NASTRAN STANDARD BY
!     G.CHAN/UNISYS, 2/1987
 
!     INPTT3  /O1,O2,O3,O4,O5/V,N,UNIT/V,N,ERRFLG/V,N,TEST  $
 
!             UNIT  = FORTRAN INPUT TAPE UNIT NO.
!                     TAPE IS REWOUND BEFORE READ IF UNIT IS NEGATIVE
!                     FORTRAN UNIT 11 (INPT) IS USED IF UNIT= 0 OR -1.
!             ERRFLG= 1, JOB TERMINATED IF DATA BLOCK ON TAPE NO FOUND
!                     0, NO TERMINATION IF DATA BLOCK NO FOUND ON TAPE
!             TEST  = 0, NO CHECK ON FILE NAMES ON TAPE AND DMAP NAMES
!                   = 1, NAMES CHECK, WILL SEARCH TAPE FOR MATCH.
 
 IMPLICIT INTEGER (a-z)
 INTEGER :: mcb(7),  NAME(2),  namx(2), subnam(2)
 DOUBLE PRECISION :: dz(1)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,     uwm,      uim
 COMMON /system/  ibuf,    nout
 COMMON /packx /  typin,   typout,   ii,      jj,      incr
 COMMON /zzzzzz/  z(1)
 COMMON /BLANK /  UNIT,    errflg,   test
 COMMON /names /  rd,      rdrew,    wrt,     wrtrew,  rew
 EQUIVALENCE      (z(1),dz(1))
 DATA             END,     head,     subnam                  /  &
     -999,    -111,     4HINPT,  4HT3           /
 
 core = korsz(z(1))
 buf1 = core - ibuf + 1
 core = buf1 - 1
 typin= 2
 typout=2
 incr = 1
 
 iu = UNIT
 IF (UNIT == 0 .OR. UNIT == -1) iu = -11
 IF (iu > 0) GO TO 10
 iu  = -iu
 irew= 0
 REWIND iu
 
 10 DO  k = 1,5
   FILE = 200 + k
   mcb(1) = FILE
   CALL rdtrl (mcb)
   IF (mcb(1) <= 0) CYCLE
   CALL gopen (FILE,z(buf1),wrtrew)
   CALL fname (FILE,NAME)
   20 READ (iu,30,ERR=160,END=180) i,namx
   30 FORMAT (i6,2A4)
   IF (i >    0) GO TO 20
   IF (i ==  END) GO TO 120
   IF (i /= head) GO TO 130
   IF (namx(1) == NAME(1) .AND. namx(2) == NAME(2)) GO TO 50
   WRITE  (nout,40) uim,namx,NAME
   40 FORMAT (a29,', DATA BLOCK ',2A4,' FOUND WHILE SEARCHING FOR ',2A4)
   IF (test == 0.0) THEN
     GO TO    70
   ELSE
     GO TO    20
   END IF
   
!     FOUND
   
   50 WRITE  (nout,60) uim,NAME
   60 FORMAT (a29,', DATA BLOCK ',2A4,' FOUND')
   70 READ   (iu,80) nr,nc,TYPE
   80 FORMAT (3I6)
   WRITE  (nout,90) NAME,nc,nr,TYPE
   90 FORMAT (/5X,'MATRIX BLOCK ',2A4,' IS OF SIZE ',i6,'(COL) BY',i5,  &
       '(ROW),  AND TYPE =',i6)
   IF (nr > core) CALL mesage (-8,nr-core,subnam)
   irew= 1
   ii  = 1
   jj  = nr
   CALL makmcb (mcb,FILE,nr,TYPE,2)
   DO  i = 1,nc
     READ (iu,100,ERR=160,END=180) (dz(j),j=1,nr)
     100 FORMAT (12X,1P,5D24.16)
     CALL pack (z,FILE,mcb)
   END DO
   CALL CLOSE (FILE,rew)
   CALL wrttrl (mcb)
   CYCLE
   
   120 IF (irew == 0) GO TO 130
   REWIND iu
   irew = 0
   GO TO 20
   130 WRITE  (nout,140) uwm,NAME
   140 FORMAT (a25,', INPTT3 FAILED TO LOCATE DATA BLOCK ',2A4,' ON ', 'TAPE')
   IF (errflg /= 0) CALL mesage (-61,0,subnam)
   REWIND iu
   irew = 0
 END DO
 RETURN
 
 160 WRITE  (nout,170) iu
 170 FORMAT ('0*** ERROR DUING READ.  TAPE UNIT',i5)
 CALL CLOSE (FILE,rew)
 CALL mesage (-61,0,subnam)
 180 WRITE  (nout,190) uwm,iu
 190 FORMAT (a25,' FROM INPTT3, EOF ENCOUNTERED ON INPUT TAPE',i4)
 CALL CLOSE  (FILE,rew)
 CALL wrttrl (mcb)
 RETURN
END SUBROUTINE inptt3
