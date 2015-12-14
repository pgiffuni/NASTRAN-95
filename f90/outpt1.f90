SUBROUTINE outpt1
     
!     COPY DATA BLOCK(S) ONTO NASTRAN USER TAPE WHICH MUST BE SET-UP.
 
!     CALL TO THIS MODULE IS
 
!     OUTPUT1   IN1,IN2,IN3,IN4,IN5//V,N,P1/V,N,P2/V,N,P3 $
 
!               P1 = 0, NO ACTION TAKEN BEFORE WRITE (DEFAULT)
!                  =+N, SKIP FORWARD N DATA BLOCKS BEFORE WRITE
!                  =-1, USER TAPE IS REWOUND BEFORE WRITE
!                  =-2, A NEW REEL IS MOUNTED BEFORE WRITE
!                  =-3, THE NAMES OF ALL DATA BLOCKS ON USER TAPE ARE
!                       PRINTED AND WRITE OCCURS AT THE END OF TAPE
!                  =-4, AN INPUT TAPE IS TO BE DISMOUNTED.
!                       A NEW OUTPUT REEL WILL THEN BE MOUNTED.
!                  =-9, WRITE EOF, REWIND AND UNLOAD.
 
!               P2 = 0, FILE NAME IS INPT (DEFAULT)
!                  = 1, FILE NAME IS INP1
!                  = 2, FILE NAME IS INP2
!                  = 3, FILE NAME IS INP3
!                  = 4, FILE NAME IS INP4
!                  = 5, FILE NAME IS INP5
!                  = 6, FILE NAME IS INP6
!                  = 7, FILE NAME IS INP7
!                  = 8, FILE NAME IS INP8
!                  = 9, FILE NAME IS INP9
 
!               P3 = TAPE ID CODE FOR USER TAPE, AN ALPHANUMERIC
!                    VARIABLE WHOSE VALUE WILL BE WRITTEN ON A USER
!                    TAPE. THE WRITTING OF THIS ITEM IS DEPENDENT ON
!                    THE VALUE OF P1 AS FOLLOWS..
!                          *P1*             *TAPE ID WRITTEN*
!                           +N                     NO
!                            0                     NO
!                           -1                    YES
!                           -2                    YES (ON NEW REEL)
!                           -3                     NO (WARNING CHECK)
!                           -4                    YES (ON NEW REEL)
!                           -9                     NO
!                    DEFAULT VALUE FOR P3 IS XXXXXXXX
 
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: tapeup,tapbit
 INTEGER :: trl(7),NAME(2),subnam(2),in(5),namex(2),ott(10),  &
     idhdr(7),idhdrx(7),p3x(2),d(3),dx(3),tapcod(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach
 COMMON /BLANK / p1,p2,p3(2) /system/ ksystm(65)  &
     /zzzzzz/ x(1)
 EQUIVALENCE     (ksystm( 1),nb  ),(ksystm( 2),nout),  &
     (ksystm( 9),nlpp),(ksystm(12),line), (ksystm(15),d(1))
 DATA    subnam/ 4HOUTP,4HT1        /
 DATA    in    / 101,102,103,104,105/
 DATA    zero  , mone,mtwo,mtre,mfor,mnin/ 0,-1,-2,-3,-4,-9/
 DATA    ott   / 4HINPT,4HINP1,4HINP2,4HINP3,4HINP4,  &
     4HINP5,4HINP6,4HINP7,4HINP8,4HINP9/
 DATA    idhdr / 4HNAST,4HRAN ,4HUSER,4H tap,4HE id,4H cod,4HE - /
 
 
 lcor = korsz(x) - 2*nb
 IF (lcor <= 0) GO TO 9908
 inbuf = lcor  + 1
 oubuf = inbuf + nb
 tapcod(1) = p3(1)
 tapcod(2) = p3(2)
 IF (p2 < 0 .OR. p2 > 9) GO TO 9904
 out = ott(p2+1)
 IF (mach >= 5) GO TO 120
 tapeup = tapbit(out)
 IF (.NOT.tapeup ) GO TO 9909
 120 IF (p1 < mnin) GO TO 9905
 IF (p1 > mnin .AND. p1 < mfor) GO TO 9905
 
 IF (p1 == mnin) GO TO 5000
 IF (p1 == mtre) GO TO 2000
 IF (p1 <= zero) GO TO 150
 
 CALL gopen (out,x(oubuf),2)
 DO  i = 1,p1
   CALL READ (*9903,*9903,out,namex,2,1,nf)
   CALL skpfil (out,1)
 END DO
 CALL CLOSE (out,2)
 GO TO 190
 
 150 IF (p1 /= mtwo .AND. p1 /= mfor) GO TO 190
 
!     P1 = -2 OR P1 = -4 IS ACCEPTABLE ONLY ON IBM OR UNIVAC
 
 IF (mach /= 2 .AND. mach /= 3) GO TO 9905
 
 iold = 3 + p1/2
 CALL gopen  (out,x(oubuf),3)
 CALL tpswit (out,iold,2,tapcod)
 
!     OPEN USER TAPE TO WRITE WITHOUT REWIND
 
 190 CALL gopen (out,x(oubuf),3)
 IF (p1 /= mone .AND. p1 /= mtwo .AND. p1 /= mfor) GO TO 195
 CALL REWIND (out)
 CALL WRITE (out,d,3,0)
 CALL WRITE (out,idhdr,7,0)
 CALL WRITE (out,p3,2,1)
 CALL eof (out)
 GO TO 195
 
 193 CALL CLOSE (out,2)
 CALL gopen (out,x(oubuf),3)
 
 195 DO  i = 1,5
   INPUT  = in(i)
   trl(1) = INPUT
   CALL rdtrl (trl)
   IF (trl(1) <= 0) CYCLE
   CALL fname (INPUT,NAME)
   
!     OPEN INPUT DATA BLOCK TO READ WITH REWIND.
   
   CALL OPEN  (*9901,INPUT,x(inbuf),0)
   CALL WRITE (out,NAME,2,0)
   CALL WRITE (out,trl(2),6,1)
   
!     LEVEL 17.5, THE ABOVE 8 WORD RECORD WAS WRITTEN OUT IN 2 RECORDS
!     2 BCD WORD NAME, AND 7 TRAILER WORDS
   
!     COPY CONTENTS OF INPUT DATA BLOCK ONTO USER TAPE.
   
   CALL cpyfil (INPUT,out,x,lcor,nf)
   
!     CLOSE INPUT DATA BLOCK WITH REWIND
   
   CALL CLOSE (INPUT,1)
   
   CALL eof (out)
   CALL page2 (-4)
   WRITE  (nout,350) uim,NAME,out,(trl(ii),ii=2,7)
   350 FORMAT (a29,' 4114', //5X,'DATA BLOCK ',2A4,  &
       ' WRITTEN ON NASTRAN FILE ',a4,', TRLR  =',6I10)
   
 END DO
 
!     CLOSE NASTRAN USER TAPE WITHOUT REWIND, BUT WITH END-OF-FILE
 
 CALL CLOSE (out,3)
 RETURN
 
!     OBTAIN LIST OF DATA BLOCKS ON USER TAPE.
 
 2000 CALL OPEN (*9902,out,x(oubuf),0)
 CALL READ (*9911,*9912,out,dx,3,0,nf)
 CALL READ (*9911,*9912,out,idhdrx,7,0,nf)
 DO  kf = 1,7
   IF (idhdrx(kf) /= idhdr(kf)) GO TO 9913
 END DO
 CALL READ (*9911,*9912,out,p3x,2,1,nf)
 IF (p3x(1) /= p3(1) .OR. p3x(2) /= p3(2)) GO TO 9914
 2006 CALL skpfil (out,1)
 kf = 0
 2007 CALL page1
 line = line + 5
 WRITE  (nout,2010) out
 2010 FORMAT (//50X,a4,14H FILE contents ,/46X,4HFILE,18X,4HNAME,//)
 2020 CALL READ (*2050,*9915,out,namex,2,1,nf)
 CALL skpfil (out,1)
 kf = kf + 1
 line = line + 1
 WRITE  (nout,2030) kf,namex
 2030 FORMAT (45X,i5,18X,2A4)
 IF (line-nlpp < 0) THEN
   GO TO  2020
 ELSE
   GO TO  2007
 END IF
 2050 CALL skpfil (out,-1)
 GO TO 193
 
 5000 CONTINUE
 CALL eof (out)
 CALL unload( out)
 RETURN
 
!     ERRORS
 
 9901 mm = -1
 GO TO 9996
 9902 WRITE  (nout,9952) sfm,out
 9952 FORMAT (a25,' 4117, SUBROUTINE OUTPT1 UNABLE TO OPEN NASTRAN FILE'  &
     ,       a4,1H.)
 line = line + 2
 GO TO 9995
 9903 WRITE  (nout,9953) ufm,p1,out,i
 9953 FORMAT (a23,' 4118, MODULE OUTPUT1 IS UNABLE TO SKIP FORWARD',i10,  &
     ' DATA BLOCKS ON PERMANENT NASTRAN FILE ',a4,1H., /5X,  &
     'NUMBER OF DATA BLOCKS SKIPPED =',i6)
 line = line + 3
 GO TO 9995
 9904 WRITE  (nout,9954) ufm,p2
 9954 FORMAT (a23,' 4119, MODULE OUTPUT1 - ILLEGAL VALUE FOR SECOND ',  &
     'PARAMETER =',i20)
 line = line + 2
 GO TO 9995
 9905 WRITE  (nout,9955) ufm,p1
 9955 FORMAT (a23,' 4120, MODULE OUTPUT1 - ILLEGAL VALUE FOR FIRST ',  &
     'PARAMETER =',i20)
 line = line + 2
 GO TO 9995
 9908 mm = -8
 INPUT = -lcor
 GO TO 9996
 9909 WRITE  (nout,9959) ufm,out
 9959 FORMAT (a23,' 4127, USER TAPE ',a4,' NOT SET UP.')
 line = line + 2
 GO TO 9995
 9911 WRITE  (nout,9961) ufm,out
 9961 FORMAT (a23,' 4128, MODULE OUTPUT1 - END-OF-FILE ENCOUNTERED ',  &
     'WHILE ATTEMPTING TO READ TAPE ID CODE ON USER TAPE ',a4)
 line = line + 2
 GO TO 9995
 9912 WRITE  (nout,9962) ufm,out
 9962 FORMAT (a23,' 4129, MODULE OUTPUT1 - END-OF-RECORD ENCOUNTERED ',  &
     'WHILE ATTEMPTING TO READ TAPE ID CODE ON USER TAPE ',a4)
 line = line + 2
 GO TO 9995
 9913 WRITE  (nout,9963) ufm,(idhdrx(kf),kf=1,7)
 9963 FORMAT (a23,' 4130, MODULE OUTPUT1 - ILLEGAL TAPE CODE HEADER = ', 7A4)
 line = line + 2
 GO TO 9995
 9914 WRITE  (nout,9964) uwm,p3x,p3
 9964 FORMAT (a25,' 4131, USER TAPE ID CODE -',2A4,'- DOES NOT MATCH ',  &
     'THIRD OUTPUT1 DMAP PARAMETER -',2A4,2H-.)
 line = line + 2
 GO TO 2006
 9915 WRITE  (nout,9965) sfm
 9965 FORMAT (a25,' 4115, MODULE OUTPUT1 - SHORT RECORD.')
 line = line + 2
 GO TO 9995
 
 9995 mm = -37
 9996 CALL mesage (mm,INPUT,subnam)
 RETURN
 
END SUBROUTINE outpt1
