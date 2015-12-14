SUBROUTINE outmsc (*,*)
     
!     COPY DATA BLOCK(S) TO FORTRAN UNIT, IN MSC/OUTPUT2 COMPATIBLE
!     RECORD FORMATS.
 
!     DMAP CALL -
!     OUTPUT2  IN1,IN2,IN3,IN4,IN5/ /V,N,P1/V,N,P2/V,N,P3/V,N,P4/V,N,P5/
!                                    V,N,P6 $
 
!     THIS ROUTINE IS CALLED ONLY BY OUTPT2
!     SEE OUTPT2 FOR PARAMETERS P1,P2,...,P6. (P6 = *MSC*)
 
!     IF P1 .NE. -9, ALTERNATE RETURN 1, OTHERWISE RETURN 2.
 
!     WRITTEN BY G.CHAN/UNISYS  3/93
 
 
 , INTENT(IN OUT)                         :: *
 , INTENT(IN OUT)                         :: *
 LOGICAL :: dp
 INTEGER :: p1,p2,p3,p4,p5,p6,endrec,endfil,out,buf1,d,  &
     inp(13),mcb(7),NAME(2),NONE(2),sub(2),tmp(2),  &
     dx(3),hdr(7),hdrx(7),tapcod(2),BLOCK(20)
 REAL :: xns(1)
 DOUBLE PRECISION :: dxns(1)
 CHARACTER (LEN=19) :: mo2
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /BLANK /  p1,p2,p3(2),p4,p5,p6(2)
 COMMON /machin/  mach
 COMMON /system/  ibuf,nout,idum4(6),nlpp,idum5(5),d(3)
 COMMON /TYPE  /  idum6(2),nwds(4)
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (xns(1),z(1))
 EQUIVALENCE      (xns(1),dxns(1))
 DATA    hdr   /  4HNAST,4HRAN ,4HFORT,4H tap,4HE id,4H cod,4HE - /
 DATA    inp   /  4HUT1 ,4HUT2 ,4HUT3 ,4HINPT,4HINP1,4HINP2,4HINP3,  &
     4HINP4,4HINP5,4HINP6,4HINP7,4HINP8,4HINP9       /
 DATA    mo2   /  '. MODULE OUTPUT2 - '      /
 DATA    NONE  ,  sub   /4H (no,4HNE) ,4HOUTP,4HUT2*              /
 
 WRITE  (nout,10) uim
 10 FORMAT (a29,'. USER REQUESTED RECORDS IN MSC/OUTPUT2 COMPATIBLE',  &
     ' RECORDS')
 
 endfil = 0
 endrec = 0
 lcor   = korsz(z(1))
 buf1   = lcor - ibuf + 1
 IF (buf1 <= 0) CALL mesage (-8,lcor,sub)
 lend   = buf1 - 1
 out    = p2
 tapcod(1) = p3(1)
 tapcod(2) = p3(2)
 IF (p1 == -9) GO TO 210
 IF (p1 == -3) GO TO 300
 IF (p1 <= -2) GO TO 620
 IF (p1 <=  0) GO TO 40
 
!     SKIP FORWARD n DATA BLOCKS, P1 = n
 
 i = 1
 20 READ (out) key
 keyx = 2
 IF (key /= keyx) GO TO 720
 READ (out) tmp
 READ (out) key
 IF (key >= 0) GO TO 740
 ASSIGN 30 TO iret
 nskip = 1
 GO TO 500
 30 i = i + 1
 IF (i <= p1) GO TO 20
 
 40 IF (p1 /= -1) GO TO 80
 REWIND out
 key = 3
 WRITE (out) key
 WRITE (out) d
 key = 7
 WRITE (out) key
 WRITE (out) hdr
 key = 2
 WRITE (out) key
 WRITE (out) p3
 endrec = endrec - 1
 WRITE (out) endrec
 WRITE (out) endfil
 endrec = 0
 WRITE  (nout,50) uim,p3
 50 FORMAT (a29,' FROM OUPUT2 MODULE.  THE LABEL IS ',2A4)
 
 80 DO  ii = 1,5
   INPUT  = 100 + ii
   mcb(1) = INPUT
   CALL rdtrl (mcb(1))
   IF (mcb(1) <= 0) CYCLE
   CALL fname (INPUT,NAME)
   IF (NAME(1) == NONE(1) .AND. NAME(2) == NONE(2)) CYCLE
   BLOCK(1) = INPUT
   nwd = nwds(mcb(5))
   dp  = mcb(5) == 2 .OR. mcb(5) == 4
   
!     OPEN INPUT DATA BLOCK TO READ WITH REWIND
   
   CALL OPEN (*600,INPUT,z(buf1),0)
   key = 2
   WRITE (out) key
   WRITE (out) NAME
   endrec = endrec - 1
   WRITE (out) endrec
   key = 7
   WRITE (out) key
   WRITE (out) mcb
   endrec = endrec - 1
   WRITE (out) endrec
   
!     COPY CONTENTS OF INPUT DATA BLOCK ONTO FILE
   
   90 CALL rectyp (INPUT,k)
   key = 1
   WRITE (out) key
   WRITE (out) k
   IF (k == 0) GO TO 130
   
!     STRING RECORD
!     BLOCK(2) = STRING TYPE, 1,2,3 OR 4
!     BLOCK(4) = FIRST (OR LAST) ROW POSITION ON A MATRIX COLUMN
!     BLOCK(5) = POINTER TO STRING, W.R.T. XNS ARRAY
!     BLOCK(6) = NO. OF TERMS IN STRING
   
   BLOCK(8) = -1
   100 CALL getstr (*170,BLOCK)
   key = BLOCK(6)*nwd
   WRITE (out) key
   
!     NEXT 3 LINES, ORIGINATED FROM MSC/OUTPUT2, DO NOT WORK FOR D.P.
!     DATA ON VAX, AND POSSIBLY SILICON-GRAPHICS. THEY ARE REPLACED BY
!     NEXT 8 LINES BELOW. BESIDE, TO WORK ON PROPER D.P. DATA BOUNDARY,
!     THE K1 IN THE FOLLOWING LINE SHOULD BE  K1 = (BLOCK(5)-1)*NWD+1
   
!     K1  = BLOCK(5)
!     K2  = K1 + KEY - 1
!     WRITE (OUT) BLOCK(4),(XNS(K),K=K1,K2)
   
   k1  = BLOCK(5)*nwd
   k2  = k1 + key -1
   IF (dp) GO TO 110
   WRITE (out) BLOCK(4),(xns(k),k=k1,k2)
   GO TO 120
   110 k1  = k1/2
   k2  = k2/2
   WRITE (out) BLOCK(4),(dxns(k),k=k1,k2)
   
   120 CALL endget (BLOCK)
   GO TO 100
   
!     NON-STRING RECORD
!     MAKE SURE EACH RECORD IS NOT LONGER THAN P4 WORDS
   
   130 CALL READ (*180,*150,INPUT,z(1),lend,0,k1)
   DO  i = 1,lend,p4
     key = lend - i + 1
     IF (key >= p4) key = p4
     k2 = i + key - 1
     WRITE (out) key
     WRITE (out) (z(k),k=i,k2)
   END DO
   GO TO 130
   150 DO  i = 1,k1,p4
     key = k1 - i + 1
     IF (key >= p4) key = p4
     k2 = i + key - 1
     WRITE (out) key
     WRITE (out) (z(k),k=i,k2)
   END DO
   
   170 endrec = endrec - 1
   WRITE (out) endrec
   GO TO 90
   
!     CLOSE INPUT DATA BLOCK WITH REWIND
   
   180 CALL CLOSE (INPUT,1)
   WRITE (out) endfil
   endrec = 0
   WRITE  (nout,190) uim,NAME,out,inp(p2-10),mcb
   190 FORMAT (a29,' 4144. DATA BLOCK ',2A4,' WRITTEN ON FORTRAN UNIT ',  &
       i3,2H (,a4,1H), /5X,'TRAILER =',6I7,i11)
   
 END DO
 
!     CLOSE FORTRAN TAPE WITHOUT END-OF-FILE AND WITHOUT REWIND
 
 RETURN 1
 
!     FINAL CALL TO OUTPUT2, P1 = -9
 
 210 WRITE (out) endfil
 RETURN 2
 
!     OBTAIN LIST OF DATA BLOCKS ON FORTRAN TAPE, P1 = -3
 
 300 REWIND out
 READ (out) key
 keyx = 3
 IF (key /= keyx) GO TO 720
 READ (out) dx
 READ (out) key
 keyx = 7
 IF (key /= keyx) GO TO 720
 READ (out) hdrx
 DO  k = 1,7
   IF (hdrx(k) /= hdr(k)) GO TO 640
 END DO
 READ (out) key
 keyx = 2
 IF (key /= keyx) GO TO 720
 READ (out) tmp
 IF (tmp(1) /= p3(1) .OR. tmp(2) /= p3(2)) GO TO 660
 320 ASSIGN 330 TO iret
 nskip = 1
 GO TO 500
 330 k = 0
 340 CALL page1
 WRITE  (nout,350) inp(p2-10),out
 350 FORMAT (//42X,'CONTENTS OF ',a4,', FORTRAN UNIT',i3, /46X,  &
     'FILE',18X,'NAME',/)
 360 READ (out) key
 IF (key < 0) THEN
   GO TO   680
 ELSE IF (key == 0) THEN
   GO TO   400
 END IF
 370 READ (out) tmp
 ASSIGN 380 TO iret
 nskip = 1
 GO TO 500
 380 k = k + 1
 WRITE  (nout,390) k,tmp
 390 FORMAT (45X,i5,18X,2A4)
 IF (MOD(k,nlpp) == 0) THEN
   GO TO   340
 ELSE
   GO TO   360
 END IF
 400 ASSIGN 80 TO iret
 nskip = k + 1
 IF (nskip > 0) REWIND out
 GO TO 500
 
!     SKIP NSKIP FILES ON FORTRAN TAPE
 
 500 IF (nskip == 0) GO TO 540
 DO  j = 1,nskip
   510 READ (out) keyx
   IF (keyx < 0) THEN
     GO TO   510
   ELSE IF (keyx == 0) THEN
     GO TO   530
   END IF
   520 IF (keyx > lcor) GO TO 700
   READ (out) (z(l),l=1,keyx)
   GO TO 510
 END DO
 540 GO TO iret, (30,80,330,380)
 
!     ERRORS
 
 600 CALL fname (INPUT,tmp)
 WRITE  (nout,610) sfm,mo2,tmp
 610 FORMAT (a25,' 4116',a19,'UNABLE TO OPEN INPUT DATA BLOCK ',2A4)
 GO TO  800
 620 WRITE  (nout,630) ufm,mo2,p1
 630 FORMAT (a23,' 4120',a19,'ILLEGAL FIRST PARAMETER ',i3)
 GO TO  800
 640 WRITE  (nout,650) ufm,mo2,hdrx
 650 FORMAT (a23,' 4130',a19,'ILLEGAL TAPE HEADER CODE ',7A4)
 GO TO  800
 660 WRITE  (nout,670) uwm,tmp,p3
 670 FORMAT (a25,' 4141. FORTRAN TAPE ID CODE - ',2A4,  &
     ' DOES NOT MATCH OUTPUT2 THIRD PARAMETER NAME - ',2A4)
 GO TO  320
 680 WRITE  (nout,690) sfm,mo2
 690 FORMAT (a25,' 4415',a19,'SHORT RECORD ENCOUNTERED')
 GO TO  800
 700 WRITE  (nout,710) ufm,lcor,key
 710 FORMAT (a23,' 2187. INSUFFICIENT WORKING CORE TO HOLD FORTRAN ',  &
     'LOGICAL RECORD.', /5X,'LENGHT OF WORKING CORE =',i11,  &
     '.   LENGTH OF FORTRAN LOGICAL RECORD =',i11)
 GO TO  800
 720 WRITE  (nout,730) sfm,key,keyx
 730 FORMAT (a25,' 2190. ILLEGAL VLUE FOR KEY =',i10,1H.,5X,  &
     'EXPECTED VALUE =',i10)
 GO TO  800
 740 WRITE  (nout,750) sfm,key
 750 FORMAT (a25,' 2190. ILLEGAL VALUE FOR KEY =',i10)
 800 CALL mesage (-61,0,sub)
 RETURN 1
 
END SUBROUTINE outmsc
