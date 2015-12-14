SUBROUTINE paraml
     
!     TO SELECT PARAMETERS FROM A GINO DATA BLOCK
 
!     PARAML  DB/ /C,N,OP/V,N,P1/V,N,P2/V,N,RSP/V,N,INTEG/V,N,RDP/
!                  V,N,BCD/V,N,SPLX/V,N,DPLX $
 
!     INPUT GINO FILE -
!       DB = TABLE  INPUT FILE IF OP='TABLEi'
!       DB = MATRIX INPUT FILE IF OP='MATRIX','NULL', etc.
!     OUTPUT GINO FILE -
!       NONE
!     INPUT PARAMETER -
!       OP    = OPERATION FLAG, ONE OF THE FOLLOWING KEY WORDS,
!               'MATRIX', 'NULL', 'PRESENCE', 'TRAILER', OR
!               'TABLE1' - ABSTRACT FROM 1 INPUT WORD TO FORM ALL OUTPUT
!                          DATA TYPE (INTEGER, S.P /D.P. REAL S.P./D.P.
!                          COMPLEX) AND 4-BYTE BCD WORD (1 WORD)
!               'TABLE2' - ABSTRACT FROM 2 INPUT WORDS TO FORM ALL
!                          OUTPUT DATA TYPE, AND 8-BYTE BCD (2 WORDS)
!               'TABLE4' - ABSTRACT FORM 4 INPUT WORDS TO FORM S.P./D.P.
!                          COMPLEX NUMBER
!               'TABLE1/2/4' OPERATES ONLY IN TABLE  DATA BLOCK, AND
!                THE OTHERS  OPERATE  ONLY IN MATRIX DATA BLOCK.
 
!                IF 'PRESENCE' IS ABBREVIATED AS 'PRES  ', THE USER
!                PARAML INFORMATION MESSAGE IS NOT ECHOED OUT.
 
!     INPUT/OUTPUT PARAMETERS -
!       P1    = RECORD NO. IF DB IS A TABLE, OR
!       P1    = ROW NO. IF DB IS A MATRIX
!               (DEFAULT=1)
!       P2    = WORD POSITION INDEX (BASED ON S.P.REAL WORD COUNT)
!               IF DB IS A TABLE, OR
!       P2    = COLUMN NUMBER, IF DB IS A MATRIX DATA BLOCK, S.P. OR
!               D.P.
!               (DEFAULT=1)
!       (ROW FIRST AND COLUMN SECOND - IN CONSISTANT WITH SCALAR MODULE)
!     OUTPUT PARAMETERS -
!       RSP   = SINGLE PRECISION REAL
!               (DATA ABSTRACTED FROM 1 OR 2 INPUT WORDS)
!       INTEG = INTEGER (DATA ABSTRACTED FROM 1 INPUT WORD)
!       RDP   = DOUBLE PREC. FLOATING NUMBERS (FROM 1 OR 2 INPUT WORDS)
!       BCD   = 8-BYTE BCD WORD, BLANK FILLED IF NECCESSARY
!       SPLX  = SINGLE PRECISION COMPLEX (FROM 1 TO 4 INPUT WORDS)
!       DPLX  = DOUBLE PRECISION COMPLEX (FROM 1 TO 4 INPUT WORDS)
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: tb1,tb2,tb4,mat,prt
 INTEGER :: mcb(7),NAME(2),ivps(1),opcd(7),fnm(2), nmvps(2),ei(3),at(2)
 REAL :: z(1),rsp,splx,sp(4),vps,x,y
 DOUBLE PRECISION :: dz(1),rdp,dplx,dp(2)
 CHARACTER (LEN=7) :: nty(4)
 CHARACTER (LEN=10) :: TYPE(4)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /xvps  /  vps(2)
 COMMON /unpakx/  ityp,ii,jj,incr
 COMMON /ilocal/  il(2),il3,il4,il5,il6,il7,il8,il9
 COMMON /system/  sysbuf,nout
 COMMON /BLANK /  op(2),p1,p2,rsp,integ,rdp,bcd(2),splx(2),dplx(2)
 COMMON /zzzzzz/  iz(1)
 EQUIVALENCE      (vps(1),ivps(1)) ,(z(1),iz(1),dz(1))
 EQUIVALENCE      (sp(1) ,  dp(1))
 DATA NAME / 4HPARA,4HML  /,BLANK/4H     /, at/ 4HAND ,4HTHRU /
 DATA opcd / 4HTABL,4HMATR,4HPRES,4HNULL,4HTRAI,4HDTI ,4HDMI  /
 DATA first/ 12 /  ,in1   / 101  /,   ei /2HE1, 2HE2, 2HE4    /
 DATA nty  / 'ZERO', 'INTEGER', 'REAL',  'BCD' /
 DATA TYPE / 'S.P. REAL ', 'D.P. REAL ', 'S.P. CMPLX', 'D.P.CMPLX'/
 
!     SUPPRESS ALL PARAML CHECKING MESSAGES IF DIAG 37 IS ON
 
 CALL sswtch (37,i)
 prt   = i == 0
 nz    = korsz(iz)
 ibuf1 = nz - sysbuf + 1
 IF (ibuf1 <= 0) GO TO 1220
 flag  = 1
 mcb(1)= in1
 CALL rdtrl (mcb)
 IF (mcb(1) > 0) GO TO 20
 
!     INPUT PURGED.  RETURN IF OP(1) IS NOT 'PRES'
 
 IF (op(1) /= opcd(3)) GO TO 1240
 flag  =-1
 CALL fndpar (-5,il5)
 IF (prt .AND. op(2) /= BLANK) WRITE (nout,40) uim,op
 10 integ = flag
 ivps(il5) = flag
 nmvps(1) = ivps(il5-3)
 nmvps(2) = ivps(il5-2)
 IF (prt .AND. op(2) /= BLANK) WRITE (nout,510) integ,nmvps
 GO TO 1240
 
 20 prec = mcb(5)
 CALL fname (in1,fnm)
 DO  j=3,9
   CALL fndpar (-j,il(j))
 END DO
 IF (op(1) == opcd(3) .AND. op(2) == BLANK) GO TO 200
 IF (op(1) == opcd(4)) GO TO 210
 IF (.NOT.prt) GO TO 45
 CALL page2 (first)
 first = 5
 WRITE  (nout,40) uim,op
 40 FORMAT (a29,' FROM PARAML MODULE  - ',2A4,2H -, /5X,  &
     '(ALL PARAML MESSAGES CAN BE SUPPRESSED BY DIAG 37)',/)
 
!     IDENTIFY OPCODE
 
 45 DO  i = 1,7
   IF (op(1) == opcd(i)) THEN
      SELECT CASE ( i )
       CASE (    1)
         GO TO 300
       CASE (    2)
         GO TO 800
       CASE (    3)
         GO TO 200
       CASE (    4)
         GO TO 210
       CASE (    5)
         GO TO 220
       CASE (    6)
         GO TO 90
       CASE (    7)
         GO TO 90
     END SELECT
   END IF
 END DO
 60 WRITE  (nout,70) ufm,op
 70 FORMAT (a23,', ILLEGAL OP REQUEST TO MODULE PARAML - ',2A4)
 80 CALL mesage (-37,0,NAME)
 
 90 IF (.NOT.prt) GO TO 60
 WRITE  (nout,100) uim
 100 FORMAT (a29,', NEW PARAMETERS USED IN PARAML MODULE:', //5X,  &
     'PARAML  DB//C,N,OP/C,N,P1/V,N,P2/V,N,RSP/V,N,INT/V,N,RDP/',  &
     'V,N,BCD/V,N,CSX/V,N,CDX  $', /13X,  &
     'OP      = OPCODE, ONE OF THE FOLLOWING KEY WORDS, BCD INPUT, N',  &
     'O DEFAULT', /23X,43H'MATRIX', 'NULL', 'PRESENCE', 'TRAILER', OR,  &
     /23X,28H'TABLE1', 'TABLE2', 'TABLE4',  &
 /13X,'P1,P2   = RECORD NO. AND WORD POSITION IF OP= TABLEi',  &
 /21X,'= ROW AND COLUMN INDEXES IF OP= MATRIX,  INTEGERS INPUT',  &
 /21X,'= P2 GIVES THE VALUE OF P1 TRAILER WORD IF OP= TRAILER',  &
     /13X,'RSP,RDP = SINGLE PRECISION AND DOUBLE PREC. REAL, OUTPUT',  &
     /23X,'(DEFAULTS ARE 0.0 AND 0.D+0,  PREVIOUS DEFAULTS WARE ONES',  &
     /13X,'INT,BCD = INTEGER AND 2-BCD WORDS OUTPUT', /23X,'INT =-1,',  &
 ' IF NULL MATRIX AND OP= NULL, OR PURGED DB AND OP= PRESENCE',  &
     /13X,'CSX,CDX = SINGLE PRECISION AND DOUBLE PRECISION COMPLEX, ',  &
     'OUTPUT', //5X,'EXAMPLE - ',  &
     'ABSTRACT THE 3RD COL. 9TH ROW ELEMENT OF KGG MATRIX, AND', /15X,  &
     'ABSTRACT THE 3RD RECORD AND 9TH WORD  OF EPT DATA BLCOK', //5X,  &
     'PARAML  KGG//*MATRIX*/C,N,9/C,N,3/V,N,R93//V,N,D93//V,N,CS93',  &
     /5X,'PARAML  EPT//*TABLE1*/C,N,3/C,N,9//V,N,I39/V,N,D39',/)
 IF (i == 6) WRITE (nout,110)
 IF (i == 7) WRITE (nout,120)
 110  FORMAT (5X,'SUGGESTION- REPLACE THE OPCODE ''DTI'' BY ''TABLE1''')
 120  FORMAT (5X,'SUGGESTION- REPLACE THE OPCODE ''DMI'' BY ''MATRIX''',  &
     /18X,'AND NOTE THAT P1 IS ROW NUMBER AND P2 IS COLUMN NO.')
 GO TO 60
 
!     OP = PRESENCE
!     TEST FOR PRESENCE OF DATA BLOCK
 
 200 GO TO 10
 
!     OP = NULL
!     TEST FOR NULL MATRIX DATA BLOCK
 
 210 IF (mcb(7) == 0) flag =-1
 GO TO 10
 
!     OP = TRAILER
!     PLACE THE (P1+1) WORD OF THE TRAILER IN P2
 
 220 IF (p1 <= 0 .OR. p1 >= 7) GO TO 230
 p2 = mcb(p1+1)
 ivps(il3) = p2
 nmvps(1) = ivps(il3-3)
 nmvps(2) = ivps(il3-2)
 IF (prt) WRITE (nout,510) p2,nmvps
 GO TO 1240
 230 WRITE  (nout,240) ufm,p1
 240 FORMAT (a23,', 2ND PARAMETER IN PARAML MODULE IS ILLEGAL',i5)
 GO TO 80
 
!     OP = TABLE
!     PROCESS TABLE TYPE DATA BLOCK
 
 300 tb1 = .false.
 tb2 = .false.
 tb4 = .false.
 IF (op(2) == ei(1)) tb1 = .true.
 IF (op(2) == ei(2)) tb2 = .true.
 IF (op(2) == ei(3)) tb4 = .true.
 IF (.NOT.tb1 .AND. .NOT.tb2 .AND. .NOT.tb4) GO TO 60
 mat = .false.
 recno = p1
 INDEX = p2
 IF (tb2) ixp1 = INDEX+1
 IF (tb4) ixp1 = INDEX+3
 atx = at(1)
 IF (tb4) atx = at(2)
 CALL OPEN (*1200,in1,iz(ibuf1),0)
 CALL skprec (in1,recno)
 CALL READ (*1210,*310,in1,iz,ibuf1-1,1,rl)
 GO TO 1220
 310 IF (INDEX > rl) GO TO 1210
 IF (il4 <= 0) GO TO 500
 
!     OUTPUT REQUEST IN S.P. REAL
 
 IF (.NOT.prt) GO TO 350
 IF (.NOT.tb1) GO TO 330
 WRITE  (nout,320) fnm,recno,INDEX
 320 FORMAT (5X,'INPUT FILE ',2A4,' RECORD',i6,' WORD',i6,13X,1H=)
 GO TO 350
 330 WRITE  (nout,340) fnm,recno,INDEX,atx,ixp1
 340 FORMAT (5X,'INPUT FILE ',2A4,' RECORD',i6,' WORDS',i6,1X,a4,i5, '  =')
 350 nmvps(1) = ivps(il4-3)
 nmvps(2) = ivps(il4-2)
 IF (tb4) GO TO 400
 IF (tb2) GO TO 355
 rsp = z(INDEX)
 IF (mat) GO TO 360
 k = numtyp(rsp)+1
 IF (k == 2 .OR. k == 4) GO TO 400
 GO TO 360
 355 k = -1
 IF (INDEX+1 > rl) GO TO 400
 sp(1) = z(INDEX  )
 sp(2) = z(INDEX+1)
!WKBI
 IF ( sp(2) == 0.0 ) dp(1) = sp(1)
!WKBR RSP = SNGL(DP(1))
 rsp = sp(1)
 k = numtyp(rsp)+1
 IF (k == 2 .OR. k == 4) GO TO 400
 360 IF (prt) WRITE (nout,370) rsp,nmvps
 370 FORMAT (1H+,70X,e15.8,'   = ',2A4)
 vps(il4) = rsp
 GO TO 500
 
 400 IF (.NOT.prt) GO TO 500
 WRITE  (nout,410) nmvps
 410 FORMAT (1H+,70X,'(INVALID REQUEST) = ',2A4)
 IF (k > 0) WRITE (nout,420) uwm,nty(k),nmvps
 IF (k == -1) WRITE (nout,430) uwm,nmvps
 420 FORMAT (a25,' - ILLEGAL OUTPUT REQUESTED. ORIG. DATA TYPE IS ',a7,  &
     ',  PARAMETER ',2A4,' NOT SAVED')
 430 FORMAT (a25,' - E-O-R ENCOUNTERED.  PARAMETER ',2A4,' NOT SAVED')
 
 500 IF (il5 <= 0 .OR. mat) GO TO 540
 
!     OUTPUT REQUEST IS INTEGER
 
 IF (.NOT.prt) GO TO 505
 IF (     tb1) WRITE (nout,320) fnm,recno,INDEX
 IF (.NOT.tb1) WRITE (nout,340) fnm,recno,INDEX,atx,ixp1
 505 nmvps(1) = ivps(il5-3)
 nmvps(2) = ivps(il5-2)
 k = 0
 IF (tb2 .OR. tb4) GO TO 520
 integ = iz(INDEX)
 k = numtyp(integ)+1
 IF (k > 2) GO TO 520
 ivps(il5) = integ
 IF (prt) WRITE (nout,510) integ,nmvps
 510 FORMAT (1H+,70X,i15,'   = ',2A4)
 GO TO 540
 
 520 IF (.NOT.prt) GO TO 540
 WRITE (nout,410) nmvps
 IF (k > 0) WRITE (nout,420) uwm,nty(k),nmvps
 IF (k == 0) WRITE (nout,530) uwm,nmvps
 530 FORMAT (a25,' - ILLEGAL INTEGER ABSTRACTION FROM 2 OR 4 DATA ',  &
     'WORDS.  OUPUT PARAMETER ',2A4,' NOT SAVED')
 GO TO 540
 
 540 IF (il6 <= 0) GO TO 600
 
!     OUTPUT REQUEST IN D.P. REAL
 
 IF (.NOT.prt) GO TO 545
 IF (     tb1) WRITE (nout,320) fnm,recno,INDEX
 IF (.NOT.tb1) WRITE (nout,340) fnm,recno,INDEX,atx,ixp1
 545 nmvps(1) = ivps(il6-3)
 nmvps(2) = ivps(il6-2)
 IF (mat) GO TO 560
 IF (tb2) GO TO 550
 IF (tb4) GO TO 590
 k = numtyp(z(INDEX))+1
 IF (k == 2 .OR. k == 4) GO TO 590
 dp(1) = DBLE(z(INDEX))
 GO TO 570
 550 k =-1
 j = 0
 IF (INDEX+1 > rl) GO TO 590
 sp(1) = z(INDEX  )
 sp(2) = z(INDEX+1)
!WKBD 9/93      X = SNGL(DP(1))
 x = sp(1)
 j = numtyp(x)+1
 IF (j == 2 .OR. j == 4) GO TO 590
 GO TO 570
 560 IF (prec == 1) dp(1) = DBLE(z(INDEX))
!WKBI
 570 IF ( sp(2) == 0.0 ) dp(1) = sp(1)
!WKBR  570 RDP = DP(1)
 rdp = dp(1)
 vps(il6  ) = sp(1)
 vps(il6+1) = sp(2)
 IF (prt) WRITE (nout,580) rdp,nmvps
 580 FORMAT (1H+,70X,d15.8,'   = ',2A4)
 GO TO 600
 
 590 IF (.NOT.prt) GO TO 600
 WRITE (nout,410) nmvps
 IF (j == 2 .OR. j == 4) k = j
 IF (k > 0) WRITE (nout,420) uwm,nty(k),nmvps
 IF (k == -1) WRITE (nout,430) uwm,nmvps
 
 600 IF (il7 <= 0 .OR. mat) GO TO 650
 
!     OUTPUT REQUEST IN BCD
 
 IF (.NOT.prt) GO TO 605
 IF (     tb1) WRITE (nout,320) fnm,recno,INDEX
 IF (.NOT.tb1) WRITE (nout,340) fnm,recno,INDEX,atx,ixp1
 605 nmvps(1) = ivps(il7-3)
 nmvps(2) = ivps(il7-2)
 k = 0
 IF (tb4) GO TO 630
 bcd(1) = iz(INDEX)
 bcd(2) = BLANK
 k = numtyp(bcd(1))+1
 IF (k /= 4) GO TO 630
 IF (tb1) GO TO 610
 k = -1
 IF (INDEX+1 > rl) GO TO 630
 bcd(2) = iz(INDEX+1)
 k = numtyp(bcd(2))+1
 IF (k /= 4) GO TO 630
 610 ivps(il7  ) = bcd(1)
 ivps(il7+1) = bcd(2)
 IF (prt) WRITE (nout,620) bcd,nmvps
 620 FORMAT (1H+,70X,2A4,'   = ',2A4)
 GO TO 650
 
 630 IF (.NOT.prt) GO TO 650
 WRITE (nout,410) nmvps
 IF (k > 0) WRITE (nout,420) uwm,nty(k),nmvps
 IF (k == 0) WRITE (nout,640) uwm,nmvps
 IF (k == -1) WRITE (nout,430) uwm,nmvps
 640 FORMAT (a25,' - ILLEGAL BCD ABSTRACTION FROM 4 DATA WORDS. ',  &
     ' PARAMETER ',2A4,'NOT SAVED')
 
 650 IF (il8 <= 0) GO TO 700
 
!     OUTPUT REQUEST IN S.P. COMPLEX
 
 IF (.NOT.prt) GO TO 655
 IF (     tb1) WRITE (nout,320) fnm,recno,INDEX
 IF (.NOT.tb1) WRITE (nout,340) fnm,recno,INDEX,atx,ixp1
 655 nmvps(1) = ivps(il8-3)
 nmvps(2) = ivps(il8-2)
 k =-1
 j = 0
 IF (tb4) GO TO 660
 splx(1) = z(INDEX)
 splx(2) = 0.0
 IF (tb1 .OR. mat) GO TO 670
 IF (INDEX+1 > rl) GO TO 690
 splx(2) = z(INDEX+1)
 GO TO 670
 660 IF (INDEX+3 > rl) GO TO 690
 sp(1)   = z(INDEX  )
 sp(2)   = z(INDEX+1)
 sp(3)   = z(INDEX+2)
 sp(4)   = z(INDEX+3)
!WKBR SPLX(1) = SNGL(DP(1))
 splx(1) = sp(1)
!WKBR SPLX(2) = SNGL(DP(2))
 splx(2) = sp(3)
 670 j = numtyp(splx(1))+1
 k = numtyp(splx(2))+1
 IF (j == 2 .OR. j == 4 .OR. k == 2 .OR. j == 4) GO TO 690
 vps(il8  ) = splx(1)
 vps(il8+1) = splx(2)
 IF (prt) WRITE (nout,680) splx,nmvps
 680 FORMAT (1H+,70X,1H(,e15.8,1H,,e15.8,1H),'  = ',2A4)
 GO TO 700
 
 690 IF (.NOT.prt) GO TO 700
 WRITE (nout,410) nmvps
 IF (j == 2 .OR. j == 4) k = j
 IF (k == 0) WRITE (nout,420) uwm,nty(k),nmvps
 IF (k == -1) WRITE (nout,430) uwm,nmvps
 
 700 IF (il9 <= 0) GO TO 1100
 
!     OUTPUT REQUEST IN D.P. COMPLEX
 
 IF (.NOT.prt) GO TO 705
 IF (     tb1) WRITE (nout,320) fnm,recno,INDEX
 IF (.NOT.tb1) WRITE (nout,340) fnm,recno,INDEX,atx,ixp1
 705 nmvps(1) = ivps(il9-3)
 nmvps(2) = ivps(il9-2)
 k =-1
 j = 0
 IF (tb4) GO TO 710
 k = numtyp(z(INDEX))+1
 IF (k == 2 .OR. k == 4) GO TO 740
 dp(1) = DBLE(z(INDEX))
 dp(2) = 0.d0
 IF (tb1 .OR. mat) GO TO 720
 IF (INDEX+1 > rl) GO TO 740
 k = numtyp(z(INDEX+1))+1
 IF (k == 2 .OR. k == 4) GO TO 740
 dp(2) = DBLE(z(INDEX+1))
 GO TO 720
 710 IF (INDEX+3 > rl) GO TO 740
 sp(1) = z(INDEX  )
 sp(2) = z(INDEX+1)
 sp(3) = z(INDEX+2)
 sp(4) = z(INDEX+3)
!WKBR X = SNGL(DP(1))
 x = sp(1)
!WKBR Y = SNGL(DP(2))
 y = sp(3)
 j = numtyp(x)+1
 k = numtyp(y)+1
 IF (j == 2 .OR. j == 4 .OR. k == 2 .OR. k == 4) GO TO 740
 dp(1)   = DBLE(z(INDEX))
 dp(2)   = 0.d0
 720 dplx(1) = dp(1)
 dplx(2) = dp(2)
 vps(il9  ) = sp(1)
 vps(il9+1) = sp(2)
 vps(il9+2) = sp(3)
 vps(il9+3) = sp(4)
 IF (prt) WRITE (nout,730) dplx,nmvps
 730 FORMAT (1H+,70X,1H(,d15.8,1H,,d15.8,1H),'  = ',2A4)
 GO TO 1100
 
 740 IF (.NOT.prt) GO TO 1100
 WRITE (nout,410) nmvps
 IF (j == 2 .OR. j == 4) k = j
 IF (k > 0) WRITE (nout,420) uwm,nty(k),nmvps
 IF (k == -1) WRITE (nout,430) uwm,nmvps
 GO TO 1100
 
!     OP = MATRIX
!     PROCESS MATRIX TYPE DATA BLOCK
 
 800 row  = p1
 col  = p2
 ityp = mcb(5)
 IF (il5 <= 0) GO TO 840
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 810 FORMAT (5X,'ELEMENT (',i5,'-ROW,',i5,'-COL) OF ',a10,' INPUT ',  &
     'FILE ',2A4,2H =)
 nmvps(1) = ivps(il5-3)
 nmvps(2) = ivps(il5-2)
 IF (.NOT.prt) GO TO 840
 WRITE  (nout,820) nmvps
 820 FORMAT (1H+,70X,'(INVALID INTEGER) = ',2A4)
 WRITE  (nout,830) uwm,nmvps
 830 FORMAT (a25,' - OUTPUT PARAMETER ',2A4,' NOT SAVED')
 840 IF (il7 <= 0) GO TO 860
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 nmvps(1) = ivps(il7-3)
 nmvps(2) = ivps(il7-2)
 IF (.NOT.prt) GO TO 860
 WRITE  (nout,850) nmvps
 850 FORMAT (1H+,70X,'(INVALID BCD WORD)= ',2A4)
 WRITE  (nout,830) uwm,nmvps
 
 860 IF (il4 <= 0 .AND. il6 <= 0 .AND. il8 <= 0 .AND. il9 <= 0) GO TO 1240
 
!     OUTPUT REQUEST - IL4 - S.P. REAL
!                      IL5 - INTEGER
!                      IL6 - D.P. REAL
!                      IL7 - BCD
!                      IL8 - S.P. COMPLEX
!                      IL9 - D.P. COMPLEX
 
 mat   = .true.
 tb1   = .false.
 tb2   = .false.
 tb4   = .false.
 recno = p2
 INDEX = p1
 rl    = 999999
 ii    = 1
 jj    = mcb(3)
 incr  = 1
 CALL gopen (in1,iz(ibuf1),0)
 CALL skprec (in1,col-1)
 CALL unpack (*1030,in1,z)
 SELECT CASE ( ityp )
   CASE (    1)
     GO TO 900
   CASE (    2)
     GO TO 910
   CASE (    3)
     GO TO 950
   CASE (    4)
     GO TO 950
 END SELECT
 
!     INPUT MATRIX PRECISION TYPE = 1, S.P. REAL
 
 900 GO TO 350
 
!     MATRIX PRECISION TYPE = 2, D.P. REAL
 
 910 IF (il4 <= 0) GO TO 920
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 rsp = dz(row)
 vps(il4) = rsp
 nmvps(1) = ivps(il4-3)
 nmvps(2) = ivps(il4-2)
 IF (prt) WRITE (nout,370) rsp,nmvps
 920 IF (il6 <= 0) GO TO 930
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 rdp   = dz(row)
 dp(1) = rdp
 vps(il6  ) = sp(1)
 vps(il6+1) = sp(2)
 nmvps(1) = ivps(il6-3)
 nmvps(2) = ivps(il6-2)
 IF (prt) WRITE (nout,580) rdp,nmvps
 930 IF (il8 <= 0) GO TO 940
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 splx(1) = dz(row)
 splx(2) = 0.0
 vps(il8  ) = splx(1)
 vps(il8+1) = splx(2)
 nmvps(1) = ivps(il8-3)
 nmvps(2) = ivps(il8-2)
 IF (prt) WRITE (nout,680) splx,nmvps
 940 IF (il9 <= 0) GO TO 1100
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 dp(1) = dz(row)
 dp(2) = 0.d0
 nmvps(1) = ivps(il9-3)
 nmvps(2) = ivps(il9-2)
 GO TO 720
 
!     INPUT MATRIX PRECISION TYPE = 3 OR 4, COMPLEX
 
 950 IF (il4 <= 0) GO TO 970
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 nmvps(1) = ivps(il4-3)
 nmvps(2) = ivps(il4-2)
 IF (.NOT.prt) GO TO 970
 WRITE  (nout,960) nmvps
 960 FORMAT (1H+,70X,' (INVALID S.P. REAL NUMBER)  = ',2A4)
 WRITE  (nout,830) uwm,nmvps
 970 IF (il6 <= 0) GO TO 990
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 nmvps(1) = ivps(il6-3)
 nmvps(2) = ivps(il6-2)
 IF (prt) WRITE (nout,980) nmvps
 980 FORMAT (1H+,70X,' (INVALID D.P.REAL NUMBER)  = ',2A4)
 990 IF (il8 <= 0 .AND. il9 <= 0) GO TO 1100
 IF (ityp == 4) GO TO 1010
 
!     INPUT MATRIX PRECISION TYPE = 3, S.P.COMPLEX
 
 IF (il8 <= 0) GO TO 1000
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 splx(1) = z(row  )
 splx(2) = z(row+1)
 vps(il8  ) = splx(1)
 vps(il8+1) = splx(2)
 nmvps(1) = ivps(il8-3)
 nmvps(2) = ivps(il8-2)
 IF (prt) WRITE (nout,680) splx,nmvps
 1000 IF (il9 <= 0) GO TO 1100
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 dp(1) = DBLE(z(row  ))
 dp(2) = DBLE(z(row+1))
 nmvps(1) = ivps(il9-3)
 nmvps(2) = ivps(il9-2)
 GO TO 720
 
!     INPUT MATRIX PRECISION TYPE = 4, D.P.COMPLEX
 
 1010 IF (il8 <= 0) GO TO 1020
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 splx(1) = SNGL(dz(row  ))
 splx(2) = SNGL(dz(row+1))
 vps(il8  ) = splx(1)
 vps(il8+1) = splx(2)
 nmvps(1) = ivps(il8-3)
 nmvps(2) = ivps(il8-2)
 IF (prt) WRITE (nout,680) splx,nmvps
 1020 IF (il9 <= 0) GO TO 1100
 IF (prt) WRITE (nout,810) row,col,TYPE(ityp),fnm
 dp(1) = dz(row  )
 dp(2) = dz(row+1)
 nmvps(1) = ivps(il9-3)
 nmvps(2) = ivps(il9-2)
 GO TO 720
 
!     NULL INPUT MATRIX ELEMENT
 
 1030 z (row  ) = 0.
 z (row+1) = 0.
 dz(row  ) = 0.d0
 dz(row+1) = 0.d0
 SELECT CASE ( ityp )
   CASE (    1)
     GO TO 900
   CASE (    2)
     GO TO 910
   CASE (    3)
     GO TO 950
   CASE (    4)
     GO TO 950
 END SELECT
 
 1100 CALL CLOSE (in1,1)
 GO TO 1240
 
!     ERRORS
 
 1200 j = -1
 GO TO 1230
 1210 j = -2
 GO TO 1230
 1220 j = -8
 1230 CALL mesage (j,in1,NAME)
 
 1240 RETURN
END SUBROUTINE paraml
