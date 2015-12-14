SUBROUTINE tabprt (iname1)
     
!     WILL PRINT TABLE - USING 1P,E13.6, I13, OR (9X,A4) FORMAT
 
!     ALL REAL NUMBERS ARE ASSUMED TO BE SINGLE PRECISION.
 
!     REVISED  3/91 BY G.CHAN/UNISYS
!     THREE PARAMETERS ARE ADDED - OP CODE (OP), RECORDD NO. (IRC), AND
!     WORD NO. (IWD)
!     THE DEFAULTS OF THESE PARAMETERS ARE - BLANK, 3, AND 3
!     OP CODE OPTIONS ARE 'PUREBCD', 'PUREFPN', AND 'PUREINT'
 
!     LAST REVISED, 12/92, BY G.CHAN/UNISYS, TO INCLUDE 3 SPECIAL TABLES
!     - KELM, MELM, BELM  - WHICH CONTAIN D.P. DATA WORDS IN 32- AND 36-
!     BIT WORD MACHINES.
 
!     IF OP CODE IS 'PUREBCD', RECORDS IRC AND THEREAFTER, AND BEGINNING
!     FROM WORD IWD OF EACH RECORD TO THE END OF THAT RECORD, ARE ALL
!     BCD  WORDS.
!     SIMILARILY FOR 'PUREINT' FOR INTEGER WORDS, AND 'PUREFPN' FOR
!     FLOATING POINT NUMBERS
 
!     THESE PARAMETER OPTIONS ARE NECESSARY BECAUSE IF THE PRINTED DATA
!     IS NOT OF STRING TYPE, SUBROUTINE NUMTYP IS CALLED TO FIND OUT
!     WHAT TYPE OF DATA IN EACH DATA WORD.  HOWEVER NUMTYP IS NOT 100
!     PERCENT FOOL-PROOF. ONCE IN A FEW THOUSANDS NUMTYP CAN NOT
!     DISTINGUISH A REAL NUMBER FROM A BCD WORD
 
!  $MIXED_FORMATS
 
 
 INTEGER, INTENT(IN)                      :: iname1
 LOGICAL :: dec
 INTEGER :: BLOCK(20),types(4),FORMAT,forms(2),jpoint,row,  &
     TYPE,flag,recf,strnbr,sysbuf,otpe,pure,bcd,fpn, op,NAME(2),icore(133)
 REAL :: xns(1),sp(3)
 DOUBLE PRECISION :: xnd(1),dcore(1)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBI
 CHARACTER (LEN=1) :: core1(2000)
 COMMON /xmssg /  ufm,uwm
 COMMON /BLANK /  op(2),irc,iwd
 COMMON /machin/  mach
 COMMON /system/  sysbuf,otpe,inx(6),nlpp,inx1(2),line,dum(42),iprc
 COMMON /output/  head1(96),head2(96)
 COMMON /zzzzzz/  core(1)
 EQUIVALENCE      (xnd(1),core(1))
 EQUIVALENCE      (xns(1),xnd(1)),  (icore(1),core(1),dcore(1)),  &
     (BLOCK(2), TYPE   ), (BLOCK(3), FORMAT),  &
     (BLOCK(4), row    ), (BLOCK(5), jpoint),  &
     (BLOCK(6), nterms ), (BLOCK(8), flag  )
!WKBI
 EQUIVALENCE      (core, core1)
 DATA    oparen,  cparen,ec,ec1,ec2,intgc,alphc,alphc1,cont,uned /  &
     4H(1X ,  4H)   ,4H,1P,,4HE13.,2H6 ,4H,i13,4H,9X,,4HA4   ,  &
     4HCONT,  4HINUE   /  d/2HD  /, NAME  / 4HTABP,4HRT      /
 DATA    BLANK ,  tabl,ebb /  1H  ,4HTABL, 1HE     /
 DATA    types /  3HRSP,3HRDP,3HCSP ,3HCDP/, forms / 3HYES ,2HNO /
 DATA    pure  ,  bcd,fpn,INT /   4HPURE, 4HBCD , 4HFPN , 4HINT  /
 DATA    nsp   ,  sp / 3, 4HKELM, 4HMELM, 4HBELM   /
 
 nz  = korsz(core) - sysbuf
 IF (nz <= 0) CALL mesage (-8,-nz,NAME)
 dec = mach == 5 .OR. mach == 6 .OR. mach == 10 .OR. mach == 21
 iname = iname1
 CALL OPEN (*190,iname,core(nz+1),0)
 DO  i = 1,96
   head2(i) = BLANK
 END DO
 head2(1) = tabl
 head2(2) = ebb
 CALL fname (iname,head2(3))
 CALL page
 head2(6) = cont
 head2(7) = uned
 head2(8) = d
 IF (iprc == 1 .OR. iname /= 101) GO TO 15
 CALL page2 (-2)
 WRITE  (otpe,13) uwm
 13 FORMAT (a25,', TABPRT MODULE ASSUMES ALL REAL DATA ARE IN S.P.,',  &
     ' D.P. DATA THEREFORE MAY BE PRINTED ERRONEOUSLY')
 15 inum     = nz/2 - 1
 inum     = MAX0(inum,133)
 ns       = inum + 1
 llen     = 0
 core(1)  = oparen
 irec     = 0
 ircd     = 999999999
 ixxx     = 999999999
 IF (op(1) /= pure .OR. op(2) == BLANK) GO TO 20
 IF (op(2) == INT) jj = 2
 IF (op(2) == fpn) jj = 3
 IF (op(2) == bcd) jj = 4
 IF (irc > 0) ircd = irc
 IF (iwd > 0) ixxx = iwd + inum
 IF (irc <= 0) ircd = 3
 IF (iwd <= 0) ixxx = 3 + inum
 20 CALL page2 (-2)
 IF (dec .AND. irec == 0) WRITE (otpe,25)
 25 FORMAT (4X,'(ALL INTEGERS EXCEEDING 16000 ARE PRINTED AS REAL ',  &
     'NUMBERS. ALL REAL NUMBERS OUTSIDE E-27 OR E+27 RANGE ',  &
     'ARE PRINTED AS INTEGERS)')
 WRITE  (otpe,30) irec
 30 FORMAT (/,' RECORD NO.',i6)
 irec = irec + 1
 DO  i  = 1,nsp
   IF (head2(3) /= sp(i)) CYCLE
   icore(1) = iname
   CALL rdtrl (icore)
   IF (icore(2) == 2) GO TO 60
 END DO
 ix   = inum
 nred = 0
 np   = inum - 1
 BLOCK(1) = iname1
 CALL rectyp (BLOCK,recf)
 IF (recf /= 0) GO TO 200
 jv   = 4
 40 ix   = ix + 1
 iout = 4
 nred = nred + 1
 np   = np + 1
 CALL READ (*170,*160,iname,core(ix),1,0,iflag)
 
 IF (irec > ircd .OR. ix > ixxx) GO TO 50
 jj = numtyp(icore(ix)) + 1
 IF (jj == 1 .AND. jv /= 4) jj = jv
 jv = jj
 50 SELECT CASE ( jj )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 140
   CASE (    3)
     GO TO 100
   CASE (    4)
     GO TO 120
 END SELECT
 
!     TABLES KELM, MELM, AND BELM - D.P. DATA ONLY
 
 60 CALL READ (*170,*170,iname,core(1),2,1,iflag)
 WRITE  (otpe,65) icore(1),icore(2)
 65 FORMAT (10X,2A4)
 70 WRITE  (otpe,30) irec
 CALL READ (*170,*80,iname,core(1),nz,1,iflag)
 CALL mesage (-8,0,NAME)
 80 np   = iflag/2
 jj   = (np+9)/10
 CALL page2 (-jj)
 irec = irec + 1
 WRITE  (otpe,90,ERR=70) (dcore(i),i=1,np)
 90 FORMAT (1X,1P,10D13.6)
 GO TO 70
 
!     REAL NUMBER  (1)
 
 100 iout = 1
 IF (llen+13 > 132) GO TO 160
 110 core(nred+1) = ec
 core(nred+2) = ec1
 core(nred+3) = ec2
 nred = nred + 2
 115 llen = llen + 13
 GO TO 40
 
!     ALPHA  (2)
 
 120 iout = 2
 IF (llen+6 > 132) GO TO 160
 130 core(nred+1) = alphc
 core(nred+2) = alphc1
 nred = nred + 1
 GO TO 115
 
!     INTEGER  (3)
 
 140 iout = 3
 IF (llen+13 > 132) GO TO 160
 150 icore(nred+1) = intgc
 GO TO 115
 
!     BUFFER FULL- END RECORD AND PRINT THE LINE
 
!     PREVIOUSLY, THE FORMAT IS IN CORE, WHICH IS DIMENSIONED TO 1.
!     THIS MAY NOT WORK IN SOME MACHINES. THE FORMAT IS NOW SPECIFIED IN
!     ICORE, WHICH IS DIMENSIONED TO 133.
!     (CORE AND ICORE ARE EQUIVALENT)
 
 160 core(nred+1) = cparen
 IF (nred >= 133) CALL mesage (-37,0,NAME)
 CALL page2 (-1)
 IF (nred == 1) GO TO 165
 IF (mach /= 2 .AND. mach /= 5 ) GO TO 162
 WRITE (otpe,icore,ERR=164) (core(i),i=ns,np)
 GO TO 164
 162 CALL wrtfmt (icore(ns),np-ns+1,core1)
 164 CONTINUE
 llen = 0
 nred = 1
 np   = inum
 
!     FINISH SEMI-PROCESSED WORD.
 
 core(inum+1) = core(ix)
 ix = inum + 1
 SELECT CASE ( iout )
   CASE (    1)
     GO TO 110
   CASE (    2)
     GO TO 130
   CASE (    3)
     GO TO 150
   CASE (    4)
     GO TO 20
 END SELECT
 
 165 WRITE  (otpe,166)
 166 FORMAT (' THIS RECORD IS NULL.')
 
!     GO TO 161 IS LOGICALLY UNSOUND. CHANG TO 164. (G.CHAN/UNISYS 1/93)
!     GO TO 161
!WKBR GO TO 164
 GO TO 162
 
 170 CALL CLOSE (iname,1)
 CALL page2 (-2)
 WRITE  (otpe,180)
180 FORMAT (//,' END OF FILE')

!     PRINT TRAILER FOR FILE

190 icore(1) = iname
CALL rdtrl (icore)
CALL page2 (-2)
WRITE  (otpe,195) (icore(i),i=2,7)
195 FORMAT ('0TRAILER WORD1 =',i8,' WORD2 =',i8,' WORD3 =',i8,  &
    ' WORD4 =',i8,' WORD5 =',i8,' WORD6 =',i8)
RETURN


!     HERE IF STRING FORMATTED RECORD

200 flag   =-1
strnbr = 1
CALL getstr (*250,BLOCK)
iform = FORMAT + 1
205 CALL page2 (-2)
WRITE  (otpe,206) strnbr,row,types(TYPE),forms(iform),nterms
206 FORMAT ('0STRING NO.',i5,'   ROW POSITION=',i5,'   STRING TYPE=',  &
    a3,'   STRING TRAILERS=',a3,'   NUMBER OF TERMS=',i5)
strnbr = strnbr + 1
SELECT CASE ( TYPE )
  CASE (    1)
    GO TO 210
  CASE (    2)
    GO TO 220
  CASE (    3)
    GO TO 230
  CASE (    4)
    GO TO 240
END SELECT

!     PRINT REAL SINGLE PRECISION STRING

210 npoint = jpoint + nterms - 1
j = jpoint
211 n = MIN0(j+7,npoint)
CALL page2 (-1)
WRITE  (otpe,212) (xns(i),i=j,n)
212 FORMAT (1X,8(1P,e15.7))
IF (n == npoint) GO TO 214
j = n + 1
GO TO 211
214 CALL endget (BLOCK)
CALL getstr (*20,BLOCK)
GO TO 205

!     PRINT STRING IN REAL DOUBLE PRECISION

220 npoint = jpoint + nterms - 1
j = jpoint
221 n = MIN0(j+7,npoint)
CALL page2 (-1)
WRITE  (otpe,222) (xnd(i),i=j,n)
222 FORMAT (1X,8(1P,d15.7))
IF (n == npoint) GO TO 224
j = n + 1
GO TO 221
224 CALL endget (BLOCK)
CALL getstr (*20,BLOCK)
GO TO 205

!     PRINT STRING IN COMPLEX SINGLE PRECISION

230 npoint = jpoint + 2*nterms - 1
j = jpoint
231 n = MIN0(j+7,npoint)
CALL page2 (-1)
WRITE  (otpe,232) (xns(i),i=j,n)
232 FORMAT (1X,4(1P,e14.7,1P,e15.7,2H//))
IF (n == npoint) GO TO 234
j = n + 1
GO TO 231
234 CALL endget (BLOCK)
CALL getstr (*20,BLOCK)
GO TO 205

!     PRINT STRING IN COMPLEX DOUBLE PRECISION

240 npoint = jpoint + 2*nterms - 1
j = jpoint
241 n = MIN0(j+7,npoint)
CALL page2 (-1)
WRITE  (otpe,242) (xnd(i),i=j,n)
242 FORMAT (1X,4(1P,d14.7,1P,d15.7,2H//))
IF (n == npoint) GO TO 244
j = n + 1
GO TO 241
244 CALL endget (BLOCK)
CALL getstr (*20,BLOCK)
GO TO 205

!     PRINT NULL COLUMN

250 CALL page2 (-1)
WRITE  (otpe,252)
252 FORMAT (5X,'NULL COLUMN')
GO TO 20

END SUBROUTINE tabprt
