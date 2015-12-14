SUBROUTINE outpt4
     
!     COPY MATRIX DATA BLOCKS ONTO A FORTRAN TAPE, BINARY OR ASCII
!     FORMATS, IN DENSE MATRIX FORM (FROM FIRST TO LAST NON-ZERO TERMS
!     OF COLUMNS), OR IN SPARSE FORM (BY STRINGS)
 
!     A LOGICAL OUTPUT RECORD, WHICH CAN BE ONE OR MORE PHYSICAL RECORES
!     BEGINS WITH 3 INTEGER WORD THEN AN ARRAY OF DATA
 
!     FIRST  INTEGER WORD = LOGICAL RECORD NUMBER, OR COLUMN NUMBER
!     SECOND INTEGER WORD = ROW POSITION OF 1ST NONZERO TERM IN COLUMN
!                         = 0, SPARSE MATRIX (BINARY ONLY)
!                       .LT.0, SPARSE MATRIX ROW POSITION (ASCII ONLY)
!     THIRD  INTEGER WORD = NW, LENGTH OF ARRAY DATA THAT FOLLOW
!                           NW IS BASED ON S.P. WORD COUNT (BINARY ONLY)
!                           NW IS DATA PRECISION TYPE DEPENDENT (ASCII)
 
!     OUTPUT4 DOES NOT HANDLE TABLE DATA BLOCK, EXECPT 6 SPECIAL TABLES
!     KELM, MELM, BELM, KDICT, MDICT, AND BDICT.
 
 
!     OUTPUT4   IN1,IN2,IN3,IN4,IN5 // V,N,P1 / V,N,P2 / V,N,P3  $
 
!     PARAMETERS P1, P2 AND P3 ARE INTEGERS
 
!     P1 = 0, NO ACTION TAKEN BEFORE WRITE (DEFAULT)
!        =-1, REWIND TAPE BEFORE WRITE
!        =-2, AT END, WRITE E-O-F MARK AND REWIND TAPE
!        =-3, BOTH -1 AND -2
!        =-9, NOT AVAILABLE
 
!     P2 = N, FORTRAN OUTPUT UNIT N (N = 11,...,24)
!        =-N, MATRIX WILL BE WRITTEN OUT IN SPARSE FORMAT ONTO UNIT N.
 
!     P3 = 1, FILE OUTPUT IN FORTRAN BINARY FORMAT (UNFORMATTED)
!        = 2, FILE OUTPUT IN BCD FORMAT (ASCII, FORMATTED)
!          .  NO MIXED INTEGERS AND REAL NUMBERS IN A FORMATTED RECORD.
!             THE RECORD LENGTH IS LESS THAN 132 BYTES.
!          .  IF INPUT MATRIX TO BE COPIED OUT IS IN S.P., INTEGERS ARE
!             WRITTEN OUT IN I13, AND S.P.REAL DATA IN 10E13.6.
!          .  IF INPUT MATRIX TO BE COPIED OUT IS IN D.P., INTEGERS ARE
!             WRITTEN OUT IN I16, AND D.P.REAL DATA IN 8D16.9.
!        = 3, FORMATS I16 AHD 8E16.9 ARE USED TO COPY INTEGERS AND S.P.
!             REAL DATA OUT TO OUTPUT TAPE. P3=3 IS USED ONLY FOR
!             MACHINE WITH LONG WORDS (60 OR MORE BITS PER WORD)
 
!     THESE OUTPUT FORMATS CAN BE CHANGED EASILY BY ALTERING FORMATS
!     40, 50, 60 AND 370. MAKE SURE AN OUTPUT LINE DOES NOT EXCEED 132
!     COLUMNS. OTHERWISE, IT WOULD BE FOLDED IN PRINTOUT OR SCREEN
!     LISTING.
 
!     WRITTEN BY G.CHAN/UNISYS  3/93
 
 LOGICAL :: sparse,bo,sp,dp,cp
 INTEGER :: p1,p2,p3,buf1,d,zero,trl(8),NAME(2),NONE(2),  &
     ix(3),BLOCK(20),inp(13),sub(2),tab1(6),tab2(6)
 REAL :: xns(1)
 DOUBLE PRECISION :: dx(1),dxns(1)
 CHARACTER (LEN=6) :: dns,spa,ds
 CHARACTER (LEN=11) :: fmd,unf,fm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBI
 CHARACTER (LEN=80) :: dsnames
!WKBI
 COMMON / dsname / dsnames(80)
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /BLANK /  p1,p2,p3
 COMMON /system/  ibuff,nout,dum1(6),nlpp,dum2(2),line,dum3(2),  &
     d(3),dum22(22),nbpw
 COMMON /machin/  mach
 COMMON /unpakx/  itu,ii,jj,incr
 COMMON /TYPE  /  prc(2),nwd(4)
 COMMON /zzzzzz/  x(1)
 EQUIVALENCE      (x(1),xns(1))
 EQUIVALENCE      (x(1),ix(1),dx(1)),(xns(1),dxns(1)),(nm1,NAME(1))
 DATA    inp   /  4HUT1 ,4HUT2 ,4HUT3 ,4HINPT,4HINP1,4HINP2,4HINP3,  &
     4HINP4,4HINP5,4HINP6,4HINP7,4HINP8,4HINP9/
 DATA    tab1  /  4HKELM,4HMELM,4HBELM,4HKDIC,4HMDIC,4HBDIC/,  &
     tab2  /  4HHKEL,4HHMEL,4HHBEL,4HHKDI,4HHMDI,4HHBDI/
 DATA    dns   ,  spa  / 'DENSE ', 'SPARSE' /  rzero,zero  / 0.,0 /
 DATA    fmd   ,  unf  / 'FORMATTED  ','UNFORMATTED' /
 DATA    NONE  ,  sub  / 4H (no,4HNE) ,4HOUTP,4HT4   /
!WKBI
 DATA    ifirst / 0 /
 
 sparse = p2 < 0
 p2   = IABS(p2)
 IF (p2 <= 10 .OR. p2 > 24) GO TO 500
 bo   = p3 /= 1
 ii   = 1
 incr = 1
 lcor = korsz(x(1))
 buf1 = lcor - ibuff
 
 fm = unf
 IF (bo) fm = fmd
!WKBNB
 IF ( bo .OR. ifirst /= 0 ) GO TO 1
 CLOSE ( UNIT=p2 )
 OPEN (UNIT=p2,STATUS='NEW',ACCESS='SEQUENTIAL',FORM=fm,ERR=500  &
     ,FILE=dsnames(p2) )
 1     CONTINUE
 ifirst = 1
!WBKNE
 IF (p1 == -1 .OR. p1 == -3) REWIND p2
 
 DO  ipt = 1,5
   ndict = 0
   INPUT = 100 + ipt
   trl(1)= INPUT
   CALL rdtrl (trl(1))
   IF (trl(1) <= 0) CYCLE
   CALL fname (INPUT,NAME)
   IF (nm1 == NONE(1) .AND. NAME(2) == NONE(2)) CYCLE
   IF (trl(7) == 0 .AND. trl(8) == 0) GO TO 250
   IF (trl(4) < 1 .OR.  trl(4) > 8) GO TO 250
   nc   = trl(2)
   nr   = trl(3)
   itu  = trl(5)
   IF (nc == 0 .OR. nr == 0 .OR. (itu < 1 .OR. itu > 4)) GO TO 250
   nwds = nwd(itu)
   IF (nr*nwds >= buf1) CALL mesage (-8,lcor,sub)
   dp   = itu == 2 .OR. itu == 4
   sp   = .NOT.dp
   cp   = sp .AND. p3 == 3 .AND. nbpw >= 60
   IF (cp) sp = .false.
   IF (bo .AND. sparse .AND. nc > 2000) WRITE (nout,10) uwm
   10 FORMAT (a25,' FROM OUTPUT4 MODULE. ON ASCII TAPE AND SPARSE ',  &
       'MATRIX OUTPUT, EACH STRING OF DATA IS WRITTEN OUT TO THE',  &
       /5X,'OUTPUT TAPE AS A FORTRAN FORMATTED REDORD. FATAL ERROR',  &
       ' COULD OCCUR WHEN NO. OF RECORDS EXCEED SYSTEM I/O LIMIT')
   
!     OPEN INPUT DATA BLOCK TO READ WITH REWIND
   
   CALL OPEN (*520,INPUT,x(buf1),0)
   CALL fwdrec (*520,INPUT)
   
   BLOCK(1) = INPUT
   
!     WRITE TRAILER RECORD ON OUTPUT TAPE.
!     SET FORM (TRL(4)) TO NEGATIVE IF ASCII RECORDS IS REQUESTED
   
   k = -trl(4)
   IF (.NOT.bo) WRITE (p2   ) nc,nr,trl(4),itu,NAME
   IF (     bo) WRITE (p2,20) nc,nr,k     ,itu,NAME
   20 FORMAT (1X,4I13,5X,2A4)
   
   IF (sparse) GO TO 100
   
!     DENSE MATRIX OUTPUT -
!     WRITE THE MATRIX COLUMNS FROM FIRST TO LAST NON-ZERO TERMS
   
   30 DO  k = 1,nc
     ii = 0
     CALL unpack (*40,INPUT,x)
     jj = (jj-ii+1)*nwds
     IF (bo) GO TO 40
     
     WRITE (p2) k,ii,jj,(x(l),l=1,jj)
     CYCLE
     
     40 m = jj/2
     IF (sp) WRITE (p2,50) k,ii,jj,( x(l),l=1,jj)
     IF (cp) WRITE (p2,60) k,ii,jj,( x(l),l=1,jj)
     IF (dp) WRITE (p2,70) k,ii,jj,(dx(l),l=1,m )
     50 FORMAT (1X,3I13,/,(1X,10E13.6))
     60 FORMAT (1X,3I16,/,(1X, 8E16.9))
     70 FORMAT (1X,3I16,/,(1X, 8D16.9))
     
   END DO
   
   GO TO 210
   
!     SPARSE MATRIX OUPUT -
!     WRITE A RECORD FOR EACH MATRIX COLUMN, IN PACKED STRINGS DATA
!     IF MATRIX IS NOT WRITTEN IN STRINGS, SEND THE MATRIX TO THE DENSE
!     MATRIX METHOD
   
   100 CALL rectyp (INPUT,k)
   IF (k /= 0) GO TO 110
   CALL REWIND (INPUT)
   CALL fwdrec (*520,INPUT)
   GO TO 30
   
!     BLOCK(2) = STRING TYPE, 1,2,3 OR 4
!     BLOCK(4) = FIRST ROW POSITION ON A MATRIX COLUMN
!     BLOCK(5) = POINTER TO STRING IN XNS ARRAY
!     BLOCK(6) = NO. OF TERMS IN STRING
   
   110 nwords   = nwd(itu)
   nword1   = nwords - 1
   DO  k = 1,nc
     BLOCK(8) = -1
     nw = 0
     120 CALL getstr (*150,BLOCK)
     IF (bo) GO TO 160
     ln = BLOCK(6)*nwords
     j1 = BLOCK(5)*nwords - nword1
     j2 = j1 + ln - 1
     
     nw = nw + 1
     ix(nw) = BLOCK(4) + 65536*BLOCK(6)
     l  = 1
     DO  j = j1,j2
       x(l+nw) = xns(j)
       l  = l + 1
     END DO
     nw = nw + ln
     IF (nw >= buf1) CALL mesage (-8,lcor,sub)
     140 CALL endget (BLOCK)
     GO TO 120
     150 IF (nw > 0) WRITE (p2) k,zero,nw,(x(j),j=1,nw)
     CYCLE
     
!     NOTE - FOR THE BCD OUTPUT RECORD, THE 2ND INTEGER WORD, ZERO
!     BEFORE, IS REPLACED BY THE NEGATIVE OF THE ROW POSITION (IN A
!     MATRIX COLUMN).
!     DOUBLE THE POINTER TO THE STRING IN XNS ARRAY, J1, IF DATA TYPE IS
!     COMPLEX, BUT NOT THE LENGTH LN.
     
     160 ln = BLOCK(6)
     j1 = BLOCK(5)
     IF (BLOCK(2) >= 3) j1 = j1*2
     j2 = j1 + ln - 1
     mrow = -BLOCK(4)
     
!                   ZERO REPLACED   EXACT LENGTH OF XNS, OR DXNS
!                               /   /
     IF (sp) WRITE (p2,50) k,mrow,ln,( xns(j),j=j1,j2)
     IF (cp) WRITE (p2,60) k,mrow,ln,( xns(j),j=j1,j2)
     IF (dp) WRITE (p2,70) k,mrow,ln,(dxns(j),j=j1,j2)
     GO TO 140
     
   END DO
   
!     WRITE AN EXTRA NCOL+1 COLUMN RECORD OUT TO P2, AND AT LEAST ONE
!     VALUE OF ZERO
   
   210 m  = 1
   k  = nc + 1
   IF (bo) GO TO 220
   WRITE (p2) k,m,m,rzero
   GO TO 230
   220 IF (sp) WRITE (p2,50) k,m,m,rzero
   IF (cp .OR. dp) WRITE (p2,60) k,m,m,rzero
   230 ds = dns
   IF (sparse) ds = spa
   WRITE  (nout,240) uim,NAME,p2,inp(p2-10),fm,ds,(trl(l),l=2,7)
   240 FORMAT (a29,' FROM OUTPUT4 MODULE. DATA BLOCK ',2A4,' WAS WRITTEN'  &
       ,       ' OUT TO FORTRAN TAPE',i3,' (',a4,')', /5X,'IN ',a11,  &
       ' RECORDS. ',a6,' MATRIX FORM.  TRAILER =',5I6,i9)
   GO TO 280
   
!     INPUT FILE IS A TABLE DATA BLOCK
!     ONLY 6 SPECIAL TABLES ARE ALLOWED
   
   250 DO  i = 1,6
     IF (nm1 == tab1(i) .OR. nm1 == tab2(i)) GO TO 290
   END DO
   IF (bo) WRITE (p2,270) uwm,INPUT,NAME,(trl(j),j=2,7)
   WRITE  (nout,270) uwm,INPUT,NAME,(trl(j),j=2,7)
   270 FORMAT (a25,'. INPUT DATA BLOCK',i5,2H, ,2A4,', IS A TABLE OR A ',  &
       'NULL MATRIX. OUTPUT4 MODULE HANDLES ONLY MATRICES', /5X, 'TRAILER =',6I6)
   280 CALL CLOSE (INPUT,1)
   CYCLE
   
!     KELM, MELM AND BELM (AND HKELM, HMELM AND HBELM) TALBES
   
   290 IF (sparse) WRITE (nout,300) uwm,NAME,p2,p3
   300 FORMAT (a25,'. PARAMETER P2 FOR SPARSE MATRIX IS MEANINGLESS FOR',  &
       ' THE ',2A8,' INPUT FILE.   P2,P3 =',2I4,/)
   CALL OPEN (*520,INPUT,x(buf1),0)
   CALL fwdrec (*520,INPUT)
   k = -trl(4)
   IF (.NOT.bo) WRITE (p2   ) nc,nr,trl(4),itu,NAME
   IF (     bo) WRITE (p2,20) nc,nr,k     ,itu,NAME
   j = 1
   k = 0
   IF (i >= 4) GO TO 310
   dp = trl(2) == 2
   sp = .NOT.dp
   cp = sp .AND. p3 == 3 .AND. nbpw >= 60
   IF (cp) sp = .false.
   310 k = k + 1
   CALL READ (*380,*320,INPUT,x,buf1-1,1,m)
   CALL mesage (-8,0,sub)
   320 IF (i >= 4) GO TO 350
   IF (bo) GO TO 330
   WRITE (p2) k,j,m,(x(l),l=1,m)
   GO TO 310
   
   330 IF (dp) GO TO 340
   IF (sp) WRITE (p2,50) k,j,m,( x(l),l=1,m)
   IF (cp) WRITE (p2,60) k,j,m,( x(l),l=1,m)
   GO TO 310
   340 m = m/2
   WRITE (p2,70) k,j,m,(dx(l),l=1,m)
   GO TO 310
   
!     KDICT, MDICT AND BDICT (AND HKDICT, HMDICT AND HBDICT) TABLES.
!     INTEGERIZE THE DAMPING CONSTANT (BY 10**8) BEFORE OUTPUT THE ARRAY
   
   350 ndict = ix(3) + 5
   DO  i = 8,m,ndict
     ix(i) = IFIX(x(i)*100000000.)
   END DO
   IF (.NOT.bo) WRITE (p2) k,j,m,(ix(l),l=1,m)
   IF (bo .AND. .NOT.cp) WRITE (p2,370) k,j,m,(ix(l),l=1,m)
   IF (bo .AND.      cp) WRITE (p2,375) k,j,m,(ix(l),l=1,m)
   370 FORMAT (1X,3I13,/,(1X,10I13))
   375 FORMAT (1X,3I13,/,(1X, 8I16))
   GO TO 310
   
   380 IF (.NOT.bo) WRITE (p2) k,j,j,zero
   IF (bo .AND. .NOT.cp) WRITE (p2,370) k,j,j,zero
   IF (bo .AND.      cp) WRITE (p2,375) k,j,j,zero
   IF (ndict /= 0) WRITE  (nout,390) uim,NAME,inp(p2-10)
   390 FORMAT (a29,'. THE DAMPING CONSTANT TERMS FROM ',2A4,' WERE ',  &
       'MULTIPLIED BY 10**8, AND INTEGERIZED', /5X,  &
       'BEFORE WRITING OUT TO ',a4,' OUTPUT FILE')
   GO TO 280
   
 END DO
 
 IF (p1 /= -2 .AND. p1 /= -3) GO TO 600
 ENDFILE p2
 REWIND  p2
 CLOSE (UNIT=p2)
 GO TO 600
 
!     ERRORS
 
 500 WRITE  (nout,510) ufm,p2
 510 FORMAT (a23,'. CANNOT OPEN OUTPUT FORTRAN FILE. UNIT =',i4)
 GO TO  540
 520 WRITE  (nout,530) uwm,INPUT
 530 FORMAT (a25,'. OUTPT4 CANNOT OPEN INPUT DATA BLOCK',i5)
 
 540 CALL mesage (-37,0,sub)
 600 RETURN
END SUBROUTINE outpt4
