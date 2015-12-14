SUBROUTINE tabpch
     
!     THE TABPCH MODULE WILL PUNCH UP TO 5 TABLES INTO DTI CARDS
 
!     DMAP CALL IS
 
!     TABPCH  IN1,IN2,IN3,IN4,IN5//P1,P2,P3,P4,P5
 
!     SINGLE FIELD CARDS WILL BE MADE UNLESS REAL NUMBERS ARE TO BE MADE
!     ALL REAL NUMBERS ARE ASSUMED TO BE SINGLE PRECISION.
 
!     LAST REVISED, 3/93, BY G.CHAN/UNISYS
!     PUNCH KELM, MELM AND BELM IN D.P. IF THESE DATA BLOCKS ARE IN D.P.
 
!  $MIXED_FORMATS
 
 INTEGER :: sysbuf    ,iz(10)    ,ifnm(5)   ,NAME(2)   ,  &
     mcb(7)    ,FILE      ,tabnm(2)  ,dti(2)    ,  &
     dtis(2)   ,idata(20) ,endrec(2) ,out       ,  &
     iform(20) ,BLANK     ,INT(2)    ,ireal(2)  ,  &
     ll(4)     ,intd(2)   ,pform(30) ,ibcd(2)   ,  &
     sp(3)     ,ibcdd(2)  ,FORM(30,2),forms(30,2)
 REAL :: rdata(20)
 DOUBLE PRECISION :: dz(1)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm       ,uwm       ,uim
 COMMON /machin/  mach
 COMMON /system/  sysbuf    ,out       ,ksystm(88),lpch
 COMMON /zzzzzz/  z(1)
 COMMON /BLANK /  n1(2,5)
 EQUIVALENCE      (z(1),iz(1),dz(1)),  (idata(1),rdata(1))
 DATA    BLANK /  1H             /
 DATA    dti   /  4HDTI , 1H     /
 DATA    dtis  /  4HDTI*, 1H     /
 DATA    endrec/  4HENDR, 4HEC   /
 DATA    forms /  4H(2A4, 26*2H  ,4H,1H+ ,4HA2,i, 4H5)  , 4H(a1,,  &
     4HA2,i ,4H5    ,24*2H  ,4H,1H+, 4HA2,i, 4H5)  /
 DATA    ibcd  /  4H,2A4, 1H     /
 DATA    ibcdd /  4H,2A4, 4H,8X  /
 DATA    ifnm  /  101, 102, 103, 104, 105/
 DATA    INT   /  4H,i8 , 1H     /
 DATA    intd  /  4H,i16, 1H     /
 DATA    iplus /  1H+            /
 DATA    ireal /  4H,e16, 4H.9   /
 DATA    istar /  1H*            /
 DATA    NAME  /  4HTABP, 4HCH   /
 DATA    ll    /  1, 1, 3, 2     /
 DATA    nsp   ,  sp  / 3, 4HKELM, 4HMELM, 4HBELM /
 
 nz    = korsz(z)
 ibuf  = nz - sysbuf + 1
 nz    = ibuf - 1
 icrq  = 10 - nz
 IF (nz <= 10) GO TO 830
 nread = nz/2  - 2
 nlist = nread + 3
 DO  j = 1,2
   DO  i = 1,30
     FORM(i,j) = forms(i,j)
   END DO
 END DO
 
!     FOR EACH  TABLE DEFINED
 
 ns = -1
 DO  i = 1,5
   mcb(1) = ifnm(i)
   CALL rdtrl (mcb)
   IF (mcb(1) <= 0) CYCLE
   
!     TABLE EXISTS SET IT UP
   
   FILE = ifnm(i)
   CALL OPEN  (*800,FILE,iz(ibuf),0)
   CALL fname (FILE,tabnm)
   io  = 0
   kmb = 4
   IF (mcb(5) == 1 .OR. mcb(5) == 3) GO TO 40
   DO  j = 1,nsp
     IF (kmb == 1 .OR. tabnm(1) /= sp(j)) CYCLE
     kmb = 1
     io  = 1
     nread = nz -1
   END DO
   IF (ns /= -1) GO TO 40
   ns = 1
   CALL page1
   WRITE  (out,30) uwm
   30 FORMAT (a25,', MODULE TABPCH ASSUMES ALL REAL DATA ARE IN S.P..',  &
       '  D.P. DATA THEREFORE MAY BE PUNCHED ERRONEOUSLY')
   IF (mach == 5 .OR. mach == 6 .OR. mach == 10 .OR. mach == 21)  &
       WRITE (out,35)
   35 FORMAT (4X,'(ALL INTEGERS EXCEEDING 16000 ARE PUNCHED AS REAL ',  &
       'NUMBERS. ALL REAL NUMBERS OUTSIDE E-27 OR E+27 RANGE ',  &
       'ARE PUNCHED AS INTEGERS)')
   
   40 CALL READ (*810,*820,FILE,iz(1),-2,0,ilen)
   irecno = 0
   ichr   = n1(1,i)
   iz(3)  = 0
   
!     SET UP FIRST RECORD
   
   iz(1) = tabnm(1)
   iz(2) = tabnm(2)
   iz(4) = mcb(2)
   iz(5) = mcb(3)
   iz(6) = mcb(4)
   iz(7) = mcb(5)
   iz(8) = mcb(6)
   iz(9) = mcb(7)
   CALL READ (*700,*50,FILE,iz(10),nread,0,ilen)
   icrq  = nread
   GO TO 830
   50 ilen  = ilen + 11
   60 iz(ilen-1) = endrec(1)
   iz(ilen  ) = endrec(2)
   GO TO 90
   
!     BRING IN NEXT RECORD
   
   70 CALL READ (*700,*80,FILE,iz(kmb),nread,io,ilen)
   icrq  = nread
   GO TO 830
   80 IF (kmb == 1) GO TO 600
   iz(3) = iz(3) + 1
   IF (ilen == 0) GO TO 70
   ilen  = ilen + 5
   GO TO 60
   
!     BUILD FORMAT VECTOR  1= INTEGER, 2 =BCD, 3=REAL
   
   90 jv = 3
   DO  k = 1,ilen
     m  = nlist + k - 1
     j  = numtyp(iz(k))
     IF (j == 0 .AND. jv /= 3) j = jv
     iz(m) = ll(j+1)
     jv = j
   END DO
   
!     MOVE DATA/FORMAT TO DATA AREA 8 FIELDS AT A TIME--SET D.F. FLAG
   
   id   = 1
   IF   = nlist
   ifrs = 1
   
!     HERE FOR EIGHT MORE WORDS
   
   110 idf = 0
   idt = 1
   ift = 1
   nf  = 1
   
!     HERE  FOR EACH FIELD
   
   120 idata(idt) = iz(id)
   iform(ift) = iz(IF)
   IF (iform(ift) == 3) idf = 1
   IF (iform(ift) /= 2) GO TO 140
   
!     BCD IS TWO WORDS
   
   idata(idt+1) = iz(id+1)
   
!     MAY BE FALSE BCD, CHECK FORMAT OF SECOND WORD ALSO
!     (SOME REAL NUMBER BIT PATTERNS LOOK LIKE BCD).
   
   IF (iz(IF+1) == 2) GO TO 130
   
!     SECOND WORD IS NOT BCD, ASSUME FIRST WORD IS REAL.
   
   idf = 1
   iform(ift) = 3
   GO TO 140
   130 idt = idt + 2
   ift = ift + 1
   id  = id  + 2
   IF  = IF  + 2
   GO TO 150
   
!     REAL OR INTEGER
   
   140 idt = idt + 1
   ift = ift + 1
   id  = id  + 1
   IF  = IF  + 1
   
!     BUMP FIELD COUNTER
   
   150 nf = nf + 1
   IF (nf >    8) GO TO 160
   IF (id < ilen) GO TO 120
   
!     FILL  WITH BLANKS
   
   idata(idt  ) = BLANK
   idata(idt+1) = BLANK
   iform(ift  ) = 2
   GO TO 130
   
!     PUNCH OUT 8 FIELDS OF DATA
   
   160 idt = 0
   IF (idf /= 0) GO TO 400
   
!     SINGLE FIELD CARD
   
   nf = 1
   170 m  = 2*nf + 2
   IF (iform(nf)-2 < 0) THEN
     GO TO   180
   ELSE IF (iform(nf)-2 == 0) THEN
     GO TO   200
   ELSE
     GO TO   210
   END IF
   
!     INTEGER
   
   180 FORM(m  ,ifrs) = INT(1)
   FORM(m+1,ifrs) = INT(2)
   
!     GET NEXT ITEM
   
   idt = idt + 1
   190 nf  = nf  + 1
   IF (nf <= 8) GO TO 170
   GO TO 220
   
!     BCD
   
   200 FORM(m  ,ifrs) = ibcd(1)
   FORM(m+1,ifrs) = ibcd(2)
   idt = idt + 2
   GO TO 190
   
!     REAL NOT LEGAL
   
   210 ip1 = -37
   GO TO 850
   
!     PUNCH OUT SINGLE CARD
   
   220 IF (ifrs /= 1) GO TO 270
   DO  j = 1,30
     pform(j) = FORM(j,1)
   END DO
   WRITE (lpch,pform,ERR=240) dti,(rdata(m),m=1,idt),ichr,irecno
   240 irecno = irecno + 1
   ifrs = 2
   DO  j = 1,30
     FORM(j,1) = forms(j,1)
   END DO
   260 IF (id >= ilen) GO TO 70
   GO TO 110
   
!     CONTINUATION CARD
   
   270 ircnm1 = irecno - 1
   DO  j = 1,30
     pform(j) = FORM(j,2)
   END DO
   WRITE (lpch,pform,ERR=290) iplus,ichr,ircnm1,(rdata(m),m=1,idt),  &
       ichr,irecno
   290 irecno = irecno + 1
   DO  j = 1,30
     FORM(j,2) = forms(j,2)
   END DO
   GO TO 260
   
!     DOUBLE FIELD CARDS
   
   400 nf = 1
   is = 1
   it = 4
   idt= 0
   m  = 2
   410 m  = m + 2
   IF (iform(nf)-2 < 0) THEN
     GO TO   420
   ELSE IF (iform(nf)-2 == 0) THEN
     GO TO   450
   ELSE
     GO TO   460
   END IF
   
!     INTEGER
   
   420 FORM(m  ,ifrs) = intd(1)
   FORM(m+1,ifrs) = intd(2)
   430 idt = idt + 1
   440 nf  = nf  + 1
   IF (m <= 8) GO TO 410
   GO TO 470
   
!     BCD
   
   450 FORM(m  ,ifrs) = ibcdd(1)
   FORM(m+1,ifrs) = ibcdd(2)
   idt = idt + 2
   GO TO 440
   
!     REAL
   
   460 FORM(m  ,ifrs) = ireal(1)
   FORM(m+1,ifrs) = ireal(2)
   GO TO 430
   
!     PUNCH OUT DOUBLE FIELD CARD
   
   470 IF (ifrs /= 1) GO TO 520
   DO  j = 1,30
     pform(j) = FORM(j,1)
   END DO
   WRITE (lpch,pform,ERR=490) dtis,(rdata(m),m=is,idt),ichr,irecno
   490 irecno = irecno + 1
   DO  j = 1,30
     FORM(j,1) = forms(j,1)
   END DO
   ifrs = 2
   510 it = 8
   m  = 2
   is = idt + 1
   GO TO 410
   
!     CONTINUATION CARD
   
   520 ircnm1 = irecno - 1
   DO  j = 1,30
     pform(j) = FORM(j,2)
   END DO
   WRITE (lpch,pform,ERR=540) istar,ichr,ircnm1,(rdata(m),m=is,idt),  &
       ichr,irecno
   540 irecno = irecno + 1
   DO  j = 1,30
     FORM(j,2) = forms(j,2)
   END DO
   IF (it == 4) GO TO 510
   GO TO 260
   
!     PUNCH KELM, MELM AND BELM IN D.P.
   
   600 IF (ilen == 0) GO TO 70
   ilen = ilen/2
   je = 0
   610 jb = je + 1
   je = je + 4
   ircnm1 = irecno
   irecno = irecno + 1
   IF (je >= ilen) GO TO 630
   WRITE  (lpch,620,ERR=840) ichr,ircnm1,(dz(j),j=jb,je),ichr,irecno
   620 FORMAT (1H*,a2,i5,1P,4D16.9,1X,a2,i5)
   GO TO 610
   630 je = ilen
   WRITE  (lpch,640,ERR=840) ichr,ircnm1,(dz(j),j=jb,je)
   640 FORMAT (1H*,a2,i5,1P,4D16.9)
   GO TO 70
   
!     CLOSE OFF FILES
   
   700 CALL CLOSE (FILE,1)
   CALL page2 (2)
   WRITE  (out,710) uim,tabnm,irecno
   710 FORMAT (a29,' 4015, TABLE ',2A4,' WAS PUNCHED OUT,',i8,' CARDS.')
 END DO
 WRITE  (lpch,730)
 730 FORMAT (1H , /,1H , /,1H )
 RETURN
 
!     ERROR MESAGES
 
 800 ip1 = -1
 GO TO 850
 810 ip1 =-2
 GO TO 850
 820 ip1 =-3
 GO TO 850
 830 ip1 = -8
 FILE = icrq
 GO TO 850
 840 ip1 = -37
 
 850 CALL mesage (ip1,FILE,NAME)
 RETURN
END SUBROUTINE tabpch
