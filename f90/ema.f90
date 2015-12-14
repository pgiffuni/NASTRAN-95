SUBROUTINE ema
     
!     DMAP SEQUENCE
 
!     EMA    GPECT,XEMD,XBLOCK/XGG/C,N,NOK4/C,N,WTMASS $
 
!            WHERE NOK4 .NE. -1 TO BUILD K4GG (USE DAMPING FACTOR),
!                       .EQ. -1 TO IGNORE DAMPING FACTOR
 
!     EMA USES TWO SCRATCH FILES
 
 EXTERNAL        lshift,rshift,andf  ,orf
 LOGICAL :: first ,last  ,piez
 INTEGER :: buf1  ,buf2  ,scrin ,scrout,scr1  ,scr2  ,gpect ,  &
     buf   ,sysbuf,rdrew ,wrtrew,rd    ,wrt   ,cls   ,  &
     clsrew,z     ,buf3  ,xemd  ,output,gpewds,gpectx,  &
     high  ,xblock,scalas,xgg   ,andf  ,prec  ,hdr   ,  &
     ppoint,union ,orf   ,op    ,rshift,col   ,oldcod,  &
     col1  ,coln  ,q     ,openr ,openw
 DOUBLE PRECISION :: zd    ,xd    ,d
!WKBI 1/95
 DOUBLE PRECISION :: xdd
 DIMENSION       buf(100)     ,scalas(32)   ,hdr(6),msg(4),mcb(7),  &
     ma1h(2)      ,zd(1) ,xd(1) ,xs(2) ,is(2) ,y(1)  , d(18) ,ihq(180)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm   ,uwm   ,uim   ,sfm
 COMMON /machin/ mach  ,ihalf ,jhalf
 COMMON /lhpwx / lhpw(5)      ,kshift
 COMMON /BLANK / nok4  ,wtmass
 COMMON /system/ ksystm(100)
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew,cls
 COMMON /zblpkx/ q(4)  ,iq
 COMMON /zzzzzz/ z(1)
 COMMON /ma1xx / buf   ,buf1  ,buf2  ,buf3  ,col   ,coln  ,col1  ,  &
     gpewds,high  ,i     ,icol  ,icolx ,idict ,ielem ,  &
     igpx  ,ilist ,imat  ,imatn ,ipvt  ,irow  ,irowp ,  &
     irowx ,jj    ,k     ,kelem ,kk    ,ii    ,k4flag,  &
     l     ,j     ,lcore ,ldict ,low   ,l1    ,l2    ,  &
     m     ,maxii ,maxipv,maxn  ,maxnpr,jnext ,n     ,  &
     nbrcol,nbrrow,nbrwds,ncol  ,ngps  ,ngrids,nhdr  ,  &
     nlist ,nlocs ,nmat  ,npvt  ,nread ,nrec  ,nrow  ,  &
     nrowp ,nsca  ,nwds  ,oldcod,op    ,prec  ,scrin , scrout, union,scalas
 EQUIVALENCE    (ksystm( 1),sysbuf), (ksystm( 2),output), (ksystm(78),ipiez )
 EQUIVALENCE    (z(1) ,zd(1)), (xs(1),xd(1),is(1)), (z(1) ,y(1) ),  &
     (buf(1),d(1)), (ipvt ,ivpt ), (buf(1),ihq(1))
 DATA   lbuf  / 100/, mcb/ 7*0 /, ma1h/ 4HEMA ,2H   /, kons/ 14 /,  &
     mdict / 6  /, hdr/ 6*0 /
 DATA   gpect , xemd  , xblock / 101, 102, 103 / ,  &
     xgg                    / 201           / ,  &
     scr1  , scr2           / 301, 302      /
 
 
!     RE-SET KONS IF HALF WORD IS LARGER THAN THAN 16 BITS
 
 IF (ihalf >= 18) kons = 16
 IF (ihalf >= 30) kons = 24
 
!     ALLOCATE BUFFERS. OPEN GPECT AND ALLOCATE A TABLE OF 4 WORDS
!     PER ELEMENT TYPE. OPEN SCRATCH FILE FOR GPECTX. OPEN XEMD.
 
 mcb(1) = 0
 mcb(2) = 0
 mcb(3) = 0
 mcb(4) = 0
 mcb(5) = 0
 mcb(6) = 0
 mcb(7) = 0
 maxii  = 2**kons - 1
 maxblk = 0
 ksft   = kshift
 mask   = lshift(jhalf,ihalf)
 buf1   = korsz(z) - sysbuf
 buf2   = buf1 - sysbuf
 buf3   = buf2 - sysbuf
 lcore  = buf1 + sysbuf - 1
 k4flag = 1
 IF (nok4 == -1) k4flag = 0
 
!     SET LOGICAL VARIABLE TRUE IF THIS IS A PIEZOELECRRIC COUPLED
!     PROBLEM AND STRUCTURAL DAMPING FLAG IS ON
 
 piez = .false.
 IF (ipiez == 1 .AND. k4flag /= 0) piez = .true.
 buf(1) = xblock
 CALL rdtrl (buf)
 IF (buf(1) < 0) GO TO 9132
 prec   = buf(2)
 buf(1) = gpect
 CALL rdtrl (buf)
 IF (buf(1) < 0) GO TO 9131
 idict = 4*buf(2) + 1
 nsil  = buf(3)
 maxel = buf(4)
 maxdof= buf(5)
 CALL gopen (gpect,z(buf1),rdrew)
 l     = (2**(ihalf+ihalf-2)-1)*2 + 1
 maxipv= rshift(l,kons)
 scrin = gpect
 scrout= scr1
 CALL gopen (xemd,z(buf3),rdrew)
 CALL OPEN (*9134,scrout,z(buf2),wrtrew)
 CALL WRITE (scrout,buf,3,1)
 
!     SET SWITCHES FOR MULTIPLICATION BY DAMPING
!     OR WEIGHT MASS FACTOR (OR BOTH)
 
 eps = ABS(wtmass-1.0)
 IF (eps < 1.e-6 .AND. k4flag == 0) ASSIGN  1370 TO kfact
 IF (eps < 1.e-6 .AND. k4flag /= 0) ASSIGN 13651 TO kfact
 IF (eps > 1.e-6 .AND. k4flag == 0) ASSIGN 13652 TO kfact
 IF (eps > 1.e-6 .AND. k4flag /= 0) ASSIGN 13653 TO kfact
 
!     FILL CORE WITH ELEMENT MATRIX DICTIONARIES. FOR EACH ELEMENT TYPE
!     STORE POINTER TO 1ST DICT AND THE NBR OF DICTS IN TABLE AT TOP OF
!     CORE ALSO STORE LENGTH OF EACH DICTIONARY AND FORMAT CODE.
 
 l = idict
 DO  i = 1,idict
   z(i)  = 0
 END DO
 maxn  = 0
 last  = .true.
 1014 CALL READ (*1026,*9135,xemd,hdr,3,0,nread)
 ielem = 4*hdr(1) - 3
 ldict = hdr(2)
 z(ielem  ) = l
 z(ielem+2) = ldict
 z(ielem+3) = hdr(3)
 1016 IF (l+ldict >= buf3) GO TO 1024
 jdict = ldict
 CALL READ (*9136,*1018,xemd,z(l),ldict,0,nread)
 l = l + ldict
 GO TO 1016
 1018 IF (nread /= 0) GO TO 9001
 high = z(l-ldict)
 z(ielem+1) = (l-z(ielem))/ldict
 GO TO 1014
 1024 last = .false.
 z(ielem+1) = (l-z(ielem))/ldict
 high = z(l-jdict)
 IF (z(ielem+1) /= 0) GO TO 1030
 z(ielem) = 0
 GO TO 1030
 1026 last = .true.
 CALL CLOSE (xemd,clsrew)
 
!     PASS GPECT (OR PARTIALLY COMPLETED GPECTX) ENTRY BY ENTRY.
!     IF ENTRY HAS BEEN COMPLETED, COPY IT OUT.  OTHERWISE, LOCATE
!     DICTIONARY FOR ELEMENT (IF IN CORE) AND ATTACH IT.
!     DETERMINE LENGTH OF LONGEST RECORD IN GPECTX.
!     DETERMINE THE MAXIMUM LENGTH OF ONE COLUMN OF AN ELEMENT MATRIX.
 
 1030 nhdr = 2
 IF (last) nhdr = 5
 maxgpe = 0
 low  = z(idict)
 ngps = 0
 
!     READ AND WRITE HEADER FOR RECORD (SIL, DOF, ETC.)
 
 1032 CALL READ (*1110,*9137,scrin,hdr,2,0,nread)
 CALL WRITE (scrout,hdr,nhdr,0)
 gpewds = nhdr
 ngps   = ngps + 1
 
!     READ FIRST WORD  OF ENTRY ON GPECT. TEST FOR DICT ALREADY ATTACHED
 
 1035 CALL READ (*9138,*1100,scrin,buf,1,0,nread)
 IF (IABS(buf(1)) > lbuf) GO TO 9002
 IF (buf(1) < 0) GO TO 1040
 
!     DICTIONARY ALREADY ATTACHED---READ REMAINDER OF ENTRY AND
!     COPY ENTRY TO GPECTX.
 
 CALL READ (*9139,*9140,scrin,buf(2),buf(1),0,nread)
 n = buf(1) + 1
 1038 CALL WRITE (scrout,buf,n,0)
 gpewds = gpewds + n
 GO TO 1035
 
!     DICTIONARY NOT ATTACHED---TRY TO LOCATE DICT IN CORE
 
 1040 m = -buf(1)
 CALL READ (*9141,*9142,scrin,buf(2),m,0,nread)
 IF (buf(2) <  low) GO TO 1044
 IF (buf(2) > high) GO TO 1042
 kelem = 4*buf(3) - 3
 IF (z(kelem) == 0) GO TO 1035
 l     = z(kelem  )
 n     = z(kelem+1)
 ldict = z(kelem+2)
 nlocs = z(kelem+3)
 CALL bisloc (*1035,buf(2),z(l),ldict,n,k)
 k = k + l - 1
 IF (k4flag /= 0 .AND. z(k+4) == 0) GO TO 1035
 GO TO 1050
 1042 IF (last) GO TO 1044
 n = m + 1
 GO TO 1038
 1044 CONTINUE
 GO TO 1035
 
!     DICTIONARY LOCATED---WRITE OUT COMPLETED ENTRY ON GPECTX
!         0             NO. OF WORDS IN ENTRY (NOT INCL THIS WORD)
!       1 - 5           ELEM ID, F, N, C, GE
!         6             LOC OF ELEMENT MATRIX COLUMNS FOR CURRENT PIVOT
!       7 - 6+NGRIDS    SIL-S OF CONNECTED GRID POINTS
 
 1050 ngrids = m - 2
 indx   = k + mdict - 1
 IF (nlocs == 1) GO TO 1056
 IF (ngrids > nlocs) GO TO 9004
 kk = 1
 DO  i = 1,ngrids
   IF (kk == 1) GO TO 10525
   
!     CHECK FOR DUPLICATE SILS - E.G. HBDY ELEMENT WITH AMBIENT PTS
   
   10522 IF (buf(kk+3) /= buf(kk+2)) GO TO 10525
   kk = kk + 1
   GO TO 10522
   10525 IF (buf(kk+3) /= hdr(1)) GO TO 10530
   
!     SIL THAT MATCHES THE PIVOT FOUND.  NOW INSURE THAT THIS SIL
!     HAS NOT BEEN ALREADY CONNECTED DUE TO A PREVIOUS ENTRY IN THIS
!     GPECT RECORD.  (CAUSED BY DUPLICATE IDS I.E. CELAS2)
   
!     GINO-LOC WILL NOW BE ZERO IF THAT IS TRUE
   
   indx = k + mdict + i - 2
   IF (z(indx) == 0.0) THEN
     GO TO 10530
   ELSE
     GO TO  1054
   END IF
   10530 kk = kk + 1
 END DO
 GO TO 1035
 1054 z(k+mdict-1) = z(indx)
 IF (z(k+1) /= 2) GO TO 1056
 buf(4) = buf(i+3)
 ngrids = 1
 1056 IF (ldict-nlocs+1 /= mdict) GO TO 9010
 n = mdict + ngrids
 CALL WRITE (scrout,n,1,0)
 CALL WRITE (scrout,z(k),mdict,0)
 maxblk = MAX0(z(indx)/ksft,maxblk)
 
!     ZERO GINO-LOC AS HAVING BEEN USED NOW.
 
 z(indx) = 0
 CALL WRITE (scrout,buf(4),ngrids,0)
 maxn   = MAX0(maxn,z(k+2))
 gpewds = gpewds + n + 1
 GO TO 1035
 
!     HERE ON END-OF-RECORD ON GPECT
 
 1100 CALL WRITE (scrout,0,0,1)
 maxgpe = MAX0(maxgpe,gpewds)
 GO TO 1032
 
!     HERE ON END-OF-FILE ON GPECT---TEST FOR COMPLETION OF GPECTX
 
 1110 CALL CLOSE (scrin ,clsrew)
 CALL CLOSE (scrout,clsrew)
 IF (ngps /= nsil) GO TO 9024
 IF (last) GO TO 1200
 
!     GPECTX NOT COMPLETE---SWITCH FILES AND MAKE ANOTHER PASS
 
 IF (scrin == gpect) scrin = scr2
 k = scrin
 scrin  = scrout
 scrout = k
 CALL gopen (scrin ,z(buf1),rdrew )
 CALL gopen (scrout,z(buf2),wrtrew)
 l = idict
 ldict = z(ielem+2)
 nlocs = z(ielem+3)
 DO  i = 1,idict
   z(i) = 0
 END DO
 last = .true.
 z(ielem  ) = idict
 z(ielem+2) = ldict
 z(ielem+3) = nlocs
 GO TO 1016
 
!     HERE WE GO NOW FOR THE ASSEMBLY PHASE---PREPARE BY ALLOCATING
!     STORAGE FOR ONE ELEMENT MATRIX COLUMN AND ITS ROW POSITIONS
 
 1200 irowp  = prec*maxn + 1
 igpx   = irowp + maxn
 first  = .true.
 gpectx = scrout
 mcb(1) = xgg
 mcb(4) = 6
 mcb(5) = prec
 mcb(6) = 0
 mcb(7) = 0
 last   = .false.
 nrec   = 0
 maxnpr = maxn*prec
 openr  = rdrew
 openw  = wrtrew
 oldcod = 0
 itab   = buf1 - maxblk
 ntab   = buf1 - 1
 IF (itab < igpx) GO TO 9011
 
!     BEGIN A PASS - OPEN GPECTX
 
 1210 ipvt = igpx
 jj   = itab - 3
 DO  indx = itab,ntab
   z(indx) = 0
 END DO
 CALL gopen (gpectx,z(buf1),openr)
 
!     READ A RECORD FROM GPECTX INTO CORE
 
 1220 IF (ivpt+maxgpe >= jj .OR. ipvt > maxipv) GO TO 1304
 CALL READ (*9144,*1222,gpectx,z(ivpt),maxgpe+1,1,nread)
 GO TO 9006
 1222 icol = ipvt + nread
 nrec = nrec + 1
 
!     MAKE A PASS THROUGH EACH ELEMENT CONNECTED TO THE PIVOT---FORM THE
!     UNION OF ALL CODE WORDS AND STORE ELEMENT POINTERS IN LIST AT THE
!     END OF OPEN CORE
 
 ppoint = lshift(ipvt,kons)
 ii     = ipvt + 5
 union  = 0
 z(ipvt+2) = 0
 z(ipvt+3) = 0
 z(ipvt+4) = ipvt
 GO TO 1234
 1231 IF (z(ii) < 0) GO TO 9007
 kk = ii - ipvt
 IF (kk > maxii) GO TO 9008
 IF (jj <= icol ) GO TO 1300
 z(jj  ) = z(ii+mdict)
 z(jj+1) = orf(ppoint,kk)
 z(jj+2) = 0
 indx    = itab + z(jj)/ksft - 1
 IF (indx > ntab) GO TO 9016
 IF (z(indx) /= 0) GO TO 1236
 z(indx) = orf(lshift(itab-jj,ihalf),itab-jj)
 GO TO 1237
 1236 jjlast = itab - andf(z(indx),jhalf)
 z(jjlast+2) = jj
 z(indx) = orf(andf(z(indx),mask),itab-jj)
 1237 jj = jj - 3
 union = orf(union,z(ii+4))
 ii = ii + z(ii) + 1
 1234 IF (ii < icol) GO TO 1231
 IF (ii /= icol) GO TO 9022
 
!     FORM THE LIST OF NON-NULL COLUMNS TO BE BUILT FOR THIS PIVOT
 
 IF (union == 0) GO TO 1280
 IF (union == oldcod) GO TO 1243
 CALL decode (union,scalas,nsca)
 oldcod = union
 1243 z(ipvt+2) = icol
 IF (icol+nsca >= jj) GO TO 1300
 ii = icol
 DO  l = 1,nsca
   z(ii) = z(ipvt) + scalas(l)
   ii = ii + 1
 END DO
 irow = ii
 
!     NOW MAKE A PASS AGAIN THROUGH EACH ELEMENT CONNECTED TO CURRENT
!     PIVOT AND FORM A LIST OF UNIQUE ROW INDICES.
 
 ii = ipvt + 5
 1252 l1 = ii + mdict + 1
 l2 = ii + z(ii)
 IF (oldcod == z(ii+4)) GO TO 1253
 icode = z(ii+4)
 CALL decode (icode,scalas,nsca)
 oldcod = z(ii+4)
 1253 CONTINUE
 kk = irow
 IF (ii /= ipvt+5) GO TO 1258
 DO  l = l1,l2
   
!     IGNORE DUPLICATE IDS AS IN SOME CELAS2 ELEMENTS ETC.
   
   IF (l > l1 .AND. z(l) == z(l-1)) CYCLE
   DO  i = 1,nsca
     z(kk) = z(l) + scalas(i)
     kk = kk + 1
     IF (kk >= jj) GO TO 1300
   END DO
 END DO
 nrow   = kk - 1
 nbrwds = kk - irow
 GO TO 1269
 1258 j = irowp
 DO  l = l1,l2
   
!     IGNORE DUPLICATE IDS AS IN SOME CELAS2 ELEMENTS ETC.
   
   IF (l > l1 .AND. z(l) == z(l-1)) CYCLE
   DO  i = 1,nsca
     z(j) = z(l) + scalas(i)
     j = j + 1
   END DO
 END DO
 IF (j > igpx) GO TO 9023
 m = j - irowp
 IF (irow+nbrwds+m >= jj) GO TO 1300
 CALL mrge (z(irow),nbrwds,z(irowp),m)
 nrow = irow + nbrwds - 1
 IF (nrow >= jj) GO TO 1300
 1269 ii = l2 + 1
 IF (ii < icol) GO TO 1252
 z(ipvt+3) = irow
 
!     NOW ALLOCATE STORAGE FOR COLUMNS OF XGG ASSOCIATED WITH THIS PIVOT
 
 imat   = nrow + 1
 nbrcol = irow - icol
 nbrrow = imat - irow
 nbrwds = prec*nbrcol*nbrrow
 nmat   = imat + nbrwds - 1
 IF (nmat >= jj) GO TO 1300
 DO  i = imat,nmat
   z(i) = 0
 END DO
 z(ipvt+4) = imat
 ii = nmat + 1
 
!     ADVANCE POINTER AND TRY TO GET ANOTHER PIVOT ALLOCATED
 
 1280 ilist = jj + 3
 npvt  = ipvt
 IF (nrec == ngps) GO TO 1310
 ipvt = ii
 GO TO 1220
 
!     HERE WHEN STORAGE EXCEEDED DURING PROCESSING OF A PIVOT.
!     IF FIRST PIVOT ON PASS, INSUFFICIENT CORE FOR MODULE.
!     OTHERWISE, BACKSPACE GPECTX AND PREPARE TO PROCESS ALL
!     PIVOTS IN CORE WHICH HAVE BEEN COMPLETELY ALLOCATED.
 
 1300 IF (ipvt == igpx) CALL mesage (-8,0,ma1h)
 CALL bckrec (gpectx)
 nrec = nrec - 1
 1304 op   = cls
 IF (ipvt == igpx) CALL mesage (-8,0,ma1h)
 GO TO 1320
 
!     HERE WHEN LAST PIVOT POINT HAS BEEN READ AND ALLOCATED
 
 1310 last = .true.
 op   = clsrew
 
!     CLOSE GPECTX. OPEN XBLOCK.
 
 1320 CALL CLOSE (gpectx,op)
 nwds = buf1 - ilist
 IF (nwds <= 0) GO TO 1402
 CALL gopen (xblock,z(buf1),rdrew)
 oldcod = 0
 
!     PASS THE LIST OF ELEMENT MATRIX POINTERS. EACH ENTRY POINTS TO THE
!     PIVOT POINT AND ELEMENT DICTIONARY IN CORE AND TO THE POSITION IN
!     THE XBLOCK FILE CONTAINING THE ASSOCIATED ELEMENT MATRIX COLUMNS.
!     WHEN PROCESSING OF ALL ENTRIES IS COMPLETE, COLUMNS OF XGG NOW IN
!     CORE ARE COMPLETE.
 
 DO  indx = itab,ntab
   IF (z(indx) == 0) CYCLE
   jj = itab - rshift(z(indx),ihalf)
   IF (jj < ilist) CYCLE
   1330 CONTINUE
   CALL filpos (xblock,z(jj))
   ipvt  = rshift(z(jj+1),kons)
   IF (ipvt > npvt) GO TO 9019
   ielem = ipvt + andf(z(jj+1),maxii)
   icol  = z(ipvt+2)
   irow  = z(ipvt+3)
   imat  = z(ipvt+4)
   
!     DECODE CODE WORD FOR ELEMENT. FORM LIST OF ROW INDICES DESCRIBING
!     TERMS IN THE ELEMENT MATRIX COLUMN. THEN CONVERT THESE INDICES TO
!     RELATIVE ADDRESSES IN XGG COLUMN IN CORE (USE LIST OF ROW INDICES
!     FOR XGG COLUMN TO DO THIS).
   
   IF (z(ielem+4) == oldcod) GO TO 1341
   icode = z(ielem+4)
   CALL decode (icode,scalas,nsca)
   oldcod = z(ielem+4)
   1341 l1 = ielem + mdict + 1
   l2 = ielem + z(ielem)
   k  = irowp
   DO  l = l1,l2
     
!     IGNORE DUPLICATE IDS AS IN SOME CELAS2 ELEMENTS ETC.
     
     IF (l > l1 .AND. z(l) == z(l-1)) CYCLE
     DO  i = 1,nsca
       z(k) = z(l) + scalas(i)
       k = k + 1
     END DO
   END DO
   nrowp = k - 1
   IF (nrowp >= igpx) GO TO 9012
   nrow = imat - 1
   id   = z(irowp)
   CALL bisloc (*9020,id,z(irow),1,(imat-irow),irowx)
   irowx = irow + irowx - 1
   DO  k = irowp,nrowp
     DO  i = irowx,nrow
       IF (z(k) == z(i)) GO TO 1347
     END DO
     GO TO 9013
     1347 z(k)  = (i-irow)*prec
     irowx = i + 1
   END DO
   nbrrow = nrowp - irowp + 1
   
!     PREPARE TO READ EACH COLUMN OF ELEMENT MATRIX
   
   ncol  = irow - 1
   icolx = icol
   nbrwds= z(ielem+3)*prec
   IF (z(ielem+2) ==  2) nbrwds = prec
   IF (nbrwds > maxnpr) GO TO 9014
   DO  i = 1,nsca
     
!     READ A COLUMN OF THE ELEMENT MATRIX AND DETERMINE ADDRESS
!     OF FIRST WORD OF ASSOCIATED COLUMN OF XGG IN CORE.
     
     CALL READ (*9145,*9146,xblock,z,nbrwds,0,nread)
     col = z(ipvt) + scalas(i)
     DO  k = icolx,ncol
       IF (col == z(k)) GO TO 1364
     END DO
     GO TO 9015
     1364 imatn = imat + (imat-irow)*(k-icol)*prec
     icolx = k + 1
     IF (z(ielem+2) /= 2) GO TO 1365
     
!     ELEMENT MATRIX IS DIAGONAL
     
     nbrrow = 1
     z(irowp) = z(irowp+i-1)
     
!     IF DAMPING OR WEIGHT MASS FACTOR (OR BOTH) PRESENT, MULTIPLY
!     EACH TERM IN THE ELEMENT MATRIX COLUMN BY THE FACTOR.
     
     1365 GO TO kfact, (1370,13651,13652,13653)
     13651 factor = y(ielem+5)
     GO TO 13654
     13652 factor = wtmass
     GO TO 13654
     13653 factor = y(ielem+5)*wtmass
     13654 CONTINUE
     IF (prec == 2) GO TO 1367
     DO  k = 1,nbrwds
       
!     FOR PIEZOELECTRIC COUPLED PROBLEMS, ANY STRUCTURAL DAMPING COEFF.
!     SHOULD MULTIPLY ONLY THE UNCOUPLED STRUCTURAL TERMS. SO, SKIP
!     EVERY 4TH TERM IN A COLUMN AND SKIP EVERY 4TH COLUMN
       
       IF (piez .AND. (i == nsca .OR. MOD(k,4) == 0)) y(k) = 0.
       y(k) = factor*y(k)
     END DO
     GO TO 1371
     1367 m = nbrwds/2
     DO  k = 1,m
       IF (piez .AND. (i == nsca.OR.MOD(k,4) == 0)) zd(k) = 0.d0
       zd(k) = factor*zd(k)
     END DO
     GO TO 1374
     
!     NOW ADD TERMS OF THE ELEMENT MATRIX INTO XGG
     
     1370 IF (prec == 2) GO TO 1374
     
!     DO ARITHMETIC IN SINGLE PRECISION
     
     1371 DO  k = 1,nbrrow
       j = imatn + z(irowp+k-1)
!WKBI 1/95
       yj   = y(j)
       y(j) = y(j) + y(k)
!WKBNB 1/95  FOLLOWING CODE WILL CAUSE A TRUE ZERO WHEN SUBTRACTING SAME NO.
       IF ( y(k) == 0.0 ) GO TO 13720
       yj   = yj / y(k)
       IF ( yj <= -.999999999998 .AND. yj >= -1.000000000001 ) y(j) = 0.0
       13720 CONTINUE
!WKBNE 1/95
     END DO
     CYCLE
     
!     DO ARITHMETIC IN DOUBLE PRECISION
     
     1374 DO  k = 1,nbrrow
       j = imatn + z(irowp+k-1)
       is(1) = z(j  )
       is(2) = z(j+1)
!WKBI 1/95
       xdd   = xd(1)
       xd(1) = xd(1) + zd(k)
!WKBNB 1/95 FOLLOWING CODE WILL CAUSE A TRUE ZERO WHEN SUBTRACTING SAME NO.
!  WITHOUT THIS CODE, A SYMMETRIC MATRIX MIGHT HAVE UNSYMMETRIC TERMS ON
!  THE HP AND ULTRIX (SEE DEMO D01011A, MATRIX KAA, COLUMN 84 ROW 70)
       IF ( zd(k) == 0.0D0 ) GO TO 1375
       xdd   = xdd / zd(k)
       IF ( xdd <= -.999999999998 .AND. xdd >= -1.000000000001 )  &
           xd(1) = 0.0D0
       1375  CONTINUE
!WKBNE 1/95
       z(j  )= is(1)
       z(j+1)= is(2)
       
     END DO
     CYCLE
     
!     END OF DO LOOPS
     
   END DO
   jj = z(jj+2)
   IF (jj >= ilist) GO TO 1330
 END DO
 
!     ALL COLUMNS OF XGG IN CORE ARE NOW COMPLETE - SEND THEM
!     OUT TO THE XGG DATA BLOCK VIA THE BLDPK ROUTINE.
 
 CALL CLOSE (xblock,clsrew)
 1402 CALL gopen (xgg,z(buf1),openw)
 ipvt = igpx
 
!     PREPARE TO PACK ALL COLUMNS FOR CURRENT PIVOT
 
 1410 col1 = z(ipvt)
 coln = col1 + z(ipvt+1) - 1
 icol = z(ipvt+2)
 irow = z(ipvt+3)
 imat = z(ipvt+4)
 ncol = irow - 1
 nrow = imat - 1
 jnext  = icol
 nxtcol = z(jnext)
 ii = imat
 DO  col = col1,coln
   
!     INITIATE PACKING BY CALLING BLDPK. TEST FOR NULL COL.
   
   CALL bldpk (prec,prec,xgg,0,0)
   IF (icol == 0) GO TO 1428
   IF (col < nxtcol) GO TO 1428
   jnext  = jnext + 1
   nxtcol = z(jnext)
   IF (jnext > ncol) nxtcol = coln + 1
   IF (prec == 2) GO TO 1426
   
!     NON-NULL COLUMN - SEND THE TERMS OUT VIA ZBLPKI
   
!     SINGLE PRECISION
   
   DO  k = irow,nrow
     iq   = z(k)
     q(1) = z(ii)
     CALL zblpki
     ii = ii + 1
   END DO
   GO TO 1428
   
!     DOUBLE PRECISION
   
   1426 DO  k = irow,nrow
     iq   = z(k)
     q(1) = z(ii  )
     q(2) = z(ii+1)
     CALL zblpki
     ii = ii + 2
   END DO
   
!     TERMINATE COLUMN BY CALLING BLDPKN
   
   1428 CALL bldpkn (xgg,0,mcb)
 END DO
 
!     LOGIC TEST TO MAKE SURE POINTERS ENDED CORRECTLY
 
 nbrwds = 5
 IF (icol == 0) GO TO 1440
 nbrwds = (imat-irow)*(irow-icol)*prec
 IF (ii-imat /= nbrwds .AND. z(ipvt+1) /= 1) GO TO 9017
 
!     TEST FOR LAST PIVOT
 
 1440 IF (ipvt >= npvt) GO TO 1450
 ipvt = imat + nbrwds
 GO TO 1410
 
!     CLOSE XGG
 
 1450 CALL CLOSE (xgg,op)
 
!     TEST FOR LAST PASS
 
 IF (last) GO TO 1490
 first = .false.
 openr = rd
 openw = wrt
 GO TO 1210
 
!     XGG NOW COMPLETE - WRITE ITS TRAILER.
 
 1490 mcb(3) = mcb(2)
 IF (mcb(2) /= z(npvt)+z(npvt+1)-1) GO TO 9018
 CALL wrttrl (mcb)
 RETURN
 
!     FATAL ERROR MESSAGES
 
 9001 msg(1) = 1016
 GO TO 9098
 9002 msg(1) = 1035
 GO TO 9098
 9004 msg(1) = 1052
 GO TO 9098
 9006 msg(1) = 1220
 GO TO 9098
 9007 msg(1) = 1231
 GO TO 9098
 9008 msg(1) = 1232
 GO TO 9098
 9010 msg(1) = 1056
 GO TO 9098
 9011 msg(1) = 1202
 GO TO 9098
 9012 msg(1) = 1344
 GO TO 9098
 9013 msg(1) = 1346
 GO TO 9098
 9014 msg(1) = 1352
 GO TO 9098
 9015 msg(1) = 1362
 GO TO 9098
 9016 msg(1) = 1238
 GO TO 9098
 9017 msg(1) = 1432
 GO TO 9098
 9018 msg(1) = 1490
 GO TO 9098
 9019 msg(1) = 1332
 GO TO 9098
 9020 msg(1) = 1345
 GO TO 9098
 9022 msg(1) = 1235
 GO TO 9098
 9023 msg(1) = 1264
 GO TO 9098
 9024 msg(1) = 1110
 GO TO 9098
 9098 WRITE  (output,9099) sfm,msg(1)
 9099 FORMAT (a25,' 3102, LOGIC ERROR EMA - ',i4)
 9097 CONTINUE
 WRITE  (output,9091)
 9091 FORMAT (/,' *** CONTENTS OF /MA1XX/')
 WRITE  (output,9092) ihq
 9092 FORMAT (5X,10I10)
 WRITE  (output,9093)
 9093 FORMAT (/,' FIRST 250 WORDS OF OPEN CORE')
 j = 250
 WRITE (output,9092) (z(i),i=1,j)
 CALL mesage (-61,0,0)
 9100 CALL fname (msg(3),msg(1))
 WRITE  (output,9101) sfm,(msg(i),i=1,4)
 9101 FORMAT (a25,' 3001, ATTEMPT TO OPEN DATA SET ',2A4,', FILE (',i4,  &
     ') IN SUBROUTINE EMA (',i4,') WHICH WAS NOT DEFINED IN FIST.')
 GO TO 9097
 9110 CALL fname (msg(3),msg(1))
 WRITE  (output,9111) sfm,(msg(i),i=1,4)
 9111 FORMAT (a25,' 3002, EOF ENCOUNTERED WHILE READING DATA SET ',2A4,  &
     ', (FILE',i5,') IN SUBROUTINE EMA (',i4,1H))
 GO TO 9097
 9120 CALL fname (msg(3),msg(1))
 WRITE  (output,9121) sfm,(msg(i),i=1,4)
9121 FORMAT (a25,' 3003, ATTEMPT TO READ PAST END OF LOGICAL RECORD IN'  &
    ,     ' DATA SET ',2A4,' (FILE',i5,') IN SUBROUTINE EMA (',i4,1H))
GO TO 9097
9131 msg(3) = gpect
msg(4) = 1002
GO TO 9100
9132 msg(3) = xblock
msg(4) = 1001
GO TO 9100
9134 msg(3) = scrout
msg(4) = 1005
GO TO 9100
9135 msg(3) = xemd
msg(4) = 1014
GO TO 9120
9136 msg(3) = xemd
msg(4) = 1017
GO TO 9110
9137 msg(3) = scrin
msg(4) = 1032
GO TO 9120
9138 msg(3) = scrin
msg(4) = 1035
GO TO 9110
9139 msg(3) = scrin
msg(4) = 1036
GO TO 9110
9140 msg(3) = scrin
msg(4) = 1036
GO TO 9120
9141 msg(3) = scrin
msg(4) = 1040
GO TO 9110
9142 msg(3) = scrin
msg(4) = 1040
GO TO 9120
9144 msg(3) = gpectx
msg(4) = 1221
GO TO 9110
9145 msg(3) = xblock
msg(4) = 1360
GO TO 9110
9146 msg(3) = xblock
msg(4) = 1360
GO TO 9120
END SUBROUTINE ema
