SUBROUTINE seemat
     
!     SUBROUTINE SEEMAT IS THE DMAP DRIVER FOR UTILITY MODULE SEEMAT
!     WHOSE DMAP CALL FOLLOWS
 
!     SEEMAT    A,B,C,D,E//C,N,PRINT(PLOT)/V,N,PFILE/C,N,FSIZE/
!                          C,N,MODIDA/C,N,MODELA/C,N,PAPERX/C,N,PAPERY
 
!     INPUT DATA BLOCKS  - A,B,C,D,E ARE MATRICES, ANY OF WHICH MAY BE
!                          PURGED.
 
!     OUTPUT DATA BLOCKS - NONE
 
!     PARAMETERS
!       1. BCD, -PRINT- MEANS USE SYSTEM PRINTER (DEFAULT).
!               -PLOT- MEANS USE SPECIFIED PLOTTER.
!       2. INTEGER, PLOT COUNTER (INPUT + OUTPUT).
!       3. INTEGER, FRAME SIZE = NUMBER OF CHARACTERS TO BE TYPED
!                   IN AN ASSUMED SQUARE FRAME (DEFAULT=100).
!       4. BCD, MODEL ID (DEFAULT=M).
!       5. INTEGER, MODEL NUMBER (DEFAULT=1).
!       6. REAL, X DIMENSION OF PLOT FRAME (DEFAULT=0.0).
!       7. REAL, Y DIMENSION OF PLOT FRAME (DEFAULT=0.0).
!      NOTE - PARAMETERS 2-7 ARE USED ONLY IF PARAMETER 1 = -PLOT-.
 
 EXTERNAL        andf,orf
 LOGICAL :: table,sq,plotit,prntit,tapbit,nobits
 INTEGER :: NAME(5),gobac,BLANK,xstar,xdolr,xdddd,seemt(2),  &
     a,b,c,it(7),sysbuf,eol,eor,iro(10),ix(1),lbl(2),  &
     andf,orf,kpp(2),plus,bcor,symbl(2),modid(2),  &
     ttl1(9),ttl2(4),ttl3(4),ttl4(3),lin(25),  &
     pp,pfile,fsize,pltter,pltype,ploter,pltbuf,two
 CHARACTER (LEN=31) :: sim
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm,sim
 COMMON /system/ sysbuf,nout,jazz1(6),nlines,jazz2(2),lnct,  &
     jazz3(26),nbpc,nbpw,ncpw
 COMMON /BLANK / pp(2),pfile,fsize,modida(2),modela,paperx,papery
 COMMON /zzzzzz/ x(1)
 COMMON /zntpkx/ z(4),iz,eol,eor
 COMMON /pltdat/ model,pltter,region(4),axmax,aymax,edge(12),  &
     skpa(9),pltype,ploter
 COMMON /xxparm/ pltbuf,kamran,nblfm,skparm(4),papsiz(2)
 COMMON /two   / two(32)
 EQUIVALENCE     (x(1),ix(1)),(iro(1),icol1),(iro(2),iblcu1),  &
     (iro(3),iblcu2),(iro(4),jblcu1),(iro(5),jblcu2),  &
     (iro(6),a,ipij1),(iro(7),b,ipij2),(iro(8),c),  &
     (it(1),nam ),(it(2),ncols),(it(3),nrows), (it(5),ityp)
 DATA    NAME  , seemt /101,102,103,104,105,4HSEEM,4HAT   /
 DATA    BLANK , ncc,xstar,xdolr,xdddd/1H ,100,1H*,1H$,1HD/
 DATA    kpp   / 4HPLOT,4H    /
 DATA    plus  / 4H+   /
 DATA    ttl1  / 4HSEEM,4HAT d,4HISPL,4HAY o,4HF ma,4HTRIX,  &
     4H dat,4HA bl,4HOCK /
 DATA    ttl2  / 4HNO. ,4HCOLU,4HMNS ,4H=   /
 DATA    ttl3  / 4HNO. ,4H  ro,4HWS  ,4H=   /
 DATA    ttl4  / 4H(tra,4HNSPO,4HSED)/
 
 ncc    = 100
 plotit =.false.
 prntit =.true.
 nlnxx  = nlines
 IF (pp(1) /= kpp(1) .OR. pp(2) /= kpp(2)) GO TO 20
 plotit =.true.
 prntit =.false.
 table  =.false.
 ncc    = fsize
 fncc   = ncc
 nlnxx  = ncc
 20 lcor   = korsz(x) - sysbuf
 
 ncc1   = ncc/4
 ncc5   = ncc - 5
 lblk   = (ncc*nlnxx-1)/32 + 3
 IF (prntit) GO TO 90
 
!     INITIALIZE PLOTTER
 
 modid(1) = modida(1)
 modid(2) = modela
 CALL fndplt (pltter,model,modid)
 papsiz(1) = paperx
 papsiz(2) = papery
 kamran = 3
 nblkfm = 0
 CALL pltset
 lcor = lcor - pltbuf
 kcor = lcor + sysbuf + 1
 IF (lcor <= 0) CALL mesage (-8,sq,seemt)
 bcor = lcor - ncc1
 IF (tapbit(ploter)) GO TO 70
 WRITE  (nout,65) uwm,ploter
 65 FORMAT (a25,' 1704, PLOT FILE -',a4,'- NOT SET UP')
 GO TO 9999
 70 IF (IABS(pltype) /= 1) table = .true.
 region(3) = AMIN1(axmax,aymax)
 region(4) = region(3)
 axmax = region(3)
 aymax = region(4)
 CALL mapset (0,0,1.01*fncc,1.01*fncc,0,0,axmax,aymax,2)
 CALL map (0.005*fncc,0.005*fncc,bllx,blly)
 CALL map (1.005*fncc,0.005*fncc,blrx,blry)
 CALL map (1.005*fncc,1.005*fncc,burx,bury)
 CALL map (0.005*fncc,1.005*fncc,bulx,buly)
 GO TO 90
 85 CALL mesage (-1,ploter,seemt)
 90 CONTINUE
 
 DO  iii = 1,5
   
   nam = NAME(iii)
   CALL rdtrl (it)
   IF (nam <= 0) CYCLE
   CALL gopen (nam,x(lcor+1),0)
   CALL fname (nam,lbl)
   sq = .true.
   IF (ncols /= nrows) sq = .false.
   nblks = 0
   ncol1 = 0
   ijmax = MAX0(ncols,nrows)
   nrows1= nrows + 1
   IF (prntit) GO TO 95
   IF (table ) GO TO 92
   pfile = pfile + 1
   CALL sopen (*85,ploter,x(kcor),pltbuf)
   CALL stplot (pfile)
   CALL map   (0.23*fncc,0.50*fncc,xxxx,yyyy)
   CALL PRINT (xxxx,yyyy,1,ttl1,9,-1)
   CALL PRINT (xxxx,yyyy,1,ttl1,9, 0)
   CALL map   (0.60*fncc,0.50*fncc,xxxx,yyyy)
   CALL PRINT (xxxx,yyyy,1,lbl,2,0)
   CALL map   (0.75*fncc,0.50*fncc,xxxx,yyyy)
   CALL PRINT (xxxx,yyyy,1,ttl4,3,0)
   CALL map   (0.40*fncc,0.40*fncc,xxxx,yyyy)
   CALL PRINT (xxxx,yyyy,1,ttl3,4,0)
   CALL map   (0.40*fncc,0.30*fncc,xxxx,yyyy)
   CALL PRINT (xxxx,yyyy,1,ttl2,4,0)
   CALL map   (0.55*fncc,0.40*fncc,xxxx,yyyy)
   CALL typint (xxxx,yyyy,1,nrows,0, 0)
   CALL map   (0.55*fncc,0.30*fncc,xxxx,yyyy)
   CALL typint (xxxx,yyyy,1,ncols,0, 0)
   CALL line (bllx,blly,bulx,buly,1,-1)
   CALL line (bllx,blly,bulx,buly,1, 0)
   CALL line (bulx,buly,burx,bury,1, 0)
   CALL line (burx,bury,blrx,blry,1, 0)
   CALL line (blrx,blry,bllx,blly,1, 0)
   CALL stplot (-1)
   92 CALL page1
   lnct = lnct + 5
   WRITE  (nout,93) lbl(1),lbl(2),ncols,nrows
   93 FORMAT (//5X,'SEEMAT PLOT FOR TRANSPOSE OF', /22X,'MATRIX DATA ',  &
       'BLOCK ',2A4,11X,'PLOT FILE ','    R','     C', /10X,  &
       'SIZE =',i6,' ROWS BY',i6,' COLUMNS')
   IF (table) GO TO 95
   WRITE  (nout,94) pfile
   94 FORMAT (1H0,62X,i5,2X,12HHEADER frame)
   95 CONTINUE
   
   
!     LOOP ON COLUMNS OF MATRIX
   
   ncol = 1
   100 CONTINUE
   CALL intpk (*2100,nam,0,ityp,0)
   
!     IF COLUMN IS NULL, RETURN FROM INTPK IS TO STATEMENT 2100
!     ITY IS TYPE OF ELEMENT STORED IN Z, NOT USED IN THIS PROGRAM
!     BLOCK IS DUMMY ENTRY NOT USED BY INTPK
   
!     LOOP ON ROWS OF MATRIX
   
   nrow = 1
   200 CONTINUE
   IF (eol /= 0) GO TO 2100
   
!     READ ELEMENT OF MATRIX INTO /ZNTPKX/
   
   CALL zntpki
   
!     COMPUTE BLOCK ID IN WHICH ELEMENT BELONGS
   
!     LOOK AT CURRENT BLOCK FIRST
   
   IF (nblks <= 0) GO TO 1045
   IF (ncol <= jblcu1 .OR. ncol > jblcu2 .OR. iz <= iblcu1 .OR.  &
       iz > iblcu2) GO TO 1020
   nblk = nblcur
   GO TO 1050
   
!     SEARCH ALL BLOCKS TO FIND OLD ONE IN WHICH ELEMENT LIES
   
   1020 DO  i2 = 1,nblks
     ip  = lblk*(i2-1) + 1
     ip1 = ip + 2
     iblcu1 = ix(ip)
     iblcu2 = iblcu1 + ncc
     jblcu1 = ix(ip+1)
     jblcu2 = jblcu1 + nlnxx
     IF (ncol <= jblcu1 .OR. ncol > jblcu2 .OR. iz <= iblcu1 .OR.  &
         iz > iblcu2) CYCLE
     nblk = i2
     GO TO 1050
   END DO
   1045 nblk = -1
   1050 IF (nblk > 0) GO TO 1100
   
!     SET UP NEW BLOCK IF THERE IS ROOM FOR IT IN CORE
   
   nblks1 = nblks + 1
   IF (lblk*nblks1 <= lcor) GO TO 1070
   WRITE  (nout,1060) swm,nblks1
   1060 FORMAT (a27,' 1701, AVAILABLE CORE EXCEEDED BY',i10,' LINE IMAGE',  &
       ' BLOCKS.')
   nblks  = -1
   GO TO 9960
   
!     SET BLOCK POINTERS AND BLANK OUT LINE IMAGE
   
   1070 ip  = lblk*nblks + 1
   ip1 = ip + 2
   ip2 = ip + lblk - 1
   DO  i = ip1,ip2
     ix(i) = 0
   END DO
   DO  ijm = 1,ijmax
     IF (ijm*ncc < iz) CYCLE
     ix(ip) = ncc*(ijm-1)
     GO TO 1075
   END DO
   kerror = 1074
   GO TO 9950
   1075 DO  ijm = 1,ijmax
     IF (ijm*nlnxx < ncol) CYCLE
     ix(ip+1) = nlnxx*(ijm-1)
     GO TO 1080
   END DO
   kerror = 1079
   GO TO 9950
   1080 iblcu1 = ix(ip)
   iblcu2 = iblcu1 + ncc
   jblcu1 = ix(ip+1)
   jblcu2 = jblcu1 + nlnxx
   nblks  = nblks1
   nblcur = nblks
   IF (nblks <= 0) GO TO 9997
   
!     INSERT BIT INTO PACKED LINE IMAGE BLOCK
   
   1100 a = ncc*(ncol-ix(ip+1)-1) + (iz-ix(ip))
   b = (a-1)/32
   c = ip1 + b
   b = a - 32*b
   ix(c) = orf(ix(c),two(b))
   
!     END OF LOOP ON ROWS
   
   nrow = nrow + 1
   IF (nrow <= nrows1) GO TO 200
   kerror = 2000
   GO TO 9950
   2100 IF (ncol-ncol1 < nlnxx) GO TO 3000
   
!     OUTPUT GROUP OF LINE IMAGE BLOCKS
   
   ASSIGN 2200 TO gobac
   GO TO 9500
   2200 nblks = 0
   ncol1 = ncol1 + nlnxx
   
!     END OF LOOP ON COLUMNS
   
   3000 ncol = ncol + 1
   IF (ncol <= ncols) GO TO 100
   
!     OUTPUT RESIDUAL LINE IMAGE BLOCKS
   
   ASSIGN 3050 TO gobac
   GO TO 9500
   3050 nblks = 0
   GO TO 9997
   
!     OUTPUT GROUP OF LINE IMAGE BLOCKS
   
   9500 CONTINUE
   IF (nblks <= 0) GO TO 9699
   DO  i = 1,nblks
     ip = lblk*(i-1) + 1
     IF (prntit) CALL page1
     i1 = ix(ip)
     j100 = i1 + ncc
     DO  ij = 1,10
       iro(ij) = i1 + 10*ij
     END DO
     IF (prntit) WRITE (nout,9520) (iro(ij),ij=1,10)
     9520 FORMAT (13H0TRANSPOSE of,9X,8HCOLUMN..,10I10)
     IF (prntit) WRITE (nout,9530) lbl(1),lbl(2)
     9530 FORMAT (8H matrix ,2A4,7X,3HROW,4X,10(9X,1H.),  &
         /23X,3H...,4X,100(1H.)/24X,1H.)
     icol1 = ix(ip+1)
     i100  = icol1 + nlnxx
     ip1   = ip - ncc1 + 1
     IF (prntit) GO TO 9535
     pfile = pfile + 1
     CALL sopen (*85,ploter,x(kcor),pltbuf)
     CALL stplot (pfile)
     CALL tipe (xxxx,yyyy,1,plus,1,-1)
     ipak  = (ncc+99)/100
     ija   = 5*ipak
     ijb   = ncc - ija
     fnccy = 1.005*fncc
     DO  ij = ija,ijb,ija
       fij = FLOAT(ij)
       CALL map (fij,fnccy,xxxx,yyyy)
       CALL tipe (xxxx,yyyy,1,plus,1,0)
     END DO
     fnccx = 1.005*fncc
     DO  ij = ija,ijb,ija
       fij = fncc - FLOAT(ij)
       CALL map (fnccx,fij,xxxx,yyyy)
       CALL tipe (xxxx,yyyy,1,plus,1,0)
     END DO
     fnccy = 0.005*fncc
     DO  ij = ija,ijb,ija
       fij = fncc - FLOAT(ij)
       CALL map (fij,fnccy,xxxx,yyyy)
       CALL tipe (xxxx,yyyy,1,plus,1,0)
     END DO
     fnccx = 0.005*fncc
     DO  ij = ija,ijb,ija
       fij = FLOAT(ij)
       CALL map (fnccx,fij,xxxx,yyyy)
       CALL tipe (xxxx,yyyy,1,plus,1,0)
     END DO
     9535 DO  ij = 1,nlnxx
       ip1   = ip1 + ncc1
       ipij1 = ip1 + 1
       ipij2 = ip1 + ncc1
       ib    = ncc*(ij-1)
       iw    = ib/32
       ib    = ib - 32*iw
       iw    = iw + ip + 2
       nobits= .true.
       IF (plotit) GO TO 9570
       DO  jj = 1,ncc1
         lin(jj) = BLANK
       END DO
       DO  jj = 1,ncc
         ib = ib + 1
         IF (ib <= 32) GO TO 9537
         ib = 1
         iw = iw + 1
         9537 IF (andf(ix(iw),two(ib)) == 0) CYCLE
         nobits = .false.
         b   = (jj-1)/4 + 1
         c   = jj - 4*(b-1)
         ixx = xstar
         IF (ix(ip+1)+ij == ncols .OR. ix(ip)+jj == nrows) ixx = xdolr
         IF (sq .AND. ix(ip+1)+ij == ix(ip)+jj) ixx = xdddd
         lin(b) = khrfn1(lin(b),c,ixx,1)
       END DO
       IF (nobits) GO TO 9560
       IF (MOD(ij,5) == 0) GO TO 9550
       WRITE  (nout,9545) (lin(jj),jj=1,ncc1)
       9545 FORMAT (28X,2H. ,25A4)
       CYCLE
       9550 icol1 = icol1 + 5
       WRITE  (nout,9555) icol1,(lin(jj),jj=1,ncc1)
       9555 FORMAT (16X,i10,4H .. ,25A4)
       CYCLE
       9560 IF (MOD(ij,5) == 0) GO TO 9565
       WRITE (nout,9545)
       CYCLE
       9565 icol1 = icol1 + 5
       WRITE (nout,9555) icol1
       CYCLE
       9570 fij = 101.0 - FLOAT(ij)
       DO  jj = 1,ncc
         ib = ib + 1
         IF (ib <= 32) GO TO 9577
         ib = 1
         iw = iw + 1
         9577 IF (andf(ix(iw),two(ib)) == 0) CYCLE
         nobits = .false.
         fjj = FLOAT(jj)
         CALL map (fjj,fij,xxxx,yyyy)
         IF (sq .AND. ix(ip+1)+ij == ix(ip)+jj) GO TO 9579
         IF (ix(ip+1)+ij == ncols .OR. ix(ip)+jj == nrows) GO TO 9578
         CALL tipe (xxxx,yyyy,1,xstar,1,0)
         CYCLE
         9578 CALL tipe (xxxx,yyyy,1,xdolr,1,0)
         CYCLE
         9579 CALL tipe (xxxx,yyyy,1,xdddd,1,0)
       END DO
     END DO
     IF (prntit) WRITE (nout,9640)
     9640 FORMAT (1H0,29X,100(1H.)/30X,10(9X,1H.))
     IF (prntit) CYCLE
     CALL stplot (-1)
     lnct = lnct + 1
     IF (lnct > nlines) CALL page1
     WRITE  (nout,9645) pfile,i100,j100
     9645 FORMAT (1H ,62X,i5,2I6)
   END DO
   
   9699 GO TO gobac, (2200,3050)
   
   9950 WRITE  (nout,9952) swm,kerror
   9952 FORMAT (a27,' 1705, LOGIC ERROR AT STATEMENT',i5,  &
       ' IN SUBROUTINE SEEMAT.')
   9960 WRITE  (nout,9962) sim,lbl
   9962 FORMAT (a31,' 1702, UTILITY MODULE SEEMAT WILL ABANDON ',  &
       'PROCESSING DATA BLOCK ',2A4 )
   9997 CALL CLOSE (nam,1)
   IF (prntit) CYCLE
   IF (table ) CYCLE
   pfile = pfile + 1
   CALL sopen (*85,ploter,x(kcor),pltbuf)
   CALL stplot (pfile)
   CALL line (bllx,blly,burx,bury,1,-1)
   CALL line (bllx,blly,burx,bury,1, 0)
   CALL line (bulx,buly,burx,bury,1, 0)
   CALL line (bulx,buly,blrx,blry,1, 0)
   CALL line (blrx,blry,bllx,blly,1, 0)
   CALL line (bllx,blly,bulx,buly,1, 0)
   CALL line (burx,bury,blrx,blry,1, 0)
   symbl(1) = 3
   symbl(2) = 6
   CALL map (0.505*fncc,0.505*fncc,xxxx,yyyy)
   CALL symbol (xxxx,yyyy,symbl,-1)
   CALL symbol (xxxx,yyyy,symbl, 0)
   CALL stplot (-1)
   lnct = lnct + 1
   WRITE  (nout,9996) pfile
   9996 FORMAT (63X,i5,2X,13HTRAILER frame)
 END DO
 9999 RETURN
 
END SUBROUTINE seemat
