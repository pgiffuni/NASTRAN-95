SUBROUTINE trlga (casecc,usetd,dlt,slt,bgpdt,sil,cstm,ap,tmldtb,  &
        itrl,iscr1,iscr2,iscr3,est,newslt,mgg,iscr4,mpt1)
     
!     THE PURPOSE OF THIS ROUTINE IS TO CONSTRUCT THE AP MATRIX
!     WHICH HAS 1 COLUMN FOR EACH FUNCTION OF TIME
!     AND TO BUILD THE TIME FUNCTION TABLE (FORMAT SHOWN IN TRLGC)
 
 
 INTEGER, INTENT(IN OUT)                  :: casecc
 INTEGER, INTENT(IN)                      :: usetd
 INTEGER, INTENT(IN)                      :: dlt
 INTEGER, INTENT(IN)                      :: slt
 INTEGER, INTENT(IN)                      :: bgpdt
 INTEGER, INTENT(IN)                      :: sil
 INTEGER, INTENT(IN)                      :: cstm
 INTEGER, INTENT(IN OUT)                  :: ap
 INTEGER, INTENT(IN)                      :: tmldtb
 INTEGER, INTENT(OUT)                     :: itrl
 INTEGER, INTENT(IN)                      :: iscr1
 INTEGER, INTENT(IN OUT)                  :: iscr2
 INTEGER, INTENT(IN)                      :: iscr3
 INTEGER, INTENT(IN)                      :: est
 INTEGER, INTENT(IN)                      :: newslt
 INTEGER, INTENT(IN)                      :: mgg
 INTEGER, INTENT(IN OUT)                  :: iscr4
 INTEGER, INTENT(IN)                      :: mpt1
 EXTERNAL        andf
 INTEGER :: sysbuf,andf,pg(7),NAME(2),slt1,bgpdt1,cstm1,sil1,  &
     mcb(7),iz(38),FILE,namt(2),gvect(30),two1,izb(4), minus(2), est1
 COMMON /BLANK / ng
 COMMON /zzzzzz/ z(1)
 COMMON /loadx / lc,slt1,bgpdt1,OLD,cstm1,sil1,isil,est1,mpt,gptt,  &
     edt,n(3),lodc,mass,nobld,idit
 COMMON /system/ ksystm(65)
 COMMON /bitpos/ isk(11),iue
 COMMON /zblpkx/ za(4),iib
 COMMON /zntpkx/ zb(4),iii,ieol,IEOR
 COMMON /two   / two1(32)
 COMMON /qvect / itran,iqvect
 EQUIVALENCE     (ksystm(1),sysbuf),(z(1),iz(1)),(zb(1),izb(1))
 DATA    NAME  / 4HTRLG,4HA   /,  namt/ 4HDLT ,4HTRLG /
 DATA    itran1, minus /4HTRAN,-1,-1  /
 
!     CORE IS ALLOCATED AS FOLLOWS -
!     . EXTERN PHASE (BUILD STATIC LOADS)
!                                                           POINTER
!     DLOAD STUFF--TLOAD ID,RECORD NO.IN DLT,SCALE FACTOR   ILLST
!     EXTERN LOAD LIST IN SLT ((NEX LENGTH)                 ISLLST
!     2  BUFFERS
!     1  G VECTOR (NG)  COMING FROM TOP
!        N.B.  EXTER WILL OPEN NEWSLT,BGPDT,CSTM,SIL
 
!     . DYNAMIC PHASE
!     DLOAD STUFF                                           ILLST
!     EXTERN LOAD LIST                                      ISLLST
!     SIL TO SILD CONVERTER (NG LENGTH)                     ISILD
!     4  BUFFERS
!     2  P SIZE VECTORS
!     COMPRESSED LIST  SILD,A,TAU                           ICLST
 
!     BRING IN DATA FROM CASECC(DLOAD ID -- TSTEP ID)
 
 nsubl = 0
 nz    = korsz(iz)
 ibuf1 = nz - sysbuf + 1
 nx    = ibuf1 - 1
 CALL gopen (casecc,iz(ibuf1),0)
 CALL fread (casecc,iz(1),166,1)
 idload = iz(13)
 itrl   = iz(38)
 CALL CLOSE (casecc,1)
 IF (idload == 0) GO TO 1020
 
!     BUILD NEW SLT
 
 CALL ssgslt (slt,newslt,est)
 
!     FIND DLOAD, TLOAD
 
 FILE = dlt
 CALL OPEN (*900,dlt,iz(ibuf1),0)
 CALL READ (*910,*10,dlt,iz(1),nx,0,iflag)
 GO TO 980
 
!     IS IT A DLOAD SET
 
 10 ndload = iz(3)
 nsimpl = iflag - 3 - ndload
 IF (ndload == 0) GO TO 100
 k = 3
 DO  i = 1,ndload
   k = k + 1
   IF (iz(k) == idload) GO TO 30
 END DO
 
!     ITS  A SIMPLE LOAD
 
 GO TO 100
 
!     PROCESS DLOAD SET
!     FORMAT OF DLOAD = SET ID,SCALE,SCALE,ID,SCALE,ID .... -1,-1
 
 30 nz1 = nx - iflag
 
!     BRING  IN  ALL  DLOADS
 
 l =  iflag + 1
 CALL READ (*910,*40,dlt,iz(l),nz1,0,i)
 GO TO 980
 
!     FIND SELECTED ID
 
 40 isel = l
 50 IF (iz(isel) == idload) GO TO 70
 60 isel = isel + 2
 IF (iz(isel+1) /= -1) GO TO 60
 isel = isel + 2
 IF (isel-l > i) GO TO 990
 GO TO 50
 
!     FOUND DLOAD SELECTED
 
 70 scale = z(isel+1)
 
!     CONVERT SCALE FACTORS TO OVERALL SCALE FACTORS
!     BUILD LIST OF TRIPLES-- TLOAD ID,RECORD NO.IN DLT, SCALE FACTOR
 
 l = isel + 2
 m = isel + i
 iflag  = m
 nsubl  = 0
 80 idload = iz(l+1)
 z(l)   = z(l)*scale
 k = ndload + 3
 DO  i = 1,nsimpl
   k = k + 1
   IF (iz(l+1) == iz(k)) GO TO 95
   
 END DO
 GO TO 990
 
!     FOUND SIMPLE ID
 
 95 iz(m  ) = iz(l+1)
 z(m+1) = z(l)
 iz(m+2) = i
 l = l + 2
 m = m + 3
 nsubl = nsubl + 1
 IF (iz(l+1) >= 0) GO TO 80
 GO TO 150
 
!     PROCESS SIMPLE LOAD REQUEST
 
 100 m = iflag + 1
 iflag  = m
 iz(m ) = idload
 z(m+1) = 1.0
 l = ndload + 3
 DO  i = 1,nsimpl
   l = l + 1
   IF (iz(l) == idload) GO TO 120
 END DO
 GO TO 990
 
!     FOUND SIMPLE LOAD
 
 120 IF (ndload /= 0) i = i + 1
 iz(m+2) = i - 1
 nsubl   = 1
 
!     MOVE STUFF TO BOTTOM OF CORE
 
 150 CALL CLOSE(dlt,1)
 illst = nz - nsubl*3 + 1
 nz    = nz - nsubl*3
 ibuf1 = nz - sysbuf + 1
 l     = iflag
 k     = illst
 DO  i = 1,nsubl
   CALL gopen  (dlt,iz(ibuf1),0)
   CALL skprec (dlt,iz(l+2))
   CALL fread  (dlt,izb,2,0)
   iz(k) = izb(2)
   CALL CLOSE (dlt,1)
   iz(k+1) = iz(l+2)
   iz(k+2) = iz(l+1)
   l = l + 3
   k = k + 3
 END DO
 
!     SET UP FOR EXTERN
 
 FILE   = newslt
 nx     = ibuf1 - 1
 isllst = illst
 noslt  = 0
 mcb(1) = slt
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) noslt = -1
 mcb(1) = sil
 mcb(3) = 0
 CALL rdtrl (mcb)
 ng     = mcb(3)
 IF (noslt /= 0) GO TO 191
 CALL OPEN (*900,newslt,iz(ibuf1),0)
 CALL READ (*910,*170,newslt,iz(1),nx,0,iflag)
 GO TO 980
 170 CALL CLOSE (newslt,1)
 m = illst
 DO  i = 1,nsubl
   DO  j = 3,iflag
     IF (iz(m) /= iz(j)) CYCLE
     
!     FOUND LOAD TO BUILD
     
     iz(j) = -IABS(iz(j))
     EXIT
   END DO
   179 m = m + 3
 END DO
 
!     ZERO LOADS NOT TO BUILD
 
 m = illst - iflag + 2
 isllst = m
 DO  j = 3,iflag
   IF (iz(j) < 0) GO TO 185
   iz(m) = 0
   GO TO 189
   185 iz(m) = IABS(iz(j))
   189 m     = m + 1
 END DO
 nex   = iflag - 2
 nz    = nz - nex
 ngrav = 0
 iharm = 0
 n1    = nex
 ibuf1 = nz - sysbuf + 1
 ibuf2 = ibuf1 - sysbuf
 
!     SET UP SCRATCH FILE FOR QLOADL
 
 itran  = itran1
 iqvect = iscr1
 CALL gopen (iscr1,iz(ibuf1),1)
 CALL makmcb (pg,iscr2,ng,2,1)
 slt1   = newslt
 bgpdt1 = bgpdt
 cstm1  = cstm
 sil1   = sil
 est1   = est
 mass   = mgg
 mpt    = mpt1
 CALL gopen (pg,iz(ibuf2),1)
 lc     = ibuf2 - 1
 CALL extern (nex,ngrav,gvect,iz(isllst),pg,n1,iharm)
 CALL CLOSE  (pg,1)
 CALL wrttrl (pg)
 CALL WRITE  (iscr1,minus,2,1)
 CALL CLOSE  (iscr1,1)
 IF (ngrav == 0) GO TO 191
 
!     DO GRAVITY LOADS
 
 mcb(1) = mgg
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) CALL mesage (-56,0,nave)
 
!     SAVE LOAD LIST IN CORE
 
 CALL gopen (iscr4,iz(ibuf2),1)
 CALL WRITE (iscr4,iz(isllst),3*nsubl+nex,1)
 CALL CLOSE (iscr4,1)
 CALL gravl1 (ngrav,gvect,iscr3,iharm)
 CALL ssg2b (mgg,iscr3,0,tmldtb,0,1,1,ap)
 CALL gravl2 (ngrav,tmldtb,pg)
 n1 = n1 + ngrav
 
!     RESTORE LOAD LIST TO CORE
 
 CALL gopen (iscr4,iz(ibuf2),0)
 CALL fread (iscr4,iz(isllst),3*nsubl+nex,1)
 CALL CLOSE (iscr4,1)
 
!     BUILD SIL  TO SILD CONVERTER
 
 191 CONTINUE
 FILE   = usetd
 CALL gopen (usetd,iz(ibuf1),0)
 mcb(1) = usetd
 CALL rdtrl (mcb)
 lusetd = mcb(2)
 CALL fread (usetd,iz(1),lusetd,1)
 CALL CLOSE (usetd,1)
 isild  = isllst - ng
 mskue  = two1(iue)
 l    = isild
 DO   i = 1,lusetd
   IF (andf(iz(i),mskue) /= 0) CYCLE
   iz(l)= i
   l    = l + 1
 END DO
 nz   = nz - ng
 
!     BEGIN LOOP ON EACH TLOAD CARD
 
 ibuf1 = nz - sysbuf + 1
 iclst = 2*lusetd + 1
 ibuf2 = ibuf1 - sysbuf
 ibuf3 = ibuf2 - sysbuf
 CALL makmcb (mcb,ap,lusetd,2,1)
 CALL gopen (ap,iz(ibuf2),1)
 iterm = 0
 CALL gopen (tmldtb,iz(ibuf3),1)
 iqvrn = 0
 ibuf4 = ibuf3 - sysbuf
 CALL gopen (iscr3,iz(ibuf4),1)
 nz    = ibuf4 - 1
 IF (nz < 5*lusetd) GO TO 980
 DO  iloop = 1,nsubl
   
!     ZERO AP AND TAU AREA
   
   k = 2*lusetd
   DO  i = 1,k
     z(i) = 0.0
   END DO
   
!     FIND APPROPRIATE STATIC LOAD
   
   k = illst + (iloop-1)*3
   scale  = z(k+2)
   idload = iz(k )
   idltr  = iz(k+1)
   IF (noslt /= 0) GO TO 300
   k = isllst - 1
   m = 0
   DO  i = 1,nex
     l = k + i
     IF (iz(l) == idload) GO TO 221
     IF (iz(l) /= 0) m = m + 1
   END DO
   GO TO 300
   
!     POSITION TO PROPER AP RECORD
   
   221 FILE = pg(1)
   CALL gopen (pg,iz(ibuf1),0)
   CALL skprec (pg,m)
   CALL intpk (*290,pg,0,1,0)
   240 IF (ieol /= 0) GO TO 290
   CALL  zntpki
   zb(1) = zb(1)*scale
   k = isild + iii - 1
   k = iz(k)
   z(k) = zb(1)
   GO TO 240
   290 CALL CLOSE (pg,1)
   
!     PROCESS DLT STUFF
   
   300 CALL gopen (dlt,iz(ibuf1),0)
   FILE = dlt
   CALL skprec (dlt,idltr)
   CALL fread (dlt,gvect,8,0)
   
!     READS AND BUILDS COMPRESSED LIST SILD,AI,TAU,FOR ALL AI.S
   
   320 CALL READ (*910,*330,dlt,izb,4,0,iflag)
   l    = izb(1)
   z(l) = zb(2) + z(l)
   z(l+lusetd) = zb(3)
   GO TO 320
   330 CALL CLOSE (dlt,1)
   iqr = 0
   ASSIGN 370 TO iretn
   m = 0
   k = iclst
   DO  i = 1,lusetd
     IF (z(i) == 0.0) CYCLE
     z(i)   = z(i)*scale
     m      = m + 1
     iz(k ) = i
     z(k+1) = z(i)
     z(k+2) = z(i+lusetd)
     k      = k + 3
   END DO
   
!     SORT ON TAU
   
   339 k = 3*m + iclst - 4
   IF (iclst > k) GO TO 1335
   DO   l = iclst,k,3
     IF (z(l+5) > z(l+2) .OR. (z(l+5) == z(l+2) .AND.  &
         iz(l+3) >= iz(l))) CYCLE
     ll      = l
     iz(k+4) = iz(l+3)
     z(k+5)  = z(l+4)
     z(k+6)  = z(l+5)
     338 iz(ll+3)= iz(ll)
     z(ll+4) = z(ll+1)
     z(ll+5) = z(ll+2)
     ll      = ll - 3
     IF (ll >= iclst .AND. (z(k+6) < z(ll+2) .OR. (z(k+6) == z(ll+2)  &
         .AND. iz(k+4) < iz(ll)))) GO TO 338
     iz(ll+3)= iz(k+4)
     z(ll+4) = z(k+5)
     z(ll+5) = z(k+6)
   END DO
   1335 CONTINUE
   
!     OUTPUT PVECTOR FOR EACH UNIQUE TAU
   
   l    = iclst
!WKBR 8/94 ALPHA  341 TAUO = Z(L+2)
   341 itauo = iz(l+2)
   CALL bldpk (1,1,ap,0,0)
   345 za(1)= z(l+1)
   iib  = iz(l)
   CALL zblpki
   l = l +3
!WKBR 8/94 ALPHA IF (L.LT.3*M+ICLST .AND. Z(L+2).EQ.TAUO) GO TO 345
   IF (l < 3*m+iclst .AND. iz(l+2) == itauo) GO TO 345
   CALL bldpkn (ap,0,mcb)
   
!     PUT OUT LINE OF TIME TABLE
   
   iterm = iterm + 1
   CALL WRITE (tmldtb,iterm,1,0)
   CALL WRITE (tmldtb,idload,1,0)
   CALL WRITE (tmldtb,gvect,1,0)
!WKBR 8/94 ALPHA CALL WRITE (TMLDTB,TAUO,1,0)
   CALL WRITE (tmldtb,itauo,1,0)
   CALL WRITE (tmldtb,gvect(3),6,0)
   CALL WRITE (tmldtb,iqr,1,0)
   IF (l >= iclst+3*m) GO TO iretn, (370,390)
   GO TO 341
   
!     FIND PROPER QVEC RECORD
   
   370 CONTINUE
   IF (noslt /= 0) CYCLE
   CALL gopen (iscr1,iz(ibuf1),0)
   FILE = iscr1
   380 CALL READ (*450,*920,iscr1,iqvid,1,0,iflag)
   IF (iqvid ==     -1) GO TO 450
   IF (iqvid == idload) GO TO 390
   CALL fwdrec (*910,iscr1)
   GO TO 380
   
!     BUILD LIST OF SILD,AI,TAU FROM QVEC STUFF
   
   390 CALL fread (iscr1,m,1,0)
   k = iclst
   IF (m == -1) GO TO 450
   DO  i = 1,m
     CALL fread (iscr1,zb,2,0)
     zb(2) = zb(2)*scale
     j = isild + izb(1) - 1
     j = iz(j)
     iz(k)  = j
     z(k+1) = zb(2)
     z(k+2) = z(j+lusetd)
     k = k + 3
   END DO
   iqvrn = iqvrn + 1
   iqr   = iqvrn
   CALL fread (iscr1,iz(k),9,0)
   CALL WRITE (iscr3,iz(k),9,0)
   ASSIGN 390 TO iretn
   GO TO 339
   
!     END OF QVECT PROCESSING
   
   450 CALL CLOSE (iscr1,1)
   
!     END OF TLOAD CARD LOOP
   
 END DO
 CALL CLOSE (ap,1)
 CALL wrttrl (mcb)
 CALL CLOSE (iscr3,1)
 
!     APPEND QVECT STUFF TO TMLDTB
 
 CALL gopen (iscr3,iz(ibuf1),0)
 FILE = iscr3
 CALL WRITE (tmldtb,0,0,1)
 CALL READ (*1010,*1010,iscr3,iz(1),nz,0,iflag)
 GO TO 980
 1010 CALL WRITE (tmldtb,iz(1),iflag,1)
 CALL CLOSE (tmldtb,1)
 mcb(1) = tmldtb
 mcb(2) = iterm
 mcb(3) = iflag
 CALL wrttrl (mcb)
 CALL CLOSE (iscr3,1)
 1020 CONTINUE
 RETURN
 
!     FATAL ERRORS
 
 900 ip1 = -1
 901 CALL mesage (ip1,FILE,NAME)
 RETURN
 903 CALL mesage (-61,0,NAME)
 RETURN
 910 ip1 = -2
 GO TO 901
 920 ip1 = -3
 GO TO 901
 980 CALL mesage (-8,0,NAME)
 GO TO 903
 990 CALL mesage (-31,idload,namt)
 RETURN
END SUBROUTINE trlga
