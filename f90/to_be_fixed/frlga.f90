SUBROUTINE frlga (dlt,frl,casecc,dit,pp,lusetd,nfreq,nload,  &
        frqset,fol,notrd)
     
!     THIS ROUTINE GENERATES LOADS INCORE AT EACH FREQUENCY
 
!     WITH ENTRY POINTS - GUST1A AND FRRD1A
!                         ======     ======
 
 
 INTEGER, INTENT(IN)                      :: dlt
 INTEGER, INTENT(IN)                      :: frl
 INTEGER, INTENT(OUT)                     :: casecc
 INTEGER, INTENT(IN OUT)                  :: dit
 INTEGER, INTENT(IN)                      :: pp
 INTEGER, INTENT(IN)                      :: lusetd
 INTEGER, INTENT(IN)                      :: nfreq
 INTEGER, INTENT(OUT)                     :: nload
 INTEGER, INTENT(OUT)                     :: frqset
 INTEGER, INTENT(IN)                      :: fol
 INTEGER, INTENT(OUT)                     :: notrd
 INTEGER :: sysbuf, icore(14), FILE,mcb(7),ihead(8),itlist(13),NAME(6)
 REAL :: fx(2)
 COMPLEX :: pow,eb,r2,r1
 DIMENSION       head(8)
 COMMON /system/ ksystm(55)
 COMMON /BLANK / xx
 COMMON /packx / it1,it2,ii,jj,incr
 COMMON /zzzzzz/ core(1)
 COMMON /condas/ pi,twophi,radeg,degra,s4pisq
 COMMON /frrdst/ ovr(152),itl(3)
 EQUIVALENCE     (core(1),icore(1)), (head(1),ihead(1),isil),  &
                 (head(2),a), (head(3),tau), (head(4),theta),  &
                 (ksystm(1),sysbuf), (ksystm(55),iprec)
 DATA    itlist/ 4,1105,11,1,1205,12,2,1305,13,3,1405,14,4 /
 DATA    NAME  / 4HDLT ,4HFRLG,4HA   ,4HGUST,4H1A  ,4HFRRD /
 DATA    ifrl  / 4HFRL /
 
!     IDENTIFICATION OF VARIABLES
 
!     NFREQ  = NUMBER OF FREQ IN SELECTED FREQ SET
!     NDONE  =  NUMBER OF FREQUENCIES CURRENTLY BUILT FOR CUR LOAD
!     LLIST  = POINTER TO START OF LOAD TABLE
!     ITABL  = POINTER TO START OF LIST OF TABLES NEEDED FOR CURRENT
!              LOAD
!     ILOAD  = POINTER TO BEGINNING OF LOADS IN CORE
!     IFL    = POINTER TO VALUES OF FREQ  FUNCTIONS
!     NBUILD = NUMBER OF FREQUENCIES WHICH CAN BE BUILT AT ONCE
!     NLOAD  = NUMBER OF LOADS FOUND IN CASE CONTROL
!     LCORE  = AMOUNT OF CORE AVAILABLE TO HOLD  LOADS + F(F)-S
!     FRQSET = SELECT FREQUENCY SET ID
!     LOADN  = SELECTED DYNAMIC LOAD
!     NDLOAD = NUMBER OF DLOAD CARDS
!     NSIMPL = NUMBER OF SIMPLE LOADS
!     NSUBL  = NUMBEL OF  SIMPLE LOADS COMPOSING PRESENT LOAD
!     NTABL  = NUBER OF TABLE ID-S IN PRESENT LOAD
!     ICDTY  = CARD TYPE CODE  1=RLOAD1,  2=RLOAD2
 
 
 GO TO 2
 
 
 ENTRY gust1a (dlt,frl,casecc,dit,pp,lusetd,nfreq,nload, frqset,fol,notrd)
!     =======================================================
 
 NAME(2) = NAME(4)
 NAME(3) = NAME(5)
 GO TO 2
 
 
 ENTRY frrd1a (dlt,frl,casecc,dit,pp,lusetd,nfreq,nload, frqset,fol,notrd)
!     =======================================================
 
 NAME(2) = NAME(6)
 NAME(3) = NAME(5)
 
 
!     INITALIZE
 
 2 it1   = 3
 it2   = 2 + iprec
 ii    = 1
 jj    = lusetd
 incr  = 1
 notrd =-1
 lcore = korsz(core(1))
 
!     PICK UP AND STORE FREQUENCY SET
 
 ibuf  = lcore - sysbuf + 1
 nz1   = ibuf  - 1
 lcore = lcore - 2*sysbuf
 nz    = lcore
 igust = 0
 IF (casecc > 0) GO TO 5
 casecc = IABS(casecc)
 igust  = 1
 5 CONTINUE
 FILE   = casecc
 CALL OPEN (*510,casecc,core(ibuf),0)
 CALL fwdrec (*530,casecc)
 CALL fread  (casecc,core,149,0)
 frqset = icore(14)
 nload  = 0
 loadn  = icore(13)
 CALL CLOSE (casecc,1)
 itl(1) = 2
 i149   = 149
 itl(2) = icore(i149)
 itl(3) = itl(2) + 1
 itld   = 1
 
!     BRING IN AND SAVE FREQ LIST -- CONVERT  W-S TO F    F = TWOPHI* W
 
 FILE  = frl
 CALL OPEN (*510,frl,core(ibuf),0)
 CALL READ (*530,*10,frl,core(1),nz1,0,iflag)
 GO TO 540
 10 DO  i = 3,iflag
   IF (icore(i) == frqset) GO TO 30
 END DO
 NAME(1) = ifrl
 CALL mesage (-31,frqset,NAME)
 30 k = i-3
 IF (k == 0) GO TO 50
 DO  i = 1,k
   CALL fwdrec (*530,frl)
 END DO
 
!     READ IN  FREQ LIST
 
 50 CALL READ (*530,*60,frl,core(1),nz1,0,nfreq)
 GO TO 540
 60 CALL CLOSE (frl,1)
 lcore = lcore - nfreq
 nz1   = nz1   - nfreq
 frqset= k + 1
 llist = nfreq + 1
 
!     CONVERT TO F
 
 DO  i = 1,nfreq
   core(i) = core(i)/twophi
 END DO
 
!     PUT HEADER ON LOAD FILE
 
 FILE = pp
 nz   = ibuf - sysbuf
 nz1  = nz1  - sysbuf
 CALL OPEN (*510,pp,core(nz),1)
 CALL fname (pp,mcb(1))
 CALL WRITE (pp,mcb(1),2,0)
 CALL WRITE (pp,core(1),nfreq,1)
 FILE = fol
 CALL OPEN (*71,fol,core(ibuf),1)
 CALL fname (fol,mcb)
 CALL WRITE (fol,mcb,2,0)
 CALL WRITE (fol,core,nfreq,1)
 CALL CLOSE (fol,1)
 mcb(1) = fol
 mcb(2) = nfreq
 mcb(3) = frqset
 CALL wrttrl (mcb)
 71 CONTINUE
 
!     SET UP MCB FOR PP
 
 mcb(1) = pp
 mcb(2) = 0
 mcb(3) = lusetd
 mcb(4) = 2
 mcb(5) = 2 + iprec
 mcb(6) = 0
 mcb(7) = 0
 
!     BEGIN LOOP ON LOADS SELECTED
 
 80 IF (nload == 0) GO TO 100
 FILE = casecc
 CALL OPEN (*510,casecc,core(ibuf),0)
 l = nload + 1
 DO  i = 1,l
   CALL fwdrec (*530,casecc)
 END DO
 CALL READ (*500,*540,casecc,core(llist),16,1,iflag)
 loadn = icore(llist+12)
 CALL CLOSE (casecc,1)
 100 nload = nload + 1
 IF (loadn == 0) GO TO 491
 ndone = 0
 lcore = nz1
 
!     FIND SELECTED LOAD IN DLT
 
 FILE = dlt
 CALL OPEN (*510,dlt,core(ibuf),0)
 CALL READ (*530,*110,dlt,core(llist),nz1,0,iflag)
 
!     IS IT A DLOAD SET
 
 110 ndload = icore(llist+2)
 nsimpl = iflag - 3 - ndload
 IF (nsimpl == 0) CALL mesage (-31,loadn,NAME)
 IF (ndload == 0) GO TO 300
 k = llist + 2
 DO  i = 1,ndload
   k = k + 1
   IF (icore(k) == loadn) GO TO 130
 END DO
 GO TO 300
 
!     PROCESS DLOAD SET
 
!     FORMAT OF DLOAD CARD = SET ID, SCALE,SCALE,ID, SCALE, ID, ...,0,-1
 
 130 nz1 = nz1 - iflag
 
!     BRING IN ALL DLOADS
 
 l = llist + iflag
 CALL READ (*530,*140,dlt,core(l),nz1,0,i)
 GO TO 540
 
!     FIND SELECTED ID
 
 140 isel  = l
 150 IF (icore(isel) == loadn) GO TO 170
 160 isel = isel + 2
 IF (icore(isel+1) /= -1) GO TO 160
 isel = isel + 2
 GO TO 150
 
!     FOUND LOAD SET  SELECTED
 
 170 scale  = core(isel+1)
 
!     CONVERT  SCALE FACTORS TO OVERALL  SCALE +ID-S TO RECORD NUMBERS-1
 
 l = isel + 2
 nsubl = 0
 180 core(l) = core(l)*scale
 k = llist + 2 + ndload
 DO  i = 1,nsimpl
   k = k + 1
   IF (icore(l+1) == icore(k)) GO TO 200
 END DO
 CALL mesage (-31,icore(l),NAME)
 
!     FOUND SIMPLE ID
 
 200 icore(l+1) = i + 1
 nsubl = nsubl + 1
 l = l + 2
 IF (icore(l+1) >= 0) GO TO 180
 
!     MOVE TO LOAD LIST AREA
 
 l = isel + 2
 k = llist
 DO  i = 1,nsubl
   icore(k)  = icore(l+1)
   core(k+1) = core(l)
   l = l + 2
   k = k + 2
 END DO
 
!     BUILD LIST OF UNIQUE TABLES NEEDED FOR NSUBL LOADS
 
 ipos  = 2
 230 ntabl = 0
 itabl = llist + 2*nsubl
 DO  i = 1,nsubl
   k = llist + (i-1)*2
   j = icore(k)
   l = j - ipos
   IF (l == 0) GO TO 250
   DO  k = 1,l
     CALL fwdrec (*530,dlt)
   END DO
   
!     READ IN DESCRIPTOR WORDS
   
   250 ipos = j + 1
   CALL READ (*530,*550,dlt,head(1),8,1,iflag)
   icdty = ihead(1)
   nt    = 4
   SELECT CASE ( icdty )
     CASE (    1)
       GO TO 251
     CASE (    2)
       GO TO 251
     CASE (    3)
       GO TO 252
     CASE (    4)
       GO TO 291
   END SELECT
   
!     TLOAD 1 CARD
   
   252 nt    = 3
   itld  = 2
   notrd = 1
   251 CONTINUE
   loop280:  DO  m = 3,nt
     IF (ihead(m) == 0) CYCLE loop280
     IF( ntabl    == 0) GO TO 270
     DO  k = 1,ntabl
       l  = itabl+k
       IF (icore(l) == ihead(m)) CYCLE loop280
     END DO
     
!     STORE NEW TABLE ID
     
     270 ntabl = ntabl + 1
     k = itabl + ntabl
     icore(k) = ihead(m)
   END DO loop280
   CYCLE
   
!     TLOAD2 CARD
   
   291 CONTINUE
   notrd = 1
 END DO
 CALL REWIND (dlt)
 lcore = lcore - ntabl - 1
 iload = itabl + ntabl + 1
 icore(itabl)  = ntabl
 GO TO 330
 
!     PROCESS SIMPLE LOAD REQUEST
 
 300 nsubl = 1
 core(llist+1) = 1.0
 l = llist + 2 + ndload
 DO  i = 1,nsimpl
   l = l + 1
   IF (icore(l) == loadn) GO TO 320
 END DO
 CALL mesage (-31,loadn,NAME)
 
!     FOUND SIMPLE LOAD  STORE RECORD NUMBER
 
 320 IF (ndload /= 0) i = i + 1
 icore(llist) = i
 ipos  = 1
 lcore = lcore - 2
 GO TO 230
 
!     ALLOCATE CORE
 
 330 lvect  = 2*lusetd
 nbuild = lcore/(lvect+ntabl*itld)
 nbuild = MIN0(nbuild,nfreq)
 IF (nbuild == 0) GO TO 540
 kk  = ntabl*nbuild
 ifl = nz - ntabl*nbuild*itld
 
!     LOOP HERE FOR FREQUENCY SPILL
 
 lcore = lcore - ntabl*nbuild
 nbuf  = lcore - sysbuf
 IF (ntabl == 0) GO TO 361
 340 CALL pretab (dit,core(iload),core(iload),core(nbuf),nbuf,l,  &
     core(itabl),itlist(1))
 DO  j = 1,ntabl
   l = itabl + j
   DO  i = 1,nbuild
     m = ndone + i
     k = ifl + nbuild*(j-1) + i - 1
     IF (itld == 2) GO TO 341
     
!                 TAB      X       F(X)
     CALL tab (core(l),core(m),core(k))
     CYCLE
     
!     TRANSFOR LOOK UP FOR TLOAD 1 CARDS
     
     341 CONTINUE
     CALL tab1 (core(l),core(m),fx(1))
     core(k   ) = fx(1)
     core(k+kk) = fx(2)
     CYCLE
   END DO
 END DO
 361 CONTINUE
 
!     READY CORE FOR BUILDING LOADS
 
 k = iload - 1
 DO  i = 1,nbuild
   DO  l = 1,lvect
     k = k + 1
     core(k) = 0.0
   END DO
 END DO
 
!     POSITION TO LOAD IN DLT
 
 ipos = 0
 DO  i = 1,nsubl
   k = llist + 2*i - 2
   l = icore(k) - ipos
   scale = core(k+1)
   IF (l  == 0) GO TO 400
   DO  j = 1,l
     CALL fwdrec (*530,dlt)
   END DO
   
!     READ IN 8 WORD LOAD ID
   
   400 ipos  = l + 1 + ipos
   CALL READ (*530,*540,dlt,head(1),8,0,iflag)
   icdty = ihead(1)
   tk1   = head(3)
   tk2   = head(4)
   nt    = 4
   SELECT CASE ( icdty )
     CASE (    1)
       GO TO 404
     CASE (    2)
       GO TO 404
     CASE (    3)
       GO TO 403
     CASE (    4)
       GO TO 435
   END SELECT
   403 nt    = 3
   
!     FIND COEFFICIENTS IN TABLE LIST
   
   404 DO  k = 3,nt
     IF (ihead(k) /= 0) GO TO 405
     ihead(k+3) = -1
     CYCLE
     405 DO  l = 1,ntabl
       m = itabl + l
       IF (icore(m) == ihead(k)) GO TO 420
     END DO
     GO TO 550
     
!     COMPUTE POINTER INTO COEF TABLE
     
     420 ihead(k+3) = ifl + (l-1)*nbuild
     IF (icdty == 3) ihead(k+4) = ifl + (l-1)*nbuild + ntabl*nbuild
   END DO
   
!     REPEATLY READ IN  4  WORDS --SIL,A,TAU,THETA
   
   435 igust1 = 0
   440 CONTINUE
   IF (igust  == 0) GO TO 442
   IF (igust1 == 1) CYCLE
   igust1 = 1
   442 CONTINUE
   CALL READ (*530,*480,dlt,ihead(1),4,0,iflag)
   IF (igust == 0) GO TO 443
   isil  = 1
   a     = 1.0
   tau   = 0.0
   theta = 0.0
   443 CONTINUE
   a     = a*scale
   theta = theta*degra
   DO  j = 1,nbuild
     IF (icdty == 4) GO TO 448
     
!     COMPUTE COEFFICIENTS
     
     c1 = 0.0
     IF (ihead(6) < 0) GO TO 445
     k  = ihead(6) + j - 1
     c1 = core(k)
     445 c2 = 0.0
     IF (ihead(7) < 0) GO TO 448
     k  = ihead(7) + j - 1
     c2 = core(k)
     448 l  = ndone + j
     m  = (j-1)*lvect + 2*isil - 2 + iload
     SELECT CASE ( icdty )
       CASE (    1)
         GO TO 450
       CASE (    2)
         GO TO 460
       CASE (    3)
         GO TO 450
       CASE (    4)
         GO TO 471
     END SELECT
     
!     RLOAD 1 CARDS OF TLOAD1 CARDS
     
     450 xlama = theta - core(l)*tau*twophi
     sinxl = SIN(xlama)
     cosxl = COS(xlama)
     core(m  ) = a*(c1*cosxl - c2*sinxl) + core(m  )
     core(m+1) = a*(c1*sinxl + c2*cosxl) + core(m+1)
     CYCLE
     
!     RLOAD2  CARDS
     
     460 xlama = theta - core(l)*tau*twophi + c2*degra
     core(m  ) = a*c1*COS(xlama) + core(m  )
     core(m+1) = a*c1*SIN(xlama) + core(m+1)
     CYCLE
     
!     TLOAD 2 CARDS
     
     471 CONTINUE
     f   = head(5)
     p   = head(6)*degra
     c   = head(7)
     ib  = head(8)  +.5
     dt  = tk2 - tk1
     rz  =-c*dt
     cz  =-dt*(f-core(l))*twophi
     
!     COMPUTE  E(B+1) (ZR2)
     
     CALL frr1a1 (rz,cz,ib+1,reb,ceb)
     eb = CMPLX(reb,ceb)
     rp =-rz
     cp = p - core(l)*twophi*tk2 + twophi*f*dt
     pow= CMPLX(rp,cp)
     r2 = CEXP(pow)*eb
     
!     COMPUTE  R1
     
     cz = -dt*(-f -core(l))*twophi
     
!     COMPUTE  E(B+1)ZR1
     
     CALL frr1a1 (rz,cz,ib+1,reb,ceb)
     eb  = CMPLX(reb,ceb)
     cp  =-p - core(l)*twophi*tk2 - twophi*f*dt
     pow = CMPLX(rp,cp)
     r1  = r2 + CEXP(pow)*eb
     
!     COMPUTE   P(W)
     r2  = CMPLX(0.,-core(l)*tau*twophi)
     pow = r1*CEXP(r2)
     cp  = (dt**(ib+1))/(2.0 *(head(8)+1.))
     rz  = REAL(pow)*a*cp
     cz  = AIMAG (pow)*a*cp
     core(m  ) = core(m  ) + rz
     core(m+1) = core(m+1) + cz
     CYCLE
   END DO
   GO TO 440
   
!     END OF STUFF IN DLT TABLE
   
 END DO
 
!     PACK OUT LOADS BUILT
 
 DO  i = 1,nbuild
   m = (i-1)*lvect + iload
   CALL pack (core(m),pp,mcb(1))
 END DO
 ndone  = ndone  + nbuild
 nbuild = MIN0(nbuild,nfreq-ndone)
 CALL REWIND (dlt)
 IF (nbuild /= 0) GO TO 340
 CALL CLOSE (dlt,1)
 GO TO 80
 
!     BUILD ZERO LOAD
 
 491 DO  i = 1,nfreq
   CALL bldpk (3,3,pp,0,0)
   CALL bldpkn (pp,0,mcb)
 END DO
 GO TO 80
 
!     EOF  ON CASECC  END OF ROUTINE
 
 500 CALL CLOSE (casecc,1)
 CALL wrttrl (mcb(1))
 CALL CLOSE (pp,1)
 RETURN
 
!     ERROR MESAGES
 
 510 ip1 = -1
 520 CALL mesage (ip1,FILE,NAME(2))
 530 ip1 = -2
 GO TO 520
 540 ip1 = -8
 GO TO 520
 550 ip1 = -7
 GO TO 520
 
END SUBROUTINE frlga
