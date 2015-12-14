SUBROUTINE fa2
     
!     THIS IS THE DMAP MODULE FA2
 
!     DMAP CALLING SEQUENCE
 
!     FA2  PHIH,CLAMA,FSAVE/PHIHL,CLAMAL,CASEYY,OVG/V,N,TSTART/C,Y,VREF/
!    1     C,Y,PRINT=YES $
 
!     ALL OUTPUTS ARE APPEND
 
!     THE PURPOSE OF THIS MODULE IS TO COPY PARTS OF PHIH, CLAMA, AND
!    1    FSAVE ONTO PHIHL, CLAMAL, CASEYY, AND OVG RESPECTIVELY
 
 EXTERNAL        lshift
 INTEGER :: sysbuf,phih,clama,fsave,phihl,clamal,caseyy,ovg,  &
     tstart,PRINT(2),mcb(7),FILE,NAME(2),fmeth,floop,  &
     mcbphl(7),mcbcl(7),mcbcc(7),mcbovg(7),buf(146),  &
     eject,iary(22),ialph(2),me(3),yes,yesb
 REAL :: xmach,kfreq,lbuf(6),iml,z(1)
 COMMON /system/ sysbuf,nout,skp(6),nlpp,mtemp,npag,nlines
 COMMON /zzzzzz/ iz(1)
 COMMON /unpakx/ itc,ii,jj,incr
 COMMON /BLANK / tstart,vref,PRINT
 EQUIVALENCE     (z(1),iz(1))
 DATA    phih  , clama,fsave,phihl,clamal,caseyy,ovg /  &
     101   , 102  ,  103,  201,   202,   203,204 /
 DATA    NAME  , no,mcbcl,mcbcc,mcbovg,iblnk         /  &
     4HFA2 , 1H , 2HNO, 21*0,4H                  /
 DATA    buf   / 146*1H                              /
 DATA    iary  / 4H poi,4HNT =,1H ,1H ,4H mac,4HH = ,1H ,1H ,  &
     4H kfr, 4HEQ= ,1H ,1H ,4H rho,4H =  ,1H ,1H ,6*1H /
 DATA    twophi/ 6.28318531          /
 DATA    me    / 1HK,  2HKE,  2HPK   /
 DATA    yes   , yesb/ 3HYES, 4HYESB /
 
!     INITIALIZE
 
 nz    = korsz(z)
 ibuf1 = nz - sysbuf + 1
 ibuf2 = ibuf1 - sysbuf
 ibuf3 = ibuf2 - sysbuf
 ibuf4 = ibuf3 - sysbuf
 nz    = ibuf4 - 1
 itc   = 3
 incr  = 1
 mcbcl(1) = clamal
 mcbcc(1) = caseyy
 mcbovg(1)= ovg
 IF (vref == 0.0) vref = 1.0
 
!     FIND PROPER METHOD
 
 FILE  = fsave
 CALL OPEN (*900,fsave,iz(ibuf1),0)
 CALL READ (*910,*920,fsave,iz(1),8,1,iflag)
 j     = 3
 fmeth = iz(j)
 meth  = me(fmeth)
 oneok = 1.e+25
 mcb(1)= fsave
 CALL rdtrl (mcb)
 floop = mcb(2)
 nloop = mcb(3)
 nvalue= mcb(7)
 j     = 6
 bref  = z(j)
 phib  = twophi*bref
 SELECT CASE ( fmeth )
   CASE (    1)
     GO TO 1000
   CASE (    2)
     GO TO 2000
   CASE (    3)
     GO TO 3000
 END SELECT
 
!     K  METHOD
 
 1000 CONTINUE
 
!     PICK UP CONSTANTS
 
 nvalue = 8
 nvalue = iz(nvalue)
 
!     COPY ONTO PHIHL
 
 IF (floop /= 1) GO TO 1010
 
!     FIRST TIME
 
 CALL gopen (phihl,iz(ibuf2),1)
 CALL CLOSE (phihl,1)
 mcbphl(1) = phih
 CALL rdtrl (mcbphl)
 mcbphl(2) = 0
 mcbphl(6) = 0
 mcbphl(7) = 0
 mcbphl(1) = phihl
 CALL wrttrl (mcbphl)
 CALL gopen  (clamal,iz(ibuf2),1)
 CALL gopen  (clama,iz(ibuf3),0)
 CALL fread  (clama,buf,146,1)
 CALL CLOSE  (clama,1)
 CALL WRITE  (clamal,buf,146,1)
 CALL WRITE  (clamal,0,0,1)
 CALL CLOSE  (clamal,1)
 CALL gopen  (caseyy,iz(ibuf2),1)
 CALL CLOSE  (caseyy,1)
 CALL gopen  (ovg,iz(ibuf2),1)
 CALL CLOSE  (ovg,1)
 
!     COPY NVALUE VECTORS TO PHIHL
 
 1010 CONTINUE
 mcb(1) = phih
 CALL rdtrl (mcb)
 ncopy = MIN0(nvalue,mcb(2))
 CALL gopen  (phih,iz(ibuf2),0)
 CALL gopen  (phihl,iz(ibuf3),0)
 CALL skpfil (phihl, 1)
 CALL skpfil (phihl,-1)
 CALL CLOSE  (phihl, 2)
 CALL gopen  (phihl,iz(ibuf3),3)
 mcbphl(1) = phihl
 CALL rdtrl  (mcbphl)
 mcbphl(7) = (2*mcbphl(7)*mcbphl(2)*mcbphl(3))/10000
 CALL cyct2b (phih,phihl,ncopy,iz,mcbphl)
 CALL CLOSE  (phih,1)
 CALL CLOSE  (phihl,1)
 CALL wrttrl (mcbphl)
 
!     PICK UP M,K,RHO FOR THIS LOOP
 
 CALL fread (fsave,iz,-3*(floop-1),0)
 CALL fread (fsave,z,3,1)
 j     = 0
 xmach = z(  1)
 kfreq = z(j+2)
 rho   = z(j+3)
 CALL fread (fsave,z,1,1)
 
!     PUT CASEYY INTO CORE
 
 CALL READ (*910,*1020,fsave,iz,nz,0,iflag)
 CALL mesage (-8,0,NAME)
 1020 CONTINUE
 CALL CLOSE (fsave,1)
 k = 39
 DO  i = 51,146
   buf(i) = iz(k)
   k = k + 1
 END DO
 
!     READY CLAMA
 
 CALL gopen  (clama,iz(ibuf1),0)
 CALL fwdrec (*910,clama)
 
!     READY CLAMAL
 
 CALL gopen  (clamal,iz(ibuf2),0)
 CALL skpfil (clamal, 1)
 CALL skpfil (clamal,-1)
 CALL bckrec (clamal)
 CALL READ   (*910,*1022,clamal,iz(iflag+1),nz,0,i)
 CALL mesage (-8,0,NAME)
 1022 CONTINUE
 CALL bckrec (clamal)
 CALL CLOSE  (clamal,2)
 CALL gopen  (clamal,iz(ibuf2),3)
 CALL WRITE  (clamal,iz(iflag+1),i,0)
 CALL rdtrl  (mcbcl)
 
!     READY CASEYY
 
 CALL gopen  (caseyy,iz(ibuf3),0)
 CALL skpfil (caseyy, 1)
 CALL skpfil (caseyy,-1)
 CALL CLOSE  (caseyy, 2)
 CALL gopen  (caseyy,iz(ibuf3),3)
 CALL rdtrl  (mcbcc)
 
!     READY OVG
 
 CALL gopen  (ovg,iz(ibuf4),0)
 CALL skpfil (ovg, 1)
 CALL skpfil (ovg,-1)
 CALL CLOSE  (ovg,2)
 CALL gopen  (ovg,iz(ibuf4),3)
 CALL rdtrl  (mcbovg)
 mcbovg(2)= mcbovg(2) + 1
 mcbcc(4) = iflag
 CALL wrttrl (mcbovg)
 mcbcc(2) = mcbcc(2) + ncopy
 CALL wrttrl (mcbcc)
 mcbcl(2) = mcbcl(2) + ncopy
 CALL wrttrl (mcbcl)
 GO TO 1042
 
!     K-E METHOD
 
 2000 CONTINUE
 
!     P - K METHOD
 
 3000 CONTINUE
 
!     READY OVG
 
 CALL gopen (ovg,iz(ibuf2),1)
 mcbovg(2) = 1
 CALL wrttrl (mcbovg)
 
!     PUT RECORD 2 OF FSAVE INTO CORE
 
 CALL READ (*910,*3010,fsave,iz(1),nz,1,iflag)
 CALL mesage (-8,0,NAME)
 3010 CONTINUE
 CALL skprec (fsave,1)
 CALL fread  (fsave,0,-51,0)
 CALL fread  (fsave,buf,96,1)
 imr   = 1
 floop = 1
 
!     COUNT RHO S
 
 nrho = 1
 IF (fmeth == 3) GO TO 3012
 irho = 1
 rho  = z(imr+2)
 imr1 = imr + 3
 3013 CONTINUE
 IF (imr1 >     iflag) GO TO 3012
 IF (rho  == z(imr1+2)) GO TO 3012
 nrho = nrho + 1
 imr1 = imr1 + 3
 GO TO 3013
 3012 CONTINUE
 3011 CONTINUE
 nv = 1
 
!     DETERMINE THE NUMBER OF M-RHO PAIRS FOR THIS GO
 
 xmach = z(imr  )
 rho   = z(imr+2)
 ncopy = 1
 imr1  = imr + 3*nrho
 3020 CONTINUE
 IF (imr1 > iflag) GO TO 1042
 IF (xmach /= z(imr1) .OR. rho /= z(imr1+2)) GO TO 1042
 ncopy = ncopy + 1
 imr1  = imr1 + 3*nrho
 GO TO 3020
 1042 CONTINUE
 
 IF (PRINT(1) == no) GO TO 1041
!     SET UP PAGE FORMATS
 
 CALL page1
 nlines = nlines + 7
 IF (PRINT(1) == yesb) WRITE (nout,1039) floop,xmach,rho,meth
 IF (PRINT(1) == yes ) WRITE (nout,1040) floop,xmach,rho,meth
 1039 FORMAT (1H0,55X,16HFLUTTER  summary, //7X,  &
     9HPOINT =  ,i3,5X,14HSIGMA value = ,f8.3,4X,  &
     16HDENSITY ratio = ,1P,e11.4,5X,9HMETHOD = ,a4, ///7X,  &
     5HKFREQ,12X, 8H1./kfreq, 9X,8HVELOCITY, 12X,7HDAMPING,  &
     9X,9HFREQUENCY,12X,20HCOMPLEX   eigenvalue)
 1040 FORMAT (1H0,55X,16HFLUTTER  summary, //7X,  &
     9HPOINT =  ,i3, 5X,14HMACH NUMBER = ,f7.4,5X,  &
     16HDENSITY ratio = ,1P,e11.4, 5X,9HMETHOD = ,a4, ///7X,  &
     5HKFREQ, 12X,8H1./kfreq, 9X,8HVELOCITY, 12X,7HDAMPING,  &
     9X,9HFREQUENCY, 12X,20HCOMPLEX   eigenvalue)
 1041 CONTINUE
 
!     SET UP FOR OVG
 
 buf(1) = 60
 buf(2) = 2002
 buf(4) = 1
 buf(5) = 10*floop
 buf(9) = 1
 buf(10)= 4
 CALL WRITE (ovg,buf,146,1)
 IF (fmeth /= 1) GO TO 1101
 DO   i = 115,146
   buf(i) = iblnk
 END DO
 CALL int2a8 (*1092,floop,ialph)
 1092 iary(3) = ialph(1)
 iary(4) = ialph(2)
 CALL re2al (xmach,ialph)
 iary(7) = ialph(1)
 iary(8) = ialph(2)
 CALL re2al (kfreq,ialph)
 iary(11) = ialph(1)
 iary(12) = ialph(2)
 CALL re2al (rho,ialph)
 iary(15) = ialph(1)
 iary(16) = ialph(2)
 k = 115
 DO  i = 1,16
   buf(k) = iary(i)
   k = k + 1
 END DO
 k = 103
 DO  i = 115,146
   iz(k) = buf(i)
   k = k + 1
 END DO
 1101 CONTINUE
 DO   i = 1,ncopy
   SELECT CASE ( fmeth )
     CASE (    1)
       GO TO 1102
     CASE (    2)
       GO TO 1150
     CASE (    3)
       GO TO 3200
   END SELECT
   
!     KE METHOD
   
   1150 CONTINUE
   IF (i /= 1 .OR. nv /= 1) GO TO 1152
   ir = iflag + 1
   j  = nvalue*2
   DO  m = 1,ncopy
     
!     READ A RECORD OF COMPLEX EIGENVALUES INTO CORE
     
     CALL fread  (fsave,iz(ir),j,1)
     CALL skprec (fsave,nrho-1)
     
!     REARRANGE THE COMPLEX EIGENVALUES IN THE RECORD IN ASCENDING
!     ORDER OF THE ABSOLUTE VALUES OF THE IMAGINARY PARTS
     
     nvalu1 = nvalue - 1
     DO  l = 1,nvalu1
       lr = ir + 2*(l-1)
       li = lr + 1
       valuer = z(lr)
       valuei = z(li)
       value  = ABS(valuei)
       INDEX  = l
       l1     = l + 1
       DO  k = l1,nvalue
         kr = ir + 2*(k-1)
         ki = kr + 1
         value1 = ABS(z(ki))
         IF (value1 >= value) CYCLE
         valuer = z(kr)
         valuei = z(ki)
         value  = value1
         INDEX  = k
       END DO
       IF (INDEX == l) CYCLE
       irr = ir  + 2*(INDEX-1)
       iri = irr + 1
       z(irr) = z(lr)
       z(iri) = z(li)
       z(lr)  = valuer
       z(li)  = valuei
     END DO
     ir = ir + j
   END DO
   
!     SELECT EACH FOR OUTPUT
   
   1152 CONTINUE
   j    = iflag + 1 + (i-1)*nvalue*2 + (nv-1)*2
   rel  = z(j)
   iml  = z(j+1)
   vout = ABS(iml)/vref
   g    = 0.0
   IF (iml /= 0.0) g = 2.*rel/iml
   kfreq= z(imr+3*i-2)
   f    = kfreq*iml/phib
   GO TO 1103
   
!     PK METHOD
   
   3200 CONTINUE
   CALL fread (fsave,lbuf,-(nv-1)*5,0)
   CALL fread (fsave,lbuf,5,1)
   rel = lbuf(1)
   iml = lbuf(2)
   kfreq = lbuf(3)
   f = lbuf(4)
   g = lbuf(5)
   vout = ABS(z(imr+3*i-2))/vref
   GO TO 1103
   
!     K METHOD
   
   1102 CONTINUE
   CALL fread (clama ,lbuf,6,0)
   CALL WRITE (clamal,lbuf,6,0)
   rel = lbuf(3)
   iml = lbuf(4)
   vout= ABS(iml)/vref
   g   = 0.0
   IF (iml /= 0.0) g = 2.0*rel/iml
   f =  kfreq*iml/(phib)
   
!     PUT OUT CASEYY
   
   CALL WRITE (caseyy,iz,iflag,1)
   1103 CONTINUE
   IF (PRINT(1) == no) GO TO 1050
   
!     PRINT OUTPUT
   
   k = eject(1)
   IF (k == 0) GO TO 1060
   IF (PRINT(1) == yesb) WRITE (nout,1039) floop,xmach,rho,meth
   IF (PRINT(1) == yes ) WRITE (nout,1040) floop,xmach,rho,meth
   nlines = nlines + 7
   1060 CONTINUE
   IF (kfreq /= 0.0) oneok = 1.0/kfreq
   WRITE  (nout,1070) kfreq,oneok,vout,g,f,rel,iml
   1070 FORMAT (1H ,5X,f8.4,5X,6(1X,1P,e14.7,3X))
   1050 CONTINUE
   
!     PUT OUT OVG PARTS
   
   lbuf(1) = vout
   lbuf(2) = 0.0
   lbuf(3) = g
   lbuf(4) = f
   CALL WRITE (ovg,lbuf,4,0)
 END DO
 floop = floop+1
 CALL WRITE (ovg,0,0,1)
 SELECT CASE ( fmeth )
   CASE (    1)
     GO TO 1031
   CASE (    2)
     GO TO 2031
   CASE (    3)
     GO TO 3331
 END SELECT
 
!     FINISH UP FOR KE METHOD
 
 2031 CONTINUE
 nv = nv + 1
 IF (nv <= nvalue) GO TO 1042
 
!     ALL MODES DONE
 
 IF (irho >= nrho) GO TO 2090
 
!     DO ANOTHER RHO
 
 irho= irho + 1
 imr = imr  + 3
 rho = z(imr+2)
 CALL skprec (fsave,ncopy*(nrho-1))
 GO TO 1042
 2090 CONTINUE
 IF (imr1 > iflag) GO TO 4000
 imr = imr1
 GO TO 3011
 
!     P-K AT POINT END
 
 3331 CONTINUE
 nv = nv + 1
 IF (nv > nvalue) GO TO 3390
 CALL skprec (fsave,-ncopy)
 GO TO 1042
 
!     ALL MODES DONE--CONSIDER MORE M-RHO VALUES
 
 3390 IF (imr1 > iflag) GO TO 4000
 imr = imr1
 GO TO 3011
 
!     DONE
 
 4000 CALL CLOSE (ovg,1)
 CALL CLOSE (fsave,1)
 RETURN
 
!     FINISH UP
 
 1031 CALL WRITE (clamal,0,0,1)
 CALL CLOSE (ovg,1)
 CALL CLOSE (clamal,1)
 CALL CLOSE (clama,1)
 CALL CLOSE (caseyy,1)
 
!     CHECK TIMES
 
 CALL klock  (now)
 CALL tmtogo (itlft)
 IF (now-tstart >= itlft .AND. floop /= nloop) GO TO 1110
 RETURN
 
!     INSUFFICIENT TIME
 
 1110 CALL mesage (45,nloop - floop,NAME)
 tstart = -1
 RETURN
 
!     ERROR MESSAGES
 
 900 ip1 = -1
 GO TO 901
 910 ip1 = -2
 GO TO 901
 920 ip1 = -3
 901 CALL mesage (ip1,FILE,NAME)
 RETURN
END SUBROUTINE fa2
