SUBROUTINE ift
     
    !     INVERSE FOURIER TRANSFORM MODULE (IFT)
 
    !     DMAP CALLING SEQ.
 
    !     IFT   UHVF,CASECC,TRL,FOL/UHVT,TOL/C,Y,IFTM
 
    INTEGER :: sysbuf,iz(1),uhvf,uhvt,casecc,trl,tol,fol,NAME(2),  &
               mcb(7),FILE,mcb1(7)

    COMMON /system/sysbuf,nout
    COMMON /packx/it1,it2,ii,jj,incr
    COMMON /unpakx/it3,ii1,jj1,incr1
    COMMON /condas/ phi,twopi
    COMMON  /zzzzzz/ z(1)
    COMMON /BLANK/ iftm

    EQUIVALENCE (z(1),iz(1))

    DATA uhvf,casecc,trl,fol,uhvt,tol/101,102,103,104,201,202/
    DATA NAME /4HIFT , 1H  /
 
    !     VARIABLE CORE
 
    !        CONTENT         LENGTH         POINTER
    !        -------         ------         -------
    !     FOL                NFREQ          IFREQ
    !     TSTEP              NGROUP*3       ITSTP
    !     UHVF               NMODES*NFREQ*2 IUHVF
    !     CK                 NBIG           ICK
    !     SK                 NBIG           ISK
    !     UDOT               NMODES*NFREW*2 IUDOT
    !     UHVT               NMODES         IUVT
 
    !     PUT FOL INTO CORE
    nz = korsz(iz)
    ibuf1 = nz-sysbuf+1
    ibuf2 = ibuf1-sysbuf
    nz = nz-2*sysbuf
    FILE = fol
    CALL OPEN(*900,fol,iz(ibuf1),0)
    CALL fread(fol ,iz,-2,0)
    CALL READ(*910,*10,fol,iz,nz,0,nfreq)
    CALL mesage(-8,0,NAME)
10  CALL CLOSE(fol ,1)
    ifreq = 1
    nz = nz-nfreq
    itstp=nfreq+1
 
    !     DEFINE BASIC SIZES
 
    mcb(1) = uhvf
    CALL rdtrl(mcb)
    nload = mcb(2)/nfreq
    nmodes = mcb(3)
    mcb(1) = uhvt
    mcb(2) = 0
    mcb(5) = 1
    mcb(6) = 0
    mcb(7) = 0
    k = nfreq + 2*(nmodes*nfreq*2) + nmodes
    IF(k > nz) CALL  mesage(-8,0,NAME)
 
    !     DETERMINE IF EQUAL FREQ - CONVERT TO W'S
 
    delw = z(ifreq+1)-z(ifreq)
    epsi = delw*1.e-6
    j = nfreq-1
    iequal = 1
    DO  i=1,j
        m = ifreq+i-1
        IF( ABS(z(m+1)-z(m)-delw) >= epsi) iequal = 0
        z(m) = z(m)*twopi
    END DO
    z(ifreq+nfreq-1) = z(ifreq+nfreq-1)*twopi
    delw = delw*twopi
 
    !     FIRST FREQUENCY MUST BE MULTIPLE OF DELW
 
    nbig = ABS(z(ifreq)/delw)+.1
    IF(ABS(FLOAT(nbig)*delw-ABS(z(ifreq))) > epsi) iequal = 0
    lll = nbig -1
 
    !     FIND TSTEP IN TRL
 
    CALL gopen(casecc,iz(ibuf1),0)
    CALL fread(casecc,0,-37,0)
    CALL fread(casecc,j,1,0)
    CALL CLOSE(casecc,1)
    FILE = trl
    CALL OPEN(*900,trl,iz(ibuf1),0)
    CALL fread(trl,mcb1,3,1)
    m = mcb1(3)
    CALL skprec(trl,m)
25 CONTINUE
   CALL fread(trl,m,1,0)
   IF(m == j) GO TO 30
   CALL fread(trl,0,0,1)
   GO TO 25
 
   !     FOUND TSTEP
 
30 CALL READ(*910,*40,trl,iz(itstp),nz,0,ngroup)
   CALL mesage(-8,0,NAME)
40 nz = nz-ngroup
   iuhvf = itstp+ngroup
   CALL CLOSE(trl,1)
   ngroup = ngroup/3
   IF(ngroup /= 1) iequal = 0
   IF( iequal == 0 ) GO TO 50
 
   !     FORCE WAT TO BE INTEGER MULTIPLE OF TWOPI/N
 
   fbig = twopi/(delw*z(itstp+1))
   nbig = fbig+.9
   z(itstp+1) = twopi/(FLOAT(nbig)*delw)
50 CONTINUE
 
   !     BUILD / WRITE TOL
 
   FILE = tol
   CALL OPEN(*900,tol,iz(ibuf1),1)
   CALL fname(tol,mcb1)
   CALL WRITE(tol,mcb1,2,0)
   delt = z(itstp+1)
   t = 0.0
   n = 0
   m = itstp
   DO  i=1,ngroup
       nstep = iz(m)
       IF(i == 1) nstep = nstep +1
       m = m+3
       DO  j=1,nstep
           CALL WRITE(tol,t,1,0)
           n = n+1
           IF(j == nstep  .AND. i /= ngroup) delt = z(m+1)
           t = t+delt
       END DO
   END DO
   CALL WRITE(tol,0,0,1)
   CALL CLOSE(tol,1)
   mcb1(1) = tol
   mcb1(2) = ngroup
   mcb1(3) = n
   mcb1(4) = 0
   mcb1(5) = 0
   mcb1(6) = 0
   mcb1(7) = 0
   CALL wrttrl(mcb1)
 
   !     BUILD TABLE OF CK, SK
 
   ick = iuhvf + 2*nmodes*nfreq
   isk = ick
   iudot = isk
   IF( iequal == 0 ) GO TO 100
   isk = ick + nbig
   iudot = isk + nbig
   m = ick
   m1 = isk
   m2 = isk
   j = iudot
   rp = COS(twopi/FLOAT(nbig))
   cp = SIN(twopi/FLOAT(nbig))
   i = m
   n = m1+1
   l = m2
   kk = j
   z(i) = 1.0
   z(l) = 0.0
65 IF(m1-i-2 < 0) THEN
       GO TO    61
   ELSE IF (m1-i-2 == 0) THEN
       GO TO    62
   ELSE
       GO TO    63
   END IF
62 cmnr = -1.
   cmnc = 0.
   GO TO 64
63 cmnr = rp*z(i) -cp*z(l)
   cmnc = cp*z(i) +rp*z(l)
64 i = i+1
   l = l+1
   m1 = m1-1
   kk = kk-1
   z(i) = cmnr
   z(l) = cmnc
   z(m1) = cmnr
   z(kk) = -cmnc
   GO TO 65
61 CONTINUE
   !     GET READY FOR OUTPUTS
 
100 CALL gopen(uhvf,iz(ibuf1),0)
   CALL gopen(uhvt,iz(ibuf2),1)
   it1 = 1
   it2 = 1
   ii = 1
   jj=nmodes
   incr = 1
   it3 = 3
   ii1 = 1
   jj1 = nmodes
   incr1 = 1
   iuvt = iudot
   IF(iftm == 2) iuvt = iuvt+2*nfreq*nmodes
   ASSIGN 235 TO ihop
 
   !     BEGIN LOOP ON LOADS
 
   DO  i=1,nload
   
       !     PUT UHVF INTO CORE
   
       DO  j=1,nfreq
           m = iuhvf+(j-1)*nmodes*2
           CALL unpack(*120,uhvf,z(m))
           CYCLE
120        CALL zeroc(z(m),2*nmodes)
       END DO
       IF(iftm /= 2) GO TO 150
       ASSIGN 236 TO ihop
   
       !     COMPUTE SPLINE FIT FOR U DOT
   
   
       !     COMPUTE A'S
   
       iap = iuvt + nmodes
       m =nfreq + iap - 1
       z(m) = 0.0
       l = nfreq-2
       IF(l <= 0) GO TO 126
       DO  j=1,l
           m = iap + nfreq -j-1
           n = ifreq + nfreq-j-1
           z(m) = (z(n) - z(n-1))/(2.*(z(n+1)-z(n-1))-(z(n+1)-z(n))*z(m+1))
       END DO
126 CONTINUE
   
    !     COMPUTE U DOT DOT
   
    DO  m1=1,nmodes
        m = iudot +(nfreq-1)*nmodes*2 +(m1-1)*2
        z(m) = 0.0
        z(m+1) = 0.0
     
        !     BEGIN BACKWARD PASS
     
        m2= iuhvf +(nfreq-1)*nmodes*2 +(m1-1)*2
        IF(l <= 0) CYCLE
        DO  j=1,l
            n2 = m
            m = m-nmodes*2
            n = ifreq + nfreq -j-1
            m2 = m2-nmodes*2
            kk = iap + nfreq-j
            ll = m2+2*nmodes
            rp = z(n+1) - z(n)
            cp = z(n) -z(n-1)
            n1 = m2-2*nmodes
            z(m) = (6. *((z(ll)-z(m2))/rp-(z(m2)-z(n1))/cp)-rp*z(kk)*z(n2)) /cp
            z(m+1) = (6.*((z(ll+1)-z(m2+1))/rp-(z(m2+1)-z(n1+1))/cp)-rp*z(kk)*  &
                z(n2+1))/cp
        END DO
    END DO
   
    !     BEGIN FORWARD PASS
   
    DO  m1=1,nmodes
        m = iudot +(m1-1)*2
        m2 = iuhvf +(m1-1)*2
        n1 = m2+2*nmodes
        ll = m+2*nmodes
        rp = z(ifreq+1) -z(ifreq)
        z(m) =(6.*(z(n1)-z(m2))/rp-rp*z(iap+1)*z(ll))/(6.*z(ifreq)+(rp)*  &
            (2.-z(iap+1)) )
        z(m+1) = 0.0
        DO  j=2,nfreq
            kk = iap+j-1
            m2 = m
            m = m+2*nmodes
            z(m) = z(kk)*(z(ll) - z(m2))
            z(m+1) = z(kk)*(z(ll+1)-z(m2+1))
            ll = ll + 2*nmodes
        END DO
    END DO
150 CONTINUE
    t = 0.0
    n = 0
    m = itstp
    delt = z(itstp+1)
   
    !     BEGIN LOOP ON TIMES
   
    DO  l=1,ngroup
        nstep = iz(m)
        IF(l == 1) nstep = nstep+1
        m = m+3
        DO  j=1,nstep
            tt = t
            CALL zeroc(z(iuvt),nmodes)
       
            !     BEGIN LOOP ON FREQUENCIES
       
            lx = lll
            DO  ll=1,nfreq
                lx = lx+1
                wn = z(ifreq+ll-1)
                IF(ll == 1) GO TO 191
                wnm1 = z(ifreq+ll-2)
191             IF(ll == nfreq) GO TO 192
                wnp1 = z(ifreq+ll)
192         CONTINUE
            IF(iequal == 0) GO TO 190
            kk = MOD(lx*n,nbig)
            ck = z(ick+kk)
            sk = z(isk+kk)
            GO TO 195
190         ck = COS(wn*tt)
            sk = SIN(wn*tt)
195     CONTINUE
         
        !     COMPUTE CMN, DMN
         
        IF(iftm /= 0) GO TO 220
         
        !     IFTM  =0
         
        cmnc = 0.0
        IF(ll == 1) GO TO 196
        IF(ll == nfreq) GO TO 197
        cmnr = (wnp1-wnm1)*.5
        GO TO 230
196 CONTINUE
    cmnr = wnp1-wn
    IF(wn == 0.0) cmnr = cmnr*.5
    GO TO 230
197 cmnr = wn -wnm1
    GO TO 230
         
!     IFTM = 1
         
220 CONTINUE
    IF(ll == 1) GO TO 221
    IF(ll > 2 .AND. iequal /= 0 .AND. ll /= nfreq) GO TO 223
    r1 = wn-wnm1
    CALL ifte2(-tt*r1,rp,cp)
    cmnr = r1*.5*rp
    cmnc = r1*.5*cp
    GO TO 222
221 cmnr = 0.
    cmnc = 0.
222 CONTINUE
    IF(ll == nfreq) GO TO 223
    r2 = wnp1-wn
    CALL ifte2(tt*r2,rp,cp)
    cmnr = cmnr+r2*.5*rp
    cmnc = cmnc+r2*.5*cp
223 IF(iftm == 2) GO TO 229
    dmnr = 0.0
    dmnc = 0.0
    GO TO 230
229 CONTINUE
         
    !     IFTM = 2
         
    im2 = iudot -2 +(ll-1)*nmodes*2
    IF(ll == 1) GO TO 224
    IF(ll > 2 .AND.  iequal /= 0 .AND. ll /= nfreq) GO TO 230
    CALL iftg(-tt*r1,rp,cp)
    r1 = -r1*r1*r1/24.
    dmnr = r1*rp
    dmnc = r1*cp
    GO TO 228
224 CONTINUE
    dmnr = 0.0
    dmnc = 0.0
228 CONTINUE
    IF(ll == nfreq) GO TO 230
    CALL iftg(tt*r2,rp,cp)
    r2 = -r2*r2*r2/24.
    dmnr = dmnr+r2*rp
    dmnc = dmnc+r2*cp
230 CONTINUE
    im1 = iuhvf-2 +(ll-1)*nmodes*2
         
    !     BEGIN LOOP ON MODES
         
    DO  kk=1,nmodes
        im = im1+2*kk
        rp = cmnr*z(im)-cmnc*z(im+1)
        cp = cmnc*z(im)+cmnr* z(im+1)
        GO TO ihop,(235,236)
236 CONTINUE
    im = im2+2*kk
    rp = rp+dmnr*z(im)-dmnc*z(im+1)
    cp = cp+dmnc*z(im)+dmnr*z(im+1)
235 CONTINUE
    z(iuvt+kk-1) = z(iuvt+kk-1) + rp*ck-cp*sk
           
!     END LOOP ON MODES
           
END DO
         
!     END LOOP ON FREQUENCIES
         
END DO
DO  kk=1,nmodes
    z(iuvt+kk-1) = z(iuvt+kk-1)/phi
END DO
CALL pack(z(iuvt),uhvt,mcb)
DO  kk =1,2
    CALL bldpk(1,1,uhvt,0,0)
    CALL bldpkn(uhvt,0,mcb)
END DO
IF(j == nstep) delt = z(m+1)
t = t + delt
n = n+1
END DO
     
!     END LOOP ON TIME
     
END DO
   
!     END LOOP ON LOADS
   
END DO
CALL CLOSE(uhvf,1)
CALL CLOSE(uhvt,1)
CALL wrttrl(mcb)
RETURN
 
!     ERROR MESSAGES
 
900 n1=-1
901 CALL mesage(n1,FILE,NAME)
CALL pexit
910 n1=-2
GO TO 901
END SUBROUTINE ift
