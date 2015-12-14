SUBROUTINE trht1c (ngroup,udvt,pd,rdd,iloop)
     
!     THIS ROUTINE  STEPS INTEGRATION PROCEDURE
 
 
 INTEGER, INTENT(IN)                      :: ngroup
 INTEGER, INTENT(IN)                      :: udvt
 INTEGER, INTENT(IN)                      :: pd
 INTEGER, INTENT(IN)                      :: rdd
 INTEGER, INTENT(IN)                      :: iloop
 INTEGER :: sysbuf, iz(1), a,        FILE,     mcb(7),   pnl1,     radlin,  &
     NAME(2),  ifn(7),   itab(4),  ll1(7)
 DOUBLE PRECISION :: dz(1)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,      uwm
 COMMON /BLANK /  beta,     tabs,     norad,    radlin,   sigma
 COMMON /system/  ksystm(63)
 COMMON /trdxx /  ktrdxx(28)
 COMMON /zzzzzz/  z(1)
 COMMON /packx /  it1,      it2,      ii,       jj,       incr
 COMMON /trhtx /  ik(7),    ib(7),    icr1,     icr2,     icr3,  &
     icr4,     icr5,     isym,     a,        icr7
 COMMON /trdd1 /  nlft1,    dit1,     nlftp1,   nout,     icount,  &
     iloop1,   moda1,    nz,       icore,    iu1,  &
     in2,      ipnl(7),  nmodes,   nstep,    pnl1,  &
     ist,      iu1dum,   deltat,   ifrst,    tabs1, sigma1,   tim1
 COMMON /unpakx/  it3,      iii,      jjj,      incr1
 COMMON /infbsx/  ill1(7),  iul1(7)
 COMMON /fbsx  /  ll1
 EQUIVALENCE      (ksystm(1),sysbuf), (ksystm(55),iprec),  &
     (ktrdxx(28),iopen), (z(1),iz(1),dz(1)), (ill1(3),mrow),     (ksystm(2),nprt)
 DATA    NAME  /  4HTRHT,  4H1C     /
 
!     SYMBOL TABLE
 
!     ICR1 IS LLL
!     ICR2 IS ULL
!     ICR5 IS INITIAL CONDITIONS
!     ICR6 IS THE  A  MATRIX
 
!     NROW     PROBLEM ORDER
!     NGROUP   NUMBER OF TRIPLES OF TIME STEPS
!     UDVT     DISPLACEMENTS AND VELOCITIES
!     PD       LOADS
!     RDD      RADIATION MATRIX
!     ILOOP    CURRENT TIME STEP GROUP
!     IBUF1    UDVT BUFFER
!     IBUF2    A    BUFFER
!     IBUF3    LLL  BUFFER
!     IBUF4    ULL  BUFFER
!     IBUF5    PD   BUFFER
!     IBUF6    PNL1 BUFFER
!     IBUF7    RDD  BUFFER
!     IBUF8    SCRATCH BUFFER(DIT,NLLOADS,SAVE STUFF ETC)
!     NZ       OPEN CORE
!     IST      OUTPUT FLAG
!     IU1,IU2  DISPLACMENT VECTOR POINTERS
!     IP1,IP2  LOAD VECTOR POINTERS
!     IN1,IN2  NON-LINEAR LOAD POINTERS
!     NOLIN    =0  MEAN NO NON-LINEAR LOADS
!     IPNT     POINTER FOR INTERNAL ZERO ROUTINE
!     FILE     FILE    FOR INTERNAL ZERO ROUTINE
!     NSTEP    NUMBER OF TIME STEPS
!     DELTAT   DELTA  T
!     NOUT     OUTPUT INCREMENT
!     H        1/ 2*DELTAT
!     ICOUNT   STEP COUNTER
!     ITLEFT   TIME LEFT
!     NORAD    =-1  NO RADIATION
!     RADLIN   =-1  NON LINEAR RADIATION
!     NLFTP1   NONLINEAR SET SELECTED BY THE USER
!     BETA,OMBETA,OPBETA  --USER BETA 1-BETA, 1+BETA
!     ISYM     0    UNSYMETRIC   1  SYMMETRIC
!     DELTA1   OLD DELTA  T
 
 iscr5  = icr5
 noload = 0
 nbust  = 0
 mcb(1) = pd
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) noload = -1
 nrow   = ik(2)
 it1    = 1
 it2    = 1
 ii     = 1
 jj     = nrow
 incr   = 1
 it3    = 1
 iii    = 1
 jjj    = nrow
 incr1  = 1
 tabs1  = tabs
 sigma1 = sigma
 nz     = korsz(z)
 igroup = nz - 3*ngroup + 1
 ibuf1  = igroup- sysbuf
 ibuf2  = ibuf1 - sysbuf
 ibuf3  = ibuf2 - sysbuf
 ibuf4  = ibuf3 - sysbuf
 ibuf5  = ibuf4 - sysbuf
 ibuf6  = ibuf5 - sysbuf
 ibuf7  = ibuf6 - sysbuf
 ibuf8  = ibuf7 - sysbuf
 nz     = ibuf8 - 1
 iloop1 = iloop
 ist    = 0
 ill1(1)= icr1
 CALL rdtrl (ill1)
 ifn(1) = icr1
 CALL rdtrl (ifn)
 iu1    = 0
 iu2    = iu1 + nrow
 ip1    = iu2 + nrow
 ip2    = ip1 + nrow
 iuk    = ip2 + nrow
 nolin  = 0
 IF (nlftp1 /= 0 .OR. norad /= -1) nolin = 1
 IF (nolin == 0) GO TO 10
 in1    = iuk + nrow
 in2    = in1 + nrow
 nz     = nz - 7*nrow
 GO TO 20
 
!     NO NON-LINEAR EFFECTS
 
 10 nz     = nz - 4*nrow
 in2    = ip2
 20 IF (nz < 0) CALL mesage (-8,0,NAME)
 icore  = in2 + nrow
 iul1(1)= icr2
 CALL rdtrl (iul1)
 ombeta = 1.0 - beta
 opbeta = 1.0 + beta
 
!     SET UP FOR CORE I/O
 
 IF (nlftp1 == 0) GO TO 21
 ifrst  = 0
 CALL trd1d
 ifrst  = 1
 21 itab(1) = a
 itab(2) = ill1(1)
 itab(3) = iul1(1)
 itab(4) = rdd
 icor = in2 + nrow + 1
 nf   = 4
 CALL gopen (a,iz(ibuf2),0)
 CALL REWIND (a)
 IF (nolin == 0 .OR. radlin /= -1 .OR. norad == -1) GO TO 30
 CALL gopen (rdd,iz(ibuf7),0)
 CALL REWIND (rdd)
 30 CONTINUE
 CALL gopen (ill1,iz(ibuf3),0)
 CALL REWIND (ill1)
 IF (isym == 1) GO TO 31
 CALL gopen (iul1,iz(ibuf4),0)
 CALL REWIND (iul1)
 31 CONTINUE
 
!     IS  THIS  A TIME  STEP CHANGE
 
 IF (iloop  /= 1) GO TO  280
 IF (noload /= 0) GO TO 33
 CALL gopen (pd,iz(ibuf5),0)
 CALL fwdrec (*440,pd)
 33 CONTINUE
 ist = -1
 CALL gopen (icr5,iz(ibuf1),0)
 
 CALL fread(icr5,iz(igroup),3*ngroup,1)
 
!     BRING IN  U0 AND UK
 
 CALL READ (*450,*35,icr5,z(iu1+1),nrow,1,nwds)
 GO TO 40
 
!     SHORT VECTOR ENCOUNTERED
 
 35 k = nwds + 1
 DO  l = k,nrow
   m = iu1 + l
   z(m) = 0.0
 END DO
 40 CONTINUE
 IF (norad == -1) GO TO 50
 CALL READ (*450,*45,icr5,z(iuk+1),nrow,1,nwds)
 GO TO 410
 
!     SHORT VECTOR ENCOUNTERED
 
 45 k = nwds + 1
 DO  l = k,nrow
   m = iuk + l
   z(m) = 0.0
 END DO
 GO TO 410
 50 CONTINUE
 CALL CLOSE (icr5,1)
 nstep  = iz(igroup) + 1
 deltat =  z(igroup+1)
 nout   = iz(igroup+2)
 h = 1.0/deltat
 CALL gopen (udvt,iz(ibuf1),1)
 CALL makmcb (mcb,udvt,nrow,2,1)
 IF (nolin == 0) GO TO 60
 CALL gopen (pnl1,iz(ibuf6),1)
 CALL makmcb (ipnl,pnl1,nrow,2,1)
 
!     LETS  GO
 
 60 icount = 1
 
!     TOP OF LOOP
 
 70 CALL tmtogo (itleft)
 IF (itleft <= 0) GO TO 230
 
!     COMPUTE  NR
 
 IF (norad  == -1) GO TO 110
 IF (radlin == -1) GO TO 90
 DO  i = 1,nrow
   l = in2 + i
   k = iuk + i
   z(l) = z(k)
 END DO
 GO TO 130
 
!     NON-CONSTANT RADIATION
 
 90 DO  i = 1,nrow
   l = iu1 + i
   k = iuk + i
   m = in2 + i
   j = iu2 + i
   
!     CHECK FOR UNSTABLE SOLUTION ABOUT TO CAUSE ARITHMETIC OVERFLOWS.
   
   IF (z(l) < 1.0E8) GO TO 98
   nbust = nbust + 1
   IF (nbust > 10) GO TO 94
   WRITE  (nprt,92) uwm,z(l),icount,i
   92 FORMAT (a25,' 3102, SUBROUTINE TRHT1C, UNSTABLE TEMP. VALUE OF',  &
       e20.8,' COMPUTED FOR TIME STEP',i5, /5X,  &
       'AT POINT NUMBER',i6,' IN THE ANALYSIS SET.')
   z(l) = 1.0E6
   GO TO 98
   94 WRITE  (nprt,96) ufm
   96 FORMAT (a23,' 3103, SUBROUTINE TRHT1C TERMINATING DUE TO ERROR ',  &
       'COUNT FOR MESSAGE 3102.')
   CALL mesage (-61,0,NAME)
   
   98 z(j) = -(z(l)+tabs)**4 + 4.0*(z(k)+tabs)**3*z(l)
   z(m) = 0.0
 END DO
 iopen  = 1
 ifn(1) = rdd
 CALL matvec (z(iu2+1),z(in2+1),ifn,iz(ibuf7))
 GO TO 130
 110 IF (nlftp1 == 0) GO TO 140
 DO  i = 1,nrow
   m = in2 + i
   z(m) = 0.0
 END DO
 130 IF (nlftp1 == 0) GO TO 140
 tim1 = tim
 CALL trd1d
 140 IF (icount /= 1 .OR. iloop /= 1) GO TO 160
 DO  i = 1,nrow
   k = ip1 + i
   z(k) = 0.0
   IF (nolin == 0) CYCLE
   l = in2 + i
   m = in1 + i
   z(m) =  z(l)
   z(k) = -z(l)
 END DO
 iopen = 0
 CALL matvec (z(iu1+1),z(ip1+1),ik,z(ibuf8))
 
!     BRING IN  NEXT P
 
 160 IF (noload /= 0) GO TO 165
 CALL unpack (*165,pd,z(ip2+1))
 GO TO 170
 165 DO  i = 1,nrow
   k = ip2 + i
   z(k) = 0.0
 END DO
 
!     ADD ALL LOAD CONTRIBUTIONS
 
 170 CONTINUE
 DO  i = 1,nrow
   l = ip1 + i
   m = ip2 + i
   z(l) = ombeta*z(l) + beta*z(m)
   IF (nolin == 0) CYCLE
   m = in1 + i
   j = in2 + i
   z(l) = z(l) + opbeta*z(j) - beta*z(m)
 END DO
 
!     MULTIPLY  IN  A MATRIX
 
 iopen  = 1
 ifn(1) = a
 CALL matvec (z(iu1+1),z(ip1+1),ifn,iz(ibuf2))
 
!     SOLVE  FOR NEXT DISPLACEMENT
 
 iopen = 1
 IF (isym == 0) CALL intfbs (z(ip1+1),z(iu2+1),iz(ibuf4))
 IF (isym /= 1) GO TO 188
 
!     ABSORBED SUBROUTINE FBSINT   SEE ALSO EQUIV.   DATA.
 
 DO  i = 1,mrow
   z(i+iu2) = z(i+ip1)
 END DO
 
!     FORWARD PASS
 
 CALL REWIND (ill1)
 CALL fwdrec (*186,ill1)
 iz(ibuf4) = ill1(1)
 ll1(1) = ill1(1)
 CALL rdtrl (ll1)
 IF (iprec /= 1) GO TO 184
 CALL fbs1 (iz(ibuf4),z(iu2+1),z(iu2+1),mrow)
 GO TO 188
 184 CALL fbs21(iz(ibuf4),z(iu2+1),z(iu2+1),mrow)
 GO TO 188
 186 CALL mesage (-2,ill1,NAME)
 
!     ABSORBED SUBROUTINE FBSINT    SEE ALSO EQUIV.   DATA.
 
 188 CONTINUE
 IF (icount == 1 .OR. icount == nstep .OR.  &
     MOD(icount+ist,nout) == 0) GO TO 200
 
!     ROTATE POINTERS
 
 190 j   = ip1
 ip1 = ip2
 ip2 = j
 j   = iu1
 iu1 = iu2
 iu2 = j
 j   = in1
 in1 = in2
 in2 = j
 tim = tim + deltat
 icount = icount + 1
 IF (icount-nstep < 0) THEN
   GO TO    70
 ELSE IF (icount-nstep == 0) THEN
   GO TO   220
 ELSE
   GO TO   230
 END IF
 
!     IT  IS OUTPUT TIME
 
 200 CALL pack (z(iu1+1),udvt,mcb)
 
!     COMPUTE  U DOT
 
 DO  i = 1,nrow
   l = ip1 + i
   m = iu1 + i
   j = iu2 + i
   z(l) = (z(j)-z(m))*h
 END DO
 CALL pack (z(ip1+1),udvt,mcb)
 
!     PUT OUT ZERO ACCERERATION VECTOR FOR LATER MODULES
 
 CALL bldpk (1,1,udvt,0,0)
 CALL bldpkn (udvt,0,mcb)
 IF (nolin == 0) GO TO 190
 CALL pack (z(in2+1),pnl1,ipnl)
 GO TO 190
 
!     END OF 1 GROUP
 
 220 IF (iloop /= ngroup) GO TO 260
 GO TO 70
 230 j = 1
 240 CALL CLOSE (udvt,j)
 CALL CLOSE (pd,j)
 CALL CLOSE (ill1,1)
 CALL CLOSE (iul1,1)
 CALL CLOSE (a,1)
 CALL wrttrl (mcb)
 IF (norad == -1) GO TO 245
 CALL CLOSE (rdd,1)
 245 IF (nolin == 0) GO TO 250
 CALL CLOSE (pnl1,j)
 CALL wrttrl (ipnl)
 250 CONTINUE
 RETURN
 
!     MORE GROUPS TO COME  SAVE STUFF
 
 260 j = 2
 CALL gopen (iscr5,iz(ibuf8),1)
 CALL WRITE (iscr5,iz(igroup),3*ngroup,1)
 IF (nolin /= 0) CALL WRITE (iscr5,iz(iuk+1),nrow,1)
 
!     SAVE   UI -1
 
 CALL WRITE (iscr5,z(iu2+1),nrow,1)
 
!     SAVE   UI
 
 CALL WRITE (iscr5,z(iu1+1),nrow,1)
 IF (nolin == 0) GO TO 270
 
!     SAVE    NI - 1
 
 CALL WRITE (iscr5,z(in2+1),nrow,1)
 
!     SAVE    NI
 
 CALL WRITE (iscr5,z(in1+1),nrow,1)
 270 CONTINUE
 CALL CLOSE (iscr5,1)
 GO TO 240
 
!     REENTRY FROM CHANGE OF TIME STEP
 
 280 CONTINUE
 CALL gopen (iscr5,iz(ibuf8),0)
 CALL fread (iscr5,iz(igroup),3*ngroup,1)
 newgrp = igroup + (iloop-1)*3
 delta1 =  z(newgrp-2)
 nstep  = iz(newgrp)
 deltat = z(newgrp+1)
 nout   = iz(newgrp+2)
 CALL gopen (pd,iz(ibuf5),2)
 h = 1.0/deltat
 CALL gopen (udvt,iz(ibuf1),3)
 mcb(1) = udvt
 CALL rdtrl (mcb(1))
 IF (nolin == 0) GO TO 290
 CALL gopen (pnl1,iz(ibuf6),3)
 ipnl(1) = pnl1
 CALL rdtrl (ipnl)
 290 CONTINUE
 
!     RESTORE  STUFF  SAVED
 
 IF (nolin /= 0) CALL fread (iscr5,z(iuk+1),nrow,1)
 CALL fread (iscr5,z(iu2+1),nrow,1)
 CALL fread (iscr5,z(iu1+1),nrow,1)
 IF (nolin == 0) GO TO 300
 CALL fread (iscr5,z(in1+1),nrow,1)
 CALL fread (iscr5,z(in2+1),nrow,1)
 300 CONTINUE
 CALL CLOSE (iscr5,1)
 
!     COMPUTE  PBAR
 
 DO  i = 1,nrow
   l = ip1 + i
   z(l) = 0.0
   IF (nolin == 0) CYCLE
   m    = in2 + i
   z(l) =-z(m)
 END DO
 iopen = 0
 CALL matvec (z(iu1+1),z(ip1+1),ik,iz(ibuf8))
 IF (ib(1) == 0) GO TO 330
 DO   i = 1,nrow
   l = iu2 + i
   m = iu1 + i
   z(l) = (z(m)-z(l))/delta1
 END DO
 iopen = 0
 CALL matvec (z(iu2+1),z(ip1+1),ib,iz(ibuf8))
 330 CONTINUE
 IF (nolin == 0) GO TO 350
 h1 = 1.0 - deltat/delta1
 h2 = deltat/delta1
 DO  i = 1,nrow
   l  = in1 + i
   m  = in2 + i
   z(l) = h2*z(l) + h1*z(m)
 END DO
 350 icount = 0
 GO TO 70
 
!     CONSTANT RADIATION
 
 410 IF(radlin == -1) GO TO 50
 DO  i = 1,nrow
   l = iuk + i
   k = in2 +i
   z(l) = -(z(l)+tabs)**4 + 4.0*(z(l)+tabs)**3*z(l)
   z(k) = 0.0
 END DO
 iopen = 1
 ifn(1) = rdd
 CALL matvec (z(iuk+1),z(in2+1),ifn, iz(ibuf7))
 DO  i = 1,nrow
   l = iuk + i
   m = in2 + i
   z(l) = z(m)
 END DO
 GO TO 50
 
!     I/O ERROR
 
 440 FILE = pd
 GO TO 460
 450 FILE = icr5
 460 CALL mesage (-2,FILE,NAME)
 RETURN
END SUBROUTINE trht1c
