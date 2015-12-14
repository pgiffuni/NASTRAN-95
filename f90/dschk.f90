SUBROUTINE dschk
     
!     MODULE TO PERFORM DIFFERENTIAL STIFFNESS CONVERGENCE TESTS
 
!     DSCHK    PGI,PGIP1,UGIP1//EPSIO,DSEPSI,NT,TOUT,TIN,DONE,SHIFT,
!                               COUNT,BETA
 
!     EPSIO    ACCEPTABLE RATIO OF ENERGY ERROR TO TOTAL ERROR(R) INPUT
!     DSEPSI   EPSI(SUB I -1)   (REAL)                            IN/OUT
!     NT       TOTAL NUMBER OF ITERATIONS ALLOWED                 INPUT
!     TOUT     START TIME FOR OUTER LOOP                          INPUT
!     TIN      START TIME FOR INNER LOOP                          INPUT
!     DONE     EXIT FLAG FOR SKIP TO SDR2                         OUTPUT
!     SHIFT    EXIT FLAG FOR SHIFT                                IN/OUT
!     COUNT    CURRENT STEP NUMBER                                IN/OUT
!     BETA     SHIFT DECISION FACTOR (REAL)                       INPUT
 
!     EXIT FLAG VALUES (IEXIT)                                    LOCAL
!          0   NOT SET
!          1   CONVERGED
!          2   DIVERGED
!          3   INSUFFICIENT TIME
!          4   ITERATION LIMIT
!          5   ZERO EPSIO
!          6   ZERO EPSI
 
 INTEGER :: pgi,pgip1,ugip1,tout,tin,done,shift,count,iz(1),  &
     sysbuf,scr1,scr2,scr3,FILE,tnow,ti,TO,tleft,beta
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /system/ sysbuf,nout,ksystm(52),iprec
 COMMON /unpakx/ ita,ii,jj,incr
 COMMON /BLANK / epsio,dsepsi,nt,tout,tin,done,shift,count,beta
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (count,ni), (z(1),iz(1))
 DATA    pgi   , pgip1,ugip1,scr1,scr2,scr3 /  &
     101   , 102,  103,  301, 302, 303  /
 
!     INITIALIZE
 
 ibuf1 = korsz(iz) - sysbuf + 1
 iexit = 0
 ifrst = shift
 shift = 1
 nf    = 1
 CALL klock (tnow)
 ti    = tnow - tin
 CALL tmtogo (tleft)
 TO    = tnow - tleft
 
!     COMPUTE DSEPSI(I)
 
 CALL ssg2b (ugip1,pgi  ,0   ,scr1,1,iprec,1,scr3)
 CALL ssg2b (ugip1,pgip1,scr1,scr2,1,iprec,2,scr3)
 
 ii   = 1
 jj   = 1
 incr = 1
 ita  = 1
 FILE = scr2
 ASSIGN 10 TO iretn
 GO TO 300
 
!     GET DENOMINATOR
 
 10 epsi = value
 FILE = scr1
 ASSIGN 20 TO iretn
 GO TO 300
 20 IF (value == 0.0) GO TO 40
 epsi  = ABS(epsi/value)
 count = count + 1
 IF (ifrst ==  -1) GO TO 30
 IF (epsi  == 0.0) GO TO 210
 xlama = ABS(dsepsi/epsi)
 IF (xlama <= 1.0) GO TO 60
 30 dsepsi = epsi
 IF (epsi > epsio) GO TO 50
 
!     CONVERGED
 
 40 iexit = 1
 done  =-1
 GO TO 220
 
!     MAKE FIRST TEST
 
 50 IF (ifrst == -1) GO TO 80
 
!     NOT FIRST TIME
 
 IF (epsio <= 0.0) GO TO 200
 nf  = ALOG(epsi/epsio)/ALOG(xlama)
 CALL klock (tnow)
 CALL tmtogo (tleft)
 ti  = tnow - tin
 TO  = tnow - tout
 GO TO 70
 
!     DIVERGED
 
 60 iexit = 2
 done  =-1
 dsepsi = epsi
 GO TO 220
 
!     CONVERGENT
 
 70 IF (nf > nt-ni) GO TO 90
 IF (ti*nf > TO+beta*ti) GO TO 100
 80 IF (tleft >= 3*ti) GO TO 120
 
!     INSUFFICIENT TIME
 
 iexit = 3
 done  =-1
 GO TO 220
 90 IF (nt-ni-beta < 0) THEN
   GO TO    80
 END IF
 
!     SET SHIFT FLAG
 
 100 shift =-1
 IF (tleft < TO+beta*ti) GO TO 80
 
!     WRAP UP FOR SHIFT
 
 done  = nf
 iexit = 0
 GO TO 220
 
!     USER LIMIT ITERATION NUMBER EXPIRED
 
 110 CONTINUE
 iexit = 4
 done  =-1
 GO TO 220
 
!     WRAP UP FOR NO SHIFT
 
 120 CONTINUE
 IF (ni >= nt) GO TO 110
 shift = 1
 done  = nf
 GO TO 220
 
!     PARAMETER ERROR, EPSIO HAS NO VALUE
 
 200 iexit = 5
 GO TO 220
 
!     AFTER SSG2B, EPSI IS ZERO DUE TO THE FIRST VAULE FROM SCR2 IS ZERO
!     WHILE VALUE FROM SCR1 IS NOT ZERO
 
 210 iexit = 6
 
!     EXIT FROM MODULE
 
 220 CALL page2 (-9)
 WRITE  (nout,230) uim,iexit,count,done,shift,dsepsi
 230 FORMAT (a29,' 7019, MODULE DSCHK IS EXITING FOR REASON',i4, /5X,  &
     'ON ITERATION NUMBER',i7,1H.,  &
     /5X,'PARAMETER VALUES ARE AS FOLLOWS',/10X,'DONE   =',i10,  &
     /10X,'SHIFT  =',i10, /10X,'DSEPSI =',1P,e14.7)
 IF (iexit >= 5) WRITE (nout,240) epsio,epsi
 240 FORMAT ( 10X,'EPSIO  =',1P,e10.3, /10X,'EPSI   =',1P,e10.3)
 RETURN
 
!     INTERNAL ROUTINE TO OBTAIN VALUE FROM MATRIX
 
 300 CALL gopen (FILE,iz(ibuf1),0)
 CALL unpack (*310,FILE,value)
 GO TO 320
 310 value = 0.0
 320 CALL CLOSE (FILE,1)
 GO TO iretn, (10,20)
END SUBROUTINE dschk
