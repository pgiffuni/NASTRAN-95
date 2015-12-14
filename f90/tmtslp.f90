SUBROUTINE tmtslp
     
!     TMTSLP TIME TESTS CPU TIMES FOR VARIOUS TYPES OF LOOPS
 
!     COMMENT FROM G.CHAN/UNISYS   5/91
!     BASICALLY THIS ROUTINE IS SAME AS TIMTS2
 
!     IF ALL TIMING CONSTANTS ARE ZEROS (OR 0.001) SYSTEM HAS A WRONG
!     CPUTIM.MDS SUBROUTINE. MOST LIKELY THE CPUTIM.MIS IS BEING USED.
 
 INTEGER :: sysbuf,buf1,buf2,END,end2,end4,TYPE,isubr(2)
 REAL :: b(1),c(1),d(1),e(16)
 COMPLEX :: ac(1),bc(1),cc(1),dc(1),adnc
 DOUBLE PRECISION :: adnd,ad(1),bd(1),cd(1),dd(1)
 COMMON /machin/  mach
 COMMON /ntime /  nitems,tgino ,tbldpk,tintpk,tpack , tunpak,tgetst,tputst,  &
     ttlrsp,ttlrdp,ttlcsp,ttlcdp, tllrsp,tllrdp,tllcsp,tllcdp,tgetsb
 COMMON /system/  sysbuf,nout,skip(74),isy77
 COMMON /zzzzzz/  a(1)
 EQUIVALENCE      (a(1),ac(1),ad(1),b(1),bc(1),bd(1),c(1),cc(1),  &
     cd(1),d(1),dc(1),dd(1)),    (e(1),tgino)
 DATA    isubr /  4HTMTS, 4HLP  /
 
!     INITIALIZE
!     DOUBLE N SIZE SINCE VAX (AND UNIX) CLOCK MAY NOT TICK FAST ENOUGH
 
 n = 50
 IF (mach >= 5) n = 100
 m = n
 
 buf1 = korsz(a) - sysbuf
 buf2 = buf1 - sysbuf
END  = n*m
IF (END >= buf1-1) CALL mesage (-8,0,isubr)

!     CPU TIME TESTS

asq  = m + n
adno = 1/(asq*asq)
adnd = adno
adnc = CMPLX(adno,adno)
end2 = END/2
end4 = END/4
DO  TYPE = 1,4
  SELECT CASE ( TYPE )
    CASE (    1)
      GO TO 10
    CASE (    2)
      GO TO 90
    CASE (    3)
      GO TO 170
    CASE (    4)
      GO TO 250
  END SELECT
  
!     REAL CPU TIME TESTS
  
  10 CONTINUE
  
IF (m > END .OR. n > END) CALL mesage (-8,0,isubr)
DO  i = 1,END
a(i) = adno
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    d(j) = a(j)*b(j) + c(j)
  END DO
END DO
CALL cputim (t2,t2,1)
ASSIGN 340 TO iret
GO TO 330
50 CONTINUE

DO  i = 1,END
a(i) = adno
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    l = i + j - 1
    d(j) = a(i)*b(l) + c(j)
  END DO
END DO
CALL cputim (t2,t2,1)
ASSIGN 350 TO iret
GO TO 330

!     DOUBLE PRECISION TESTS

90 CONTINUE

IF (m > end2 .OR. n > end2) CALL mesage (-8,0,isubr)
DO  i = 1,end2
  ad(i) = adnd
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    dd(j) = ad(j)*bd(j) + cd(j)
  END DO
END DO
CALL cputim (t2,t2,1)
ASSIGN 360 TO iret
GO TO 330
130 CONTINUE

DO  i = 1,end2
  ad(i) = adnd
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    l = i + j - 1
    dd(j) = ad(i)*bd(l) + cd(j)
  END DO
END DO
CALL cputim (t2,t2,1)
ASSIGN 370 TO iret
GO TO 330

!     COMPLEX SINGLE PRECISION TESTS

170 CONTINUE

IF (m > end2 .OR. n > end2) CALL mesage (-8,0,isubr)
DO  i = 1,end2
  ac(i) = adnc
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    dc(j) = ac(j)*bc(j) + cc(j)
  END DO
END DO
CALL cputim (t2,t2,1)
ASSIGN 380 TO iret
GO TO 330
210 CONTINUE

DO  i = 1,end2
  ac(i) = adnc
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    l = i + j - 1
    dc(j) = ac(i)*bc(l) + cc(j)
  END DO
END DO
CALL cputim (t2,t2,1)
ASSIGN 390 TO iret
GO TO 330

!     DOUBLE PRECISION COMPLEX TESTS

250 CONTINUE

IF (m > end4 .OR. n > end4) CALL mesage (-8,0,isubr)
DO  i = 1,end2
  ad(i) = adnd
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    
!     D(J) AND D(J+1) CALCULATIONS WERE REVERSED
!     IN ORDER TO COUNTERACT THE ITERATIVE BUILD UP
    
    dd(j+1) = ad(j)*bd(j  ) - ad(j+1)*bd(j+1) + cd(j  )
    dd(j  ) = ad(j)*bd(j+1) + ad(j+1)*bd(j  ) + cd(j+1)
  END DO
END DO
CALL cputim (t2,t2,1)
ASSIGN 400 TO iret
GO TO 330
290 CONTINUE

DO  i = 1,end2
  ad(i) = adnd
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    l = i + j - 1
    dd(j  ) = ad(i)*bd(l  ) - ad(i+1)*bd(l+1) + cd(j  )
    dd(j+1) = ad(i)*bd(l+1) + ad(i+1)*bd(l  ) + cd(j+1)
  END DO
END DO
CALL cputim (t2,t2,1)
ASSIGN 410 TO iret


!     INTERNAL ROUTINE TO STORE TIMING DATA IN /NTIME/ COMMON BLOCK

330 time = t2 - t1
itot = m*n
tperop = 1.0E6*time/itot
GO TO iret, (340,350,360,370,380,390,400,410)
340 ttlrsp = tperop
GO TO 50
350 tllrsp = tperop
CYCLE
360 ttlrdp = tperop
GO TO 130
370 tllrdp = tperop
CYCLE
380 ttlcsp = tperop
GO TO 210
390 tllcsp = tperop
CYCLE
400 ttlcdp = tperop
GO TO 290
410 tllcdp = tperop
END DO

!     MAKE SURE ALL TIME CONTSTANTS ARE OK

DO  i = 1,nitems
  IF (isy77 == -3 .AND. e(i) < 0.001) e(i) = 0.001
  IF (isy77 /= -3 .AND. e(i) < 1.e-7) e(i) = 1.e-7
END DO
IF (isy77 /= -3) GO TO 460
WRITE  (nout,440) nitems,nitems,e
440 FORMAT ('0*** NASTRAN SYSTEM MESSAGE. IF THESE',i4,' NEW TIMING',  &
    ' CONSTANTS ARE HARD-CODED INTO THE LABEL COMMON /NTIME/ OF',  &
    /5X, 'SUBROUTINE SEMDBD, COMPILE, AND RE-LINKE LINK 1, THE ',  &
    'COMPUTATIONS OF THESE CONSTANTS IN ALL NASTRAN JOBS WILL',/5X,  &
    'BE ELIMINATED.',  /5X,'OR TO ACCOMPLISH THE SAME RESULT, ',  &
    'EDIT THE TIM-LINE IN THE NASINFO FILE TO INCLUDE THESE',i4,  &
    ' NEW',/5X,'TIMING CONSTANTS', //5X,9F8.3, /5X,7F8.3,//)
CALL pexit
460 CALL sswtch (35,j)
IF (j /= 0) CALL tmtsot

RETURN
END SUBROUTINE tmtslp
