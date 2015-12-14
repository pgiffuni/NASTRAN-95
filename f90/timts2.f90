SUBROUTINE timts2
     
!     TIMTS2 TIME TESTS CPU TIMES FOR VARIOUS TYPES OF LOOPS
 
 
 INTEGER :: sysbuf, output, buf1, buf2, END, end2, end4,  &
     p, opt1, opt2, TYPE, NAME(4), tig(16), med(16), los(16), isubr(2)
 REAL :: b(1), c(1), d(1)
 DOUBLE PRECISION :: adnd, ad(1), bd(1), cd(1), dd(1)
 COMPLEX :: ac(1), bc(1), cc(1), dc(1), adnc
 COMMON /BLANK / n, m, TYPE, opt1, opt2
 COMMON /system/ ksystm(65)
 COMMON /zzzzzz/ a(1)
 EQUIVALENCE     (ksystm(1),sysbuf), (ksystm(2),output),  &
     (a(1),ac(1),ad(1), b(1),bc(1),bd(1), c(1),cc(1),cd(1), d(1),dc(1),dd(1))
 DATA    tig   / 1H ,4HTIGH, 4HT( r,4HSP ) , 1H ,4HTIGH, 4HT( r,4HDP ) ,  &
     1H ,4HTIGH, 4HT( c,4HSP ) , 1H ,4HTIGH, 4HT( c,4HDP ) /
 DATA    med   / 1H ,4HMEDI, 4HUM(r,4HSP ) , 1H ,4HMEDI, 4HUM(r,4HDP ) ,  &
     1H ,4HMEDI, 4HUM(c,4HSP ) , 1H ,4HMEDI, 4HUM(c,4HDP ) /
 DATA    los   / 1H ,4HLOOS, 4HE (r,4HSP ) , 1H ,4HLOOS, 4HE (r,4HDP ) ,  &
     1H ,4HLOOS, 4HE (c,4HSP ) , 1H ,4HLOOS, 4HE (c,4HDP ) /
 DATA    isubr / 4HTIMT, 4HS2  /, m8/-8/
 
!     INITIALIZE
 
 CALL page1
 WRITE  (output,11) n,m,TYPE,opt1
 11 FORMAT (1H  , 20X, 25HNASTRAN time test d   n =, i4, 5H, m =, i4 ,  &
     8H, TYPE =,i4, 8H, opt1 =,i4)
 buf1 = korsz(a) - sysbuf
 buf2 = buf1 - sysbuf
END = n*m
IF (END >= buf1-1) CALL mesage (m8,0,isubr)

!     CPU TIME TESTS

p = 4*(TYPE-1) + 1
asq  = m + n
adno = 1/(asq*asq)
adnd = adno
adnc = CMPLX(adno,adno)
end2 = END/2
end4 = END/4
SELECT CASE ( TYPE )
  CASE (    1)
    GO TO 105
  CASE (    2)
    GO TO 106
  CASE (    3)
    GO TO 107
  CASE (    4)
    GO TO 108
END SELECT

!     REAL CPU TIME TESTS

105 CONTINUE

IF (m > END .OR. n > END) CALL mesage (m8,0,isubr)
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
iret = 1
NAME(1) = tig(p  )
NAME(2) = tig(p+1)
NAME(3) = tig(p+2)
NAME(4) = tig(p+3)
GO TO 500
501 CONTINUE

DO  i = 1,END
a(i) = adno
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    d(j) = a(i)*b(j) + c(j)
  END DO
END DO
CALL cputim (t2,t2,1)
iret = 2
NAME(1) = med(p  )
NAME(2) = med(p+1)
NAME(3) = med(p+2)
NAME(4) = med(p+3)
GO TO 500
502 CONTINUE

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
iret = 3
NAME(1) = los(p  )
NAME(2) = los(p+1)
NAME(3) = los(p+2)
NAME(4) = los(p+3)
GO TO 500

!     DOUBLE PRECISION TESTS

106 CONTINUE

IF (m > end2 .OR. n > end2) CALL mesage (m8,0,isubr)
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
iret = 4
NAME(1) = tig(p  )
NAME(2) = tig(p+1)
NAME(3) = tig(p+2)
NAME(4) = tig(p+3)
GO TO 500
504 CONTINUE

DO  i = 1,end2
  ad(i) = adnd
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    dd(j) = ad(i)*bd(j) + cd(j)
  END DO
END DO
CALL cputim (t2,t2,1)
iret = 5
NAME(1) = med(p  )
NAME(2) = med(p+1)
NAME(3) = med(p+2)
NAME(4) = med(p+3)
GO TO 500
505 CONTINUE

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
iret = 6
NAME(1) = los(p  )
NAME(2) = los(p+1)
NAME(3) = los(p+2)
NAME(4) = los(p+3)
GO TO 500

!     COMPLEX SINGLE PRECISION TESTS

107 CONTINUE

IF (m > end2 .OR. n > end2) CALL mesage (m8,0,isubr)
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
iret = 7
NAME(1) = tig(p  )
NAME(2) = tig(p+1)
NAME(3) = tig(p+2)
NAME(4) = tig(p+3)
GO TO 500
507 CONTINUE

DO  i = 1,end2
  ac(i) = adnc
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    dc(j) = ac(i)*bc(j) + cc(j)
  END DO
END DO
CALL cputim (t2,t2,1)
iret = 8
NAME(1) = med(p  )
NAME(2) = med(p+1)
NAME(3) = med(p+2)
NAME(4) = med(p+3)
GO TO 500
508 CONTINUE

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
iret = 9
NAME(1) = los(p  )
NAME(2) = los(p+1)
NAME(3) = los(p+2)
NAME(4) = los(p+3)
GO TO 500

!     DOUBLE PRECISION COMPLEX TESTS

108 CONTINUE

IF (m > end4 .OR. n > end4) CALL mesage (m8,0,isubr)
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
iret = 10
NAME(1) = tig(p  )
NAME(2) = tig(p+1)
NAME(3) = tig(p+2)
NAME(4) = tig(p+3)
GO TO 500
510 CONTINUE

DO  i = 1,end2
  ad(i) = adnd
END DO
CALL cputim (t1,t1,1)
DO  i = 1,n
  DO  j = 1,m
    dd(j  ) = ad(i)*bd(j  ) - ad(i+1)*bd(j+1) + cd(j  )
    dd(j+1) = ad(i)*bd(j+1) + ad(i+1)*bd(j  ) + cd(j+1)
  END DO
END DO
CALL cputim (t2,t2,1)
iret = 11
NAME(1) = med(p  )
NAME(2) = med(p+1)
NAME(3) = med(p+2)
NAME(4) = med(p+3)
GO TO 500
511 CONTINUE

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
iret = 12
NAME(1) = los(p  )
NAME(2) = los(p+1)
NAME(3) = los(p+2)
NAME(4) = los(p+3)
GO TO 500
600 CONTINUE
RETURN


!     INTERNAL ROUTINE TO WRITE OUTPUT ONTO THE OUTPUT FILE

500 time = t2 - t1
itot = m*n
tperop = 1.0E6*time/itot
IF (iret == 2 .OR. iret == 5 .OR. iret == 8 .OR. iret == 11)  &
    WRITE (output,998) NAME,itot,time,tperop

IF (iret /= 2 .AND. iret /= 5 .AND. iret /= 8 .AND. iret /= 11)  &
    WRITE (output,999) NAME,itot,time,tperop

998 FORMAT (1H0, 4A4, ' CPU TIME FOR ', i9,  &
    ' OPERATIONS = ', e12.5, ' SECONDS'/  &
    1X , 16X, ' CPU TIME FOR ', '      ONE',  &
    ' OPERATION  = ', e12.5, ' MICROSECONDS')

999 FORMAT (1H0, 4A4, ' CPU TIME FOR ', i9,  &
    ' OPERATIONS = ', e12.5, ' SECONDS'/  &
    1X , 16X, ' CPU TIME FOR ', '      ONE',  &
    ' OPERATION  = ', e12.5, ' MICROSECONDS',  &
    '  ---  DATA FOR USE IN COMMON /NTIME/')

SELECT CASE ( iret )
  CASE (    1)
    GO TO 501
  CASE (    2)
    GO TO 502
  CASE (    3)
    GO TO 600
  CASE (    4)
    GO TO 504
  CASE (    5)
    GO TO 505
  CASE (    6)
    GO TO 600
  CASE (    7)
    GO TO 507
  CASE (    8)
    GO TO 508
  CASE (    9)
    GO TO 600
  CASE (   10)
    GO TO 510
  CASE (   11)
    GO TO 511
  CASE (   12)
    GO TO 600
END SELECT
END SUBROUTINE timts2
