SUBROUTINE frd2c (a,b,x,scr1,scr2,scr3,scr4,scr5,nload,nfreq)
     
!     SOLVE A X = B
!     USE INCORE DECOMP IF POSSIBLE
 
 
 INTEGER, INTENT(IN)                      :: a
 INTEGER, INTENT(IN)                      :: b
 INTEGER, INTENT(IN OUT)                  :: x
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN OUT)                  :: scr2
 INTEGER, INTENT(IN OUT)                  :: scr3
 INTEGER, INTENT(IN OUT)                  :: scr4
 INTEGER, INTENT(IN OUT)                  :: scr5
 INTEGER, INTENT(IN)                      :: nload
 INTEGER, INTENT(IN OUT)                  :: nfreq
 INTEGER :: sysbuf,out,ta(7), tb(7),tx(7)
 DIMENSION       zz(1)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /system/ sysbuf,out,dum(52),iprec
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /unpakx/ iout,inn,nnn,incr1
 COMMON /frd2bc/ ih,ip
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (zz(1),z(1))
 
 icore= korsz(z)
 incr = 1
 ii   = 1
 inn  = 1
 incr1= 1
 iout = 3
 IF (ih == 0 .AND. iprec == 2) iout = 4
 
!     IH IN /FRD2BC/ IS INITIALIZED BY ROUTINE FRRD2.
!     (COMPLEX D.P. ARITHMETIC IS USED IF IH=0)
 
 ito = iout
 iti = ito
 
!     DECIDE IF INCORE IS POSSIBLE
 
 ta(1) = a
 CALL rdtrl (ta)
 tb(1) = b
 CALL rdtrl (tb)
 na    = ta(2)
 nb    = tb(3)*nload
 ibuf1 = icore - sysbuf
 ncore = na*na*2 + nb*2 + nb*2 + sysbuf
 
!     IF IH=0, COMPLEX D.P. COMPUTATION WILL BE USED.  NOTICE THAT THE
!     ROUTINE INCORE IS WRITTEN ONLY FOR COMPLEX S.P. OPERATION.
 
 IF (ih == 0) GO TO 102
 IF (ncore > icore) GO TO 100
 
!     DO INCORE
 
 ia = 1
 CALL gopen (a,z(ibuf1),0)
 nnn = ta(3)
 incr1 = nnn
 n = na + na
 DO  i = 1,n,2
   CALL unpack (*11,a,z(i))
   CYCLE
   11 DO  k = 1,n,2
     l = (k-1)*nnn
     z(i+l  ) = 0.0
     z(i+l+1) = 0.0
   END DO
 END DO
 CALL CLOSE (a,1)
 
!     GET FREQ FROM B
 
 ib   = nnn*nnn*2 + 1
 nnn  = tb(3)
 incr1= nload
 n1   = nnn + nnn
 j    = tb(2)/nload - 1
 m    = 0
 CALL gopen (b,z(ibuf1),0)
 CALL skprec (b,nfreq-1)
 DO  i = 1,nload
   CALL unpack (*31,b,z(ib+m))
   GO TO 33
   31 DO  k = 1,n1,2
     l = (k-1)*nload + ib + m
     z(l  ) = 0.0
     z(l+1) = 0.0
   END DO
   33 IF (i /= nload) CALL skprec (b,j)
   m = m+2
 END DO
 CALL CLOSE (b,1)
 ix = nload*nnn*2 + ib
 CALL incore (z(ia),na,z(ib),z(ix),nload)
 nn = na
 CALL gopen (x,z(ibuf1),1)
 CALL makmcb (tx,x,nn,tb(4),ito)
 incr = nload
 j = ix
 DO  i = 1,nload
   CALL pack (z(j),x,tx)
   j = j + 2
 END DO
 CALL CLOSE (x,1)
 CALL wrttrl (tx)
 GO TO 1000
 
!     USE FILE SOLVE
 
 100 IF (ip /= 0) GO TO 102
 ip = ncore - icore
 WRITE  (out,101) uim,ip
 101 FORMAT (a29,' 2437, ADDITIONAL CORE NEEDED FOR IN-CORE ',  &
     'DECOMPOSITION IN FRRD2 MODULE IS',i8,' WORDS.')
 102 CALL cfactr (a,scr1,scr2,scr3,scr4,scr5,iopt)
 icore = korsz(zz)
 ibuf1 = icore - sysbuf
 ibuf2 = ibuf1 - sysbuf
 CALL gopen (b,zz(ibuf1),0)
 CALL gopen (scr3,zz(ibuf2),1)
 iout = 3
 IF (ih == 0 .AND. iprec == 2) iout = 4
 incr1 = 1
 j = tb(2)/nload - 1
 nn = tb(3)
 CALL makmcb (tx,scr3,nn,tb(4),ito)
 CALL skprec (b,nfreq-1)
 DO  i = 1,nload
   CALL cyct2b (b,scr3,1,zz,tx)
   IF (i /= nload) CALL skprec (b,j)
 END DO
 CALL CLOSE (scr3,1)
 CALL CLOSE (b,1)
 CALL wrttrl (tx)
 CALL cfbsor (scr1,scr2,scr3,x,iopt)
 1000 RETURN
END SUBROUTINE frd2c
