SUBROUTINE exfort (rw,u,f,buf,nwds,prec,dbuf)
!*****
 
!         *** IBM 360/370, VAX/780 VERSION ***
 
!     EXFORT PERFORMS FORTRAN FORMATTED IO FOR MODULE EXIO
 
!*****
 
 INTEGER, INTENT(IN OUT)                  :: rw
 INTEGER, INTENT(IN OUT)                  :: u
 INTEGER, INTENT(IN OUT)                  :: f
 INTEGER, INTENT(IN OUT)                  :: buf(nwds)
 INTEGER, INTENT(IN)                      :: nwds
 INTEGER, INTENT(IN OUT)                  :: prec
 DOUBLE PRECISION, INTENT(OUT)            :: dbuf(1)
 INTEGER :: fp,FMT,frmt(10)
 
 COMMON /BLANK /  x1(26),lbuf
 COMMON /exio2p/  nf,fp(5,1)
 COMMON /exio2f/  FMT(1)
!     COMMON /EXIO2X/  ==> /ZZEXO2/ UNIVAC ONLY
 
 DATA    leof  / 4H&eof/
 
 IF (nwds <= 0) RETURN
 IF (f    <= 0) GO TO 8
 ifmt = fp(1,f)-1
 DO  i = 1,10
   frmt(i) = FMT(ifmt+i)
 END DO
 8 SELECT CASE ( rw )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 20
   CASE (    3)
     GO TO 80
   CASE (    4)
     GO TO 150
 END SELECT
 10 SELECT CASE ( prec )
   CASE (    1)
     GO TO 30
   CASE (    2)
     GO TO 50
 END SELECT
 20 SELECT CASE ( prec )
   CASE (    1)
     GO TO 40
   CASE (    2)
     GO TO 60
 END SELECT
 
!     READ -- SINGLE PRECISION
 
 30 READ (u,frmt,ERR=35) buf
 35 IF (buf(1) == leof) GO TO 70
 RETURN
 
!     WRITE -- SINGLE PRECISION
 
 40 WRITE (u,frmt,ERR=45) buf
 45 RETURN
 
!     READ -- DOUBLE PRECISION
 
 50 n = nwds/3
 READ (u,frmt) (buf(4*i-3),dbuf(2*i),i=1,n)
 RETURN
 
!     WRITE -- DOUBLE PRECISION
 
 60 n = nwds/3
 WRITE (u,frmt) (buf(4*i-3),dbuf(2*i),i=1,n)
 RETURN
 
!     END OF FILE
 
 70 buf(3) = -1
 RETURN
 
!     POSITION THE FILE
 
 80 SELECT CASE ( nwds )
   CASE (    1)
     GO TO 90
   CASE (    2)
     GO TO 100
   CASE (    3)
     GO TO 100
 END SELECT
 90 REWIND u
 RETURN
 100 n = lbuf/33+1
 DO  i = 1,n
   BACKSPACE u
 END DO
 120 READ (u,160) n
 IF (n /= leof) GO TO 120
 BACKSPACE u
 RETURN
 
!     WRITE LOGICAL EOF
 
 150 n = lbuf/33
 DO  i = 1,n
   WRITE  (u,160) leof
   160 FORMAT (a4,128X)
 END DO
 RETURN
END SUBROUTINE exfort
