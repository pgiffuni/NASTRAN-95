SUBROUTINE ddamat
     
!     DDAMAT  A,B/C/C,Y,GG=1.  $
 
!     DDAMAT  TAKES THE OUTER PRODUCT OF MATRICES A AND B, AND MULTIPLES
!     BY GG TO GET C, I.E.  CIJ=GG*(AIJ*BIJ).  ALSO, IF B HAS ONLY ONE
!     COLUMN, AND NUMBER OF COLUMNS OF A .GT. 1, THEN USE THAT COLUMN
!     ON EACH COLUMN OF A.
 
 INTEGER :: a,b,c,buf1,buf2,buf3
 DOUBLE PRECISION :: dz(1), dgg, ggdz
 DIMENSION        nam(2),mcb(7)
 COMMON /unpakx/  jout,iii,nnn,jncr
 COMMON /packx /  iin,iout,ii,nn,incr
 COMMON /system/  ibuf(80)
 COMMON /BLANK /  gg
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (iprec,ibuf(55)),(z(1),dz(1))
 DATA    a,b,c /  101,102,201  /
 DATA    nam   /  4HDDAM,4HAT  /
 
!     SET PACK AND UNPACK PARAMETER
 
 jout = iprec
 iin  = iprec
 iout = iprec
 incr = 1
 jncr = 1
 ii   = 1
 iii  = 1
 
!     SET OPEN CORE
 
 lcore = korsz(z)
 buf1  = lcore - ibuf(1) + 1
 buf2  = buf1 - ibuf(1)
 buf3  = buf2 - ibuf(1)
 lcore = buf3 - 1
 IF (lcore <= 0) GO TO 1008
 
 mcb(1) = a
 CALL rdtrl (mcb)
 ncola  = mcb(2)
 nrowa  = mcb(3)
 mcb(1) = b
 CALL rdtrl (mcb)
 ncolb  = mcb(2)
 nrowb  = mcb(3)
 IF (nrowa /= nrowb) GO TO 1007
 IF (lcore < 2*nrowa*iprec) GO TO 1008
 
!     NO. OF COLUMNS OF A AND B MUST BE EQUAL OR
!     NO. OF COLUMNS OF B MUST BE 1
 
 IF (ncola == ncolb) GO TO 5
 IF (ncolb ==     1) GO TO 5
 GO TO 1007
 
 5 nn  = nrowa
 nnn = nrowa
 mcb(1) = c
 mcb(2) = 0
 mcb(3) = nrowa
 mcb(6) = 0
 mcb(7) = 0
 IF (iprec == 2) dgg = gg
 
 CALL gopen (a,z(buf1),0)
 CALL gopen (b,z(buf2),0)
 CALL gopen (c,z(buf3),1)
 
!     UNPACK A COLUMN OF A AND B, COMPUTE PRODUCTS, AND PACK TO C.
!     IF I.GT.1 AND B=1, USE THE ONE COLUMN OF B OVER AGAIN.
 
 DO  i = 1,ncola
   
   inull = 0
   SELECT CASE ( iprec )
     CASE (    1)
       GO TO 10
     CASE (    2)
       GO TO 40
   END SELECT
   10 ggz = gg
   CALL unpack (*11,a,z(1))
   GO TO 12
   11 inull = 1
   12 IF (i > 1 .AND. ncolb == 1) GO TO 20
   CALL unpack (*15,b,z(nrowa+1))
   GO TO 20
   15 inull = 1
   DO  j = 1,nrowa
     z(nrowa+j) = 0.
   END DO
   20 IF (inull == 1) ggz = 0.
   DO  j = 1,nrowa
     z(j) = ggz*z(j)*z(nrowa+j)
   END DO
   CALL pack (z(1),c,mcb)
   CYCLE
   40 ggdz = dgg
   CALL unpack (*41,a,dz(1))
   GO TO 42
   41 inull = 1
   42 IF (i > 1 .AND. ncolb == 1) GO TO 50
   CALL unpack (*45,b,dz(nrowa+1))
   GO TO 50
   45 inull = 1
   DO  j = 1,nrowa
     dz(nrowa+j) = 0.d0
   END DO
   50 IF (inull == 1) ggdz = 0.d0
   DO  j = 1,nrowa
     isub  = nrowa + j
     dz(j) = ggdz*dz(j)*dz(isub)
   END DO
   CALL pack (dz(1),c,mcb)
   
!     DO ANOTHER COLUMN
   
 END DO
 
 CALL wrttrl (mcb)
 CALL CLOSE (a,1)
 CALL CLOSE (b,1)
 CALL CLOSE (c,1)
 RETURN
 
!     FATAL ERROR MESSAGE
 
 1007 k = -7
 GO TO 1010
 1008 k = -8
 1010 CALL mesage (k,0,nam)
 RETURN
END SUBROUTINE ddamat
