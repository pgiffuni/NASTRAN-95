SUBROUTINE invtr (*,x,dx)
     
!     INVTR WILL INVERT A LOWER OR UPPER TRIANGULAR MATRIX
 
 
!     FILEA    =  MATRIX CONTROL BLOCK FOR THE INPUT FILE A
!     FILEB    =  MATRIX CONTROL BLOCK FOR THE OUTPUT MATRIX B
!     SCRFIL   =  SCRATCH FILE (NEEDED ONLY FOR AN UPPER TRIANGLE)
!     NX       =  NUMBER OF CELLS OF CORE AVAILABLE AT X
!     PREC     =  DESIRED PRECISION OF ARITHMETIC OPERATIONS
!     X        =  BLOCK OF AVAILABLE CORE
!     DX       =  SAME BLOCK AS X, BUT TYPED DOUBLE PRECISION
 
 REAL, INTENT(OUT)                        :: x(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dx(1)
 INTEGER :: rd,rdrew,wrt,wrtrew,rew,norew,eofnrw,rc,tra2,  &
     core,eol,outbuf,typea,typeb,prec,typear,  &
     forma,sysbuf,bakskp,forskp,CMPLX,tra,tra1,scrfil, filea,fileb,t
 DOUBLE PRECISION :: da(2),dtemp
 DIMENSION  NAME(2),t(7)
 COMMON /TYPE  / prc(2),nwds(4),rc(10)
 COMMON /system/ sysbuf
!     COMMON /DESCRP/ LENGTH,MAJOR(1)
 COMMON /zntpkx/ a(4),ii,eol
 COMMON /invtrx/ filea(7),fileb(7),scrfil,nx,prec
 COMMON /names / rd,rdrew,wrt,wrtrew,rew,norew,eofnrw
 COMMON /packx / it1,it2,iy,jy ,incry
 COMMON /unpakx/ itx1,ix,jx ,incrx
 EQUIVALENCE     (a(1),da(1)),(filea(3),nrow),(filea(4),forma),  &
     (filea(5),typea),(fileb(5),typeb)
 DATA    NAME  / 4HINVT,4HR     /, t /7*0/
 
!     INITIALIZE
 
 typear = rc(typea)
 typear = rc(typear) + prec - 1
 incr   = nwds(typear)
 it1    = typear
 it2    = typeb
 itx1   = typear
 incrx  = 1
 incry  = 1
 fileb(2) = 0
 fileb(6) = 0
 fileb(7) = 0
 iobuf = nx - sysbuf
 CMPLX = rc(typear)
 core  = iobuf - 1
 CALL gopen (fileb,x(iobuf),1)
 CALL CLOSE (fileb,norew)
 IF (forma == 5) GO TO 500
 IF (forma /= 4) GO TO 1000
 
!     INVERT A LOWER TRIANGULAR MATRIX
 
 bakskp = nrow
 forskp = 1
 SELECT CASE ( typear )
   CASE (    1)
     GO TO 1
   CASE (    2)
     GO TO 2
   CASE (    3)
     GO TO 3
   CASE (    4)
     GO TO 4
 END SELECT
 1 ASSIGN 50  TO tra
 ASSIGN 110 TO tra1
 GO TO  5
 2 ASSIGN 60  TO tra
 ASSIGN 120 TO tra1
 GO TO  5
 3 ASSIGN 70  TO tra
 ASSIGN 130 TO tra1
 GO TO  5
 4 ASSIGN 80  TO tra
 ASSIGN 140 TO tra1
 5 CONTINUE
 
!     ALLOCATE CORE STORAGE
 
 CALL gopen (filea,x(iobuf),0)
 j = 1
 
!     SOLVE QUADRATIC FOR K
 
 10 m  = nrow - j + 1
 l  = 2*m + 1
 k  = m
 IF (l*l <= 8/incr*core) GO TO 20
 a1 = l*l - 8/incr*core
 k  = SQRT(a1)
 k  = k + 1
 k  = (l-k)/2
 IF (k <= 0) GO TO 1040
 
!     GENERATE COLUMNS J THROUGH J+K OF THE IDENTITY MATRIX (STORE
!     ONLY THE LOWER TRIANGLE IN CORE)
 
 20 l = (m*k-(k*(k-1))/2)*incr
 DO  i = 1,l
   x(i) = 0.
 END DO
 l = 1
 IF (prec == 2) GO TO 41
 DO  i = 1,k
   x(l) = 1.
   l = l + (m-i+1)*incr
 END DO
 GO TO 44
 41 DO  i = 1,k
   dx(l) = 1.d0
   l = l + (m-i+1)*CMPLX
 END DO
 44 CONTINUE
 
!     READ MATRIX A ONE ELEMENT AT A TIME, ADDING IN TERMS TO THE
!     IDENTITY MATRIX
 
 l  = 1
 ll = 1
 
!     II =  COLUMN INDEX
!     M  =  HEIGTH OF TRAPAZOID
!     K  =  LENGTH OF TRAPAZOID
 
 DO  i = j,nrow
   CALL intpk (*1050,filea,0,typear,0)
   45 IF (eol == 0.0) THEN
     GO TO    46
   ELSE
     GO TO  1050
   END IF
   46 CALL zntpki
   IF (i /= ii) GO TO 45
   l1 = 0
   DO  i1 = 1,ll
     in1 = (l-1)*CMPLX + 1 + l1
     GO TO tra, (50,60,70,80)
     50 x(in1) = x(in1)/a(1)
     GO TO 90
     60 dx(in1) = dx(in1)/da(1)
     GO TO 90
     70 temp     = (a(1)*x(in1  ) + a(2)*x(in1+1))/(a(1)*a(1) + a(2)*a(2))
     x(in1+1) = (a(1)*x(in1+1) - a(2)*x(in1  ))/(a(1)*a(1) + a(2)*a(2))
     x(in1)   = temp
     GO TO 90
     80 dtemp    = (da(1)*dx(in1  ) +da(2)*dx(in1+1))/(da(1)**2 +da(2)**2)
     dx(in1+1)= (da(1)*dx(in1+1) -da(2)*dx(in1  ))/(da(1)**2 +da(2)**2)
     dx(in1)  = dtemp
     90 l1 = l1 + (m-i1)*CMPLX
   END DO
   100 IF (eol == 1) GO TO 190
   CALL zntpki
   l1 = 0
   DO  i1 = 1,ll
     in2 = (l-1)*CMPLX + 1+ l1
     in1 = in2 + (ii-i)*CMPLX
     GO TO tra1, (110,120,130,140)
     110 x(in1) = x(in1) - a(1)*x(in2)
     GO TO 150
     120 dx(in1) = dx(in1) - da(1)*dx(in2)
     GO TO 150
     130 x(in1  ) = x(in1  ) - a(1)*x(in2  ) + a(2)*x(in2+1)
     x(in1+1) = x(in1+1) - a(1)*x(in2+1) - a(2)*x(in2  )
     GO TO 150
     140 dx(in1  ) = dx(in1  ) - da(1)*dx(in2  ) + da(2)*dx(in2+1)
     dx(in1+1) = dx(in1+1) - da(1)*dx(in2+1) - da(2)*dx(in2  )
     150 l1 = l1 + (m-i1)*CMPLX
   END DO
   GO TO 100
   190 ll = ll + 1
   IF (ll > k) ll = k
   l = l + 1
 END DO
 forskp = forskp + k
 bakskp = bakskp - k
 i1 = rew
 IF (bakskp < forskp) i1 = norew
 CALL CLOSE (filea ,i1)
 CALL gopen (fileb,x(iobuf),wrt)
 l  = 1
 iy = j
 jy = nrow
 DO  i = 1,k
   CALL pack (x(l),fileb,fileb)
   iy = iy + 1
   l  = l  + (m-i+1)*incr
 END DO
 CALL CLOSE (fileb,norew)
 j  = j + k
 IF (j <= nrow) GO TO 206
 CALL gopen (fileb,x(iobuf),wrt)
 CALL CLOSE (fileb,rew)
 RETURN
 
 206 CONTINUE
 CALL gopen (filea,x(iobuf),rd)
 IF (forskp > bakskp) GO TO 220
 CALL skprec (filea,forskp)
 GO TO 10
 220 CALL skprec (filea,-bakskp)
 GO TO 10
 
!     INVERT UPPER TRIANGULAR MATRIX
 
 500 SELECT CASE ( typear )
   CASE (    1)
     GO TO 510
   CASE (    2)
     GO TO 520
   CASE (    3)
     GO TO 530
   CASE (    4)
     GO TO 540
 END SELECT
 510 ASSIGN 600 TO tra
 ASSIGN 700 TO tra1
 ASSIGN 770 TO tra2
 GO TO  550
 520 ASSIGN 610 TO tra
 ASSIGN 710 TO tra1
 ASSIGN 780 TO tra2
 GO TO  550
 530 ASSIGN 610 TO tra
 ASSIGN 720 TO tra1
 ASSIGN 790 TO tra2
 GO TO  550
 540 ASSIGN 630 TO tra
 ASSIGN 730 TO tra1
 ASSIGN 800 TO tra2
 
!     REWRITE UPPER TRIANGULAR MATRIX ON SCRATCH FILE
 
 550 inbuf  = iobuf
 forskp = nrow + 1
 bakskp = 0
 outbuf = inbuf - sysbuf
 IF (outbuf < nrow+1) GO TO 1040
 CALL gopen (filea,x(iobuf),0)
 
!     POSITION FILE AT LAST RECORD
 
 CALL skprec (filea,nrow)
 
!     REWRITE THE INPUT MATRIX ON A SCRATCH FILE WITH THE RECORDS
!     WRITTEN IN THE REVERSE ORDER AND THE COLUMNS INVERTED
 
 CALL gopen (scrfil,x(outbuf),1)
 it2 = typear
 DO  i = 1,nrow
   ix = 1
   jx = 0
   CALL bckrec (filea)
   CALL unpack (*1050,filea,x)
   CALL bckrec (filea)
   kk = jx - ix + 1
   k  = kk/2
   IF (k == 0) GO TO 641
   kk = kk + 1
   DO  j = 1,k
     l  = kk - j
     GO TO tra, (600,610,630)
     600 temp = x(j)
     x(j) = x(l)
     x(l) = temp
     CYCLE
     610 dtemp = dx(j)
     dx(j) = dx(l)
     dx(l) = dtemp
     CYCLE
     630 dtemp = dx(j)
     dx(j) = dx(l)
     dx(l) = dtemp
     dtemp = dx(j+1)
     dx(j+1) = dx(l+1)
     dx(l+1) = dtemp
   END DO
   641 CONTINUE
   iy = nrow - jx + 1
   jy = nrow - ix + 1
   CALL pack (x,scrfil,t)
 END DO
 it1 = typear
 it2 = typeb
 CALL CLOSE (filea,rew)
 CALL CLOSE (scrfil,eofnrw)
 CALL gopen (scrfil,x(iobuf),0)
 CALL skprec (scrfil,nrow)
 CALL CLOSE (scrfil,norew)
 
!     ALLOCATE CORE
 
 j  = 0
 650 m  = j + 1
 CALL gopen (scrfil,x(iobuf),rd)
 k  = nrow - j
 IF (k*m+k*(k-1)/2 < core/incr) GO TO 652
 a1 = (2*m-1)**2 + 8*core/incr
 k  = SQRT(a1)
 k  = (-(2*m-1)+k)/2
 IF (k <= 0) GO TO 1040
 652 bakskp = bakskp + k
 forskp = forskp - k
 
!     POSITION SCRATCH FILE
 
 IF (forskp > bakskp) GO TO 660
 CALL REWIND (scrfil)
 CALL skprec (scrfil,forskp)
 GO TO 665
 660 CALL skprec (scrfil,-bakskp)
 665 CONTINUE
 
!     GENERATE UPPER TRIANGLE OF THE IDENTITY MATRIX
 
 lend = (m*k+k*(k-1)/2)*incr
 DO  i = 1,lend
   x(i) = 0.
 END DO
 l = m
 IF (prec == 2) GO TO 676
 DO  i = 1,k
   x(l) = 1.
   l = l + (i+m)*incr
 END DO
 GO TO 680
 676 DO  i = 1,k
   dx(l) = 1.d0
   l = l + (i+m)*CMPLX
 END DO
 680 CONTINUE
 
!     READ UPPER TRIANGLE ONE ELEMENT AT A TIME, ADDING IN
!     APPROPIATE TERMS TO THE IDENTITY MATRIX
 
 IF (prec == 2) lend = lend/2
 j = j + k
 l = 1
 DO  jj = 1,j
   CALL intpk (*1050,scrfil,0,typear,0)
   CALL zntpki
   i = nrow - ii + 1
   IF (i /= j-jj+1) GO TO 1050
   l1 = 0
   DO  i1 = 1,l
     in1 = lend - l*CMPLX - l1 + 1
     GO TO tra1, (700,710,720,730)
     700 x(in1) = x(in1)/a(1)
     GO TO 740
     710 dx(in1) = dx(in1)/da(1)
     GO TO 740
     720 temp   = (a(1)*x(in1  ) + a(2)*x(in1+1))/(a(1)*a(1) + a(2)*a(2))
     in2    = in1 + 1
     x(in2) = (a(1)*x(in1+1) - a(2)*x(in1  ))/(a(1)*a(1) + a(2)*a(2))
     x(in1) = temp
     GO TO 740
     730 dtemp  = (da(1)*dx(in1  ) + da(2)*dx(in1+1))/(da(1)**2 + da(2)**2)
     in2    = in1 + 1
     dx(in2)= (da(1)*dx(in1+1) - da(2)*dx(in1  ))/(da(1)**2 + da(2)**2)
     dx(in1)= dtemp
     740 CONTINUE
     l1 = l1 + (m+k-1-i1)*CMPLX
   END DO
   760 IF (eol == 1) GO TO 901
   CALL zntpki
   l1 = 0
   i  = j - jj - nrow + ii
   DO  i1 = 1,l
     in2 = lend - l*CMPLX - l1 + 1
     in1 = in2 - i*CMPLX
     GO TO tra2, (770,780,790,800)
     770 x(in1) = x(in1) - a(1)*x(in2)
     GO TO 810
     780 dx(in1) = dx(in1) - da(1)*dx(in2)
     GO TO 810
     790 x(in1  ) = x(in1  ) - a(1)*x(in2  ) + a(2)*x(in2+1)
     x(in1+1) = x(in1+1) - a(1)*x(in2+1) - a(2)*x(in2  )
     GO TO 810
     800 dx(in1  ) = dx(in1  ) - da(1)*dx(in2  ) + da(2)*dx(in2+1)
     dx(in1+1) = dx(in1+1) - da(1)*dx(in2+1) - da(2)*dx(in2  )
     810 l1 = l1 + (m+k-1-i1)*CMPLX
   END DO
   GO TO 760
   901 l = l + 1
 END DO
 CALL CLOSE (scrfil,norew)
 CALL gopen (fileb,x(iobuf),wrt)
 l = j - k + 1
 ll= 1
 DO  i = 1,k
   iy = 1
   jy = l
   CALL pack (x(ll),fileb,fileb)
   l = l + 1
   ll = ll + (m+i-1)*incr
 END DO
 CALL CLOSE (fileb,norew)
 IF (j < nrow) GO TO 650
 CALL gopen (fileb,x(iobuf),wrt)
 CALL CLOSE (fileb,rew)
 CALL gopen (scrfil,x(iobuf),rd)
 CALL CLOSE (scrfil,rew)
 GO TO 2000
 
 1000 no = -7
 GO TO 1100
 1040 no = -8
 GO TO 1100
 1050 RETURN 1
 1100 CALL mesage (no,0,NAME)
 
 2000 RETURN
END SUBROUTINE invtr
