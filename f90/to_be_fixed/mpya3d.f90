SUBROUTINE mpya3d (aa,bb,nrow,band,cc)
     
!     WITH ENTRY MPYA3S (A,B,NROW,BAND,C)
 
!     WAS NAMED DATBAD/DATBAS IN UAI CODE
 
!     THESE ROUTINES PERFORM TRIPLE MATRIX MULTIPLY OF THE FORM
 
!                          T
!                 C = C + A * B * A
 
!     ON TWO INCOMING ROW-LOADED MATRICES A AND B, AND ADD THEM TO
!     MATRIX C
 
!     THE INCOMING MATRICES MUST BE SQUARE (AND OBVIOUSLY OF THE SAME
!     SIZE, NROW.) AND
!     SYMMETRICAL (SINCE WE OPERATE ONLY ON LOWER TRIANGULAR MATRICES)
 
!     MATRIX A CAN BE A PSUEDO-DIAGONAL MATRIX, I.E. A MATRIX HAVING
!     SQUARE PARTITIONS OF NON-ZERO TERMS ALONG ITS DIAGONAL.
!     THESE PARTITIONS ARE OF THE SIZE  BAND X BAND.
!     NOTE THAT NROW MUST BE AN INTEGER MULTIPLE OF BAND.
 
!     THIS ALGORITHM IS SUITABLE FOR TRIPLE MULTIPLIES INVOLVING GLOBAL
!     TRANSFORMATIONS.
 
 
 
 DOUBLE PRECISION, INTENT(IN)             :: aa(1)
 DOUBLE PRECISION, INTENT(IN)             :: bb(1)
 INTEGER, INTENT(IN)                      :: nrow
 INTEGER, INTENT(IN)                      :: band
 DOUBLE PRECISION, INTENT(OUT)            :: cc(1)
 
 REAL :: a(1) ,b(1) ,c(1)
 DOUBLE PRECISION :: dd
 
 
!     DOUBLE PRECISION VERSION
 
 ii = 0
 loop50:  DO  ib = 1,nrow
   ia1 = ((ib-1)/band+1)*band
   
   DO  id1 = 1,nrow,band
     id2 = id1 + band - 1
     IF (id1 > ia1) GO TO 50
     
     id11n = (id1-1)*nrow
     DO  id = id1,id2
       jj = id11n
       dd = 0.0D0
       
       DO  ic = id1,id2
         ibic = ii + ic
         icid = jj + id
         IF (aa(icid) == 0.0D0) GO TO 10
         dd = dd + bb(ibic)*aa(icid)
         10 jj = jj + nrow
       END DO
       
       IF (dd == 0.0D0) CYCLE
       kk = (id-1)*nrow
       
       DO  ia = id,ia1
         ibia = ii + ia
         IF (aa(ibia) == 0.0D0) GO TO 20
         iaid = kk + id
         cc(iaid) = cc(iaid) + dd*aa(ibia)
         20 kk = kk + nrow
       END DO
       
     END DO
   END DO
   ii = ii + nrow
 END DO loop50
 
!     COPY THE LOWER TRIANGLE TO THE UPPER
 
 kk = nrow - 1
 ii = 0
 DO  i = 1,kk
   ib = i + 1
   jj = i*nrow
   DO  j = ib,nrow
     cc(ii+j) = cc(jj+i)
     jj = jj + nrow
   END DO
   ii = ii + nrow
 END DO
 
 RETURN
 
 
 ENTRY mpya3s (a,b,nrow,band,c)
!     ==============================
 
!     SINGLE PRECISION VERSION
 
 ii = 0
 loop150:  DO  ib = 1,nrow
   ia1 = ((ib-1)/band+1)*band
   
   DO  id1 = 1,nrow,band
     id2 = id1 + band - 1
     IF (id1 > ia1) GO TO 150
     
     id11n = (id1-1)*nrow
     DO  id = id1,id2
       jj = id11n
       dd = 0.0D0
       
       DO  ic = id1,id2
         ibic = ii + ic
         icid = jj + id
         IF (a(icid) == 0.0) GO TO 110
         dd = dd + DBLE(b(ibic))*DBLE(a(icid))
         110 jj = jj + nrow
       END DO
       IF (dd == 0.0D0) CYCLE
       kk = (id-1)*nrow
       
       DO  ia = id,ia1
         ibia = ii + ia
         IF (a(ibia) == 0.0) GO TO 120
         iaid = kk + id
         c(iaid) = SNGL(DBLE(c(iaid)) + dd*DBLE(a(ibia)))
         120 kk = kk + nrow
       END DO
       
     END DO
   END DO
   ii = ii + nrow
 END DO loop150
 
!     COPY THE LOWER TRIANGLE TO THE UPPER
 
 kk = nrow - 1
 ii = 0
 DO  i = 1,kk
   ib = i + 1
   jj = i*nrow
   DO  j = ib,nrow
     c(ii+j) = c(jj+i)
     jj = jj + nrow
   END DO
   ii = ii + nrow
 END DO
 
 RETURN
END SUBROUTINE mpya3d
