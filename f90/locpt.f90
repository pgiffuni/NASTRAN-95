SUBROUTINE locpt ( n,p, m,s,k,ks, eps, loc)
     
!     DETERMINES POSITION OF EACH OF N POINTS (P) RELATIVE TO SURFACE
!        BOUNDED BY M POINTS (S)
!        ALL POINTS IN THE SAME COORDINATE SYSTEM
!     KS IS THE (UNIT) VECTOR NORMAL TO THE SURFACE
!     LOC(I) IS FLAG INDICATING POSITION OF POINT I RELATIVE TO SURFACE:
!        LOC= 1    WHEN POINT WITHIN SURFACE BOUNDRY
!        LOC= 0    WHEN POINT IS ON SURFACE BOUNDRY
!        LOC= -1   WHEN POINT OUTSIDE SURFACE BOUNDRY
 
!     MAXIMUM OF 4 POINTS MAY BE LOCATED RELATIVE
!        TO SURFACE WITH MAXIMUM OF 4 SIDES
!        WHOSE ENDPOINTS ARE IN K(2,M)
 
 
 INTEGER, INTENT(IN)                      :: n
 DOUBLE PRECISION, INTENT(IN)             :: p(3,4)
 INTEGER, INTENT(IN)                      :: m
 DOUBLE PRECISION, INTENT(IN)             :: s(3,4)
 INTEGER, INTENT(IN)                      :: k(2,1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: ks(3)
 DOUBLE PRECISION, INTENT(IN)             :: eps(2)
 INTEGER, INTENT(OUT)                     :: loc(1)
 
 DOUBLE PRECISION :: ve(3,4), vp(3), v(3),vemag(4)
 DOUBLE PRECISION :: vmag, vpmag, vdotk, edotp
 DOUBLE PRECISION :: dvmag, dadotb
 
 
 
 
!     EPS ARRAY FOR SIGNIFICANCE TESTING
!        EPS(1) IS AREA, ANGLE LIMIT
!        EPS(2) IS LENGTH LIMIT
 
 
 
!     SET UP VECTORS ALONG EACH SURFACE EDGE
 
 DO  NE=1,m
   k1= k(1,NE)
   k2= k(2,NE)
   
   DO  i=1,2
!     VE IS VECTOR ALONG SURFACE EDGE
     ve(i,NE)= s(i,k2) -s(i,k1)
   END DO
   ve(3,NE)=0.d0
   vemag(NE)= dvmag( ve(1,NE), eps(2) )
 END DO
 
!      DETERMINE LOCATION OF POINT RELATIVE TO SURFACE
 
 DO  np= 1,n
!        (PRESET POINT FLAG TO INTERIOR CODE)
   loc(np)= 1
   
   DO  NE= 1,m
     k1= k(1,NE)
     
     DO   i=1,2
!     VP IS VECTOR FROM FIRST END OF EDGE VECTOR TO POINT
       vp(i)= p(i,np) -s(i,k1)
     END DO
     vp(3)= 0.d0
     vpmag= dvmag(vp,eps(2))
     
!     V= VE CROSS VP
     CALL daxb(ve(1,NE), vp,v)
     vmag= dvmag(v,eps(1))
!     VDOTK= (VE CROSS VP) DOT K,  K NORMAL TO PLANE OF SURFACE
     vdotk= dadotb(v, ks)
!     EDOTP IS VE DOT VP
     edotp= dadotb( ve(1,NE), vp)
     IF (vpmag <= eps(2) )    GO TO 37
     IF (vdotk > eps(1))           CYCLE
!                                      INSIDE THIS EDGE
     IF (         (vdotk < -eps(1) )         .OR.  &
         (edotp <= eps(1))           .OR.  &
         (vemag(NE) +eps(2) < vpmag) )        GO TO 35
!                                                          OUTSIDE
     GO TO 37
!                  ON THIS EDGE
   END DO
!     POINT IS WITHIN SURFACE BOUNDRY IF NOT OUTSIDE ANY EDGE
!        AND NOT ON SURFACE BOUNDRY
   CYCLE
   
!     POINT IS OUTSIDE SURFACE BOUNDRY IF OUTSIDE ANY EDGE
   35 loc(np)= -1
   CYCLE
   
   
!     POINT IS ON BOUNDRY WHEN ANGLE IS EFFECTIVELY ZERO
!        OR (EFFECTIVELY) COINCIDENT WITH EDGE POINT
   37 loc(np)= 0
   
 END DO
 RETURN
END SUBROUTINE locpt
