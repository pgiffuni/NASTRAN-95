SUBROUTINE ateig(m,a,rr,ri,iana,ia,b,rra,rri)
     
!     ..................................................................
 
!        SUBROUTINE ATEIG
 
!        PURPOSE
!           COMPUTE THE EIGENVALUES OF A REAL ALMOST TRIANGULAR MATRIX
 
!        USAGE
!           CALL ATEIG(M,A,RR,RI,IANA,IA)
 
!        DESCRIPTION OF THE PARAMETERS
!           M      ORDER OF THE MATRIX
!           A      THE INPUT MATRIX, M BY M
!           RR     VECTOR CONTAINING THE REAL PARTS OF THE EIGENVALUES
!                  ON RETURN
!           RI     VECTOR CONTAINING THE IMAGINARY PARTS OF THE EIGEN-
!                  VALUES ON RETURN
!           IANA   VECTOR WHOSE DIMENSION MUST BE GREATER THAN OR EQUAL
!                  TO M, CONTAINING ON RETURN INDICATIONS ABOUT THE WAY
!                  THE EIGENVALUES APPEARED (SEE MATH. DESCRIPTION)
!           IA     SIZE OF THE FIRST DIMENSION ASSIGNED TO THE ARRAY A
!                  IN THE CALLING PROGRAM WHEN THE MATRIX IS IN DOUBLE
!                  SUBSCRIPTED DATA STORAGE MODE.
!                  IA=M WHEN THE MATRIX IS IN SSP VECTOR STORAGE MODE.
 
!        REMARKS
!           THE ORIGINAL MATRIX IS DESTROYED
!           THE DIMENSION OF RR AND RI MUST BE GREATER OR EQUAL TO M
 
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!           NONE
 
!        METHOD
!           QR DOUBLE ITERATION
 
!        REFERENCES
!           J.G.F. FRANCIS - THE QR TRANSFORMATION---THE COMPUTER
!           JOURNAL, VOL. 4, NO. 3, OCTOBER 1961, VOL. 4, NO. 4, JANUARY
!           1962.  J. H. WILKINSON - THE ALGEBRAIC EIGENVALUE PROBLEM -
!           CLARENDON PRESS, OXFORD, 1965.
 
!     ..................................................................
 
 
 INTEGER, INTENT(IN)                      :: m
 DOUBLE PRECISION, INTENT(IN OUT)         :: a(1)
 DOUBLE PRECISION, INTENT(OUT)            :: rr(1)
 DOUBLE PRECISION, INTENT(OUT)            :: ri(1)
 INTEGER, INTENT(OUT)                     :: iana(1)
 INTEGER, INTENT(IN)                      :: ia
 REAL, INTENT(OUT)                        :: b(1)
 REAL, INTENT(OUT)                        :: rra(1)
 REAL, INTENT(OUT)                        :: rri(1)
 
 
 DOUBLE PRECISION :: prr(2),pri(2)
 DOUBLE PRECISION :: alpha,cap,d,delta,eps,eta,e10,e6,e7,g1,g2,g3
 DOUBLE PRECISION :: pan,pan1,psi1,psi2,r,rmod,s,t,u,v
 INTEGER :: p,p1,q
 
 e7=1.0E-8
 e6=1.0E-6
 e10=1.0E-10
 delta=0.5
 maxit=30
 
!        INITIALIZATION
 
 n=m
 20 n1=n-1
 in=n1*ia
 nn=in+n
 IF(n1 == 0) THEN
   GO TO  1300
 END IF
 30 np=n+1
 
!        ITERATION COUNTER
 
 it=0
 
!        ROOTS OF THE 2ND ORDER MAIN SUBMATRIX AT THE PREVIOUS
!        ITERATION
 
 DO  i=1,2
   prr(i)=0.0
   pri(i)=0.0
 END DO
 
!        LAST TWO SUBDIAGONAL ELEMENTS AT THE PREVIOUS ITERATION
 
 pan=0.0
 pan1=0.0
 
!        ORIGIN SHIFT
 
 r=0.0
 s=0.0
 
!        ROOTS OF THE LOWER MAIN 2 BY 2 SUBMATRIX
 
 n2=n1-1
 in1=in-ia
 nn1=in1+n
 n1n=in+n1
 n1n1=in1+n1
 60 t=a(n1n1)-a(nn)
 u=t*t
 v=4.0D0*a(n1n)*a(nn1)
 IF (DABS(v)-u*e7 > 0.0) THEN
   GO TO    65
 ELSE
   GO TO   100
 END IF
 65 t=u+v
 IF (DABS(t)-DMAX1(u,DABS(v))*e6 > 0.0) THEN
   GO TO    68
 END IF
 67 t=0.0
 68    u=(a(n1n1)+a(nn))/2.0D0
 v=DSQRT(DABS(t))/2.0D0
 IF(t < 0.0) THEN
   GO TO   140
 END IF
 70 IF(u < 0.0) THEN
   GO TO    80
 END IF
 75 rr(n1)=u+v
 rr(n)=u-v
 GO TO 130
 80 rr(n1)=u-v
 rr(n)=u+v
 GO TO 130
 100 IF(t < 0.0) THEN
   GO TO   120
 END IF
 110 rr(n1)=a(n1n1)
 rr(n)=a(nn)
 GO TO 130
 120 rr(n1)=a(nn)
 rr(n)=a(n1n1)
 130 ri(n)=0.0
 ri(n1)=0.0
 GO TO 160
 140 rr(n1)=u
 rr(n)=u
 ri(n1)=v
 ri(n)=-v
 160 IF(n2 > 0) THEN
   GO TO   180
 ELSE
   GO TO  1280
 END IF
 
!        TESTS OF CONVERGENCE
 
 180 n1n2=n1n1-ia
 rmod=rr(n1)*rr(n1)+ri(n1)*ri(n1)
 eps=e10*DSQRT(rmod)
 IF(DABS(a(n1n2))-eps > 0.0) THEN
   GO TO   240
 ELSE
   GO TO  1280
 END IF
 240   IF(DABS(a(nn1))-e10*DABS(a(nn)) > 0.0) THEN
   GO TO   250
 ELSE
   GO TO  1300
 END IF
 250   IF(DABS(pan1-a(n1n2))-DABS(a(n1n2))*e6 > 0.0) THEN
   GO TO   260
 ELSE
   GO TO  1240
 END IF
 260   IF(DABS(pan-a(nn1))-DABS(a(nn1))*e6 > 0.0) THEN
   GO TO   300
 ELSE
   GO TO  1240
 END IF
 300 IF(it-maxit < 0) THEN
   GO TO   320
 ELSE
   GO TO  1240
 END IF
 
!        COMPUTE THE SHIFT
 
 320 j=1
 DO  i=1,2
   k=np-i
   IF(DABS(rr(k)-prr(i))+DABS(ri(k)-pri(i))-delta*(DABS(rr(k))  &
       +DABS(ri(k)))) 340,360,360
   340 j=j+i
   360 CONTINUE
 END DO
 SELECT CASE ( j )
   CASE (    1)
     GO TO 440
   CASE (    2)
     GO TO 460
   CASE (    3)
     GO TO 460
   CASE (    4)
     GO TO 480
 END SELECT
 440 r=0.0
 s=0.0
 GO TO 500
 460 j=n+2-j
 r=rr(j)*rr(j)
 s=rr(j)+rr(j)
 GO TO 500
 480 r=rr(n)*rr(n1)-ri(n)*ri(n1)
 s=rr(n)+rr(n1)
 
!        SAVE THE LAST TWO SUBDIAGONAL TERMS AND THE ROOTS OF THE
!        SUBMATRIX BEFORE ITERATION
 
 500 pan=a(nn1)
 pan1=a(n1n2)
 DO  i=1,2
   k=np-i
   prr(i)=rr(k)
   pri(i)=ri(k)
 END DO
 
!        SEARCH FOR A PARTITION OF THE MATRIX, DEFINED BY P AND Q
 
 p=n2
 ipi = n1n2
 DO  j=2,n2
   ipi = ipi - ia - 1
   IF(DABS(a(ipi))-eps > 0.0) THEN
     GO TO   530
   ELSE
     GO TO   600
   END IF
   530 ipip=ipi+ia
   ipip2=ipip+ia
   d=a(ipip)*(a(ipip)-s)+a(ipip2)*a(ipip+1)+r
   IF(d == 0.0) THEN
     GO TO   560
   END IF
   540   IF(DABS(a(ipi)*a(ipip+1))*(DABS(a(ipip)+a(ipip2+1)-s)+  &
       DABS(a(ipip2+2)))-DABS(d)*eps) 620,620,560
   560 p=n1-j
 END DO
 600 q=p
 GO TO 680
 620 p1=p-1
 q = p
 DO  i = 1,p1
   ipi = ipi - ia - 1
   IF(DABS(a(ipi))-eps > 0.0) THEN
     GO TO   660
   ELSE
     GO TO   680
   END IF
   660 q=q-1
 END DO
 
!        QR DOUBLE ITERATION
 
 680 ii=(p-1)*ia+p
 DO  i=p,n1
   ii1=ii-ia
   iip=ii+ia
   IF(i-p == 0) THEN
     GO TO   700
   ELSE
     GO TO   720
   END IF
   700 ipi=ii+1
   ipip=iip+1
   
!        INITIALIZATION OF THE TRANSFORMATION
   
   g1=a(ii)*(a(ii)-s)+a(iip)*a(ipi)+r
   g2=a(ipi)*(a(ipip)+a(ii)-s)
   g3=a(ipi)*a(ipip+1)
   a(ipi+1)=0.0
   GO TO 780
   720 g1=a(ii1)
   g2=a(ii1+1)
   IF(i-n2 > 0) THEN
     GO TO   760
   END IF
   740 g3=a(ii1+2)
   GO TO 780
   760 g3=0.0
   780   cap=DSQRT(g1*g1+g2*g2+g3*g3)
   IF(cap == 0.0) THEN
     GO TO   860
   END IF
   800 IF(g1 < 0.0) THEN
     GO TO   820
   ELSE
     GO TO   840
   END IF
   820 cap=-cap
   840 t=g1+cap
   psi1=g2/t
   psi2=g3/t
   alpha=2.0D0/(1.0D0+psi1*psi1+psi2*psi2)
   GO TO 880
   860 alpha=2.0
   psi1=0.0
   psi2=0.0
   880 IF(i-q == 0) THEN
     GO TO   960
   END IF
   900 IF(i-p == 0) THEN
     GO TO   940
   END IF
   920 a(ii1)=-cap
   GO TO 960
   940 a(ii1)=-a(ii1)
   
!        ROW OPERATION
   
   960 ij=ii
   DO  j=i,n
     t=psi1*a(ij+1)
     IF(i-n1 < 0) THEN
       GO TO   980
     ELSE
       GO TO  1000
     END IF
     980 ip2j=ij+2
     t=t+psi2*a(ip2j)
     1000 eta=alpha*(t+a(ij))
     a(ij)=a(ij)-eta
     a(ij+1)=a(ij+1)-psi1*eta
     IF(i-n1 < 0) THEN
       GO TO  1020
     ELSE
       GO TO  1040
     END IF
     1020 a(ip2j)=a(ip2j)-psi2*eta
     1040 ij=ij+ia
   END DO
   
!        COLUMN OPERATION
   
   IF(i-n1 < 0) THEN
     GO TO  1080
   END IF
   1060 k=n
   GO TO 1100
   1080 k=i+2
   1100 ip=iip-i
   DO  j=q,k
     jip=ip+j
     ji=jip-ia
     t=psi1*a(jip)
     IF(i-n1 < 0) THEN
       GO TO  1120
     ELSE
       GO TO  1140
     END IF
     1120 jip2=jip+ia
     t=t+psi2*a(jip2)
     1140 eta=alpha*(t+a(ji))
     a(ji)=a(ji)-eta
     a(jip)=a(jip)-eta*psi1
     IF(i-n1 < 0) THEN
       GO TO  1160
     ELSE
       GO TO  1180
     END IF
     1160 a(jip2)=a(jip2)-eta*psi2
     1180 CONTINUE
   END DO
   IF(i-n2 < 0) THEN
     GO TO  1200
   ELSE
     GO TO  1220
   END IF
   1200 ji=ii+3
   jip=ji+ia
   jip2=jip+ia
   eta=alpha*psi2*a(jip2)
   a(ji)=-eta
   a(jip)=-eta*psi1
   a(jip2)=a(jip2)-eta*psi2
   1220 ii=iip+1
 END DO
 it=it+1
 GO TO 60
 
!        END OF ITERATION
 
 1240  IF(DABS(a(nn1))-DABS(a(n1n2)) < 0.0) THEN
   GO TO  1300
 END IF
 
!        TWO EIGENVALUES HAVE BEEN FOUND
 
 1280 iana(n)=0
 iana(n1)=2
 n=n2
 IF(n2 > 0) THEN
   GO TO    20
 ELSE
   GO TO  1400
 END IF
 
!        ONE EIGENVALUE HAS BEEN FOUND
 
 1300 rr(n)=a(nn)
 ri(n)=0.0
 iana(n)=1
 IF(n1 > 0) THEN
   GO TO  1320
 ELSE
   GO TO  1400
 END IF
 1320 n=n1
 GO TO 20
 1400  CONTINUE
 k=0
 DO  i=1,m
   rra(i)=rr(i)
   rri(i)=ri(i)
   DO  j=1,m
     k=k+1
     b(k)=a(k)
   END DO
 END DO
 RETURN
END SUBROUTINE ateig
