SUBROUTINE mbgae(ajjl,in17,a,f,df,f1,df1,f2,df2,q,q1,q2,mood)
     
!     MULTIPLY SUM OBTAINED PREVIOUSLY BY SCRIPT A FACTOR
 
 
 INTEGER, INTENT(IN OUT)                  :: ajjl
 INTEGER, INTENT(IN OUT)                  :: in17
 COMPLEX, INTENT(OUT)                     :: a(1)
 REAL, INTENT(IN)                         :: f(1)
 REAL, INTENT(IN)                         :: df(1)
 REAL, INTENT(IN)                         :: f1(1)
 REAL, INTENT(IN)                         :: df1(1)
 REAL, INTENT(IN)                         :: f2(1)
 REAL, INTENT(IN)                         :: df2(1)
 COMPLEX, INTENT(IN)                      :: q(1)
 COMPLEX, INTENT(IN)                      :: q1(1)
 COMPLEX, INTENT(IN)                      :: q2(1)
 INTEGER, INTENT(OUT)                     :: mood
 LOGICAL :: cntrl2 , cntrl1 , crank1 , crank2 , asym , debug
 
 REAL :: mach
 
 
 COMMON /mboxc/ njj ,crank1,crank2,cntrl1,cntrl2,nbox,  &
     npts0,npts1,npts2,asym,gc,cr,mach,beta,ek,ekbar,ekm,  &
     boxl,boxw,boxa ,ncb,nsb,nsbd,ntote,kc,kc1,kc2,kct,kc1t,kc2t
 COMMON /system/ sysbuf,n6
 COMMON /amgmn  / mcb(7)
 DATA    debug /.false./
 
 gck  =  gc * boxw
 DO  i=1,njj
   a(i) = (0.0,0.0)
 END DO
 DO  i=1,npts0
   CALL fread(in17,f  ,kct ,0)
   CALL fread(in17,df ,kct ,0)
   DO  j=1,kc
     a(i)      = a(i)      + CMPLX( df(j),-ek*f(j))*q(j)
   END DO
   IF( kc == kct )  CYCLE
   kcc  =  kc + 1
   DO      j = kcc,kct
     a(i)       = a(i)      + f(j)*q(j)
   END DO
 END DO
 IF( .NOT. cntrl1 ) GO TO 1660
 jj = npts0
 DO  i=1,npts1
   CALL fread(in17,f1 ,kc1t,0)
   CALL fread(in17,df1,kc1t,0)
   DO  j=1,kc1
     a(i+jj)   = a(i+jj)   + CMPLX( df1(j),-ek*f1(j))*q1(j)
   END DO
   IF( kc1 == kc1t )  CYCLE
   kcc1  =  kc1 + 1
   DO      j = kcc1,kc1t
     a(i+jj)    =  a(i+jj)   + f1(j)*q1(j)
   END DO
 END DO
 1660 IF( .NOT. cntrl2 ) GO TO 1700
 jj = jj + npts1
 DO  i=1,npts2
   CALL fread(in17,f2 ,kc2t,0)
   CALL fread(in17,df2,kc2t,0)
   DO  j=1,kc2
     a(i+jj)   = a(i+jj)   + CMPLX( df2(j),-ek*f2(j))*q2(j)
   END DO
   IF( kc2 == kc2t )  CYCLE
   kcc2  = kc2 + 1
   DO      j = kcc2,kc2t
     a(i+jj)    = a(i+jj)   + f2(j)*q2(j)
   END DO
 END DO
 1700 CONTINUE
 CALL bckrec(in17)
 DO  i=1,njj
   a(i) = a(i) * gck
 END DO
 CALL pack(a,ajjl,mcb)
 
!     PRINT OUT GENERALIZED AERODYNAMIC FORCE COEFFICIENTS
 
 IF(.NOT.debug) RETURN
 IF(mood > 1) GO TO 2100
 WRITE  (n6 , 1900 )  mach ,  boxl , ek , boxw
 1900 FORMAT  ( 1H1 , 31X , 30HGENERALIZED aerodynamic force  &
     , 12HCOEFFICIENTS / 1H0 , 9X , 11HMACH NUMBER , f9.3 ,  &
     40X , 10HBOX length , f12.6 / 1H0  &
     , 9X , 33HREDUCED frequency  ( root chord ) , f10.5 , 17X  &
     , 9HBOX width , f13.6 / 1H0 , 42X , 21H- -  a ( i , j )  - -  &
     / 6H-  row , 9X , 4HREAL , 10X , 4HIMAG , 14X , 4HREAL , 10X  &
     , 4HIMAG , 14X , 4HREAL , 10X , 4HIMAG )
 2100 WRITE(n6,2000) mood, (a(j),j=1,njj)
 2000 FORMAT  ( 1H0 , i4 , 3 ( e18.4 , e14.4 ) / ( 1H0 , 4X , 3 ( e18.4  &
     , e14.4 ) ) )
 RETURN
END SUBROUTINE mbgae
