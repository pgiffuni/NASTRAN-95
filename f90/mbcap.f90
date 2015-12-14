SUBROUTINE mbcap(nphi,capphi)
     
 
 INTEGER, INTENT(IN OUT)                  :: nphi
 COMPLEX, INTENT(OUT)                     :: capphi(1)
 REAL :: km , kbar , mach , w(10), p(10)
 
 COMMON /mboxc/ njj ,crank1,crank2,cntrl1,cntrl2,nbox,  &
     npts0,npts1,npts2,asym,gc,cr,mach,beta,ek,ekbar,ekm,  &
     boxl,boxw,boxa ,ncb,nsb,nsbd,ntote,kc,kc1,kc2,kct,kc1t,kc2t
 EQUIVALENCE  ( km , ekm ) , ( kbar , ekbar )
 DATA  w / 0.0506143,0.1111905,0.1568533,0.1813419,0.1813419,  &
     0.1568533,0.1111905,0.0506143,0.0,0.0/,  &
     p / 0.0198551,0.1016667,0.2372338,0.4082826,0.5917174,  &
     0.7627662,0.8983333,0.9801449,0.0,0.0/
 
 DO      i = 1 , nphi
   capphi(i)  =  ( 0.0 , 0.0 )
 END DO
 
!     COMPUTE CAPPHI FOR RECEIVING BOX
 
 IF ( kbar <= 0.0 )   GO TO  400
 DO      i = 1 , 8
   j  =  9 - i
   arg  =  kbar * p(j) / 2.0
   arg1  =  w(i) * zj ( arg / mach ) / 2.0
   capphi(1)  =  capphi(1) + CMPLX ( -COS ( arg ) * arg1 ,  &
       SIN ( arg ) * arg1 )
 END DO
 GO TO 500
 
 400  capphi(1)  =  ( -0.5 , 0.0 )
 
!     COMPUTE REMAINING CAPPHI
 
 500  nphi  =  1
 xb  =  0.5
 xu  =  xb + 1.0
 DO      i = 2 , ncb
   xl  =  -0.5
   xr  =  xl + 1.0
   DO      j = 1 , i
     nphi  =  nphi + 1
     DO      l = 1 , 8
       x  =  xb + p(l)
       arg  =  kbar * x
       arg1  =  w(l) * GO ( x , xr , xl , km ) / 3.14159265
       capphi(nphi)  =  capphi(nphi) - CMPLX ( COS ( arg ) * arg1 ,  &
           -SIN ( arg ) * arg1 )
     END DO
     xl  =  xr
     xr  =  xr + 1.0
   END DO
   
   xb  =  xu
   xu  =  xb + 1.0
 END DO
 
 DO      i = 1 , nphi
   capphi(i)  =  boxw * capphi(i)
 END DO
 RETURN
END SUBROUTINE mbcap
