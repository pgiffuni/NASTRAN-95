SUBROUTINE mbprit(aw,ac,at)
     
!     SUBROUTINE TO PRINT GEOMETRY DATA
 
 
 REAL, INTENT(IN OUT)                     :: aw
 REAL, INTENT(IN OUT)                     :: ac
 REAL, INTENT(IN OUT)                     :: at
 LOGICAL :: cntrl2 , cntrl1 , crank1 , crank2 , asym
 COMMON /system/ sys,n6
 COMMON /mboxa/ x(12),y(12),tang(10),ang(10),cotang(10)
 COMMON /mboxc/ njj ,crank1,crank2,cntrl1,cntrl2,nbox,  &
     npts0,npts1,npts2,asym,gc,cr,mach,beta,ek,ekbar,ekm,  &
     boxl,boxw,boxa ,ncb,nsb,nsbd,ntote,kc,kc1,kc2,kct,kc1t,kc2t
 
 WRITE  (n6 , 200 )  cntrl2 , cntrl1 , crank1 , crank2 , asym
 200  FORMAT  ( 1H1 , 35X , 27HSUPERSONIC mach box PROGRAM / 1H0 , 43X  &
     , 12HCONTROL DATA / l20 , 9X , 6HCNTRL2 / l20 , 9X  &
     , 6HCNTRL1 / l20 , 9X , 21HCRANK  (leading edge)  &
     / l20 , 9X , 22HCRANK  (trailing edge) / l20 , 9X , 14HANTI-symmetric / l20 )
 
 WRITE  (n6 , 300 )   ( i , x(i) , y(i) , tang(i) , ang(i) , i=1,7)
 300  FORMAT  (1H- , 42X , 13HGEOMETRY DATA / 1H0 , 8X , 1HN , 11X , 1HX  &
     , 17X , 1HY , 16X , 4HTANG , 14X , 3HANG / ( i10 , 4E18.6 ) )
 
 WRITE  (n6 , 400 )  ( i , x(i) , y(i) , tang(i) , i = 8 , 10)  &
     , ( i , x(i) , y(i) , i = 11 , 12 )
 400  FORMAT(i10,3E18.6/i10,3E18.6/i10,3E18.6/(i10,2E18.6))
 
 WRITE  (n6 , 500 )   aw , ac , at
 500  FORMAT  ( 1H0 , 5X , 23HAREA of main (semispan) , 11X  &
     , 15HAREA of cntrl1 , 18X , 14HAREA of cntrl2 / e22.6,e34.6,e29.6)
 RETURN
END SUBROUTINE mbprit
