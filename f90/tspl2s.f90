SUBROUTINE tspl2s (ts7)
     
!    TRANSVERSE SHEAR ROUTINE2 FOR CTRPLT1 - SINGLE PRECISION VERSION
 
 
 REAL, INTENT(OUT)                        :: ts7(60)
 
 COMMON /sma1io/ x,y
 
 DO  i=1,60
   ts7(i)=0.0
 END DO
 x2=x*x
 xy=x*y
 y2=y*y
 x3=x2*x
 x2y=x2*y
 xy2=x*y2
 y3=y2*y
 ts7(   4)=2.0
 ts7(   7)=6.0*x
 ts7(   8)=2.0*y
 ts7(  11)=12.0*x2
 ts7(  12)=6.0*xy
 ts7(  13)=2.0*y2
 ts7(  16)=20.0*x3
 ts7(  17)=6.0*xy2
 ts7(  18)=2.0*y3
 ts7(  26)=2.0
 ts7(  29)=2.0*x
 ts7(  30)=6.0*y
 ts7(  33)=2.0*x2
 ts7(  34)=ts7(12)
 ts7(  35)=12.0*y2
 ts7(  37)=2.0*x3
 ts7(  38)=6.0*x2y
 ts7(  39)=12.0*xy2
 ts7(  40)=20.0*y3
 ts7(  45)=2.0
 ts7(  48)=4.0*x
 ts7(  49)=4.0*y
 ts7(  52)=6.0*x2
 ts7(  53)=8.0*xy
 ts7(  54)=6.0*y2
 ts7(  57)=12.0*x2y
 ts7(  58)=ts7(39)
 ts7(  59)=8.0*y3
 RETURN
END SUBROUTINE tspl2s
