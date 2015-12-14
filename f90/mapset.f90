SUBROUTINE mapset(x1,y1,x2,y2,ki1,kj1,ki2,kj2,l)
     
!     POINT 1 IS LOWER LEFT CORNER OF FRAME
!     POINT 2 IS UPPER RIGHT CORNER OF FRAME
!     I,J ARE IN PLOTTER UNITS
!     X,Y ARE IN PHYSICAL UNITS
!     L IS OUTPUT FLAG, 1=I,J ARE INTEGER, 2=I,J ARE REAL
 
 EQUIVALENCE (i1,zi1),(j1,zj1),(i2,zi2),(j2,zj2)
 EQUIVALENCE (i,zi),(j,zj)
 
 i1=ki1
 j1=kj1
 i2=ki2
 j2=kj2
 ll=l
 
 IF(l == 2) GO TO 100
 a=FLOAT(i2-i1)/(x2-x1)
 b=FLOAT(i1)-a*x1
 c=FLOAT(j2-j1)/(y2-y1)
 d=FLOAT(j1)-c*y1
 RETURN
 100 a=(zi2-zi1)/(x2-x1)
 b=zi1-a*x1
 c=(zj2-zj1)/(y2-y1)
 d=zj1-c*y1
 RETURN
 
 
!***********************************************************************
 
 ENTRY map(x,y,ki,kj)
 IF(ll == 2) GO TO 200
 i=a*x+b+0.5
 j=c*y+d+0.5
 GO TO 300
 200 zi=a*x+b
 zj=c*y+d
 300 ki=i
 kj=j
 RETURN
 
END SUBROUTINE mapset
