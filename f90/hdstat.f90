SUBROUTINE hdstat(mt,nit,ixr,x21,y21,z21,iia,iv,a,b,c,  &
        ik,xa,ya,za,ccc,xxx,lz)
     
 
!     THIS SUBROUTINE TAKES THE PTS OF INTERSECTION DETERMINED BY
!     SUBROUTINE SOLVE AND PICKS THE COORDINATES WITH THE MAX AND
!     MIN X COORDINATES PROVIDED THEY LIE ON THE INTERIOR/BOUNDARY
!     OF BOTH ELEMENTS.
 
 
 
 INTEGER, INTENT(IN)                      :: mt
 INTEGER, INTENT(OUT)                     :: nit
 INTEGER, INTENT(IN)                      :: ixr
 REAL, INTENT(OUT)                        :: x21(1)
 REAL, INTENT(OUT)                        :: y21(1)
 REAL, INTENT(OUT)                        :: z21(1)
 INTEGER, INTENT(OUT)                     :: iia(1)
 INTEGER, INTENT(IN)                      :: iv(1)
 REAL, INTENT(IN)                         :: a
 REAL, INTENT(IN)                         :: b
 REAL, INTENT(IN)                         :: c
 INTEGER, INTENT(IN OUT)                  :: ik
 REAL, INTENT(IN OUT)                     :: xa(1)
 REAL, INTENT(IN OUT)                     :: ya(1)
 REAL, INTENT(OUT)                        :: za(1)
 REAL, INTENT(IN)                         :: ccc(1)
 REAL, INTENT(IN)                         :: xxx(1)
 INTEGER, INTENT(IN)                      :: lz
 INTEGER :: xcc
 
 COMMON/hdptrs/xdum,xcc
 COMMON/zzzzzz/rz(1)
 COMMON/go3/l0,l1,l00,l01,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13
 
 exx=.015
 nx=0
 IF(mt == 0)GO TO 160
 loop50:  DO  jx=1,mt
   ei=0
   10 ei=ei+.1
   IF(ei >= 1.) GO TO 160
   d=ei*xa(jx)-ya(jx)
   loop40:  DO  jo=1,2
     m=iv(jo)
     jc=l13+(m-1)*lz
     jxc=l12+(m-1)*5
     nk=xxx(5+jxc)
     i=0
     ib=nk*5
     
     
!     DETERMINE IF THE PROJECTION OF THE POINT OF INTERSECTION
!     BELONGS TO THE INTERIOR OF BOTH PLANES.
     
     
     DO  j=1,ib,5
       exx=.015
       nsub=j+1+jc
       IF(ABS(ccc(nsub)) >= 100.)exx=ALOG10(ABS(ccc(nsub)))
       ve=xa(jx)
       IF(ccc(j+jc) == 0.)ve=ya(jx)
       s=ve-ccc(j+3+jc)
       s1=ve-ccc(j+4+jc)
       t=ccc(j+jc)*ya(jx)+ccc(j+1+jc)*xa(jx)+ccc(j+2+jc)
       IF((ABS(t) < exx).AND.(s*s1 <= 0.))CYCLE loop40
       t=-ccc(j+2+jc)+ccc(j+jc)*d
       r=ei*ccc(j+jc)+ccc(j+1+jc)
       IF(r == 0.)CYCLE
       t=t/r
       IF(t < xa(jx))CYCLE
       IF(ccc(j+jc) /= 0.)GO TO 20
       t=ei*t-d
       20 CONTINUE
       IF((t == ccc(j+3+jc)).OR.(t == ccc(j+4+jc)))GO TO 10
       s=t-ccc(j+3+jc)
       s1=t-ccc(j+4+jc)
       IF(s*s1 > 0.)CYCLE
       i=i+1
     END DO
     IF(MOD(i,2) == 0)CYCLE loop50
   END DO loop40
   nx=nx+1
   xa(nx)=xa(jx)
   ya(nx)=ya(jx)
   za(nx)=za(jx)
 END DO loop50
 IF(nx == 0)GO TO 160
 
 
 
!     THIS CODE FINDS THE MAX/MIN X-COORDINATES(Y-COORDINATES) AND
!     STORES THEM. FUTHERMORE BOTH THE EQUATION OF LINE AND POINTS(2)
!     ARE TREATED LIKE ADDITIONAL EDGES. IN THIS WAY, THE ALGORITHM NEED
!     NOT BE DISTURBED. ESSENTIALLY,THEN,THIS TRICK IS TRANSPARENT TO
!     THE REST OF THE PROGRAM.
 
 
 amaxx=-(10**6)
 aminx=-amaxx
 amaxy=amaxx
 aminy=aminx
 is=5+(ik-1)*5+l12
 is=xxx(is)
 DO  ji=1,nx
   IF(a == 0.)GO TO 80
   IF(xa(ji) >= aminx)GO TO 60
   aminx=xa(ji)
   yi=ya(ji)
   zi=za(ji)
   60 IF(xa(ji) <= amaxx)GO TO 70
   amaxx=xa(ji)
   yii=ya(ji)
   zii=za(ji)
   70 CONTINUE
   CYCLE
   80 CONTINUE
   IF(ya(ji) >= aminy)GO TO 90
   aminy=ya(ji)
   xi=xa(ji)
   zi=za(ji)
   90 CONTINUE
   IF(ya(ji) <= amaxy)GO TO 100
   xii=xa(ji)
   amaxy=ya(ji)
   zii=za(ji)
   100 CONTINUE
 END DO
 nit=nit+1
 k=5*(nit-1+is)+1
 rz(xcc+k-1)=a
 rz(xcc+k  )=b
 rz(xcc+k+1)=c
 IF (a == 0.) GO TO 120
 rz(xcc+k+2)=aminx
 rz(xcc+k+3)=amaxx
 amin=aminx
 amax=amaxx
 ye=yii
 ze=zii
 GO TO 130
 120 CONTINUE
 rz(xcc+k+2)=aminy
 rz(xcc+k+3)=amaxy
 amin=xi
 amax=xii
 yi=aminy
 ye=amaxy
 ze=zii
 130 CONTINUE
 ig=ixr+nit*3
 x21(ig-2)=amin
 y21(ig-2)=yi
 z21(ig-2)=zi
 DO  jk=1,2
   ie=ig-jk+1
   x21(ie)=amax
   y21(ie)=ye
   z21(ie)=ze
 END DO
 DO  jk=1,2
   iia(ig-jk)=0
 END DO
 iia(ig)=1
 tx=(amax-amin)**2
 ty=(ye-yi)**2
 dx=(tx+ty)**.5
 IF(dx < .001)nit=nit-1
 160 RETURN
END SUBROUTINE hdstat
