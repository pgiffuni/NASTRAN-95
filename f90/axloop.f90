SUBROUTINE axloop (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
     
 
    REAL, INTENT(IN)                         :: buf(50)
    INTEGER, INTENT(IN)                      :: ibuf(50)
    REAL, INTENT(IN)                         :: xx
    REAL, INTENT(IN)                         :: yy
    REAL, INTENT(IN)                         :: zz
    REAL, INTENT(OUT)                        :: hc1
    REAL, INTENT(OUT)                        :: hc2
    REAL, INTENT(OUT)                        :: hc3
    INTEGER :: otpe
 
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm,uwm
    COMMON /system/ sysbuf,otpe
    COMMON /BLANK / idum(3),epse
 
    pi    = 3.1415926536
    piby2 = 1.5707963268
    fpi   = 12.56637062
    c     = 1.
 
    xj  = buf(1)
    iaxi= ibuf(2)
    x1  = buf(3)
    y1  = buf(4)
    z1  = buf(5)
    x2  = buf(6)
    y2  = buf(7)
    z2  = buf(8)
    xc  = buf(9)
    yc  = buf(10)
    zc  = buf(11)
 
    !     FOR NOW, ICID = 0
 
    icid = ibuf(12)
 
    !     CHECK FOR AXISYMMETRIC PROBLEM
 
    IF (iaxi /= 1) GO TO 10
    xc = 0.
    yc = 0.
    zc = z1
    x2 = 0.
    y2 = x1
    z2 = z1
10 CONTINUE
 
   !     DETERMINE THE DIRECTION OF THE CURRENT LOOP AXIS
 
   cx = x1 - xc
   cy = y1 - yc
   cz = z1 - zc
   bx = x2 - xc
   by = y2 - yc
   bz = z2 - zc
 
   !     THE VECTOR AN IS NORMAL TO THE PLANE OF THE LOOP
 
   anx = cy*bz - cz*by
   any = cz*bx - cx*bz
   anz = cx*by - cy*bx
   at1 = SQRT(anx*anx + any*any + anz*anz)
   at2  = bx*bx + by*by + bz*bz
   rad2 = cx*cx + cy*cy + cz*cz
   radius = SQRT(rad2)
   xiacpi = (xj*rad2*pi)/c
 
   anx = anx/at1
   any = any/at1
   anz = anz/at1
 
   !     THE VECTOR R IS FROM THE CENTER OF LOOP TO THE FIELD POINT
 
   rx = xx - xc
   ry = yy - yc
   rz = zz - zc
 
   r2 = rx*rx + ry*ry + rz*rz
   r  = SQRT(r2)
 
   !     AT (OR NEAR) CENTER OF LOOP TEST
 
   IF (r >= .001) GO TO 218
   costhe = 1.
   sinthe = 0.
   sqar2s = SQRT(rad2+r2)
   rx  = anx
   ry  = any
   rz  = anz
   rpx = 0.
   rpy = 0.
   rpz = 0.
   GO TO 220
218 CONTINUE
 
    rx = rx/r
    ry = ry/r
    rz = rz/r
    costhe = anx*rx + any*ry + anz*rz
    sinthe = SQRT(1. - costhe*costhe)
 
    !     ON (OR VERY NEAR) AXIS OF LOOP TEST
 
    IF (sinthe >= .000001) GO TO 219
    costhe = 1.
    sinthe = 0.
    sqar2s = SQRT(rad2+r2)
    rx  = anx
    ry  = any
    rz  = anz
    rpx = 0.
    rpy = 0.
    rpz = 0.
    GO TO 220
219 CONTINUE
 
    sqar2s = SQRT(rad2 + r2 + (2.*radius*r*sinthe))
    realk2 = (4.*radius*r*sinthe)/(rad2+r2+(2.*radius*r*sinthe))
    realk  = SQRT(realk2)
    xiacr  = (xj*radius)/(c*r)
 
    !     A CROSS R, NORMAL TO THE PLANE OF A AND R
 
    tx = any*rz - anz*ry
    ty = anz*rx - anx*rz
    tz = anx*ry - any*rx
 
    !     (A CROSS R) CROSS R, NORMAL TO THE PLANE OF R AND (A AND R)
 
    trpx = ty*rz - tz*ry
    trpy = tz*rx - tx*rz
    trpz = tx*ry - ty*rx
    at3  = SQRT(trpx*trpx + trpy*trpy + trpz*trpz)
 
    !     RPERP, PERPENDICULAR TO THE VECTOR FROM THE CENTER TO THE FIELD PT
 
    rpx = trpx/at3
    rpy = trpy/at3
    rpz = trpz/at3
 
    !     FOR SMALL POLAR ANGLE OR SMALL RADIUS USE ALTERNATIVE APPROX.
 
    IF (realk2 < .0001) GO TO 220
 
    !     COMPUTE ELLIPTIC INTEGRAL OF FIRST KIND
 
    f = 1.
    deltf1 = 1.
    DO  n = 1,15000
        xn2  = 2.*FLOAT(n)
        xn21 = xn2 - 1.
        deltf1 = deltf1*(xn21/xn2)*realk
        deltf2 = deltf1*deltf1
        f = f + deltf2
        IF (ABS(deltf2/f) <= epse) GO TO 250
    END DO
    delf = ABS(deltf2/f)
    WRITE (otpe,245) uwm,xx,yy,zz,xc,yc,zc,x1,y1,z1,x2,y2,z2,delf,epse
245 FORMAT (a25,', CONVERGENCE OF ELLIPTIC INTEGRAL IS UNCERTAIN. ',  &
        'GRID OR INTEGRATION POINT AT COORDINATES', /5X,  &
        1P,3E15.6,'  IS TOO CLOSE TO CURRENT LOOP WITH CENTER AT',  &
        /5X,1P,3E15.6,' AND 2 POINTS AT ',1P,3E15.6, /5X,4HAND ,1P,  &
        3E15.6,' COMPUTATIONS WILL CONTINUE WITH LAST VALUES', /5X,  &
        'CONVERGENCE VALUE WAS ',1P,e15.6, ' CONVERGENCE CRITERION IS ',1P,e15.6)
250 f = piby2*f
 
    !     COMPUTE ELLIPTIC INTEGRAL OF SECOND KIND
 
    e = 1.
    delte1 = 1.
    DO  n = 1,15000
        xn2  = 2.*FLOAT(n)
        xn21 = xn2-1.
        delte1 = delte1*(xn21/xn2)*realk
        delte2 = (delte1*delte1)/xn21
        e = e - delte2
        IF (ABS(delte2/e) <= .000001) GO TO 270
    END DO
    dele = ABS(delte2/e)
    WRITE (otpe,245) uwm,xx,yy,zz,xc,yc,zc,x1,y1,z1,x2,y2,z2,dele
270 e = piby2*e
 
    !     COMPUTE THE RADIAL COMPONENT OF THE MAGNETIC FIELD
 
    br = xiacr*(costhe/sinthe)*(e/sqar2s)*(realk2/(1.-realk2))
 
    !     COMPUTE THE POLAR COMPONENT OF THE MAGNETIC FIELD
 
    bthe = xiacr*(1./(sqar2s*radius*r*sinthe))*  &
        (((((2.*r2)-((r2+(radius*r*sinthe))*realk2))/ (1.-realk2))*e)-(2.*r2*f))
 
    !     GO TO THE RESOLUTION OF FIELD COMPONENTS
 
    GO TO 230
 
!     ALTERNATIVE APPROXIMATION FOR SMALL K**2
 
!     COMPUTE THE RADIAL COMPONENT OF THE MAGNETIC FIELD
 
220 CONTINUE
    br = xiacpi*costhe*(((2.*rad2)+(2.*r2)+(radius*r*sinthe))/ ((sqar2s)**5))
 
    !     COMPUTE THE POLAR COMPONENT OF THE MAGNETIC FIELD
 
    bthe = -xiacpi*sinthe* (((2.*rad2)-r2+(radius*r*sinthe))/((sqar2s)**5))
 
!     RESOLVE MAGNETIC FIELD COMPONENTS INTO RECTANGULAR COMPONENTS
 
230 CONTINUE
    hcx = rx*br + rpx*bthe
    hcy = ry*br + rpy*bthe
    hcz = rz*br + rpz*bthe
    hc1 = hcx/fpi
    hc2 = hcy/fpi
    hc3 = hcz/fpi

    RETURN
END SUBROUTINE axloop
