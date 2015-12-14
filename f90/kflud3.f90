SUBROUTINE kflud3
     
!     THIS ROUTINE GENERATES THE PSUEDO STIFFNESS MATRIX TERMS
!     FOR THE TRIANGULAR FLUID ELEMENT
 
!     THE ECPT DATA IS THE FOLLOWING
 
!         FIELD         SYMBOL
!           1             ID
!           2             SIL1
!           3             SIL2
!           4             SIL3
!           5             RHO
!           6             BULK
!           7             N
!           8             CSF
!           9             R1
!           10            Z1
!           11            -
!           12            CSF
!           13            R2
!           14            Z2
!           15            -
!           16            CSF
!           17            R3
!           18            Z3
!           19            -
!           20            -
 
 LOGICAL :: nogo
 INTEGER :: out      ,np(3)    ,necpt(100)
 DOUBLE PRECISION :: constd   ,dpi      ,  &
     r        ,r1       ,r2       ,r3       ,  &
     z1       ,z2       ,z3       ,deth     ,  &
     h        ,ra       ,rb       ,za       ,  &
     zb       ,dr       ,dz       ,beta     ,  &
     blog     ,dzr      ,dzr2     ,bet2     ,  &
     r12      ,r22      ,g00      ,g10      ,  &
     g20      ,g01      ,g11      ,g02      ,  &
     rn       ,piro     ,prn2     ,kq       , kvec     ,kg
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /condad/  constd(5)
 COMMON /system/  sysbuf   ,out      ,nogo
 COMMON /sma1io/  dum1(10) ,ifkgg
 COMMON /sma1cl/  iopt4    ,k4ggsw   ,npvt
 COMMON /sma1et/  ecpt(100)
 COMMON /sma1dp/  r        ,r1       ,r2       ,r3       ,  &
     z1       ,z2       ,z3       ,deth     ,  &
     h(9)     ,ra       ,rb       ,za       ,  &
     zb       ,dr       ,dz       ,beta     ,  &
     blog     ,dzr      ,dzr2     ,bet2     ,  &
     r12      ,r22      ,g00      ,g10      ,  &
     g20      ,g01      ,g11      ,g02      ,  &
     rn       ,piro     ,prn2     ,kq(9)    ,  &
     kvec(3)  ,kg(3)    ,iret     ,ipt      , jc       ,ir
 EQUIVALENCE      (constd(1),dpi)    ,(ecpt(1),necpt(1))
 
!     SELECT POINTS FOR COUNTERCLOCKWISE ORDER
 
 np(1) = necpt(2)
 np(2) = necpt(3)
 np(3) = necpt(4)
 r1    = ecpt( 9)
 z1    = ecpt(10)
 r2    = ecpt(13)
 z2    = ecpt(14)
 r3    = ecpt(17)
 z3    = ecpt(18)
 r     = (r2-r1)*(z3-z1) - (r3-r1)*(z2-z1)
 IF (r < 0.0) THEN
   GO TO    10
 ELSE IF (r == 0.0) THEN
   GO TO  2000
 ELSE
   GO TO    20
 END IF
 10 np(2) = np(3)
 np(3) = necpt(3)
 r2    = ecpt(17)
 r3    = ecpt(13)
 z2    = ecpt(18)
 z3    = ecpt(14)
 20 IF (r1 <= 0.0D0 .OR. r2 <= 0.0D0 .OR. r3 <= 0.0D0) GO TO 1000
 deth  = DABS(r)
 h(1)  = (r2*z3-r3*z2)/deth
 h(4)  = (r3*z1-r1*z3)/deth
 h(7)  = (r1*z2-r2*z1)/deth
 h(2)  = (z2-z3)/deth
 h(5)  = (z3-z1)/deth
 h(8)  = (z1-z2)/deth
 h(3)  = (r3-r2)/deth
 h(6)  = (r1-r3)/deth
 h(9)  = (r2-r1)/deth
 
!     THE INTEGRAL PARAMETERS ARE THE SUM DUE TO SIDES 1-2,2-3,3-1.
 
 g00   = 0.0
 g01   = 0.0
 g02   = 0.0
 g10   = 0.0
 g11   = 0.0
 g20   = 0.0
 iret  = 1
 ra    = r1
 rb    = r2
 za    = z1
 zb    = z2
 GO TO 500
 100 iret  = 2
 ra    = r2
 rb    = r3
 za    = z2
 zb    = z3
 GO TO 500
 110 iret  = 3
 ra    = r3
 rb    = r1
 za    = z3
 zb    = z1
 
!     THE INTEGRAL PARAMETERS ARE CALCULATED BELOW
 
 500 dr    = rb - ra
 dz    = zb - za
 IF (dr**2/deth <= 1.0D-6) GO TO 140
 beta  = za - ra*dz/dr
 bet2  = beta**2
 blog  = beta*DLOG(ra/rb)
 dzr   = dz/dr
 dzr2  = dzr**2
 r12   = (ra**2-rb**2)/2.0D0
 r22   = (ra**3-rb**3)/3.0D0
 g00   = g00 + blog - dz
 g10   = g10 - beta*dr  + r12*dzr
 g20   = g20 + beta*r12 + dzr*r22
 g01   = g01 + blog*beta/2.0D0 - beta*dz + dzr2*r12/2.0D0
 g11   = g11 - bet2*dr/2.0D0 + beta*dzr*r12 + dzr2*r22/2.0D0
 g02   = g02 + bet2*blog/3.0D0 - bet2*dz + beta*dzr2*r12 + dzr*dzr2*r22/3.0D0
 140 CONTINUE
 SELECT CASE ( iret )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 110
   CASE (    3)
     GO TO 120
 END SELECT
 120 CONTINUE
 
!     FORM THE PSUEDO STIFFNESS MATRIX USING THE PARAMETERS
 
 rn    = necpt(7)
 IF (ecpt(5) <= 0.0) RETURN
 piro  = dpi/DBLE(ecpt(5))
 IF(necpt(7) == 0) piro=piro*2.0D0
 prn2  = piro*rn**2
 kq(1) = prn2*g00
 kq(2) = prn2*g10
 kq(3) = prn2*g01
 kq(4) = kq(2)
 kq(5) = (piro + prn2)*g20
 kq(6) = prn2*g11
 kq(7) = kq(3)
 kq(8) = kq(6)
 kq(9) = piro*g20 + prn2*g02
 DO  i = 1,3
   ipt   = i - 1
   IF (npvt == np(i)) GO TO 210
 END DO
 RETURN
 
 210 ipt   = 3*ipt + 1
 CALL gmmatd (h(ipt),1,3,0,kq,3,3,0,kvec)
 CALL gmmatd (kvec,1,3,0,h(1),3,3,1,kg)
 jc    = npvt
 DO  i = 1,3
   ir    = np(i)
   CALL sma1b (kg(i),ir,jc,ifkgg,0.0D0)
 END DO
 2000 RETURN
 
 1000 ir    = necpt(1)/1000
 WRITE  (out,5001) ufm,ir
 5001 FORMAT (a23,' 5001, NEG. OR ZERO RADIUS DETECTED FOR CFLUID3 OR',  &
     ' CFLUID4 ELEMENT',i12)
 nogo = .true.
 RETURN
END SUBROUTINE kflud3
