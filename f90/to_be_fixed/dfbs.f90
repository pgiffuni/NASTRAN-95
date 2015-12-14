SUBROUTINE dfbs
     
!     FBS   L,U,B/X/V,Y,ISYM=0/V,Y,KSIGN=1/V,Y,IPREC=0/V,Y,ITYPE=0 $
 
!     ISYM  =  1  USE FBS
!           = -1  USE GFBS
!           =  0  CHOOSE WHICH BASED ON SUPPLIED INPUT
!     KSIGN =  1, SOLVE LUX= B
!             -1,       LUX=-B
!     IPREC = REQUESTED PRECISION - DEFAULT BASED ON INPUT OR SYSTEM(55)
!     ITYPE = REQUESTED TYPE OF X - DEFAULT IS LOGICAL CHOICE ON INPUT
 
!     REVISED  12/91 BY G.CHAN/UNISYS
!     FATAL ERROR IN FBS (NOT GFBS) IF INPUT MATRIX IS NOT A LOWER
!     TRIANGULAR FACTOR
 
 INTEGER :: l,u,b,x,sbnm(2),dosi(3),refus(3),outpt,scr
 DIMENSION       zz(1)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /BLANK / isym, ksign, iprec, itype
 COMMON /system/ ksystm(65)
 COMMON /fbsx  / il(7),iu(7),ib(7),ix(7),inx,ip1,is1,iscr
 COMMON /gfbsx / jl(7),ju(7),jb(7),jx(7),jnx,jp1,js1
!ZZ   COMMON /ZZDFB1/ Z(1)
 COMMON /zzzzzz/ z(20000)
!ZZ   COMMON /ZZDFB2/ ZZ(1)
 EQUIVALENCE     (zz(1),z(1))
 EQUIVALENCE     (ksystm(55),kprec),(ksystm(2),outpt)
 DATA    l, u, b, x, scr   / 101,102,103,201,301 /
 DATA    sbnm  / 4HDFBS,1H /
 DATA    dosi  / 4HSING, 4HDOUB, 4HMLTP/,  refus / 2*3H   ,3HREF/
 
 
 ju(1) = u
 CALL rdtrl (ju)
 10 IF (isym < 0) THEN
   GO TO   150
 ELSE IF (isym == 0) THEN
   GO TO    20
 ELSE
   GO TO    30
 END IF
 20 isym  = -1
 IF (ju(1) < 0) isym = 1
 GO TO 10
 
!     SET UP CALL TO FBS
 
 30 nogo  = 0
 il(1) = l
 CALL rdtrl (il)
 IF (il(1) > 0) GO TO 40
 CALL mesage (30,198,l)
 nogo  = 1
 40 CONTINUE
 IF (il(4) /= 4) GO TO 100
 n     = il(2)
 ib(1) = b
 CALL rdtrl (ib)
 IF (nogo == 0) GO TO 50
 CALL mesage (-30,199,sbnm)
 50 CONTINUE
 inx   = korsz(z)
 iprec1= MAX0(il(5),ib(5),iu(5))
 IF (iprec1 > 2) iprec1 = iprec1 - 2
 IF (iprec1 < 1 .OR. iprec1 > 2) iprec1 = kprec
 IF (iprec == iprec1 .OR. iprec == 0) GO TO 70
 IF (iprec < 1 .OR. iprec > 2) iprec = 3
 WRITE  (outpt,60) swm,dosi(iprec),refus(iprec),sbnm,dosi(iprec1)
 60 FORMAT (a27,' 2163, REQUESTED ',a4,'LE PRECISION ',a3,' USED BY ',  &
     2A4,2H. ,a4,'LE PRECISION IS LOGICAL CHOICE')
 IF (iprec /= 3) iprec1 = iprec
 70 iprec = iprec1
 ip1   = iprec1
 is1   = ksign
 ltype = iprec1
 IF (il(5) == 3 .OR. il(5) == 4 .OR. iu(5) == 3 .OR. iu(5) == 4 .Or
!WKBR spr 93014  1   .IL(5).EQ.3 .OR. IL(5).EQ.4)  LTYPE = IPREC1 + 2  &
 .ib(5) == 3 .OR. ib(5) == 4)  ltype = iprec1 + 2
 IF (itype == 0 .OR. itype == ltype) GO TO 90
 jj    = 1
 IF (itype < 1 .OR. itype > 4) jj = 3
 WRITE  (outpt,80) swm,itype,refus(jj),sbnm,ltype
 80 FORMAT (a27,' 2164, REQUESTED TYPE ',i4,2H, ,a3,' USED BY ',2A4,  &
     '. TYPE ',i4,' IS LOGICAL CHOICE.')
 IF (jj /= 3) ltype = itype
 90 itype = ltype
 ix(5) = itype
 ix(1) = x
 iscr  = scr
 CALL fbs (z,z)
 ix(3) = n
 ix(4) = 2
 IF (ix(3) == ix(2)) ix(4) = 1
 CALL wrttrl (ix)
 GO TO 200
 
 100 CALL fname (il(1),il(2))
 WRITE  (outpt,110) il(2),il(3),il(4)
 110 FORMAT ('0*** INPUT MATRIX ',2A4,' TO FBS MODULE IS NOT A LOWER ',  &
     'TRIANGULAR FACTOR.  FORM =',i4)
 CALL errtrc ('DFBS    ',110)
 GO TO 200
 
!     SET UP CALL TO GFBS
 
 150 jl(1) = l
 CALL rdtrl (jl)
 n     = jl(2)
 jb(1) = b
 CALL rdtrl (jb)
 jnx   = korsz(zz)
 iprec1= MAX0(jl(5),jb(5),ju(5))
 IF (iprec1 > 2) iprec1 = iprec1 - 2
 IF (iprec1 < 1 .OR. iprec1 > 2) iprec1 = kprec
 IF (iprec == iprec1 .OR. iprec == 0) GO TO 160
 IF (iprec < 1 .OR. iprec > 2) iprec = 3
 WRITE (outpt,60) swm,dosi(iprec),refus(iprec),sbnm,dosi(iprec1)
 IF (iprec /= 3) iprec1 = iprec
 160 iprec = iprec1
 jp1   = iprec1
 js1   = ksign
 jx(1) = x
 ltype = iprec1
 IF (jl(5) == 3 .OR. jl(5) == 4 .OR. ju(5) == 3 .OR. ju(5) == 4 .Or  &
     .jl(5) == 3 .OR. jl(5) == 4) ltype = iprec1 + 2
 IF (itype == 0 .OR. itype == ltype) GO TO 170
 jj    = 1
 IF (itype < 1 .OR. itype > 4) jj = 3
 WRITE (outpt,80) swm,itype,refus(jj),sbnm,ltype
 IF (jj /= 3) ltype = itype
 170 itype = ltype
 jx(5) = itype
 CALL gfbs (zz,zz)
 jx(3) = n
 jx(4) = 2
 IF (jx(3) == jx(2)) jx(4) =  1
 CALL wrttrl (jx)
 
 200 RETURN
END SUBROUTINE dfbs
