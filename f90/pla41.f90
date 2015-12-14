SUBROUTINE pla41
     
!     THIS ROUTINE APPENDS DISPLACEMENT VECTOR INFORMATION TO THE
!     ECPTNL DATA BLOCK AND CREATES A SCRATCH DATA BLOCK, ECPTS, OF
!     THIS MERGED INFORMATION.  ECPTS IS PROCESSED BY SUBROUTINE PLA41.
 
 INTEGER :: sysbuf,clsrw,buffr1,buffr2,ugv,ecptnl,ecpts, eor,FILE,eltype
 DIMENSION       NAME(2),mcbugv(7),nwords(40),ngpts(40), xecpt(100),iecpt(100)
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ z(1)
 COMMON /unpakx/ itypeb,iunpk,junpk,incupk
 EQUIVALENCE     (xecpt(1),iecpt(1))
 DATA    ugv   , ecptnl,ecpts / 106,103,301   /
 DATA    NAME  / 4HPLA4,4H1   /
 DATA    eor   , neor,clsrw   / 1,0,1         /
 
!    1        ROD       BEAM      TUBE      SHEAR     TWIST
!    2        TRIA1     TRBSC     TRPLT     TRMEM     CONROD
!    3        ELAS1     ELAS2     ELAS3     ELAS4     QDPLT
!    4        QDMEM     TRIA2     QUAD2     QUAD1     DAMP1
!    5        DAMP2     DAMP3     DAMP4     VISC      MASS1
!    6        MASS2     MASS3     MASS4     CONM1     CONM2
!    7        PLOTEL    REACT     QUAD3     BAR       CONE
!    8          X         X         X         X         X
 
 DATA   nwords/ 20,        0,       19,        0,        0,  &
     33,        0,        0,       27,       20,  &
     0,        0,        0,        0,        0,  &
     32,       27,       32,       38,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,       45,        0,  &
     0,        0,        0,        0,        0 /
 DATA   ngpts / 2,        2,        2,        4,        4,  &
     3,        3,        3,        3,        2,  &
     2,        2,        2,        2,        4,  &
     4,        3,        4,        4,        2,  &
     2,        2,        2,        2,        2,  &
     2,        2,        2,        2,        2,  &
     2,        0,        0,        2,        2,  &
     0,        0,        0,        0,        0 /
 
!     INITIALIZE
 
 izmax  = korsz(z)
 buffr1 = izmax  - sysbuf
 buffr2 = buffr1 - sysbuf
 left   = buffr2 - 1
 
!     READ THE DISPLACEMENT VECTOR INTO OPEN CORE.
 
 FILE = ugv
 CALL gopen (ugv,z(buffr1),0)
 mcbugv(1) = ugv
 CALL rdtrl (mcbugv)
 IF (left < mcbugv(3)) CALL mesage (-8,0,NAME)
 itypeb = 1
 iunpk  = 1
 junpk  = mcbugv(3)
 incupk = 1
 CALL unpack (*9050,ugv,z(1))
 CALL CLOSE  (ugv,clsrw)
 
!     OPEN THE ECPTNL AND ECPTS FILES.
 
 CALL gopen (ecpts ,z(buffr1),1)
 CALL gopen (ecptnl,z(buffr2),0)
 
!     READ AND WRITE THE PIVOT POINT
 
 10 CALL READ  (*60,*9030,ecptnl,npvt,1,neor,iflag)
 CALL WRITE (ecpts,npvt,1,neor)
 
!     READ ELEMENT TYPE
 
 20 CALL READ (*9020,*50,ecptnl,eltype,1,neor,iflag)
 j = nwords(eltype)
 IF (j <= 0) CALL mesage (-30,114,iecpt(1))
 
!     READ THE ECPT ENTRY FOR THIS ELEMENT.
 
 CALL fread (ecptnl,xecpt,j,0)
 
!     APPEND DISPLACEMENT VECTOR TO THE ECPT ENTRY
 
 j = j + 1
 nwds = 3
 IF (eltype == 34) nwds = 6
 nogpts = ngpts(eltype)
 DO  i = 1,nogpts
   INDEX = iecpt(i+1)
   DO  k = 1,nwds
     xecpt(j) = z(INDEX)
     INDEX = INDEX + 1
     j = j + 1
   END DO
 END DO
 
!     THE ECPT ENTRY IS NOW COMPLETE.  WRITE IT OUT.
 
 CALL WRITE (ecpts,eltype, 1,neor)
 CALL WRITE (ecpts,xecpt,j-1,neor)
 GO TO 20
 
!     AN EOR HAS BEEN READ ON ECPTNL.  WRITE EOR ON ECPTS.
 
 50 CALL WRITE (ecpts,0,0,eor)
 GO TO 10
 
!     PROCESSING IS COMPLETE.  CLOSE FILES.
 
 60 CALL CLOSE (ecptnl,clsrw)
 CALL CLOSE (ecpts,clsrw)
 RETURN
 
!     FATAL ERRORS
 
 9020 CALL mesage (-2,FILE,NAME)
 9030 CALL mesage (-3,FILE,NAME)
 9050 CALL mesage (-30,83,NAME)
 RETURN
END SUBROUTINE pla41
