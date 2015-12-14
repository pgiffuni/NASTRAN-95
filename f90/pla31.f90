SUBROUTINE pla31
     
!     THIS ROUTINE READS THE INCREMENTAL DISPLACEMENT VECTOR INTO CORE
!     AND APPENDS THE PROPER DISPLACEMENT VALUES TO THE ESTNL ENTRY FOR
!     EACH ELEMENT, THEREBY CREATING THE ESTNLS, THE ESTNL SCRATCH FILE,
!     WHICH IS PROCESSED BY SUBROUTINE PLA32.
 
 INTEGER :: bufsz,bufr1,bufr2,delugv,estnl,estnls,FILE,eor,  &
     clsrw,iz(1),iestbk(100),estwds(40),eltype
 DIMENSION       NAME(2),ngpts(40),mcbugv(7),estbk(100)
 COMMON /BLANK / icom
 COMMON /system/ bufsz
 COMMON /zzzzzz/ z(1)
 COMMON /unpakx/ itypeb,iunpk,junpk,incupk
 EQUIVALENCE     (z(1),iz(1)),(estbk(1),iestbk(1))
 DATA    NAME  / 4HPLA3,4H1   /
 DATA    delugv, estnl,estnls / 104,105,301/
 DATA    eor   , neor,clsrw   / 1,0,1      /
 
!    1        ROD       BEAM      TUBE      SHEAR     TWIST
!    2        TRIA1     TRBSC     TRPLT     TRMEM     CONROD
!    3        ELAS1     ELAS2     ELAS3     ELAS4     QDPLT
!    4        QDMEM     TRIA2     QUAD2     QUAD1     DAMP1
!    5        DAMP2     DAMP3     DAMP4     VISC      MASS1
!    6        MASS2     MASS3     MASS4     CONM1     CONM2
!    7        PLOTEL    REACT     QUAD3     BAR       CONE
!    8        TRIARG    TRAPRG    TORDRG    CORE      CAP
 
 DATA    estwds/ 21,        0,       20,        0,        0,  &
     38,        0,        0,       27,       21,  &
     0,        0,        0,        0,        0,  &
     32,       32,       37,       43,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,       50,        0,  &
     0,        0,        0,        0,        0 /
 DATA    ngpts / 2,        2,        2,        4,        4,  &
     3,        3,        3,        3,        2,  &
     2,        2,        2,        2,        4,  &
     4,        3,        4,        4,        2,  &
     2,        2,        2,        2,        2,  &
     2,        2,        2,        2,        2,  &
     2,        0,        0,        2,        2,  &
     3,        4,        2,        4,        2 /
 
!     DETERMINE SIZE OF CORE, DEFINE BUFFERS AND INITIALIZE CORE
!     POINTERS AND COUNTERS
 
 izmax = korsz (z)
 bufr1 = izmax - bufsz
 bufr2 = bufr1 - bufsz
 left  = bufr2 - 1
 idisp = 0
 
!     OPEN THE DISPLACEMENT VECTOR FILE AND READ THE DISPLACEMENT VECTOR
!     INTO OPEN CORE.
 
 FILE = delugv
 CALL gopen (delugv,z(bufr1),0)
 mcbugv(1) = delugv
 CALL rdtrl (mcbugv)
 IF (left < mcbugv(3)) CALL mesage (-8,0,NAME)
 itypeb = 1
 iunpk  = 1
 junpk  = mcbugv(3)
 incupk = 1
 CALL unpack (*9040,delugv,z(idisp+1))
 CALL CLOSE (delugv,clsrw)
 
!     BUILD THE SCRATCH FILE ESTNLS
 
 CALL gopen (estnl,z(bufr1),0)
 CALL gopen (estnls,z(bufr2),1)
 
!     READ AN ELEMENT TYPE FROM ESTNL AND WRITE IT ON ESTNLS.
 
 10 CALL READ (*60,*9030,estnl,eltype,1,neor,iflag)
 nwdsrd = estwds(eltype)
 IF (nwdsrd <= 0) CALL mesage (-30,91,eltype)
 CALL WRITE (estnls,eltype,1,neor)
 
!     READ AN ESTNL ENTRY
 
 20 j = nwdsrd
 CALL READ (*9020,*50,estnl,estbk,j,neor,iflag)
 nogpts = ngpts(eltype)
 IF (nogpts <= 0) CALL mesage (-30,92,eltype)
 
!     APPEND THE DISPLACEMENT VECTORS ONTO THE ESTBK.
 
 nwds = 3
 j = j + 1
 IF (eltype == 1 .OR. eltype == 3 .OR. eltype == 10 .OR.  &
     eltype == 6 .OR. eltype == 17 .OR. eltype == 18 .OR.  &
     eltype == 19 .OR. eltype == 34) nwds = 6
 DO  i = 1,nogpts
   INDEX = idisp + iestbk(i+1)
   DO  k = 1,nwds
     estbk(j) = z(INDEX)
     INDEX = INDEX + 1
     j = j + 1
   END DO
 END DO
 
!     THE APPENDED ESTNL ENTRY, WHICH IS AT ESTBK IS NOW COMPLETE.
 
 CALL WRITE (estnls,estbk,j-1,neor)
 GO TO 20
 
!     WRITE AN EOR ON THE ESTNLS FILE.
 
 50 CALL WRITE (estnls,0,0,eor)
 GO TO 10
 
!     PROCESSING IS NOW COMPLETE
 
 60 CALL CLOSE (estnl,clsrw)
 CALL CLOSE (estnls,clsrw)
 RETURN
 
!     FATAL ERRORS
 
 9020 CALL mesage (-2,FILE,NAME)
 9030 CALL mesage (-3,FILE,NAME)
 9040 CALL mesage (-5,delugv,NAME)
 RETURN
END SUBROUTINE pla31
