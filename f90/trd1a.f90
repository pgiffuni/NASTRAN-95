SUBROUTINE trd1a (casexx,trl,ic,nlftp,ngroup,moda1)
     
!     THIS ROUTINE BUILDS THE INITIAL CONDITIONS TABLE, PUTS TSTEP STUFF
!     IN CORE AND EXTRACTS THE NLFTP POINTER
 
!     THIS ROUTINE IS SUITABLE FOR SINGLE PRECISION OPERATION
 
 
 INTEGER, INTENT(IN OUT)                  :: casexx
 INTEGER, INTENT(IN)                      :: trl
 INTEGER, INTENT(IN OUT)                  :: ic
 INTEGER, INTENT(OUT)                     :: nlftp
 INTEGER, INTENT(OUT)                     :: ngroup
 INTEGER, INTENT(IN)                      :: moda1
 INTEGER :: sysbuf    , intrl(2) , iz(160)   ,mcb(7)   ,FILE     ,NAME(2)
 
 COMMON   /system/  sysbuf
 COMMON   /packx /  it1       ,it2      ,ii       ,jj       , incr
 COMMON   /zzzzzz / z(1)
 
 EQUIVALENCE        (z(1)     ,iz(1))
 
 DATA               NAME                ,intrl              /  &
     4HTRD1    ,4HA      ,4HTRL    ,4HTRD    /
 
!     IDENTIFICATION VARIABLES
 
!     NGROUP        NUMBER OF CHANGES OF TIME STEP
 
!     ITSTEP        SELECTED TSTEP ID
 
!     NLFTP         SELECTED NON-LINEAR LOAD ID
 
!     ICP           SELECTED INITIAL CONDITION ID
 
!     LUD           LENGTH OF INITIAL CONDITION--D SET
 
!     IGROUP        POINTER TO TSTEP STUFF
 
 
 
!     INITIALIZE
 
 it1 = 1
 it2 = 1
 ii  = 1
 incr= 1
 nz  = korsz (z)
 nx  = nz
 
!     PICK UP AND STORE CASECC POINTERS
 
 ibuf1 = nz -sysbuf +1
 CALL gopen (casexx,iz(ibuf1),0)
 CALL fread (casexx,iz,166,1)
 CALL CLOSE (casexx,1)
 itstep = iz(38)
 icp    = iz(9)
 nlftp  = iz(160)
 IF (icp /= 0 .AND. moda1 == 1) GO TO 920
 
!     BUILD INITIAL CONDITION FILE
 
 CALL gopen (ic,iz(ibuf1),1)
 ibuf2 = ibuf1-sysbuf
 nz    = nz - 2*sysbuf
 icrq  =-nz
 IF (icrq > 0) GO TO 980
 FILE = trl
 CALL OPEN (*900,trl,iz(ibuf2),0)
 CALL READ (*910,*10,trl,iz(1),nz,0,iflag)
 icrq = nz
 GO TO 980
 10 lud  = iz(iflag)
 jj   = lud
 icrq = 2*lud - nz
 IF (icrq > 0) GO TO 980
 l    = iz(3)
 itrl = l
 
!     ZERO I. C.
 
 ivel  = ibuf2- lud-1
 idisp = ivel -lud
 DO  i = 1,lud
   k = ivel +i
   z(k) = 0.0
   k = idisp +i
   z(k) = 0.0
 END DO
 CALL makmcb (mcb,ic,lud,2,1)
 IF (icp   == 0) GO TO 80
 IF (iz(3) == 0) GO TO 40
 iflag = iflag-1
 DO  i = 4,iflag
   IF (iz(i) == icp) GO TO 50
 END DO
 40 itstep = icp
 GO TO 940
 50 k = i-4
 l = iflag -i
 CALL skprec (trl,k)
 70 CALL READ (*910,*80,trl,iz(1),3,0,iflag)
 k    = iz(1) +idisp
 i2   = 2
 z(k) = z(k) + z(i2)
 k    = iz(1) + ivel
 z(k) = z(k) + z(i2+1)
 GO TO 70
 80 CALL pack (z(idisp+1),ic,mcb)
 CALL pack (z(ivel +1),ic,mcb)
 CALL CLOSE (ic,1)
 CALL wrttrl (mcb)
 CALL skprec (trl,l)
 
!     BRING TSTEP STUFF INTO CORE
 
 100 itrl = itrl +1
 CALL READ (*940,*110,trl,iz(1),nz,0,iflag)
 icrq = nz
 GO TO 980
 110 IF (iz(1) /= itstep) GO TO 100
 
!     TSTEP CARD FOUND
 
 CALL CLOSE (trl,1)
 ngroup = (iflag-1)/3
 
!     MOVE TSTEP STUFF TO BOTTOM OF CORE
 
 nz = nx - iflag +1
 igroup = nz+1
 DO  i = 2,iflag
   k = igroup +i-2
   iz(k) = iz(i)
 END DO
 RETURN
 
!     ERROR MESSAGES
 
 900 ip1 = -1
 901 CALL mesage (ip1,FILE,NAME)
 RETURN
 910 ip1 = -2
 GO TO 901
 920 ip1  = -51
 FILE = icp
 GO TO 901
 940 CALL mesage (-31,itstep,intrl)
 RETURN
 980 ip1 = -8
 FILE = icrq
 GO TO 901
END SUBROUTINE trd1a
