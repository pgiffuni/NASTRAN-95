SUBROUTINE modac2(     nv,inp1,iout)
     
!     MODAC2  REDUCES THE SIZE OF INP1 (BY REMOVING SELECTED COLUMNS)
 
!     CORE IS LAIDED OUT AS FOLLOWS
 
!         CONTENTS            LENGTH  TYPE   POINTER
!         --------            ------  ----   -------
 
!         NEW TIMES           NFN      R     IFN
!         KEEP/REMOVE         NFO      I     IKR
!         COPIED COLUMN       MCB(3)   R     ICOL
 
!         2  BUFFERS          SYSBUF   I     IBUF1
!                             SYSBUF   I     IBUF2
 
!     VARIABLES
 
!     NV       NUMBER OF COLUMS TO PROCESS TOGETHER (MINUS SAYS ADD HEAD
!     INP1     COPY FROM THIS FILE
!     IOUT     COPY TO  THIS  FILE
 
 
 
 
 INTEGER, INTENT(IN)                      :: nv
 INTEGER, INTENT(IN)                      :: inp1
 INTEGER, INTENT(IN)                      :: iout
 INTEGER :: iz,sysbuf,NAME(2),ihd(2),mcb(7),FILE
 REAL :: z(1)
 COMMON /unpakx/itc,ii,jj,incr
 COMMON /system/ sysbuf
 COMMON /modac3/ nfo,nfn,nz
 COMMON /zzzzzz/ iz(1)
 EQUIVALENCE (z(1),iz(1))
 DATA  NAME /4HMODA,4HC2  /
 
!     ALLOCATE CORE
 
 mcb(1) =iout
 CALL rdtrl(mcb)
 IF ( mcb(1) <= 0) RETURN
 mcb(1) =  inp1
 CALL rdtrl(mcb)
 IF (mcb(1)  <= 0) RETURN
 nload = mcb(2)/(nfo*IABS(nv))
 ifn =1
 ikr = ifn + nfn
 icol = ikr + nfo
 ibuf1 = nz -sysbuf+1
 ibuf2 = ibuf1- sysbuf
 IF ( icol + mcb(3) + 2*sysbuf > nz) CALL mesage(-8,0,NAME)
 
!     OPEN  FILES
 
 FILE = inp1
 CALL gopen(inp1,iz(ibuf1),0)
 FILE = iout
 CALL OPEN(*900,iout,iz(ibuf2),1)
 CALL fname(iout,ihd)
 CALL WRITE(iout,ihd,2,0)
 IF ( nv  > 0) GO  TO  10
 CALL  WRITE(iout,z,nfn,0)
 10 CALL  WRITE(iout,0,0,1)
 
!     SET UP MATRIX TRAILER
 
 FILE = inp1
 mcb(2) =0
 mcb(6) =0
 mcb(7) =0
 mcb(1) = iout
 itc = mcb(5)
 incr = 1
 inv = IABS(nv)
 DO  m = 1,nload
   k = ikr -1
   DO    i =1,nfo
     k  =k+1
     IF( iz(k) == 0)  GO TO 20
     
!     KEEP COLUMN
     
     CALL cyct2b(inp1,iout,inv,iz(icol),mcb)
     CYCLE
     
!     SKIP COLUMN
     
     20 DO  j = 1,inv
       CALL fwdrec(*910,inp1)
     END DO
   END DO
 END DO
 
!     CLOSE  UP
 
 CALL CLOSE(inp1,1)
 CALL CLOSE(iout,1)
 CALL wrttrl(mcb)
 RETURN
 
!     ERROR MESSAGES
 
 900 ip1= -1
 901 CALL mesage(ip1,FILE,NAME)
 910 ip1 = -2
 GO TO 901
END SUBROUTINE modac2
