SUBROUTINE rand6(xycb,buffer,npoint,iz,INPUT,lcore)
     
!     ANALYSIS OF REQUESTS AND BUILDS LIST
 
 
 INTEGER, INTENT(IN)                      :: xycb
 INTEGER, INTENT(IN OUT)                  :: buffer(1)
 INTEGER, INTENT(OUT)                     :: npoint
 INTEGER, INTENT(IN OUT)                  :: iz(1)
 INTEGER, INTENT(IN)                      :: INPUT
 INTEGER, INTENT(IN OUT)                  :: lcore
 INTEGER :: FILE,NAME(2),ilist(6),psdf,auto, itype(13,5)
 DATA NAME,psdf,auto/4HRAND,4H6   ,2,3/
 DATA itype / 13,4HDISP,1,4HVELO,2,4HACCE,3,4HDISP,8,4HVELO,9,4HACCE,10,  &
     3,4HLOAD,5,           10*0, 3,4HSPCF,4,           10*0,  &
     3,4HSTRE,6,           10*0, 3,4HELFO,7,           10*0 /
! *****
!     XYCB     XY OUTPUT REQUESTS
!     BUFFER   SYSTEM BUFFER
!     NPOINT   NUMBER OF POINTS REQUESTED FOR THIS FILE
!     IZ       LIST OF REQUESTS
!     INPUT    CURRENT FILE
!     ILIST    LIST OF REQUEST FROM XYCB   6  WORDS PER
!     SUBC,FILE,ID,COMP,OPER,DEST
!     PSDF     KEY FOR POWER SPECTRAL DENSITY FUNCTION
!     AUTO     KEY FOR AUTOCORRELATION FUNCTION
!     ITYPE    LIST OF DATA TYPES ON EACH INPUT FILE
!     IREQ     PSDF =1 , AUTO =2  BOTH = 3
!     IP       POINTER INTO  IZ  FOR LAST POINT(SAME POINT MAY OCCUR
!                MANY TIMES IN XYCB
 
!     LIST FORMAT
!     FILE,ID,COMP,IREQ,DEST
 
 
 
 
 
 
!     FIND  ACCEPTABLE MNEUMONICS
 
 k = INPUT -103
 ntype =  itype(1,k)
 ip =-4
 npoint = 0
 
!     OPEN XYCB
 
 FILE =xycb
 CALL OPEN(*90,xycb,buffer(1),0)
 CALL fwdrec(*910,xycb)
 
!     SKIP PROSE RECORD
 
 CALL fwdrec(*40,xycb)
 
!     READ DATA RECORD 6 WORDS AT A TIME
 
 5 CALL READ(*40,*40,xycb,ilist(1),6,0,i)
 
!     IS DATA BLOCK PROPER
 
 DO  i=2, ntype,2
   IF(ilist(2) == itype(i+1,k)) GO TO 20
 END DO
 
!     GO TO NEXT REQUEST
 
 GO TO 5
 
!     CHECK FOR RANDOM REQUEST
 
 20 IF(ilist(5) == psdf) GO TO 25
 IF(ilist(5) == auto) GO TO 30
 GO TO 5
!     PSDF REQUEST
 
 25 ireq =1
 GO TO 31
 
!     AUTOCORRELATION REQUEST
 
 30 ireq =2
 
!     STORE  IN LIST
 
 31 IF(npoint == 0) GO TO 35
 
!     IS THIS A NEW POINT
 
 IF(iz(ip) /= itype(i,k)) GO TO 35
 IF(iz(ip+1) /= ilist(3) .OR. iz(ip+2) /= ilist(4)) GO TO 35
 
!     ANOTHER REQUEST FOR SAME POINT
 
 IF( iz(ip+3) == 3 .OR. iz(ip+3) == ireq) GO TO 32
 iz(ip+3) = iz(ip+3) + ireq
 32 IF(iz(ip+4) ==  3  .OR. iz(ip+4) == ilist(6))GO TO 5
 iz(ip+4) = iz(ip+4) + ilist(6)
 GO TO 5
 
!     ADD POINT TO LIST
 
 35 npoint = npoint +1
 ip = ip +5
 IF (ip+5 > lcore) GO TO 905
 iz(ip) = itype(i,k)
 iz(ip+1) = ilist(3)
 iz(ip+2) = ilist(4)
 iz(ip+3) = ireq
 iz(ip+4) = ilist(6)
 GO TO 5
 
!     GET OUT
 
 40 CALL CLOSE(xycb,1)
 
!     SAVE ORIGINAL COMPONENT IN THE FIFTH LIST LORD
 
 IF(npoint == 0) GO TO 90
 DO  i = 1,npoint
   l = (i-1)*5+1
   iz(l+4) = iz(l+2)
   IF(k < 4) iz(l+2) = iz(l+2) -2
 END DO
 90 RETURN
 
!     FILE ERRORS
 
 905 npoint = npoint+9
 GO TO 90
 910 ip1= -2
 911 CALL mesage(ip1,FILE,NAME(1))
 GO TO 911
END SUBROUTINE rand6
