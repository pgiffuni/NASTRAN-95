SUBROUTINE optp1d (elop,pr,pl)
     
!     PROPERTY OPTIMIZER   SET POINTERS TO PLIMIT
 
 
 INTEGER, INTENT(IN)                      :: elop(2,1)
 INTEGER, INTENT(IN OUT)                  :: pr(1)
 REAL, INTENT(OUT)                        :: pl(1)
 INTEGER :: count,ycor,b1p1,scrth1, sysbuf,outtap,plp,pid,NAME(2),nkl(2)
 REAL :: kl
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / skp1(2),count,skp2(2),ycor,b1p1,npow,  &
     skp3(2),nprw,nwdsp,nklw,skp4(6),scrth1
 COMMON /optpw1/ klwds,kl(4)
 COMMON /system/ sysbuf,outtap
 COMMON /names / nrd,noeor,nwrt,nweor
 EQUIVALENCE     (nkl(1),kl(1))
 DATA    NAME  / 4H opt,4HPID  /
 
 nogo = 0
 plp  = 1
 GO TO 15
 
!     READ A NEW ELEMENT TYPE
 
 10 CALL fread (scrth1,0,0,nweor)
 15 l   = 0
 npl = 0
 CALL READ (*150,*180,scrth1,itp,1,noeor,i)
 IF (itp <= npow) GO TO 40
 
 20 CALL page2 (-2)
 WRITE  (outtap,30) sfm,NAME,itp,l
 30 FORMAT (a25,' 2301,',2A4,' FILE OPTIMIZATION PARAMETER INCORRECT',  &
     ' AS',2I8)
 nogo = nogo + 1
 GO TO 140
 
 40 ip1 = elop(2,itp)
 ip2 = elop(2,itp+1) - 1
 npr = ip2 - ip1
 IF (npr <= 0) GO TO 10
 CALL fread (scrth1,l,1,noeor)
 IF (l <= 0) GO TO 20
 
 CALL fread (scrth1,nkl(1),4,noeor)
 l = l - 1
 
!     SEQUENTIAL SEARCH ON PLIMIT AND PROPERTY DATA
!     LPL -- LAST PLIMIT POINTED TO (BY ILL).
!     NPL -- NUMBER OF PLIMIT FOR THIS ELEMENT TYPE IN CORE.
!     PLP -- POINTER FIRST PLIMIT  --    --     -- .
 
 lpl = -9877
 
 DO  ipr = ip1,ip2,nwdsp
   pid = pr(ipr)
   
   50 IF (pid-nkl(1) < 0.0) THEN
     GO TO    70
   ELSE IF (pid-nkl(1) == 0.0) THEN
     GO TO    80
   END IF
   
!     CHECK UPPER RANGE PLIMIT
   
   60 CONTINUE
   IF (pid-nkl(2) > 0.0) THEN
     GO TO    70
   ELSE
     GO TO    80
   END IF
   
!     READ NEXT PLIMIT INTO CORE
   
   70 IF (l <= 0) GO TO 140
   CALL fread (scrth1,nkl(1),4,noeor)
   l = l - 1
   GO TO 50
   
!     PLIMIT EXISTS - SEE IF MATCHES LAST
   
   80 IF (lpl == l) GO TO 120
   
!     DOESNOT - CHECK IF PREVIOUS ENTRY
   
   IF (npl == 0) GO TO 100
   DO  lpl = plp,loc,2
     IF (pl(lpl)   /= kl(3)) CYCLE
     IF (pl(lpl+1) == kl(4)) GO TO 110
   END DO
   
!     NEW PLIMIT
   
   100 IF (npl+plp+1 > ycor) GO TO 190
   npl = npl + 2
   loc = npl + plp - 2
   pl(loc  ) = kl(3)
   pl(loc+1) = kl(4)
   lpl = l
   ill = loc
   GO TO 120
   
!     PREVIOUS MATCH
   
   110 ill = lpl
   lpl = l
   
!     LOAD POINTER
   
   120 pr(ipr+5) = ill
   
 END DO
 
 140 plp = plp + npl
 GO TO 10
 
!     END-OF-FILE
 
 150 nklw = plp + npl - 1
 160 IF (nogo > 0) count = -1
 RETURN
 
!     ILLEGAL EOR
 
 180 CALL mesage (-3,scrth1,NAME)
 
!     INSUFFICIENT COREINTERNAL ELEMENT NUMBER PRINTED
 
 190 CALL page2 (-2)
 WRITE  (outtap,200) ufm,NAME,b1p1,itp
 200 FORMAT (a23,' 2298, INSUFFICIENT CORE ',2A4,1H(,i10, ' ), PROPERTY',i9)
 nklw = -plp
 GO TO 160
END SUBROUTINE optp1d
