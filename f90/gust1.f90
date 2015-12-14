SUBROUTINE gust1(casecc,dit,dlt,frl,pp,fol,gustl,nfreq,nload,  &
        xo,v,nogust,casnew)
     
!     THE PURPOSE OF THI ROUTINE IS TO GERATE PP,GUSTL,FOL.
 
!     THE ROUTINE PROCEEDS AS FOLLOWS
 
!         FIND  GUST CARD(NO-CARDS--SET NOGUST=1 AND RETURN)
!         PUT GUST CARDS IN CORE
!         READ CASECC -- BUILD GUSTL
!           SUPPLU DLOAD =   FROM GUST =
 
!         CALL GUST1A WITH NEW CASECC
 
 
 INTEGER, INTENT(IN)                      :: casecc
 INTEGER, INTENT(IN)                      :: dit
 INTEGER, INTENT(IN OUT)                  :: dlt
 INTEGER, INTENT(IN OUT)                  :: frl
 INTEGER, INTENT(IN OUT)                  :: pp
 INTEGER, INTENT(IN OUT)                  :: fol
 INTEGER, INTENT(IN OUT)                  :: gustl
 INTEGER, INTENT(IN OUT)                  :: nfreq
 INTEGER, INTENT(IN OUT)                  :: nload
 REAL, INTENT(OUT)                        :: xo
 REAL, INTENT(OUT)                        :: v
 INTEGER, INTENT(OUT)                     :: nogust
 INTEGER, INTENT(IN OUT)                  :: casnew
 INTEGER :: sysbuf,NAME(2), FILE,igust(2),lgust(5)
 REAL :: z(1),rgust(5)
 COMMON /system/sysbuf
 COMMON /zzzzzz/ iz(1)
 EQUIVALENCE (iz(1),z(1)),(rgust(1),lgust(1))
 DATA  NAME /4HGUST,1H1 /,igust /1005,10 /
 DATA igst /178/
 
!     INITIALIZE
 
 nz = korsz(iz)
 ibuf1 = nz-sysbuf
 ibuf2 = ibuf1-sysbuf
 ibuf3 = ibuf2-sysbuf
 nz = ibuf3-1
 nogust =-1
 nogo =0
 CALL preloc(*1000,iz(ibuf1),dit)
 CALL locate(*1000,iz(ibuf1),igust,idx)
 
!     PUT  GUST CARDS IN CORE
 
 FILE =dit
 CALL READ(*910,*10,dit,iz,nz,0,nlgust)
 CALL mesage(-8,0,NAME)
 10 CONTINUE
 CALL CLOSE(dit,1)
 icc = nlgust+1
 CALL gopen(casecc,iz(ibuf1),0)
 CALL gopen(casnew,iz(ibuf2),1)
 CALL gopen(gustl,iz(ibuf3),1)
 nz = nz - nlgust
 
!     BLAST READ A CASE CONTROL RECORD INTO CORE
 
 20 CONTINUE
 FILE = casecc
 CALL READ(*100,*30,casecc,iz(icc),nz,0,lcc)
 CALL mesage(-8,0,NAME)
 30 CONTINUE
 igsid = iz(icc+igst)
 iz(icc+12) = igsid
 CALL zeroc(rgust,5)
 IF( igsid ==  0)  GO TO 90
 
!     FIND GUST ID AMONG GUST CARDS
 
 DO  i = 1 ,nlgust,5
   IF( iz(i) == igsid) GO TO 50
 END DO
 CALL mesage(31,igsid,NAME)
 nogo =1
 GO TO 90
 
!     FOUND GUST CARD
 
 50 CONTINUE
 iz(icc+12) = iz(i+1)
 igust(1) =igsid
 lgust(2) = iz(i+1)
 rgust(3) = z(i+2)
 rgust(4) = z(i+3)
 rgust(5) = z(i+4)
 xo =  rgust(4)
 v  =  rgust(5)
 nogust = 1
 
!     PUT OUT GUSTL /CASNEW
 
 90 CALL WRITE(casnew,iz(icc),lcc,1)
 CALL WRITE(gustl,lgust,5,1)
 GO TO 20
 
!     END OF FILE ON CASECC
 100 CONTINUE
 IF( nogo == 1) CALL mesage(-61,0,NAME)
 CALL CLOSE(casecc,1)
 CALL CLOSE(gustl,1)
 CALL CLOSE(casnew,1)
 
!     CALL GUST1A FOR LOADS(W)
 
 CALL gust1a (dlt, frl, -casnew, dit, pp, 1, nfreq, nload, frqset, fol, notrd)
 CALL dmpfil(-pp,iz,nz)
 1000 CALL CLOSE(dit,1)
 RETURN
 
!     FILE  ERRORS
 
 910 ip1 = -2
 CALL mesage (ip1, FILE, NAME)
 RETURN
END SUBROUTINE gust1
