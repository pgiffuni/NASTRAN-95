SUBROUTINE dmpmat (ifile,z      ,lz)
     
!     DUMPS A FILE ON DIAG 20 SETTING.
 
 
 INTEGER, INTENT(IN OUT)                  :: ifile
 REAL, INTENT(IN OUT)                     :: z(2)
 INTEGER, INTENT(IN)                      :: lz
 INTEGER :: sysbuf,outpt,buf,FILE,NAME(2)
 INTEGER :: itrail(7)
 
 COMMON /system/ sysbuf,outpt
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 COMMON /unpakx/ iout,irow,nrow,incr
 DATA    NAME  / 4HDMPF,2HIL  /
 
 1 FORMAT (1H0,i5,10(1X,i10,1X))
 2 FORMAT (1H ,5X,10(1X,1P,e11.4))
 3 FORMAT (1H ,5X,10(4X,a4,4X))
 
!      CALL SSWTCH (20,L)
!      IF (L .EQ. 0) RETURN
 
 FILE = IABS(ifile)
 buf  = lz - sysbuf + 1
 IF (buf < 6) GO TO 91
 lcore = (buf-1)/5
 lcore = lcore*5
 CALL OPEN (*90,FILE,z(buf),rdrew)
 WRITE  (outpt,10) FILE
 10 FORMAT (14H1DUMP of FILE ,i3)
 GO TO 100
 
 70 itrail(1) = FILE
 CALL CLOSE (FILE,clsrew)
 CALL rdtrl (itrail)
 WRITE  (outpt,80) (itrail(i),i=1,7)
 80 FORMAT (4H0EOF ,//,8H0TRAILER  ,/,7(1X,i12 /))
 90 RETURN
 
 91 CALL mesage (8,0,NAME)
 GO TO 90
 
 100 CALL READ (*70,*101,FILE,z,2,1,iwords)
 101 WRITE  (outpt,102) z(1),z(2)
 102 FORMAT (14H0HEADER record  ,/1H0,2A4)
 itrail(1) = FILE
 CALL rdtrl (itrail)
 ncols = itrail(2)
 IF (ncols > 300) ncols = 100
 iout = 1
 incr = 1
 IF (ncols > 0) THEN
   GO TO   110
 ELSE
   GO TO    70
 END IF
 110 DO  j = 1,ncols
   WRITE  (outpt,115) j
   115 FORMAT (7H0COLUMN  ,i5)
   irow = 0
   nrow = 0
   CALL unpack (*140,FILE,z)
   WRITE  (outpt,118) irow,nrow
   118 FORMAT (1H+,20X,3HROW  ,i4,11H   thru row   ,i5)
   IF (nrow > 300) nrow = 100
   nels = nrow - irow + 1
   IF (nels <= 0) CYCLE
   WRITE  (outpt,119) (z(k),k=1,nels)
   WRITE  (outpt,1119) (z(k),k=1,nels)
   1119 FORMAT( 10Z9)
   119 FORMAT (1P,10E13.4)
   CYCLE
   140 WRITE  (outpt,141)
   141 FORMAT (13H null column  )
 END DO
 GO TO 70
 
END SUBROUTINE dmpmat
