SUBROUTINE emsg(nchar,no,isys,iwf,itext)
     
 INTEGER, INTENT(IN OUT)                  :: nchar
 INTEGER, INTENT(IN OUT)                  :: no
 INTEGER, INTENT(IN OUT)                  :: isys
 INTEGER, INTENT(IN)                      :: iwf
 INTEGER, INTENT(IN OUT)                  :: itext(1)
 INTEGER :: imsg(2,4)
 COMMON /system/sysbuf(41)
!     ISYS = 1   USER           IWF  =     1   WARNING
!          = 2   SYSTEM             =     2   FATAL
 
 EQUIVALENCE (ncpw,sysbuf(41)),(nout,sysbuf(2)),(imach,sysbuf(22))
 DATA imsg /4HUSER,1H ,4HSYST,4HEM  ,4HWARN,4HING ,4HFATA,1HL  /
 
 nword = (nchar + ncpw-1)/ncpw
 nline = (nchar + 9 + 131)/ 132   +2
 CALL page2(-nline)
 no1=IABS(no)
 k = iwf +2
 WRITE(nout,10)(imsg(i,isys),i=1,2),(imsg(m,k),m=1,2),no1
 10 FORMAT(1H0,4H*** ,4A4, i4,1H,)
 IF(nchar == 0) RETURN
 SELECT CASE ( imach )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 20
   CASE (    4)
     GO TO 50
   CASE (    5)
     GO TO 30
 END SELECT
 
!     7094
 
 20 WRITE (nout,25)(itext(i),i=1,nword)
 25 FORMAT(10X, 20A6,a2)
 GO TO 60
 
!     360/370
 
 30 WRITE(nout,35) (itext(i),i=1,nword)
 35 FORMAT(10X, 30A4,a2)
 GO TO 60
 
!     CDC
 
 50 WRITE(nout,55) (itext(i),i=1,nword)
 55 FORMAT(10X,12A10,a2)
 60 IF (no < 0) CALL mesage(-61,0,0)
 RETURN
END SUBROUTINE emsg
