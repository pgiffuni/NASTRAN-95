SUBROUTINE cpyfil (infile,oufile,area,n,count)
     
!     CPYFIL COPIES RECORDS FROM INFILE TO OUFILE UNTIL AN END-OF-FILE
!     ON INFILE IS ENCOUNTERED. AN END-OF-FILE IS NOT WRITTEN ON OUFILE.
!     RECTYPE IS CALLED PRIOR TO COPYING EACH RECORD. STRING RECORDS ARE
!     COPIED USING CPYSTR. NORMAL RECORDS ARE COPIED USING READ/WRITE.
!     UPON EXIT, INFILE IS POSITIONED IMMEDIATELY AFTER THE END-OF-FILE
!     AND OUFILE IS POSITIONED AFTER THE LAST RECORD WRITTEN.
 
!     THIS ROUTINE DOES NOT OPEN NOR CLOSE ANY FILE, NOR WRITE ANY
!     MATRIX TRAILER
 
 
 INTEGER, INTENT(IN)                      :: infile
 INTEGER, INTENT(IN)                      :: oufile
 INTEGER, INTENT(IN)                      :: area(2)
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(OUT)                     :: count
 INTEGER :: inblk(15),oublk(15), TYPE,eor
 COMMON /BLANK/ iparam
 
!     INITIALIZE STRING COMMUNICATION BLOCKS AND DETERMINE RECORD TYPE
 
 inblk(1) = infile
 oublk(1) = oufile
 count    = 0
 10 CALL rectyp (infile,TYPE)
 IF (TYPE /= 0) GO TO 50
 
!     COPY NORMAL RECORD
 
 20 eor  = 1
 CALL READ (*60,*30,infile,area,n,0,nwds)
 eor  = 0
 nwds = n
 30 IF (count /= 0 .OR. iparam /= -1111) GO TO 40
 CALL fname (infile,area(nwds+1))
 IF (area(nwds+1) == area(1) .AND. area(nwds+2) == area(2))  &
     CALL fname (oufile,area)
 40 CALL WRITE (oufile,area,nwds,eor)
 count = count + nwds
 IF (eor == 0.0) THEN
   GO TO    20
 ELSE
   GO TO    10
 END IF
 
!     COPY STRING RECORDS
 
 50 CALL cpystr (inblk,oublk,0,0)
 count = count + oublk(13)
 GO TO 10
 
!     RETURN WHEN END-OF-FILE IS ENCOUNTERED
 
 60 RETURN
END SUBROUTINE cpyfil
