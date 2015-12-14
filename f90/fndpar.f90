SUBROUTINE fndpar (np2,INDEX)
     
!     FNDPAR FINDS THE INDEX INTO THE  VPS FOR PARAMETER NUMBER NP
!     IN THE CURRENT OSCAR (THIS PARAMETER MUST BE VARIABLE)
 
 
 INTEGER, INTENT(IN OUT)                  :: np2
 INTEGER, INTENT(OUT)                     :: INDEX
 EXTERNAL        andf
 INTEGER :: oscar,NAME(2),andf
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /oscent/ oscar(7)
 COMMON /system/ sysbuf,nout
 COMMON /sem   / mask,mask2,mask3
 DATA    NAME  / 4HFNDP,4HAR  /
 
 nip  = oscar(7)
 itype= andf(oscar(3),7)
 i    = 8 + 3*nip
 IF (itype == 2) GO TO 100
 nop  = oscar(i)
 i    = i + 3*nop + 1
 100 CONTINUE
 i    = i + 1
 np1  = oscar(i)
 np   = IABS(np2)
 IF (np <= np1) GO TO 120
 IF (np2 <=  0) GO TO 200
 WRITE  (nout,110) ufm,np
 110 FORMAT (a23,' 3123, PARAMETER NUMBER',i6,' NOT IN DMAP CALL.')
 CALL mesage (-61,0,NAME)
 120 CONTINUE
 np1 = np - 1
 k   = i + 1
 IF (np1 == 0) GO TO 170
 DO  i = 1,np1
   m   = oscar(k)
   IF (m < 0) THEN
     GO TO   140
   ELSE
     GO TO   150
   END IF
   
!     VARTABLE
   
   140 k = k + 1
   CYCLE
   
!     CONSTANT
   
   150 k = k + 1 + m
 END DO
 
!     K POINTS  TO WANTED OSCAR WORD
 
 170 IF (oscar(k) < 0) GO TO 190
 IF (np2 <= 0) GO TO 200
 WRITE  (nout,180) ufm,np
 180 FORMAT (a23,' 3124, PARAMETER NUMBER',i6,' IS NOT A VARIABLE.')
 CALL mesage (-61,0,NAME)
 190 INDEX = andf(oscar(k),mask3)
 RETURN
 
!     PARAMETER SPORT NOT SUPPLIES
 
 200 INDEX = -1
 RETURN
END SUBROUTINE fndpar
