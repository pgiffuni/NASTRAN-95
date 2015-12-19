SUBROUTINE sgena (TYPE,buf,mcb,ifile,icode,iextra,ofile,ocode,  &
        oextra)
     
!     THIS ROUTINE READS SUBSTRUCTURING CONSTRAINT AND DYNAMIC PROPERTY
!     CARDS AND CONVERTS THEM TO NASTRAN FORMAT
 
!     INPUTS -
 
!     TYPE   - BCD CARD NAME
!     BUF    - GINO BUFFER FOR INPUT FILE
!     MCB    - MATRIX CONTROL BLOCK FOR INPUT FILE
!     IFILE  - INPUT FILE NAME
!     ICODE  - LOCATE CODE FOR INPUT CARD TYPE
!     IEXTRA - NUMBER OF EXTRA WORDS (AFTER GRID) TO BE READ
!     OFILE  - OUTPUT FILE NAME
!     OCODE  - LOCATE CODE FOR OUTPUT CARD TYPE
!     OEXTRA - NUMBER OF EXTRA WORDS (AFTER GRID) TO BE WRITTEN
 
 
 INTEGER, INTENT(IN)                      :: TYPE(2)
 INTEGER, INTENT(IN OUT)                  :: buf(1)
 INTEGER, INTENT(OUT)                     :: mcb(7)
 INTEGER, INTENT(IN OUT)                  :: ifile
 INTEGER, INTENT(OUT)                     :: icode(4)
 INTEGER, INTENT(IN OUT)                  :: iextra
 INTEGER, INTENT(IN OUT)                  :: ofile
 INTEGER, INTENT(IN)                      :: ocode(4)
 INTEGER, INTENT(IN)                      :: oextra
 EXTERNAL        andf,complf,orf
 INTEGER :: z,sysbuf,outt,two,subnam(2),card(20),comp,  &
     cin(6),code,cexist(6),andf,complf,orf
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / idry,NAME(2)
 COMMON /sgencm/ nono,nss,iptr
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,outt
 COMMON /two   / two(32)
 DATA    subnam/ 4HSGEN,4HA    /
 
!     LOCATE CARDS ON FILE
 
 CALL locate (*200,buf(1),icode(1),icd)
 icode(4) = 1
 
!     WRITE HEADER RECORD ON OUTPUT FILE
 
 CALL WRITE (ofile,icode(1),3,0)
 
!     READ SID AND SUBSTRUCTURING NAME FROM CARD
 
 10 CALL READ (*1002,*150,ifile,card,3,0,nwds)
 card(4) = card(1)
 n = 6 + oextra
 DO  i = 5,n
   card(i) = 0
 END DO
 
!     FIND SUBSTRUCTURE
 
 DO  i = 1,nss
   inam = 2*i + 3
   IF (z(inam) == card(2) .AND. z(inam+1) == card(3)) GO TO 50
 END DO
 
!     SUBSTRUCTURE NOT FOUND - SKIP OVER DATA
 
 CALL page2 (-4)
 WRITE (outt,63290) uwm,(card(j),j=2,3),TYPE,NAME
 40 CALL fread (ifile,card,2+iextra,0)
 IF (card(1) < 0.0) THEN
   GO TO    10
 ELSE
   GO TO    40
 END IF
 
!     FOUND SUBSTRUCTURE NAME
 
 50 ipt  = iptr + i - 1
 igrd = z(ipt)
 ngrd = (z(ipt+1) - z(ipt))/3
 
!     PROCESS GRID-COMPONENT PAIRS
 
 60 CALL fread (ifile,card(5),2+iextra,0)
 igrid = card(5)
 IF (igrid == -1) GO TO 10
 IF (igrid ==  0) GO TO 60
 comp = card(6)
 IF (comp == 0) comp = 1
 card(6) = 0
 CALL bisloc (*80,igrid,z(igrd),3,ngrd,igr)
 ig = igr + igrd - 1
 npro = 0
 70 IF (z(ig-3) /= z(ig)) GO TO 90
 IF (ig <= igrd) GO TO 90
 ig = ig - 3
 GO TO 70
 
!     BAD GRID
 
 80 nono = 1
 CALL page2 (-3)
 WRITE (outt,60220) ufm,(card(j),j=2,3),igrid,comp,TYPE,NAME
 GO TO 60
 
!     SPLIT COMPONENTS
 
 90 CALL splt10 (comp,cin,ncin)
 100 code = z(ig+2)
 IF (code == 0) code = 1
 isil = z(ig+1)
 CALL decode (code,cexist,nc)
 
!     FIND ACTUAL REMAINING COMPONENTS AND WRITE CONVERTED DATA TO
!     OUTPUT FILE
 
 DO  j  = 1,nc
   DO  jg = 1,ncin
     IF (cin(jg)-cexist(j)-1 == 0.0) THEN
       GO TO   110
     ELSE
       GO TO   120
     END IF
     110 npro = npro + 1
     card(5) = isil + j - 1
     CALL WRITE (ofile,card(4),3+oextra,0)
     120 CONTINUE
   END DO
 END DO
 IF (npro >= ncin) GO TO 60
 IF (z(ig+3) /= z(ig)) GO TO 80
 IF ((ig+3) >= (igrd+3*ngrd)) GO TO 80
 ig = ig + 3
 GO TO 100
 
!     FINISH PROCESSING CARDS BY CLOSING OUTPUT FILE RECORD
 
 150 CALL WRITE (ofile,0,0,1)
 
!     TURN OFF TRAILER FOR INPUT CARD TYPE
 
 j = (icode(2) - 1)/16
 i = icode(2) - 16*j
 mcb(j+2) = andf(complf(two(i+16)),mcb(j+2))
 
!     TURN ON TRAILER FOR OUTPUT CARD TYPE
 
 j = (ocode(2) - 1)/16
 i = ocode(2) - 16*j
 mcb(j+2) = orf(two(i+16),mcb(j+2))
 
!     RETURN
 
 200 RETURN
 
!     ERRORS
 
 1002 CALL mesage (-2,ifile,subnam)
 RETURN
 60220 FORMAT (a23,' 6022, SUBSTRUCTURE ',2A4,', GRID POINT',i9,  &
     ', COMPONENTS',i9,1H, /30X,'REFERENCED ON ',2A4,  &
 ' CARD, DO NOT EXIST ON SOLUTION STRUCTURE ',2A4)
   63290 FORMAT (a25,' 6329, SUBSTRUCTURE ',2A4,' REFERENCED ON ',2A4,  &
       ' CARD', /30X,'IS NOT A COMPONENT BASIC SUBSTRUCTURE OF ',  &
       'SOLUTION STRUCTURE ',2A4,/30X,'THIS CARD WILL BE IGNORED')
 END SUBROUTINE sgena
