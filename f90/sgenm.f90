SUBROUTINE sgenm (ntype,ifile,sfile,ofile,icode,ocode,ctypes,  &
        ctypeo)
     
!     THIS SUBROUTINE MERGES CONVERTED SUBSTRUCTURING DATA WITH EXISTING
!     NASTRAN DATA
 
!     INPUTS
!     NTYPE  - NUMBER OF DIFFERENT SUBSTRUCTURING CARDS
!     IFILE  - INPUT FILE NAME
!     SFILE  - SCRATCH FILE NAME
!     OFILE  - OUTPUT FILE NAME
!     ICODE  - LOCATE CODES FOR INPUT CARD TYPES
!     OCODE  - LOCATE CODES FOR OUTPUT CARD TYPES
!     CTYPES - BCD NAMES OF SUBSTRUCTURING CARDS
!     CTYPEO - BCD NAMES OF CORRESPONDING NASTRAN CARDS
 
 
 INTEGER, INTENT(IN)                      :: ntype
 INTEGER, INTENT(IN)                      :: ifile
 INTEGER, INTENT(IN)                      :: sfile
 INTEGER, INTENT(IN OUT)                  :: ofile
 INTEGER, INTENT(IN)                      :: icode(4,1)
 INTEGER, INTENT(IN)                      :: ocode(4,1)
 INTEGER, INTENT(IN OUT)                  :: ctypes(2,8)
 INTEGER, INTENT(IN OUT)                  :: ctypeo(2,8)
 INTEGER :: buf1,buf2,buf3,z,sysbuf,outt,card(3), nlimit(3),subnam(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / idry,NAME(2)
 COMMON /sgencm/ nono,nss,iptr,buf1,buf2,buf3,nz
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,outt
 DATA    nlimit/ 3*2147483647 /
 DATA    subnam/ 4HSGEN,4HM    /
 
!     OPEN FILES
 
 CALL gopen (ifile,z(buf1),0)
 CALL gopen (sfile,z(buf2),0)
 CALL gopen (ofile,z(buf3),1)
 
!     READ HEADER FROM IFILE - DETERMINE IF SUBSTRUCTURING OR NASTRAN
!     CARD
 
 FILE = ifile
 10 CALL READ (*1002,*1003,ifile,card,3,0,idx)
 IF (card(1) == nlimit(1)) GO TO 70
 DO  i = 1,ntype
   IF (icode(1,i) /= card(1) .OR. icode(2,i) /= card(2)) CYCLE
   
!     SKIP RECORD IF SUBSTRUCTURING CARD
   
   CALL fwdrec (*70,ifile)
   GO TO 10
 END DO
 DO  i = 1,ntype
   IF (ocode(1,i) /= card(1) .OR. ocode(2,i) /= card(2)) CYCLE
   
!     FATAL ERROR IF BOTH SUBSTRUCTURING AND NASTRAN CARDS
   
   IF (icode(4,i) == 0) GO TO 40
   nono = 1
   j = ocode(4,i)
   WRITE (outt,6330) ufm,NAME,(ctypes(k,j),k=1,2),(ctypeo(k,j),k=1,2)
   CALL fwdrec (*70,ifile)
   GO TO 10
 END DO
 
!     COPY RECORD FROM IFILE TO OUTPUT
 
 40 CALL WRITE (ofile,card,3,0)
 50 CALL READ  (*1002,*60,ifile,z,nz,0,nwds)
 CALL WRITE (ofile,z,nz,0)
 GO TO 50
 60 CALL WRITE (ofile,z,nwds,1)
 GO TO 10
 
!     COPY RECORD FROM SFILE TO OUTPUT
 
 70 i1 = 1
 80 DO  i = i1,ntype
   IF (icode(4,i) == 1) GO TO 100
 END DO
 GO TO 150
 100 CALL fread (sfile,card,3,0)
 CALL WRITE (ofile,ocode(1,i),3,0)
 FILE = sfile
 110 CALL READ  (*1002,*120,sfile,z,nz,0,nwds)
 CALL WRITE (ofile,z,nz,0)
 GO TO 110
 120 CALL WRITE (ofile,z,nwds,1)
 i1 = i + 1
 GO TO 80
 
!     CLOSE FILES
 
 150 CALL WRITE (ofile,nlimit,3,1)
 CALL CLOSE (ifile,1)
 CALL CLOSE (sfile,1)
 CALL CLOSE (ofile,1)
 RETURN
 
!     ERRORS
 
 1002 m = -2
 GO TO 2000
 1003 m = -3
 2000 CALL mesage (m,FILE,subnam)
 RETURN
 
 6330 FORMAT (a23,' 6330, SOLUTION SUBSTRUCTURE ',2A4,3H - ,2A4,' AND ',  &
     2A4,' CARDS CANNOT BE USED TOGETHER.', /30X, 'USE EITHER ONE, BUT NOT BOTH.')
END SUBROUTINE sgenm
