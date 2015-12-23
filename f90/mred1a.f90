SUBROUTINE mred1a (mode)
     
!     THIS SUBROUTINE PROCESSES THE BDYC DATA FOR THE FIXED
!     IDENTIFICATION SET (FIXSET) AND THE BOUNDARY IDENTIFICATION SET
!     (BNDSET) FOR THE MRED1 MODULE.
 
!     INPUT DATA
!     GINO - GEOM4    - BDYC DATA
!     MODE            - PROCESSING OPERATION FLAG
!                     = 1, PROCESS FIXED ID SET
!                     = 2, PROCESS BOUNDARY ID SET
 
!     OUTPUT DATA
!     GINO - USETX    - S,R,B DEGREES OF FREEDOM
 
!     PARAMETERS
!     INPUT  - GBUF1  - GINO BUFFER
!              KORLEN - CORE LENGTH
!              BNDSET - BOUNDARY SET IDENTIFICATION NUMBER
!              FIXSET - FIXED SET IDENTIFICATION NUMBER
!              IO     - OUTPUT OPTION FLAG
!              KORUST - STARTING ADDRESS OF USET ARRAY
!              NCSUBS - NUMBER OF CONTRIBUTING SUBSTRUCTURES
!              NAMEBS - BEGINNING ADDRESS OF BASIC SUBSTRUCTURES NAMES
!              KBDYC  - BEGINNING ADDRESS OF BDYC DATA
!              NBDYCC - NUMBER OF BDYC CARDS
!              USETL  - NUMBER OF WORDS IN USET ARRAY
!     OUTPUT - NOUS   - FIXED POINTS FLAG
!                       .GE.  0, FIXED POINTS DEFINED
!                       .EQ. -1, NO FIXED POINTS DEFINED
!              DRY    - MODULE OPERATION FLAG
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(IN OUT)                  :: mode
 EXTERNAL        rshift,andf
 LOGICAL :: ponly
 DIMENSION       array(3),bdyc(2),modnam(2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / oldnam(2),dry,idum1,nous,skipm,idum2(3),gbuf1,  &
     idum3(4),korlen,idum4(2),bndset,fixset,idum5,io,  &
     idum6(6),ncsubs,namebs,idum7(3),kbdyc,nbdycc, usetl,korust,idum14(5),ponly
 COMMON /zzzzzz/ z(1)
 COMMON /system/ idum8,iprntr,idum9(6),nlpp,idum10(2),line
 COMMON /two   / itwo(32)
 COMMON /bitpos/ idum11(5),ul,ua,uf,idum12,un,idum13(11),ui
 DATA    geom4 , bdyc  /102,910,9/
 DATA    modnam/ 4HMRED,4H1A     /
 
!     TEST PROCESSING MODE FLAG
 
 IF (mode == 2) GO TO 10
 
!     TEST FIXED SET ID FLAG AND SET FIXED INDEX
 
 IF (fixset == 0 .OR. skipm == -1) GO TO 260
 setid  = fixset
 ishift = 10
 GO TO 20
 
!     SET BOUNDARY INDEX
 
 10 IF (bndset == 0) GO TO 240
 setid  = bndset
 ishift = 1
 
!     ALLOCATE USET ARRAY AND TEST OPEN CORE LENGTH
 
 IF (nous == 1) GO TO 40
 20 kbdyc = korust + usetl
 IF (kbdyc >= korlen) GO TO 200
 
!     TURN UL, UA, UF, UN, AND UI BITS ON IN USET ARRAY
 
 ibits = itwo(ul) + itwo(ua) + itwo(uf) + itwo(un) + itwo(ui)
 DO  i = 1,usetl
   z(korust+i-1) = ibits
 END DO
 
!     READ BOUNDARY SET (BDYC) BULK DATA FOR REQUESTED FIXED SET
!     ID (FIXSET) OR BOUNDARY SET ID (BNDSET)
 
 40 ifile = geom4
 CALL preloc (*170,z(gbuf1),geom4)
 CALL locate (*230,z(gbuf1),bdyc,iflag)
 50 CALL READ (*180,*230,geom4,array,1,0,iflag)
 IF (array(1) == setid) GO TO 70
 60 CALL READ (*180,*190,geom4,array,3,0,iflag)
 IF (array(3) == -1) GO TO 50
 GO TO 60
 
!     SET ID FOUND, STORE AT Z(KBDYC+NWDS)
 
 70 nwds = 0
 80 CALL READ (*180,*190,geom4,z(kbdyc+nwds),3,0,iflag)
 IF (z(kbdyc+nwds+2) == -1) GO TO 110
 
!     CHECK THAT SUBSTRUCTURE IS A COMPONENT OF STRUCTURE BEING
!     REDUCED
 
 DO  i = 1,ncsubs
   j = 2*(i-1)
   IF ((z(namebs+j) == z(kbdyc+nwds)) .AND. (z(namebs+j+1) ==  &
       z(kbdyc+nwds+1))) GO TO 100
 END DO
 
!     SUBSTRUCTURE IS NOT A COMPONENT
 
 IF (mode == 1) WRITE (iprntr,900) uwm,z(kbdyc+nwds), z(kbdyc+nwds+1)
 IF (mode == 2) WRITE (iprntr,901) uwm,z(kbdyc+nwds), z(kbdyc+nwds+1)
 dry = -2
 GO TO 80
 
!     SAVE BASIC SUBSTRUCTURE INDEX
 
 100 z(kbdyc+nwds+3) = i
 nwds = nwds + 4
 IF (kbdyc+nwds >= korlen) GO TO 200
 GO TO 80
 
!     CHECK FOR DUPLICATE BDYC SUBSTRUCTURE NAMES
 
 110 nwds = nwds/4
 IF (nwds <= 1) GO TO 125
 i = nwds - 1
 DO  j = 1,i
   k  = j + 1
   ii = 4*(j-1)
   DO  l = k,nwds
     ll = 4*(l-1)
     IF (z(kbdyc+ii  ) /= z(kbdyc+ll  )) CYCLE
     IF (z(kbdyc+ii+1) /= z(kbdyc+ll+1)) CYCLE
     WRITE (iprntr,902) ufm,oldnam,array(1)
     dry = -2
   END DO
 END DO
 
!     TEST OUTPUT OPTION
 
 125 CONTINUE
 IF (andf(rshift(io,ishift),1) == 0) GO TO 150
 IF (nwds == 0) GO TO 150
 line = nlpp + 1
 DO  i = 1,nwds
   IF (line <= nlpp) GO TO 130
   CALL page1
   IF (mode == 1) WRITE (iprntr,903) fixset
   IF (mode == 2) WRITE (iprntr,904) bndset
   line = line + 7
   130 j = 4*(i-1)
   IF (mode == 1) WRITE (iprntr,905) z(kbdyc+j),z(kbdyc+j+1), z(kbdyc+j+2)
   IF (mode == 2) WRITE (iprntr,906) z(kbdyc+j),z(kbdyc+j+1), z(kbdyc+j+2)
   line = line + 1
 END DO
 
!     SORT BDYC DATA ON SET ID
 
 150 nbdycc = nwds
 IF (nbdycc <= 1) GO TO 270
 nwds = 4*nbdycc
 CALL sort (0,0,4,3,z(kbdyc),nwds)
 GO TO 270
 
!     PROCESS SYSTEM FATAL ERRORS
 
 170 imsg = -1
 GO TO 220
 180 imsg = -2
 GO TO 210
 190 imsg = -3
 GO TO 210
 200 imsg = -8
 ifile = 0
 210 CALL CLOSE (geom4,1)
 220 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
 
!     PROCESS MODULE FATAL ERRORS
 
 230 IF (mode == 1) WRITE (iprntr,907) uwm,fixset
 IF (mode == 2) WRITE (iprntr,908) uwm,bndset
 dry = -1
 GO TO 250
 240 IF (ponly) GO TO 280
 WRITE (iprntr,909) ufm
 dry = -2
 250 CALL sofcls
 CALL CLOSE (geom4,1)
 RETURN
 
!     NO FIXED ID SET DATA
 
 260 nous = -1
 
!     END OF PROCESSING
 
 270 CALL CLOSE (geom4,1)
 280 CONTINUE
 
 900 FORMAT (a25,' 6622, A FIXED SET HAS BEEN SPECIFIED FOR ',2A4,  &
     ', BUT IT IS NOT A COMPONENT OF',/32X,'THE PSEUDOSTRUCTURE'  &
     ,      ' BEING PROCESSED.  THE FIXED SET WILL BE IGNORED.')
 901 FORMAT (a25,' 6604, A BOUNDARY SET HAS BEEN SPECIFIED FOR ',2A4,  &
     ', BUT IT IS NOT A COMPONENT OF',/32X,'THE PSEUDOSTRUCTURE'  &
     ,      ' BEING PROCESSED.  THE BOUNDARY SET WILL BE IGNORED.')
 902 FORMAT (a23,' 6623, SUBSTRUCTURE ',2A4,  &
     ' HAS DUPLICATE NAMES IN BDYC DATA SET ',i8,1H.)
 903 FORMAT (1H0,43X,'SUMMARY OF COMBINED FIXED SET NUMBER ',i8, //57X,  &
     'BASIC      FIXED', /54X,'SUBSTRUCTURE  SET ID', /58X, 'NAME      NUMBER',/)
 904 FORMAT (1H0,43X,'SUMMARY OF COMBINED BOUNDARY SET NUMBER ',i8,  &
     //57X,'BASIC      BOUNDARY', /54X,'SUBSTRUCTURE   SET ID',  &
     /58X,'NAME       NUMBER',/)
 905 FORMAT (56X,2A4,3X,i8)
 906 FORMAT (56X,2A4,4X,i8)
 907 FORMAT (a25,' 6621, FIXED SET',i9,' SPECIFIED IN CASE CONTROL ',  &
     'HAS NOT BEEN DEFINED BY BULK DATA.')
 908 FORMAT (a25,' 6606, BOUNDARY SET',i9,' SPECIFIED IN CASE CONTROL',  &
     ' HAS NOT BEEN DEFINED BY BULK DATA.')
 909 FORMAT (a23,' 6603, A BOUNDARY SET MUST BE SPECIFIED FOR A REDUCE'  &
     ,      ' OPERATION.')
 
 RETURN
END SUBROUTINE mred1a
