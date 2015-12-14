SUBROUTINE bdat01
     
!     THIS SUBROUTINE PROCESSES CONCT1 BULK DATA GENERATING
!     CONNECTION ENTRIES IN TERMS OF GRID POINT ID NUMBERS
!     CODED TO THE PSEUDO-STRUCTURE ID NUMBER.
!     THESE ARE THEN WRITTEN ON SCR1.
 
 EXTERNAL        rshift,andf
 LOGICAL :: tdat,PRINT
 INTEGER :: io(9),id(14),is(7),ic(7),scr1,conset,geom4,aaa(2),  &
     flag,buf1,conct1(2),buf2,outt,andf,rshift,combo
 DIMENSION       ibits(32),jbits(32),NAME(14),ihd(16)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc, geom4,casecc
 COMMON /zzzzzz/ z(1)
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
     mcon,restct(7,7),isort,origin(7,3),iprint
 COMMON /cmb004/ tdat(6)
 COMMON /cmbfnd/ inam(2),ierr
 COMMON /output/ ititl(96),ihead(96)
 COMMON /BLANK / step,idry
 DATA    aaa   / 4HBDAT,4H01     / , conct1 / 110,41 /
 DATA    ihd   / 4H  su , 4HMMAR , 4HY of , 4H con , 4HNECT ,  &
     4HION  , 4HENTR , 4HIES  , 4HSPEC , 4HIFIE ,  &
     4HD by , 4H con , 4HCT1  , 4HBULK , 4H dat , 4HA    /
 DATA    iblnk / 4H     /
 
 DO  i = 1,96
   ihead(i) = iblnk
 END DO
 j = 1
 DO  i = 73,88
   ihead(i) = ihd(j)
   j = j + 1
 END DO
 PRINT = .false.
 IF (andf(rshift(iprint,2),1) == 1) PRINT = .true.
 np2 = 2*npsub
 DO  i = 1,np2,2
   j = i/2 + 1
   NAME(i  ) = combo(j,1)
   NAME(i+1) = combo(j,2)
 END DO
 ifile = scr1
 CALL OPEN (*320,scr1,z(buf2),1)
 CALL locate (*400,z(buf1),conct1,flag)
 ifile = geom4
 30 CALL READ (*300,*210,geom4,id,2,0,n)
 nss   = id(1)
 nssp1 = nss + 1
 IF (id(2) == conset) GO TO 50
 40 CALL READ (*300,*310,geom4,id,1,0,nnn)
 IF (id(1) /= -1) GO TO 40
 GO TO 30
 50 nwd = 2*nss
 IF (.NOT.PRINT) GO TO 70
 CALL page
 CALL page2 (6)
 WRITE  (outt,60) (NAME(kdh),kdh=1,np2)
 60 FORMAT (/24X,74HNOTE  grid point id numbers have been coded TO the  &
     component substructure, /30X,75HWITHIN a given pseudostructure by  &
     - 1000000*component no. + actual grid id.,//15X,22HCONNECTED   co  &
     nnection,23X,33HGRID point id for pseudostructure, /18X,3HDOF,9X,  &
     4HCODE,3X,7(3X,2A4)/)
 70 CONTINUE
 
!     MAKING IT TO 50 IMPLIES THAT CONCT1 DATA EXISTS
 
 tdat(1) = .true.
 CALL READ (*300,*310,geom4,id,nwd,0,nnn)
 DO  i = 1,nss
   j = 2*(i-1)
   CALL finder (id(1+j),is(i),ic(i))
   IF (ierr /= 1) CYCLE
   WRITE  (outt,80) ufm,id(1+j),id(2+j)
   80 FORMAT (a23,' 6522, THE BASIC SUBSTRUCTURE ',2A4, /30X,  &
       'REFERED TO BY A CONCT1 BULK DATA CARD CAN NOT BE FOUND ',  &
       'IN THE PROBLEM TABLE OF CONTENTS.')
   idry = -2
 END DO
 100 DO  i = 1,9
   io(i) = 0
 END DO
 DO  i = 1,nssp1
   CALL READ (*300,*310,geom4,id(i),1,0,nnn)
   IF (id(i) == -1) GO TO 30
 END DO
 DO  i = 1,nss
   DO  j = 1,nss
     IF (i == j) CYCLE
     IF (is(i) == is(j) .AND. id(i+1) /= 0 .AND. id(j+1) /= 0) GO TO 150
   END DO
 END DO
 GO TO 170
 150 kk = 2*is(i) - 1
 WRITE  (outt,160) ufm,id(i+1),id(j+1),NAME(kk),NAME(kk+1)
 160 FORMAT (a23,' 6536, MANUAL CONNECTION DATA IS ATTEMPTING TO ',  &
     'CONNECT', /31X,'GRID POINTS',i9,5X,4HAND ,i8, /31X,  &
     'WHICH ARE BOTH CONTAINED IN PSEUDOSTRUCTURE ',2A4)
 idry = -2
 170 CALL encode (id(1))
 io(1) = id(1)
 isum  = 0
 DO  i = 1,nss
   IF (id(i+1) == 0) CYCLE
   IF (id(i+1) /= 0) isum = isum + 2**(is(i)-1)
   m = 2 + is(i)
   io(m) = ic(i)*1000000 + id(i+1)
 END DO
 io(2) = -1*isum
 nwd   = 2 + npsub
 CALL WRITE (scr1,io,nwd,1)
 IF (.NOT.PRINT .OR. idry == -2) GO TO 200
 CALL bitpat (io(1),ibits)
 CALL bitpat (IABS(io(2)),jbits)
 CALL page2 (1)
 WRITE (outt,190) (ibits(kdh),kdh=1,2),(jbits(kdh),kdh=1,2),  &
     (io(kdh+2),kdh=1,npsub)
 190 FORMAT (16X,a4,a2,6X,a4,a3,2X,7(3X,i8))
 200 CONTINUE
 GO TO 100
 210 CONTINUE
 GO TO 400
 
 300 imsg = -2
 GO TO 330
 310 imsg = -3
 GO TO 330
 320 imsg = -1
 330 CALL mesage (imsg,ifile,aaa)
 400 CALL CLOSE  (scr1,2)
 RETURN
END SUBROUTINE bdat01
