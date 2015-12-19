SUBROUTINE bdat02
     
!     THIS SUBROUTINE PROCESSES CONCT BULK DATA AND WRITES CONNECTION
!     ENTRIES IN TERMS OF CODED GRID POINT ID NUMBERS ON SCR1
 
 EXTERNAL        rshift,andf
 LOGICAL :: tdat,PRINT
 INTEGER :: scr1,outt,buf1,buf2,conct(2),flag,geom4,id(2),  &
     comp,nams(4),io(9),aaa(2),conset,andf,rshift,combo
 DIMENSION       ibits(32),jbits(32),NAME(14),ihd(16)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc, geom4,casecc
 COMMON /zzzzzz/ z(1)
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
     mcon,restct(7,7),isort,origin(7,3),iprint
 COMMON /cmb004/ tdat(6)
 COMMON /output/ ititl(96),ihead(96)
 COMMON /cmbfnd/ inam(2),ierr
 COMMON /BLANK / step,idry
 DATA    aaa   / 4HBDAT,4H02       / , conct  / 210,2  /
 DATA    ihd   / 4H  su , 4HMMAR , 4HY of , 4H con , 4HNECT ,  &
     4HION  , 4HENTR , 4HIES  , 4HSPEC , 4HIFIE ,  &
     4HD by , 4H con , 4HCT   , 4HBULK , 4H dat , 4HA    /
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
 IF (andf(rshift(iprint,3),1) == 1) PRINT = .true.
 np2 = 2*npsub
 DO  i = 1,np2,2
   j = i/2 + 1
   NAME(i  ) = combo(j,1)
   NAME(i+1) = combo(j,2)
 END DO
 ifile = scr1
 CALL OPEN (*220,scr1,z(buf2),3)
 CALL locate (*180,z(buf1),conct,flag)
 ifile = geom4
 40 CALL READ (*200,*170,geom4,id,1,0,n)
 IF (id(1) == conset) GO TO 60
 CALL READ (*200,*210,geom4,id,-5,0,n)
 50 CALL READ (*200,*210,geom4,id,2,0,n)
 IF (id(1)+id(2) /= -2) GO TO 50
 GO TO 40
 60 CALL READ (*200,*210,geom4,comp,1,0,n)
 IF (.NOT.PRINT) GO TO 80
 CALL page
 CALL page2 (6)
 WRITE  (outt,70) (NAME(kdh),kdh=1,np2)
 70 FORMAT (/24X,74HNOTE  grid point id numbers have been coded TO the &
     &component substructure ,/30X,75HWITHIN a given pseudostructure by &
     &- 1000000*component no. + actual grid id., //15X,22HCONNECTED   c&
     &onnection,23X,33HGRID point id for pseudostructure/18X,3HDOF,9X, &
     4HCODE,3X,7(3X,2A4)/)
 80 CONTINUE
 tdat(2) = .true.
 CALL encode (comp)
 CALL READ (*200,*210,geom4,nams,4,0,n)
 CALL finder (nams(1),is1,ic1)
 IF (ierr /= 1) GO TO 90
 WRITE (outt,100) ufm,nams(1),nams(2)
 idry = -2
 90 CONTINUE
 CALL finder (nams(3),is2,ic2)
 IF (ierr /= 1) GO TO 110
 WRITE (outt,100) ufm,nams(3),nams(4)
 idry = -2
 100 FORMAT (a23,' 6523, THE BASIC SUBSTRUCTURE ',2A4, /30X,  &
     'REFERED TO BY A CONCT  BULK DATA CARD CAN NOT BE FOUND ',  &
     'IN THE PROBLEM TABLE OF CONTENTS.')
 110 CALL READ (*200,*210,geom4,id,2,0,n)
 
 IF (id(1)+id(2) == -2) GO TO 40
 IF (is1 /= is2) GO TO 130
 kk = 2*is1 - 1
 WRITE  (outt,120) ufm,id(1),id(2),NAME(kk),NAME(kk+1)
 120 FORMAT (a23,' 6536, MANUAL CONNECTION DATA IS ATTEMPTING TO ',  &
     'CONNECT', /31X,'GRID POINTS',i9,5X,4HAND ,i8, /31X,  &
     'WHICH ARE BOTH CONTAINED IN PSEUDOSTRUCTURE ',2A4)
 idry = -2
 130 CONTINUE
 DO  i = 1,9
   io(i) = 0
 END DO
 io(1) = comp
 io(2) = 2**(is1-1) + 2**(is2-1)
 io(2+is1) = ic1*1000000 + id(1)
 io(2+is2) = ic2*1000000 + id(2)
 nwd = 2 + npsub
 CALL WRITE (scr1,io,nwd,1)
 IF (.NOT.PRINT .OR. idry == -2) GO TO 160
 CALL bitpat (io(1),ibits)
 CALL bitpat (io(2),jbits)
 CALL page2 (1)
 WRITE (outt,150) (ibits(kdh),kdh=1,2),(jbits(kdh),kdh=1,2),  &
     (io(kdh+2),kdh=1,npsub)
 150 FORMAT (16X,a4,a2,6X,a4,a3,2X,7(3X,i8))
 160 CONTINUE
 GO TO 110
 170 CONTINUE
 180 CALL CLOSE (scr1,1)
 RETURN
 
 200 imsg = -2
 GO TO 230
 210 imsg = -3
 GO TO 230
 220 imsg = -1
 230 CALL mesage (imsg,ifile,aaa)
 RETURN
END SUBROUTINE bdat02
