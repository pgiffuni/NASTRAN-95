SUBROUTINE bdat04
     
!     THIS SUBROUTINE PROCESSES THE RELES BULK DATA.
 
 EXTERNAL        rshift,andf

 LOGICAL :: name,tdat,print,pager

 INTEGER :: scr2,buf2,buf1,reles(2),flag,geom4,id(2),aaa(2),  &
            conset,ip(6),icc(6),andf,scbdat,rshift,ihd(96),   &
            outt,ibas(2)

 DIMENSION  ibits(32),jbits(32),kbits(32)

 CHARACTER (LEN=23) :: ufm

 COMMON /xmssg / ufm
 COMMON /zzzzzz/ z(1)
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon, sctoc,geom4,casecc
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
     mcon,restct(7,7),isort,origin(7,3),iprint
 COMMON /cmb004/ tdat(6)
 COMMON /cmbfnd/ inam(2),ierr
 COMMON /output/ ititl(96),ihead(96)
 COMMON /system/ xxx,iot,junk(6),ipage,line,itline,maxlin,idat(3)
 COMMON /blank / step,idry

 DATA    ihd   / 11*4H    ,4H  SU,4HMMAR,4HY OF,4H PRO,4HCESS,	   &
                 4HED R,4HELES,4H BUL,4HK DA,4HTA  ,18*4H    ,	   &
                 4H   B,4HASIC,2*4H    ,4H GRI,4HD   ,4H     ,	   &
                 4HREQU,4HESTE,4HD   ,4H  IN,4HTERN,  4HAL   ,	   &
                 4H   C,4HURRE,4HNT  ,4H  DO,4HF TO,  4H BE  ,	   &
                 13*4H    ,4HSUBS,4HTRUC,4HTURE,4H   P,4HOINT,	   &
                 4H ID ,4H    ,4H REL,4HEASE,4H    ,4H  PO,4HINT , &
                 4HNO. ,4H    ,4H DOF,4H    ,4H   R,4HELEA,4HSED , &
                 6*4H       /
 DATA    reles / 410,4      / , aaa / 4HBDAT,4H04    /
 
 DO  i = 1,96
   ihead(i) = ihd(i)
 END DO
 pager = .true.
 PRINT = .false.
 IF (andf(rshift(iprint,8),1) == 1) PRINT = .true.
 ifile = scbdat
 CALL OPEN (*200,scbdat,z(buf2),0)
 CALL skpfil (scbdat,3)
 CALL CLOSE  (scbdat,2)
 CALL OPEN (*200,scbdat,z(buf2),3)
 ifile = scr2
 CALL locate (*170,z(buf1),reles,flag)
 ifile = geom4
 20 CALL READ (*210,*160,geom4,id,1,0,n)
 IF (id(1) == conset) GO TO 40
 30 CALL READ (*210,*220,geom4,id,2,0,n)
 IF (id(1)+id(2) /= -2) GO TO 30
 GO TO 20
 40 NAME = .true.
 IF (pager .AND. PRINT) CALL page
 pager = .false.
 tdat(4) = .true.
 50 CALL READ (*210,*220,geom4,id,2,0,n)
 IF (id(1)+id(2) /= -2) GO TO 60
 CALL WRITE (scbdat,id,0,1)
 GO TO 20
 60 IF (.NOT.NAME) GO TO 100
 CALL finder (id,is,ic)
 ibas(1) = id(1)
 ibas(2) = id(2)
 IF (ierr /= 1) GO TO 90
 WRITE  (outt,70) ufm,(id(k),k=1,2)
 70 FORMAT (a23,' 6517, THE BASIC SUBSTRUCTURE  ',2A4, /30X,  &
     'REFERED TO BY A RELES  BULK DATA CARD CAN NOT BE FOUND ',  &
     'IN THE PROBLEM TABLE OF CONTENTS.')
 idry = -2
 80 CALL READ (*210,*220,geom4,id,2,0,n)
 IF (id(1)+id(2) /= -2) GO TO 80
 GO TO 20
 90 CONTINUE
 CALL WRITE (scbdat,is,1,0)
 NAME = .NOT.NAMe
 GO TO 50
 100 CALL fndgrd (is,ic,id(1),ip,icc,n)
 IF (ierr /= 1) GO TO 120
 WRITE  (outt,110) ufm,id(1),inam
 110 FORMAT (a23,' 6515, GRID POINT',i10,' BASIC SUBSTRUCTURE ',2A4,  &
     ' DOES NOT EXIST.')
 idry = -2
 GO TO 50
 120 CALL encode (id(2))
 CALL bitpat (id(2),ibits)
 DO  i = 1,n
   iccc = andf(id(2),icc(i))
   CALL bitpat (iccc,jbits)
   icc(i) = andf(icc(i),63)
   CALL bitpat (icc(i),kbits)
   IF (iccc == 0) CYCLE
   IF (.NOT.PRINT ) GO TO 140
   WRITE (outt,130) ibas,id(1),ibits(1),ibits(2),ip(i),kbits(1),  &
       kbits(2),jbits(1),jbits(2)
   130 FORMAT (35X,2A4,5X,i8,7X,a4,a2,6X,i8,6X,a4,a2,6X,a4,a2)
   140 CONTINUE
   CALL WRITE (scbdat,ip(i),1,0)
   CALL WRITE (scbdat,iccc, 1,0)
 END DO
 GO TO 50
 160 CONTINUE
 170 CALL CLOSE (scbdat,1)
 RETURN
 
 200 imsg = -1
 GO TO 230
 210 imsg = -2
 GO TO 230
 220 imsg = -3
 230 CALL mesage (imsg,ifile,aaa)
 
 RETURN
END SUBROUTINE bdat04
