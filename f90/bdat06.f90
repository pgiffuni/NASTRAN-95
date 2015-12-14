SUBROUTINE bdat06
     
!     THIS SUBROUTINE PROCESSES THE GTRAN BULK DATA
 
 EXTERNAL        rshift,andf
 LOGICAL :: PRINT,tdat
!WKBI 8/94 ALPHA-VMS
 INTEGER :: geom4, scr1
 INTEGER :: scr2,buf3,scbdat,buf2,buf1,gtran(2),flag,id(5),  &
     combo,score,z,aaa(2),outt,buf4,andf,rshift,ihd(96)
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
 COMMON /BLANK / step,idry
 DATA   gtran  / 1510,15 / , aaa/ 4HBDAT,4H06   /
 DATA   ihd  / 11*4H    ,4H  su,4HMMAR,4HY of,4H pro,4HCESS,4HED g,  &
     4HTRAN,4H bul,4HK da,4HTA  ,19*4H    ,4H pse,4HUDO-,4H    ,  &
     4H    ,4H com,4HPONE,4HNT  ,4H    ,4H   t,4HRANS,2*4H     ,  &
     4HGRID,4H    ,4HREFE,4HRENC,4HE   ,14*4H    ,4H  st,4HRUCT,  &
     4HURE ,4HNO. ,4H   s,4HTRUC,4HTURE,4H no.,4H    ,4H  se   ,  &
     4HT id,2*4H    ,4H id ,4H    ,4HTRAN,4HS  i,4HD   ,7*2H   /
 
 ifile = scr1
 kk    = 0
 PRINT = .false.
 IF (andf(rshift(iprint,5),1) == 1) PRINT = .true.
 DO  i = 1,96
   ihead(i) = ihd(i)
 END DO
 CALL OPEN (*100,scr2,z(buf3),1)
 ifile = scbdat
 CALL locate (*80,z(buf1),gtran,flag)
 IF (PRINT) CALL page
 ifile = geom4
 20 CALL READ (*110,*60,geom4,id,1,0,n)
 DO  i = 1,npsub
   IF (id(1) == combo(i,3)) GO TO 40
 END DO
 CALL READ (*110,*120,geom4,id,-4,0,n)
 GO TO 20
 40 tdat(6) = .true.
 kk = kk + 1
 CALL READ (*110,*120,geom4,id(2),4,0,n)
 CALL finder (id(2),is,ic)
 IF (ierr /= 1) GO TO 50
 WRITE (outt,210) ufm,id(2),id(3)
 idry = -2
 50 CONTINUE
 IF (PRINT) CALL page2 (1)
 IF (PRINT) WRITE (outt,200) is,ic,id(1),id(4),id(5)
 id(3) = id(1)
 id(1) = is
 id(2) = ic
 id(4) = ic*1000000 + id(4)
 z(buf4+kk) = id(5)
 CALL WRITE (scr2,id,5,0)
 GO TO 20
 60 CALL WRITE (scr2,id,0,1)
 CALL CLOSE (scr2,1)
 IF (.NOT.tdat(6)) GO TO 80
 ifile = scr2
 CALL OPEN (*100,scr2,z(buf3),2)
 CALL READ (*110,*70,scr2,z(score),lcore,0,nn)
 GO TO 130
 70 CALL sort (0,0,5,1,z(score),nn)
 CALL WRITE (scbdat,z(score),nn,1)
 80 CALL eof (scbdat)
 z(buf4) = kk
 CALL CLOSE (scr2,1)
 IF (PRINT) CALL page2 (3)
 IF (PRINT) WRITE (outt,220)
 RETURN
 
 100 imsg = -1
 GO TO 140
 110 imsg = -2
 GO TO 140
 120 imsg = -3
 GO TO 140
 130 imsg = -8
 140 CALL mesage (imsg,ifile,aaa)
 RETURN
 
 200 FORMAT (36X,i1,14X,i5,8X,i8,4X,i8,4X,i8)
 210 FORMAT (a23,' 6530, THE BASIC SUBSTRUCTURE ',2A4, /30X,  &
     'REFERED TO BY A GTRAN BULK DATA CARD WHICH CANNOT BE ',  &
     'FOUNDD IN THE PROBLEM TABLE OF CONTENTS.')
 220 FORMAT (/5X,'NOTE - THE PSEUDOSTRUCTURE AND COMPONENT NUMBERS RE',  &
     'FER TO THEIR POSITIONS IN THE PROBLEM TABLE OF CONTENTS.')
END SUBROUTINE bdat06
