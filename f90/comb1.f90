SUBROUTINE comb1
     
!     THIS IS THE MODULE FOR THE COMBINATION OF SUBSTRUCTURES.
 
!     IT IS PRIMARILY AN INITIALIZER AND DRIVER CALLING THE ROUTINES
!     NECESSARY TO PROCESS THE COMBINE. THE SUBROUTINES ARE
 
!          CMCASE  -  READS THE CASECC DATA BLOCK AND INITIALIZES
!                     PARAMETERS FOR THE COMBINE OPERATION.
!          CMTOC   -  GENERATES THE TABLE OF CONTENTS OF PSEUDO-
!                     STRUCTURES BEING COMBINED AND THEIR COMPONENT
!                     BASIC SUBSTRUCTURES.
!          BDAT01  -  PROCESSES THE CONCT1 BULK DATA.
!          BDAT02  -  PROCESSES THE CONCT  BULK DATA.
!          BDAT03  -  PROCESSES THE TRANS  BULK DATA.
!          BDAT04  -  PROCESSES THE RELES  BULK DATA.
!          BDAT05  -  PROCESSES THE GNEW   BULK DATA.
!          BDAT06  -  PROCESSES THE GTRAN  BULK DATA.
!          CMSFIL  -  GENERATES SUBFIL - THE BASIC FILE USED TO STORE
!                     THE DATA NECESSARY TO AFFECT THE COMBINATION.
!          CMCONT  -  GENERATES THE CONNECTION ENTRIES TO BE USED.
!          CMCKCD  -  CHECKS VALIDITY OF MANUALLY-SPECIFIED CONNECTIONS
!          CMAUTO  -  PROCESSES USERS REQUEST FOR AUTOMATIC
!                     COMBINATION OF SUBSTRUCTURES.
!          CMRELS  -  APPLIES ANY MANUAL RELEASE DATA TO THE SYSTEM.
!          CMCOMB  -  PROCESSES MULTIPLY CONNECTED POINTS.
!          CMDISC  -  PROCESSES GRID POINTS NOT TO BE CONNECTED.
!          CMSOFO  -  GENERATES NEW SOF ITEMS FOR THE RESULTANT
!                     COMBINED STRUCTURE.
!          CMHGEN  -  GENERATES THE DOF TRANSFORMATION MATRIX FOR
!                     EACH COMPONENT TO THE COMBINATION
 
 LOGICAL :: tdat,conect,iauto,tran,mcon,tocopn,lonly
 INTEGER :: scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc,geom4,  &
     casecc,buf1,buf2,buf3,buf4,buf5,score,dry,step,  &
     sys,outt,aaa(2),restct,sccstm,scr3
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc,  &
     geom4,casecc,sccstm,scr3
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,intp,outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
     mcon,restct(7,7),isort,origin(7,3),iprint,tocopn
 COMMON /cmb004/ tdat(6),nipnew,cnam(2),lonly
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sys(69)
 COMMON /BLANK / step,dry
 DATA    aaa   / 4HCOMB,4H1    /
 
 IF (dry == 0 .OR. dry == -2) GO TO 210
 scr1   = 301
 scr2   = 302
 scbdat = 303
 scsfil = 304
 scconn = 305
 scmcon = 306
 sctoc  = 307
 sccstm = 308
 scr3   = 309
 geom4  = 102
 casecc = 101
 DO  i = 1,6
   tdat(i) = .false.
 END DO
 DO  i = 1,7
   DO  j = 1,3
     origin(i,j) = 0.0
   END DO
 END DO
 lonly = .false.
 intp  = sys(4)
 outt  = sys(2)
 ibuf  = sys(1)
 
 nz    = korsz(z(1))
 buf1  = nz - ibuf - 2
 buf2  = buf1 - ibuf
 buf3  = buf2 - ibuf
 buf4  = buf3 - ibuf
 buf5  = buf4 - ibuf
 ib1   = buf5 - ibuf
 ib2   = ib1  - ibuf
 ib3   = ib2  - ibuf
 score = 1
 lcore = ib3 - 1
 IF (lcore > 0) GO TO 30
 CALL mesage (8,0,aaa)
 dry   = -2
 GO TO 210
 
 30 CALL OPEN (*120,scconn,z(buf2),1)
 CALL CLOSE (scconn,2)
 CALL sofopn (z(ib1),z(ib2),z(ib3))
 
 CALL cmcase
 IF (dry == -2) GO TO 130
 CALL cmtoc
 IF (.NOT.lonly) GO TO 40
 CALL cmsofo
 GO TO 70
 
 40 ifile = geom4
 CALL preloc (*90,z(buf1),geom4)
 IF (.NOT.conect) GO TO 50
 CALL bdat01
 CALL bdat02
 50 ifile = scbdat
 CALL OPEN (*120,scbdat,z(buf2),1)
 CALL bdat05
 CALL bdat06
 CALL bdat03
 CALL CLOSE (geom4,1)
 
 CALL cmsfil
 CALL preloc (*100,z(buf1),geom4)
 CALL bdat04
 CALL CLOSE (geom4,1)
 IF (dry == -2) GO TO 150
 60 IF (tdat(1) .OR. tdat(2)) CALL cmcont
 IF (dry == -2) GO TO 170
 CALL cmauto
 IF (tdat(1) .OR. tdat(2)) CALL cmckcd
 IF (dry == -2) GO TO 170
 IF (tdat(4)) CALL cmrels
 CALL cmmcon (nce)
 nps  = npsub + 1
 ndof = 6
 IF (mcon) CALL cmcomb (nps,nce,ndof,z)
 IF (dry == -2) GO TO 170
 
 CALL cmckdf
 IF (dry == -2) GO TO 170
 CALL cmdisc
 CALL cmsofo
 CALL cmhgen
 
 70 CALL sofcls
 IF (tocopn) CALL CLOSE (sctoc,1)
 WRITE  (outt,80) uim
 80 FORMAT (a29,' 6521, MODULE COMB1 SUCCESSFULLY COMPLETED.')
 GO TO 210
 
 90 IF (conect .OR. tran) GO TO 100
 ifile = scbdat
 CALL OPEN (*120,scbdat,z(buf2),1)
 CALL eof (scbdat)
 CALL CLOSE (scbdat,1)
 CALL cmsfil
 IF (.NOT.conect) GO TO 60
 
!     ERRORS
 
 100 WRITE  (outt,110) ufm
 110 FORMAT (a23,' 6510, THE REQUESTED COMBINE OPERATION REQUIRES ',  &
     'SUBSTRUCTURE BULK DATA WHICH HAS NOT BEEN GIVEN.')
 GO TO  190
 120 CALL mesage (1,scbdat,aaa)
 GO TO  170
 130 WRITE  (outt,140) ufm
 140 FORMAT (a23,' 6535, MODULE COMB1 TERMINATING DUE TO ABOVE ',  &
     'SUBSTRUCTURE CONTROL ERRORS.')
 GO TO  200
 150 WRITE  (outt,160) ufm
 160 FORMAT (a23,' 6536, MODULE COMB1 TERMINATING DUE TO ABOVE ERRORS',  &
     ' IN BULK DATA.')
 GO TO  190
 170 WRITE  (outt,180) ufm
 180 FORMAT (a23,' 6537, MODULE COMB1 TERMINATING DUE TO ABOVE ERRORS')
 190 IF (tocopn) CALL CLOSE (sctoc,1)
 200 dry = -2
 CALL sofcls
 210 RETURN
END SUBROUTINE comb1
