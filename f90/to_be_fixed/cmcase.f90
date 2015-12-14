SUBROUTINE cmcase
     
!     THIS SUBROUTINE PROCESSES THE CASE CONTROL DATA BLOCK
 
 EXTERNAL        orf
 LOGICAL :: iauto,tran,conect,lf(3),lonly,srch
 INTEGER :: casecc,buf2,step,z,cnam,restct,combo,outt,auto,  &
     orf,ncnam(2),ihd(96),ibits(32),conset,mnem(11),  &
     snam(7,2),idir(3),comp(7,2),symt(7),trans(7), isym(15,2),aaa(2),pora,papp
 DIMENSION       az(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /cmb001/ junk(8),casecc
 COMMON /cmb002/ buf1,buf2,junk1(6),outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
     mcon,restct(7,7),isort,origin(7,3),iprint
 COMMON /cmb004/ tdat(6),nipnew,cnam(2),lonly
 COMMON /zzzzzz/ z(1)
 COMMON /output/ ititl(96),ihead(96)
 COMMON /system/ xxx,iot,munk(6),nlpp,junk3(2),line,junk2(2), idat(3)
 COMMON /BLANK / step,idry,pora
 EQUIVALENCE     (z(1),az(1))
 DATA    nmnem / 11 /, idir/ 1HX, 1HY, 1HZ /,  auto/ 4HAUTO / ,  &
     aaa   / 4HCMCA, 4HSE   /
 DATA    mnem  / 4HOPTS, 4HSORT, 4HNAMC, 4HNAMS, 4HTOLE, 4HCONN,  &
     4HCOMP, 4HTRAN, 4HSYMT, 4HSEAR, 4HOUTP/
 DATA    isym  / 4,2,1,6,6,5,5,3,3,6*7,1HX,1HY,1HZ,2HXY,2HYX,2HXZ,  &
     2HZX,2HYZ,2HZY,3HXYZ,3HXZY,3HYXZ,3HYZX,3HZXY, 3HZYX /
 DATA    ihd   / 74*4H    , 4H sum , 4HMARY , 4H of  , 4HCASE ,  &
     4H con ,   4HTROL , 4H for , 4H com , 4HBINE ,  &
     4H ope ,   4HRATI , 4HON   , 10*4H           /
 DATA    nheqss/ 4HEQSS /
 DATA    papp  , loap,lods/ 4HPAPP , 4HLOAP , 4HLODS /
 
!     OPEN CASECC DATA BLOCK AND READ INTO OPEN CORE
 
 srch = .false.
 ierr = 0
 DO  i = 1,96
   ihead(i) = ihd(i)
 END DO
 ifile = casecc
 CALL OPEN (*580,casecc,z(buf2),0)
 nrec = step
 IF (nrec == 0) GO TO 30
 DO  i = 1,nrec
   CALL fwdrec (*580,casecc)
 END DO
 30 CALL READ (*570,*590,casecc,z(1),5,0,nnn)
 i = 2
 nwdscc = z(i  )
 npsub  = z(i+1)
 CALL READ (*570,*40,casecc,z(1),nwdscc,1,nnn)
 40 jj = 0
 kk = 0
 iprint = 0
 
!     INITIALIZE COMBO AND RESTCT ARRAYS
 
 DO  i = 1,7
   DO  j = 1,5
     combo(i,j) = 0
   END DO
   DO  j = 1,7
     restct(i,j) = 0
   END DO
 END DO
 
!     INITIALIZE COMP,TRANS,AND SYMT ARRAYS
 
 conect = .false.
 tran   = .false.
 DO  i = 1,7
   symt(i) = 0
   trans(i)= 0
   DO  j = 1,2
     comp(i,j) = 0
   END DO
 END DO
 DO  i = 1,3
   lf(i) = .false.
 END DO
 cnam(1) = 0
 cnam(2) = 0
 
!     PROCESS CASE CONTROL MNEMONICS
 
 loop350:  DO  i = 1,nwdscc,3
   DO  j = 1,nmnem
     IF (z(i) /= mnem(j)) CYCLE
     SELECT CASE ( j )
       CASE (    1)
         GO TO 120
       CASE (    2)
         GO TO 130
       CASE (    3)
         GO TO 160
       CASE (    4)
         GO TO 170
       CASE (    5)
         GO TO 180
       CASE (    6)
         GO TO 190
       CASE (    7)
         GO TO 200
       CASE (    8)
         GO TO 220
       CASE (    9)
         GO TO 230
       CASE (   10)
         GO TO 260
       CASE (   11)
         GO TO 290
     END SELECT
   END DO
   CYCLE loop350
   120 iauto = .false.
   IF (z(i+1) == auto) iauto = .true.
   CYCLE loop350
   
   130 DO  l = 1,3
     IF (z(i+1) == idir(l)) GO TO 150
   END DO
   isort = 1
   CYCLE loop350
   150 isort = l
   CYCLE loop350
   
   160 IF (lf(1)) GO TO 300
   lf(1)   = .true.
   cnam(1) = z(i+1)
   cnam(2) = z(i+2)
   CYCLE loop350
   
   170 jj = jj + 1
   snam(jj,1) = z(i+1)
   snam(jj,2) = z(i+2)
   CYCLE loop350
   
   180 IF (lf(2)) GO TO 300
   lf(2) = .true.
   toler = az(i+2)
   CYCLE loop350
   
   190 IF (lf(3)) GO TO 300
   lf(3)  = .true.
   conset = z(i+2)
   conect = .true.
   CYCLE loop350
   
   200 kk = kk + 1
   comp(kk,1) = z(i+1)
   comp(kk,2) = z(i+2)
   DO  lindx = 1,npsub
     IF (z(i+1) == snam(lindx,1) .AND. z(i+2) == snam(lindx,2)) CYCLE loop350
   END DO
   WRITE (outt,630) ufm,z(i+1),z(i+2)
   ierr = 1
   CYCLE loop350
   
   220 trans(kk) = z(i+2)
   tran = .true.
   CYCLE loop350
   
   230 DO  l = 1,15
     IF (z(i+1) == isym(l,2)) GO TO 250
   END DO
   ierr = 1
   WRITE (outt,620) ufm,z(i+1),comp(kk,1),comp(kk,2)
   CYCLE loop350
   250 symt(kk) = isym(l,1)
   CYCLE loop350
   
   260 DO  l = 1,npsub
     IF (z(i+1) == snam(l,1) .AND. z(i+2) == snam(l,2)) GO TO 280
   END DO
   WRITE (outt,630) ufm,z(i+1),z(i+2)
   ierr = 1
   CYCLE loop350
   280 srch = .true.
   restct(lindx,l) = 1
   restct(l,lindx) = 1
   CYCLE loop350
   
   290 iprint = orf(iprint,z(i+2))
   CYCLE loop350
   
   300 SELECT CASE ( j )
     CASE (    1)
       GO TO 350
     CASE (    2)
       GO TO 350
     CASE (    3)
       GO TO 310
     CASE (    4)
       GO TO 350
     CASE (    5)
       GO TO 320
     CASE (    6)
       GO TO 330
   END SELECT
   310 WRITE (outt,740) ufm
   GO TO 340
   320 WRITE (outt,750) ufm
   GO TO 340
   330 WRITE (outt,760) ufm
   340 ierr = 1
 END DO loop350
 
!     IF NO SEARCH OPTIONS SPECIFIED - SEARCH ALL POSSIBLE CONNECTIONS
 
 IF (srch) GO TO 370
 DO  i = 1,7
   DO  j = 1,7
     restct(i,j) = 1
   END DO
 END DO
 370 CONTINUE
 DO  i = 1,npsub
   DO  j = 1,npsub
     IF (snam(i,1) == comp(j,1) .AND. snam(i,2) == comp(j,2)) GO TO 390
   END DO
   combo(i,1) = snam(i,1)
   combo(i,2) = snam(i,2)
   combo(i,3) = 0
   combo(i,4) = 0
   CYCLE
   390 combo(i,1) = snam(i,1)
   combo(i,2) = snam(i,2)
   combo(i,3) = trans(j)
   combo(i,4) = symt(j)
 END DO
 CALL CLOSE (casecc,1)
 CALL page
 WRITE (outt,690) npsub
 IF (iauto) WRITE (outt,700)
 IF (.NOT. iauto) WRITE (outt,710)
 IF (.NOT.(iauto .OR. conect)) GO TO 550
 410 IF (conect) WRITE (outt,720) conset
 IF (cnam(1) == 0 .AND. cnam(2) == 0) GO TO 560
 WRITE (outt,640) cnam
 CALL fdsub (cnam,itest)
 IF (itest /= -1) GO TO 500
 IF (pora == papp) GO TO 540
 420 IF (.NOT.lf(2)) GO TO 530
 WRITE (outt,670) toler
 CALL decode (iprint,ibits,nflg)
 IF (nflg == 0) ibits(1) = 0
 IF (nflg == 0) GO TO 440
 DO  i = 1,nflg
   ibits(i) = ibits(i) + 1
 END DO
 440 CONTINUE
 WRITE (outt,810) (ibits(kdh),kdh=1,nflg )
 450 DO  i = 1,npsub
   WRITE (outt,770) i,combo(i,1),combo(i,2)
   ncnam(1) = combo(i,1)
   ncnam(2) = combo(i,2)
   CALL sfetch (ncnam,nheqss,3,itest)
   IF (itest == 4) WRITE (outt,780) ufm,ncnam
   IF (itest == 4) idry = -2
   IF (combo(i,3) /= 0) WRITE (outt,790) combo(i,3)
   IF (combo(i,4) == 0) CYCLE
   DO  mj = 1,15
     IF (combo(i,4) == isym(mj,1)) EXIT
   END DO
   470 WRITE (outt,800) isym(mj,2)
 END DO
 490 IF (ierr == 1) idry = -2
 GO TO 610
 500 litm = lods
 IF (pora == papp) litm = loap
 CALL sfetch (cnam,litm,3,itest)
 lonly = .false.
 IF (itest == 3) GO TO 520
 IF (pora == papp) GO TO 510
 WRITE (outt,650) ufm
 ierr = 1
 GO TO 420
 
!     OPTIONS PA YET LOAP ITEM ALREADY EXISTS
 
 510 WRITE (outt,820) ufm,cnam
 ierr = 1
 GO TO 490
 
!     NEW LODS ONLY DEFINED
 
 520 lonly = .true.
 RETURN
 
 530 WRITE (outt,660) ufm
 ierr = 1
 GO TO 450
 
!     OPTIONS PA YET SUBSTRUCTURE DOES NOT EXIST
 
 540 WRITE (outt,830) ufm,cnam
 ierr = 1
 GO TO 490
 550 WRITE (outt,680) ufm
 ierr = 1
 GO TO 410
 560 WRITE (outt,730) ufm
 ierr = 1
 GO TO 490
 570 imsg = -2
 GO TO 600
 580 imsg = -1
 GO TO 600
 590 imsg = -3
 600 CALL mesage (imsg,ifile,aaa)
 610 CONTINUE
 RETURN
 
 620 FORMAT (a23,' 6505, THE SYMMETRY OPTION ',a4,  &
     ' CONTAINS AN INVALID SYMBOL.')
 630 FORMAT (a23,' 6506, THE COMPONENT SUBSTRUCTURE ',2A4,  &
     ' IS NOT ONE OF THOSE ON THE COMBINE CARD.')
 640 FORMAT (/10X,38HTHE resultant pseudostructure NAME is ,2A4)
 650 FORMAT (a23,' 6508, THE NAME SPECIFIED FOR THE RESULTANT ',  &
     'PSEUDOSTRUCTURE', /32X,'ALREADY EXISTS ON THE SOF.')
 660 FORMAT (a23,' 6504, A TOLERANCE MUST BE SPECIFIED FOR A COMBINE ',  &
     'OPERATION.')
 670 FORMAT (/10X,32HTHE tolerance on connections is ,e15.6)
 680 FORMAT (a23,' 6501, THE MANUAL COMBINE OPTION HAS BEEN SPECIFIED',  &
     ', BUT NO CONNECTION SET WAS GIVEN.')
 690 FORMAT (/10X,'THIS JOB STEP WILL COMBINE ',i1,' PSEUDOSTRUCTURES')
 700 FORMAT (/10X,40HCONNECTIONS are generated automatically. )
 710 FORMAT (/10X,35HCONNECTIONS are specified manually. )
 720 FORMAT (/10X,25HTHE connection set id is ,i8)
 730 FORMAT (a23,' 6502, NO NAME HAS BEEN SPECIFIED FOR THE RESULTANT',  &
     ' COMBINED PSEUDOSTRUCTURE.')
 740 FORMAT (a23,' 6519, REDUNDANT NAMES FOR RESULTANT PSEUDOSTRUCTURE'  &
     ,      ' HAVE BEEN SPECIFIED.')
 750 FORMAT (a23,' 6520, REDUNDANT VALUES FOR TOLER HAVE BEEN ', 'SPECIFIED.')
 760 FORMAT (a23,' 6512, REDUNDANT CONNECTION SET ID S HAVE BEEN ',  &
     'SPECIFIED.')
 770 FORMAT (/10X, 27HCOMPONENT substructure no. ,i1,8H NAME = ,2A4)
 780 FORMAT (a23,' 6507, THE SUBSTRUCTURE ',2A4,' DOES NOT EXIST ON ',  &
     'THE SOF FILE')
 790 FORMAT (/15X, 15HTRANS set id = ,i8)
 800 FORMAT (15X,22HSYMMETRY directions = ,a4)
 810 FORMAT (/10X,30HTHE PRINT control options are ,25I3)
 820 FORMAT (a23,' 6533, OPTIONS PA HAS BEEN SPECIFIED BUT THE LOAP ',  &
     'ITEM ALREADY EXISTS FOR SUBSTRUCTURE ',2A4)
 830 FORMAT (a23,' 6534, OPTIONS PA HAS BEEN SPECIFIED BUT THE ',  &
     'SUBSTRUCTURE ',2A4,' DOES NOT EXIST.', /30X,  &
'YOU CANNOT APPEND SOMETHING TO NOTHING.')
END SUBROUTINE cmcase
