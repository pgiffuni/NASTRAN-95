SUBROUTINE cmrd2c (iter)
     
!     THIS SUBROUTINE PERFORMS THE GUYAN REDUCTION ON THE STRUCTURE
!     POINTS FOR THE CMRED2 MODULE.
 
!     INPUT  DATA
!     GINO - KII    - KII PARTITION MATRIX
!            KIB    - KIB KIB PARTITION MATRIX
!     SOF  - GIMS   - G TRANSFORMATION MATRIX FOR BOUNDARY POINTS OF
!                      ORIGINAL SUBSTRUCTURE
 
!     OUTPUT DATA
!     SOF  - LMTX   - LII PARTITION MATRIX
!            GIMS   - G TRANSFORMATION MATRIX FOR BOUNDARY POINTS OF
!                    ORIGINAL SUBSTRUCTURE
 
!     PARAMETERS
!     INPUT- GBUF   - GINO BUFFER
!            ISCR   - SCRATCH FILE NUMBER ARRAY
!            KORLEN - LENGTH OF OPEN CORE
!            KORBGN - BEGINNING ADDRESS OF OPEN CORE
!            OLDNAM - NAME OF SUBSTRUCTURE BEGING REDUCED
!            RSAVE  - DECOMPOSITION SAVE FLAG
!     OTHERS-KII    - KII PARTITION MATRIX FILE NUMBER
!            LII    - LII PARTITION MATRIX FILE NUMBER
!            SYMTRY - KII SYMMETRY FLAG
 
 
 INTEGER, INTENT(IN OUT)                  :: iter
 LOGICAL :: rsave,restor,symtry
 INTEGER :: dry,sbuf1,sbuf2,sbuf3,oldnam,z,power,chlsky,uiitc,  &
     scr,powerc,b,bbar,u,gibt,prec,SIGN,ugfbs,gibfbs,  &
     prec1,atrlr,attrlr,gib,uii,dblkor,upper,him
 DOUBLE PRECISION :: detr,deti,mindia,det,mindc,dz
 DIMENSION       itrlr(7),modnam(2),itmlst(3),dz(1)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / idum1,dry,idum4(4),sbuf1,sbuf2,sbuf3,idum3(11),  &
     otfile(6),iscr(11),korlen,korbgn,oldnam(2), idum2(8),rsave,idum6(4),lstzwd
 COMMON /zzzzzz/ z(1)
 COMMON /sfact / kiit(7),liit(7),iscrq(7),iscra,iscrb,nzsf,  &
     detr,deti,power,iscrc,mindia,chlsky
 COMMON /cdcmpx/ kiitc(7),liitc(7),uiitc(7),scr(3),det(2),powerc,  &
     nx,mindc,b,bbar
 COMMON /fbsx  / liifbs(7),u(7),kibt(7),gibt(7),nzfbs,prec,SIGN
 COMMON /gfbsx / ligfbs(7),ugfbs(7),kigfbs(7),gibfbs(7),nzgfbs, prec1,ISIGN
 COMMON /trnspx/ atrlr(7),attrlr(7),lcore,nscrth,iscrth(8)
 COMMON /system/ idum5,iprntr
 EQUIVALENCE     (kib,iscr(2)),(kbi,iscr(3)),(kii,iscr(4)),  &
     (lii,iscr(8)),(uii,iscr(9)),(him,iscr(10)), (gib,iscr(11)),(dz(1),z(1))
 DATA    modnam/ 4HCMRD,4H2C  /
 DATA    lower , upper /4,5   /
 DATA    itmlst/ 4HLMTX,4HGIMS,4HHORG/
 
!     PREFORM GUYAN REDUCTION
 
 IF (dry == -2) GO TO 130
 restor = .false.
 
!     TRANSPOSE KII, KBI
 
 IF (iter == 1) GO TO 8
 IF (symtry) GO TO 37
 1 dblkor = (korbgn/2) + 1
 lcore  = lstzwd - ((2*dblkor)-1)
 nscrth = 5
 DO  i = 1,nscrth
   iscrth(i) = iscr(4+i)
 END DO
 DO  i = 1,2
   itrlr(1) = kii
   IF (i == 2) itrlr(1) = kbi
   CALL rdtrl (itrlr)
   DO  j = 1,7
     atrlr(j) = itrlr(j)
     attrlr(j) = itrlr(j)
   END DO
   attrlr(2) = itrlr(3)
   attrlr(3) = itrlr(2)
   CALL trnsp (dz(dblkor))
   CALL wrttrl (attrlr)
 END DO
 IF (restor) CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 IF (restor) GO TO 39
 restor = .true.
 CALL sofcls
 GO TO 16
 
!     DECOMPOSE INTERIOR STIFFNESS MATRIX
!        (SYMMETRIC)
 
!                                 T
!        **   **   **   ** **   **
!        *     *   *     * *     *
!        * KII * = * LII * * LII *
!        *     *   *     * *     *
!        **   **   **   ** **   **
 
 8 CALL sofcls
 kiit(1) = kii
 CALL rdtrl (kiit)
 IF (kiit(4) /= 6) GO TO 12
 symtry = .true.
 iprc = 1
 ityp = 0
 IF (kiit(5) == 2 .OR. kiit(5) == 4) iprc = 2
 IF (kiit(5) >= 3) ityp = 2
 itype = iprc + ityp
 CALL makmcb (liit,lii,kiit(3),lower,itype)
 iscrq(1) = iscr(5)
 iscra  = iscr(6)
 iscrb  = iscr(7)
 iscrc  = iscr(9)
 chlsky = 0
 power  = 1
 dblkor = (korbgn/2) + 1
 nzsf   = lstzwd - ((2*dblkor)-1)
 CALL sdcomp (*40,dz(dblkor),dz(dblkor),dz(dblkor))
 CALL wrttrl (liit)
 GO TO 18
 
!     DECOMPOSE INTERIOR STIFFNESS MATRIX
!        (UNSYMMETRIC)
 
!        **   **   **   ** **   **
!        *     *   *     * *     *
!        * KII * = * LII * * UII *
!        *     *   *     * *     *
!        **   **   **   ** **   **
 
 12 symtry = .false.
 16 kiitc(1) = kii
 CALL rdtrl (kiitc)
 ityp = 0
 iprc = 1
 IF (kiitc(5) == 2 .OR. kiitc(5) == 4) iprc = 2
 IF (kiitc(5) >= 3) ityp = 2
 itype = iprc + ityp
 CALL makmcb (liitc,lii,kiitc(3),lower,itype)
 CALL makmcb (uiitc,uii,kiitc(3),upper,itype)
 scr(1) = iscr(5)
 scr(2) = iscr(6)
 scr(3) = iscr(7)
 b      = 0
 bbar   = 0
 dblkor = (korbgn/2) + 1
 nx = lstzwd - ((2*dblkor)-1)
 CALL cdcomp (*42,dz(dblkor),dz(dblkor),dz(dblkor))
 CALL wrttrl (liitc)
 CALL wrttrl (uiitc)
 
!     SAVE LII AS LMTX ON SOF
 
 18 IF (iter == 2 .OR. .NOT.rsave) GO TO 20
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 ifile = lii
 CALL mtrxo (lii,oldnam,itmlst(1),0,itest)
 item = itmlst(1)
 IF (itest /= 3) GO TO 70
 CALL sofcls
 
!     SOLVE STRUCTURE REDUCTION TRANSFORMATION MATRIX
!        (SYMMETRIC)
 
!                       T
!        **   ** **   ** **   **    **   **
!        *     * *     * *     *    *     *
!        * LII * * LII * * GIB * = -* KIB *
!        *     * *     * *     *    *     *
!        **   ** **   ** **   **    **   **
 
 20 IF (.NOT.symtry) GO TO 32
 kibt(1) = kib
 IF (iter == 2) kibt(1) = kbi
 CALL rdtrl (kibt)
 DO  i = 1,7
   liifbs(i) = liit(i)
 END DO
 iprc = 1
 ityp = 0
 IF (kibt(5) == 2 .OR. kibt(5) == 4) iprc = 2
 IF (liit(5) == 2 .OR. liit(5) == 4) iprc = 2
 IF (kibt(5) >= 3) ityp = 2
 IF (liit(5) >= 3) ityp = 2
 itype = iprc + ityp
 CALL makmcb (gibt,gib,kibt(3),kibt(4),itype)
 nzfbs = lstzwd - ((2*dblkor)-1)
 prec  = kibt(5) - 2
 SIGN  = -1
 CALL fbs (dz(dblkor),dz(dblkor))
 CALL wrttrl (gibt)
 GO TO 36
 
!     SOLVE STRUCTURE REDUCTION TRANSFORMATION MATRIX
!        (UNSYMMETRIC)
 
!        **   ** **   ** **   **    **   **
!        *     * *     * *     *    *     *
!        * LII * * UII * * GIB * = -* KIB *
!        *     * *     * *     *    *     *
!        **   ** **   ** **   **    **   **
 
 32 kigfbs(1) = kib
 IF (iter == 2) kigfbs(1) = kbi
 CALL rdtrl (kigfbs)
 DO  i = 1,7
   ligfbs(i) = liitc(i)
   ugfbs(i)  = uiitc(i)
 END DO
 iprc = 1
 ityp = 0
 IF (kigfbs(5) == 2 .OR. kigfbs(5) == 4) iprc = 2
 IF (liitc(5) == 2 .OR. liitc(5) == 4) iprc = 2
 IF (uiitc(5) == 2 .OR. uiitc(5) == 4) iprc = 2
 IF (kigfbs(5) >= 3) ityp = 2
 IF (liitc(5)  >= 3) ityp = 2
 IF (uiitc(5)  >= 3) ityp = 2
 itype = iprc + ityp
 CALL makmcb (gibfbs,gib,kigfbs(3),kigfbs(4),itype)
 nzgfbs = lstzwd - ((2*dblkor)-1)
 prec1  = iprc
 ISIGN  = -1
 CALL gfbs (dz(dblkor),dz(dblkor))
 CALL wrttrl (gibfbs)
 
!     SAVE GIB AS GIMS ON SOF
 
 36 IF (restor) GO TO 1
 IF (iter == 2) GO TO 39
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 ifile = gib
 CALL mtrxo (gib,oldnam,itmlst(2),0,itest)
 item = itmlst(2)
 IF (itest /= 3) GO TO 70
 GO TO 39
 
!     KII SYMMETRIC, GIBBAR = GIB
 
 37 item = itmlst(2)
 CALL mtrxi (gibbar,oldnam,item,0,itest)
 IF (itest /= 1) GO TO 70
 39 CONTINUE
 GO TO 130
 
!     PROCESS SYSTEM FATAL ERRORS
 
 40 WRITE (iprntr,903) uwm,oldnam
 GO TO 44
 42 WRITE (iprntr,904) uwm,oldnam
 44 imsg  = -37
 ifile = 0
 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 GO TO 130
 
!     PROCESS MODULE FATAL ERRORS
 
 70 SELECT CASE ( itest )
   CASE (    1)
     GO TO 80
   CASE (    2)
     GO TO 82
   CASE (    3)
     GO TO 84
   CASE (    4)
     GO TO 90
   CASE (    5)
     GO TO 100
   CASE (    6)
     GO TO 120
 END SELECT
 80 WRITE (iprntr,900) ufm,modnam,item,oldnam
 dry = -2
 GO TO 130
 
 82 WRITE (iprntr,902) ufm,modnam,item,oldnam
 dry = -2
 GO TO 130
 
 84 imsg = -1
 GO TO 110
 90 imsg = -2
 GO TO 110
 100 imsg = -3
 110 CALL smsg (imsg,item,oldnam)
 GO TO 130
 
 120 WRITE (iprntr,901) ufm,modnam,item,oldnam
 dry = -2
 130 RETURN
 
 900 FORMAT (a23,' 6211, MODULE ',2A4,' - ITEM ',a4,  &
     ' OF SUBSTRUCTURE ',2A4,' HAS ALREADY BEEN WRITTEN.')
 901 FORMAT (a23,' 6632, MODULE ',2A4,' - NASTRAN MATRIX FILE FOR I/O',  &
     ' OF SOF ITEM ',a4,', SUBSTRUCTURE ',2A4,', IS PURGED.')
 902 FORMAT (a23,' 6215, MODULE ',2A4,' - ITEM ',a4,  &
     ' OF SUBSTRUCTURE ',2A4,' PSEUDO-EXISTS ONLY.')
 903 FORMAT (a25,' 6311, SDCOMP DECOMPOSITION FAILED ON KII MATRIX ',  &
     'FOR SUBSTRUCTURE ',2A4)
 904 FORMAT (a23,' 6635, CDCOMP DECOMPOSITION FAILED ON KII MATRIX ',  &
     'FOR SUBSTRUCTURE ',2A4)
 
END SUBROUTINE cmrd2c
