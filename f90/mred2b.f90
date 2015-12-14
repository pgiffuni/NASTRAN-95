SUBROUTINE mred2b
     
!     THIS SUBROUTINE PERFORMS THE GUYAN REDUCTION ON THE STRUCTURE
!     POINTS FOR THE MRED2 MODULE.
 
!     INPUT DATA
!     GINO   - KII    - KII PARTITION MATRIX
!     SOF    - GIMS   - G TRANSFORMATION MATRIX FOR BOUNDARY POINTS OF
!                       ORIGINAL SUBSTRUCTURE
 
!     OUTPUT DATA
!     GINO   - LII    - LII PARTITION MATRIX
!     SOF    - LMTX   - LII PARTITION MATRIX
!              GIMS   - G TRANSFORMATION MATRIX FOR BOUNDARY POINTS OF
!                       ORIGINAL SUBSTRUCTURE
 
!     PARAMETERS
!     INPUT  - GBUF   - GINO BUFFER
!              ISCR   - SCRATCH FILE NUMBER ARRAY
!              KORLEN - LENGTH OF OPEN CORE
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              OLDNAM - NAME OF SUBSTRUCTURE BEGING REDUCED
!              BOUNDS - OLDBOUNDS OPTION FLAG
!              RSAVE  - DECOMPOSITION SAVE FLAG
!     OTHERS - KIB    - KIB PARTITION MATRIX FILE NUMBER
!              KII    - KII PARTITION MATRIX FILE NUMBER
!              LII    - LII PARTITION MATRIX FILE NUMBER (ISCR11)
 
 LOGICAL :: bounds,rsave
 INTEGER :: dry,sbuf1,sbuf2,sbuf3,oldnam,z,power,chlsky,u,  &
     gibt,prec,SIGN,gib,dblkor,dmr
 DOUBLE          PRECISION detr,deti,mindia,dz
 DIMENSION       itrlr(7),modnam(2),itmlst(2),dz(1)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /BLANK / idum1,dry,idum4(4),sbuf1,sbuf2,sbuf3,infile(12),  &
     idum6(6),iscr(10),korlen,korbgn,oldnam(2),  &
     idum2(8),bounds,idum3,rsave,idum7(4),lstzwd,iscr11
 COMMON /zzzzzz/ z(1)
 COMMON /sfact / kiit(7),liit(7),iscrq(7),iscra,iscrb,nzsf,  &
     detr,deti,power,iscrc,mindia,chlsky
 COMMON /fbsx  / liifbs(7),u(7),kibt(7),gibt(7),nzfbs,prec,SIGN
 COMMON /system/ idum5,iprntr
 EQUIVALENCE     (dmr,infile(11)),(gib,iscr(6)),(dz(1),z(1)),  &
     (kib,iscr(2)),(kii,iscr(3)),(lii,iscr11)
 DATA    modnam/ 4HMRED,4H2B  /
 DATA    lower / 4 /
 DATA    itmlst/ 4HLMTX,4HGIMS/
 
!     TEST FOR GUYAN REDUCTION
 
 IF (dry == -2) GO TO 140
 IF (.NOT.bounds) GO TO 10
 itrlr(1) = dmr
 CALL rdtrl (itrlr)
 IF (itrlr(1) < 0) GO TO 35
 item = itmlst(1)
 CALL softrl (oldnam,item,itrlr)
 IF (itrlr(1) == 1) GO TO 35
 
!     DECOMPOSE INTERIOR STIFFNESS MATRIX
 
!                                 T
!        **   **   **   ** **   **
!        *     *   *     * *     *
!        * KII * = * LII * * LII *
!        *     *   *     * *     *
!        **   **   **   ** **   **
 
 10 CALL sofcls
 kiit(1) = kii
 CALL rdtrl (kiit)
 CALL makmcb (liit,lii,kiit(3),lower,kiit(5))
 iscrq(1) = iscr(6)
 iscra  = iscr(7)
 iscrb  = iscr(8)
 iscrc  = iscr(9)
 power  = 1
 chlsky = 0
 dblkor = 1 + korbgn/2
 nzsf   = lstzwd - 2*dblkor - 1
 CALL sdcomp (*40,dz(dblkor),dz(dblkor),dz(dblkor))
 CALL wrttrl (liit)
 
!     SAVE LII AS LMTX ON SOF
 
 IF (.NOT. rsave) GO TO 20
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 ifile = lii
 item  = itmlst(1)
 CALL mtrxo (lii,oldnam,item,0,itest)
 IF (itest /= 3) GO TO 70
 IF (bounds) GO TO 35
 CALL sofcls
 
!     SOLVE STRUCTURE REDUCTION TRANSFORMATION MATRIX
 
!                       T
!        **   ** **   ** **   **    **   **
!        *     * *     * *     *    *     *
!        * LII * * LII * * GIB * = -* KIB *
!        *     * *     * *     *    *     *
!        **   ** **   ** **   **    **   **
 
 20 IF (bounds) GO TO 32
 kibt(1) = kib
 CALL rdtrl (kibt)
 DO  i = 1,7
   liifbs(i) = liit(i)
 END DO
 CALL makmcb (gibt,gib,kibt(3),kibt(4),kibt(5))
 nzfbs = lstzwd - 2*dblkor
 prec  = kibt(5)
 SIGN  = -1
 CALL fbs (dz(dblkor),dz(dblkor))
 CALL wrttrl (gibt)
 
!     SAVE GIB AS GIMS ON SOF
 
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 ifile = gib
 item  = itmlst(2)
 CALL mtrxo (gib,oldnam,item,0,itest)
 IF (itest /= 3) GO TO 70
 GO TO 35
 32 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 35 CONTINUE
 GO TO 140
 
!     PROCESS SYSTEM FATAL ERRORS
 
 40 WRITE  (iprntr,45) swm,oldnam
 45 FORMAT (a27,' 6311, SDCOMP DECOMPOSITION FAILED ON KII MATRIX ',  &
     'FOR SUBSTRUCTURE ',2A4)
 imsg  = -37
 ifile = 0
 CALL mesage (imsg,ifile,modnam)
 GO TO 140
 
!     PROCESS MODULE FATAL ERRORS
 
 70 SELECT CASE ( itest )
   CASE (    1)
     GO TO 80
   CASE (    2)
     GO TO 80
   CASE (    3)
     GO TO 80
   CASE (    4)
     GO TO 90
   CASE (    5)
     GO TO 100
   CASE (    6)
     GO TO 120
 END SELECT
 80 imsg = -9
 GO TO 130
 90 imsg = -2
 GO TO 110
 100 imsg = -3
 110 CALL smsg (imsg,item,oldnam)
 GO TO 140
 120 imsg = -10
 130 dry = -2
 CALL smsg1 (imsg,item,oldnam,modsam)
 140 RETURN
 
END SUBROUTINE mred2b
