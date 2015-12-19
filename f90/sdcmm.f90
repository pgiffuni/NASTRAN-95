SUBROUTINE sdcmm (z,mset,msze,matrix,uset,gpl,sil,subnam)
     
!     THIS ROUTINE WRITES THE EXTERNAL ID AND COMPONENT ID FOR VARIOUS
!     MATRIX ERROR CONDITIONS.
!     SCRATCH1 CONTAINS 3 WORDS/ERROR, EACH MESSAGE BEING 1 RECORD
!         WORD 1 = COLUMN * 10 + ERROR CODE
!         WORD 2 = INPUT DIAGONAL
!         WORD 3 = OUTPUT DIAGONAL
!     SUBROUTINE -MXCID- (NON-SUBSTRUCTURING) IS CALLED TO SUPPLY IDENT.
!     DATA FOR EACH COLUMN.  FOR SUBSTRUCTURING -MXCIDS- IS CALLED - IT
!     RETURNS TWO WORDS/COLUMN PLUS THE BCD NAME OF THE SUBSTRUCTURES AT
!     THE START OF CORE.  IN EITHER CASE, THE 1ST WORD IS 10*ID +
!         COMPONENT.
!     THE SCRATCH FILE IS READ AND THE EXTERNAL ID INDEXED DIRECTLY.E
!     NOTE - THAT EACH COLUMN MAY GENERATE MORE THAN 1 MESSAGE.
!     OPEN CORE IS Z(1) TO Z(BUF-1).  TWO BUFFERS FOLLOW Z(BUF)
 
 
 INTEGER, INTENT(IN)                      :: z(1)
 INTEGER, INTENT(IN OUT)                  :: mset
 INTEGER, INTENT(IN OUT)                  :: msze
 INTEGER, INTENT(IN OUT)                  :: matrix
 INTEGER, INTENT(IN OUT)                  :: uset
 INTEGER, INTENT(IN OUT)                  :: gpl
 INTEGER, INTENT(IN OUT)                  :: sil
 INTEGER, INTENT(IN OUT)                  :: subnam(2)
 INTEGER :: buf,buf2, EXIT(8),ERR(14),filmsg,gpid(4),  &
            in(3),iner(4),n(7),NAME(2), typ(6)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBNB 8/94
 CHARACTER (LEN=4) :: ctyp(6)
 REAL :: xgpid(4), xin(3)
!WKBNE 8/94
!WKBI  8/94
 EQUIVALENCE     (ctyp,typ), (xgpid,gpid), (xin,in)
 COMMON /xmssg / ufm,uwm
 COMMON /sdcq  / nerr(2),noglev,buf,filmsg
 COMMON /names / krd2,krr0,skpn(3), kcl2
 COMMON /system/ ksystm(69)
 EQUIVALENCE     (ksystm(1),nbufsz),(ksystm(2),iout), (ksystm(69),isubst)
 DATA    ERR   / 4HNULL, 4HCOL., 4HZERO, 4HDIAG, 4HNEG., 4HDIAG,  &
     4HSING, 4HTEST, 4HBAD , 4HCOL., 4HNON-, 4HCONS, 4HZERO, 4HDIAG/
 DATA    iner  / 4HINPU, 2HT   , 4HDECM, 2HP   /
 DATA    EXIT  / 4HCONT, 4HINUE, 4HAT e, 4HND  , 4HAT s, 4HTART,  &
     4HIN d, 4HECMP/
 DATA    NAME  / 4HSDCM, 2HM   /
 DATA    iblk  / 4H    /
 
 buf2 = buf + nbufsz
 n(1) = 0
 n(2) = 0
 n(3) = 0
 n(4) = 0
 n(5) = 0
 n(6) = 0
 n(7) = 0
 IF (buf <= 0) GO TO 50
 
!     GENERATE EXTERNAL ID
 
 IF (isubst == 0) GO TO 5
 
!     SUBSTRUCTURING - READ EQSS FILE ON THE SOF
 
!     4 BUFFERS NEEDED
 
 i = buf - 2*nbufsz
 IF (i <= 3*msze) GO TO 50
 nwds = 2
 CALL mxcids (*50,z,mset,msze,nwds,uset,i,subnam)
 nstart = i - 1
 GO TO  7
 
 5 nstart = 0
 nwds   = 1
 
!     2 BUFFERS NEEDED
 
 CALL mxcid (*50,z,mset,msze,nwds,uset,gpl,sil,buf)
 
 7 CALL OPEN (*110,filmsg,z(buf2),krr0)
 CALL page2 (3)
 WRITE  (iout,10) uwm
 10 FORMAT (a25,' 2377A, MATRIX CONDITIONING ERRORS GIVEN WITH ',  &
     'EXTERNAL ID', /5X,'GID - C  INPUT-DIAG.   DECOMP-DIAG.',  &
     6X,'TYPE',17X,'SUBSTRUCTURE')
 
 ASSIGN 30 TO iret
 typ(5) = iblk
 typ(6) = iblk
 IF (isubst /= 0) ASSIGN 27 TO iret
 
!     LOOP ON MESSAGES - 0 COLUMN IS FLAG TO QUIT
 
 20 CALL fread (filmsg,in,3,1)
 IF (in(1) == 0) GO TO 200
 i = in(1)/10
 j = in(1) - i*10
 l = nstart + i*nwds
 gpid(1) = z(l)/10
 gpid(2) = z(l) - gpid(1)*10
 gpid(3) = in(2)
 gpid(4) = in(3)
 
!     INTERNAL FUNCTION
 
 25 CONTINUE
 IF (j <= 0 .OR. j > 7) GO TO 100
 k = 2*j - 1
 typ(1) = ERR(k)
 typ(2) = ERR(k+1)
 k = 1
 IF (j > 1 .AND. j < 7) k = 3
 typ(3) = iner(k  )
 typ(4) = iner(k+1)
 n(j) = n(j) + 1
 CALL page2 (2)
 GO TO iret, (27,30,80)
 
 27 typ(5) = z(2*l-1)
 typ(6) = z(2*l  )
 30 CONTINUE
!WKBR 8/94      WRITE  (IOUT,40) GPID,TYP
 WRITE ( iout, 40 ) gpid(1), gpid(2), xgpid(3), xgpid(4), ctyp
 40 FORMAT (1H0,i9,2H -,i2,1P,2E14.6,3X,2A5,2H/ ,a4,a2,6HMATRIX,2X, 2A4)
 GO TO 20
 
!     INSUFFICIENT CORE IN -MATCID-
 
 50 CALL page2 (3)
 WRITE  (iout,60) uwm
 60 FORMAT (a25,' 2377B, MATRIX CONDITIONING ERRORS GIVEN WITH ',  &
     'INTERNAL ID', /,5X,'COLUMN  INPUT DIAG.   DECOMP-DIAG.', 6X,'TYPE')
 
 CALL OPEN (*110,filmsg,z(buf2),krr0)
 ASSIGN 80 TO iret
 
!     LOOP
 
 70 CONTINUE
 CALL fread (filmsg,in,3,1)
 IF (in(1) == 0) GO TO 200
 i = in(1)/10
 j = in(1) - i*10
 in(1) = i
 GO TO 25
 
 80 CONTINUE
!WKBR 8/94  WRITE  (IOUT,90) IN,TYP
 WRITE  (iout,90) in(1), xin(2), xin(3), ctyp
 90 FORMAT (1H0,i8,1P,2E14.6,3X,2A5,2H/ ,a4,a2,6HMATRIX,2X,2A4)
 GO TO 70
 
!     ILLEGAL DATA
 
 100 CALL mesage (7,filmsg,NAME)
 GO TO 200
 
!     SCRATCH FILE NOT AVAILABLE
 
 110 CALL mesage (1,filmsg,NAME)
 
!     ALL DONE, SUMMARIZE
 
 200 CALL page2 (11)
 WRITE  (iout,210) matrix,msze,n
 210 FORMAT (1H0,3X,10HFOR matrix,i4,6H, size,i8,/i9,13H null columns,  &
     /i9,15H zero diagonals, /i9,19H negative diagonals, /i9,  &
     31H singularity tolerance exceeded, /i9,12H bad columns, /i9,  &
     24H nonconservative columns, /i9,23H zero diagonals (INPUT))
 
!     CHECK FOR EXIT CONDITIONS
 
 i = 2*noglev + 1
 
!     NOTE - NOGLEV OF 4 ALSO HAS NEGATIVE PARM(1)
 
 IF (noglev == 4) i = 7
 j = i + 1
 WRITE  (iout,220) EXIT(i),EXIT(j)
 220 FORMAT (1H0,3X,13HABORT code = ,2A4)
 CALL CLOSE (filmsg,kcl2)
 RETURN
END SUBROUTINE sdcmm
