SUBROUTINE scalar
     
!     CONVERTS MATRIX ELEMENT TO PARAMETER
 
!     SCALAR   MTX//C,N,ROW/C,N,COL/V,N,RSP/V,N,RDP/V,N,SPLX/V,N,DPLX $
 
!     INPUT GINO FILE
!       MTX = ANY MATRIX, S.P. OR D.P.; REAL OR COMPLEX
!     OUTPUT GINO FILE
!       NONE
!     INPUT PARAMETERS
!       ROW, COL = ROW AND COLUMN OF MTX (DEFAULT ARE 1,1)
!     OUTPUT PARAMETERS
!       RSP  = VALUE OF MTX(ROW,COL), REAL SINGLE PRECISION
!       RDP  = VALUE OF MTX(ROW,COL), REAL DOUBLE PRECISION
!       SPLX = VALUE OF MTX(ROW,COL), S.P. COMPLEX
!       DPLX = VALUE OF MTX(ROW,COL), D.P. COMPLEX
 
!     ORIGINALY WRITTEN BY R. MITCHELL, GSFC, NOV. 1972
 
!     COMPLETELY REWRITTEN BY G.CHAN/UNISYS IN JUNE 1988, SUCH THAT THE
!     OUTPUT PARAMETERS ARE SAVED CORRECTLY ACCORDING TO THEIR PRECISION
!     TYPES. (THE PRTPARM MODULE WILL BE ABLE TO PRINT THEM OUT
!     CORRECTLY.) PLUS IMPROVED MESSAGES (WHICH CAN BE SUPPRESSED BY
!     DIAG 37)
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: noprt
 INTEGER :: NAME(2),ia(7),fnm(2),pnm(2)
 REAL :: rsp,splx,a,vps(1),sp(4)
 DOUBLE PRECISION :: da(2),rdp,dplx(2),dp(2)
 CHARACTER (LEN=1) :: ufm*23,uwm*25,uim*29,TYPE(4)*10
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /xvps  /  ivps(1)
 COMMON /system/  sysbuf,nout
 COMMON /zntpkx/  a(4),ii,eol,eor
 COMMON /zzzzzz/  core(1)
 COMMON /BLANK /  bk(1),row,col,rsp,r2(2),splx(2),d4(4)
 EQUIVALENCE      (r2(1),rdp),(d4(1),dplx(1)),(ia(2),ncol),  &
     (ia(3),nrow),(ia(4),FORM),(ia(5),prec),  &
     (da(1),a(1)),(dp(1),sp(1)),(vps(1),ivps(1))
 DATA  in1,NAME/  101, 4HSCAL,4HAR  / , first / 12 /
 DATA  TYPE    /  'S.P. REAL ' ,    'D.P. REAL '   ,  &
     'S.P. CMPLX' ,    'D.P. CMPLX'   /
 
!     SUPPRESS ALL SCALAR MESSAGES IF DIAG 37 IS ON
 
 CALL sswtch (37,i)
 noprt = i == 1
 
!     MOVE VARIALBES IN /BLANK/ BY ONE WORD TO GET BY WORD BOUNDARY
!     ALIGNMENT SITUATION
 
 j = 12
 DO  i = 1,11
   bk(j) = bk(j-1)
   j = j - 1
 END DO
 
!     INITIALIZATION
 
 lcore = korsz(core)
 ibuf  = lcore - sysbuf + 1
 IF (ibuf < 1) GO TO 400
 rsp     = 0.
 splx(1) = 0.
 splx(2) = 0.
 rdp     = 0.d0
 dplx(1) = 0.d0
 dplx(2) = 0.d0
 dp(1)   = 0.d0
 dp(2)   = 0.d0
 CALL fname (in1,fnm)
 CALL page2 (first)
 first   = 3
 
!     GET STATUS OF INPUT MATRIX
!     CHECK FOR PURGED INPUT OR OUT OF RANGE INPUT PARAMETERS
 
 ia(1) = in1
 CALL rdtrl (ia)
 IF (ia(1) <   0) GO TO 410
 IF (row > nrow) GO TO 420
 
 SELECT CASE ( FORM )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 20
   CASE (    3)
     GO TO 40
   CASE (    4)
     GO TO 60
   CASE (    5)
     GO TO 70
   CASE (    6)
     GO TO 20
   CASE (    7)
     GO TO 50
   CASE (    8)
     GO TO 30
 END SELECT
!     SQUARE, RECTANGULAR OR SYMMETRIC MATRIX
 
 20 IF (col > ncol) GO TO 420
 GO TO 100
 
!     IDENTITY MATRIX
 
 30 IF (row /= col) GO TO 200
 rsp = 1.0
 rdp = 1.d0
 splx(1) = 1.
 dplx(1) = 1.d0
 GO TO 200
 
!     DIAGONAL MATRIX
 
 40 IF (row /=  col) GO TO 200
 IF (col > nrow) GO TO 420
!     SET COL TO 1 FOR SPECIAL DIAGONAL FORMAT
 col = 1
 GO TO 100
 
!     ROW VECTOR
!     SWITCH ROW AND COLUMN FOR PROPER INDEXING
 
 50 row = col
 col = 1
 GO TO 100
 
!     LOWER TRIANGULAR MATRIX (UPPER HALF= 0)
 
 60 IF (col-row > 0.0) THEN
   GO TO   200
 ELSE
   GO TO   100
 END IF
 
!     UPPER TRIANGULAR MATRIX (LOWER HALF= 0)
 
 70 IF (row-col > 0.0) THEN
   GO TO   200
 END IF
 
!     OPEN INPUT FILE AND SKIP HEADER RECORD AND UNINTERSTING COLUMNS
 
 100 CALL OPEN (*410,in1,core(ibuf),0)
 CALL skprec (in1,col)
 
!     READ AND SEARCH COLUMN CONTAINING DESIRED ELEMENT.
!     RECALL THAT DEFAULT VALUE WAS SET TO ZERO
 
 CALL intpk (*200,in1,0,prec,0)
 
!     FETCH ONE ELEMENT
!     CHECK FOR DESIRED ELEMENT
!     IF INDEX HIGHER, IT MEANS ELEMENT WAS 0.
 
 110 CALL zntpki
 IF (ii-row < 0) THEN
   GO TO   120
 ELSE IF (ii-row == 0) THEN
   GO TO   130
 ELSE
   GO TO   200
 END IF
 
!     CHECK FOR LAST NON-ZERO ELEMENT IN COLUMN.
 
 120 IF (eol > 0.0) THEN
   GO TO   200
 ELSE
   GO TO   110
 END IF
 
!     MOVE VALUES TO OUTPUT PARAMETER AREA.
!     CHECK PRECISION OF INPUT VALUE.
 
 130 SELECT CASE ( prec )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 150
   CASE (    3)
     GO TO 170
   CASE (    4)
     GO TO 180
 END SELECT
 
 140 rsp = a(1)
 rdp = DBLE(rsp)
 GO TO 160
 
 150 rdp = da(1)
 rsp = SNGL(rdp)
 160 splx(1) = rsp
 dplx(1) = rdp
 GO TO 200
 
 170 splx(1) = a(1)
 splx(2) = a(2)
 dplx(1) = DBLE(splx(1))
 dplx(2) = DBLE(splx(2))
 GO TO 190
 
 180 dplx(1) = da(1)
 dplx(2) = da(2)
 splx(1) = SNGL(dplx(1))
 splx(2) = SNGL(dplx(2))
 190 rsp = 0.0
 rdp = 0.d0
 
!     MOVE VALUES TO OUTPUT PARAMETERS AS REQUESTED BY USER, AND
!     SAVE PARAMETERS
 
 200 IF (noprt) GO TO 215
 CALL page2 (3)
 WRITE  (nout,210) uim
 210 FORMAT (a29,' FROM SCALAR MODULE -', /5X,  &
     '(ALL SCALAR MESSAGES CAN BE SUPPRESSED BY DIAG 37)')
 215 CALL fndpar (-3,j)
 IF (j <= 0) GO TO 260
 pnm(1) = ivps(j-3)
 pnm(2) = ivps(j-2)
 IF (prec >= 3) GO TO 240
 vps(j) = rsp
 IF (noprt) GO TO 260
 WRITE (nout,220) rsp,pnm
 220 FORMAT (73X,e15.8,4H  = ,2A4)
 WRITE (nout,230) row,col,TYPE(prec),fnm
 230 FORMAT (1H+,4X,'ELEMENT (',i5,'-ROW,',i5,'-COL) OF ',a10,' INPUT',  &
     ' FILE ',2A4,2H =)
 GO TO 260
 240 WRITE  (nout,250) uwm,pnm
 250 FORMAT (a25,' - INVALID OUTPUT REQUEST.', /5X,'ORIG. ELEM. IN ',  &
     'COMPLEX FORM. OUTPUT PARAMETER ',2A4,' NOT SAVED)',/)
 260 CALL fndpar (-4,j)
 IF (j <= 0) GO TO 290
 pnm(1)   = ivps(j-3)
 pnm(2)   = ivps(j-2)
 IF (prec >= 3) GO TO 280
 dp(1)    = rdp
 vps(j  ) = sp(1)
 vps(j+1) = sp(2)
 IF (noprt) GO TO 290
 WRITE (nout,270) rdp,pnm
 270 FORMAT (73X,d15.8,4H  = ,2A4)
 WRITE (nout,230) row,col,TYPE(prec),fnm
 GO TO 290
 280 WRITE (nout,250) uwm,pnm
 290 CALL fndpar (-5,j)
 IF (j <= 0) GO TO 310
 vps(j  ) = splx(1)
 vps(j+1) = splx(2)
 pnm(1)   = ivps(j-3)
 pnm(2)   = ivps(j-2)
 IF (noprt) GO TO 310
 WRITE (nout,300) splx,pnm
 300 FORMAT (73X,1H(,e15.8,1H,,e15.8,1H),4H  = ,2A4)
 WRITE (nout,230) row,col,TYPE(prec),fnm
 310 CALL fndpar (-6,j)
 IF (j <= 0) GO TO 330
 dp(1)    = dplx(1)
 dp(2)    = dplx(2)
 vps(j  ) = sp(1)
 vps(j+1) = sp(2)
 vps(j+2) = sp(3)
 vps(j+3) = sp(4)
 pnm(1)   = ivps(j-3)
 pnm(2)   = ivps(j-2)
 IF (noprt) GO TO 330
 WRITE (nout,320) dplx,pnm
 320 FORMAT (73X,1H(,d15.8,1H,,d15.8,1H),4H  = ,2A4)
 WRITE (nout,230) row,col,TYPE(prec),fnm
 
!     CLOSE INPUT UNIT AND RETURN
 
 330 CALL CLOSE (in1,1)
 RETURN
 
!     ERROR MESSAGES, SET THEM ALL TO NON-FATAL
 
!     NOT ENOUGH CORE FOR GINO BUFFER
 
 400 j = 8
 GO TO 430
 
!     INPUT FILE ERROR
 
 410 j = 1
 GO TO 430
 
!     INVALID ROW OR COLUMN NUMBER
 
 420 j = 7
 
 430 CALL mesage (j,in1,NAME)
 RETURN
 
END SUBROUTINE scalar
