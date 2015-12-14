SUBROUTINE trnsp (core)
     
!     OUT-OF-CORE MATRIX TRANSPOSE USING 1 TO 8 SCRATCH FILES - NASTRAN
!     ORIGINAL ROUTINE.
 
!     (SEE TRANSP FOR IN-CORE MATRIX TRANSPOSE FOR UPPER TRIAG. MATRIX,
!      AND TRNSPS FOR OUT-OF-CORE MATRIX TRANSPOSE WITH 1 SCRATCH FILE,
!      A NASTRAN NEW ROUTINE)
 
!     REVERT TO NASTRAN ORIGINAL TRNSP IF DIAG 41 IS ON, OR 94TH WORD OF
!     /SYSTEM/ IS 1000. OTHERWISE SEND THE TRANSPOSE JOB TO THE NEW
!     TRNSPS ROUTINE, EXECPT LOWER AND UPPER TRIANGULAR MATRICES
 
 
 REAL, INTENT(OUT)                        :: core(1)
 INTEGER :: scrth,otpe,sysbuf,trb1
 DIMENSION  trb1(7,8),a(2),iparm(2),NAME(2),zero(4)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /trnspx/ namea ,ncola ,nrowa ,iforma,itypa ,ia(2),  &
     nameat,ncolat,nrowat,iforat,itypat,iat(2), lcare,nscrh,scrth(8)
 COMMON /machin/ mach
 COMMON /system/ sysbuf,otpe,skip(91),ksys94
 COMMON /packx / iotyp,iotypa,ii,jj,incr
 COMMON /unpakx/ iotyp1,is1,nrow1,incr1
 DATA    iparm / 4HTRAN,4HPOSE/, zero / 4*0.0   /
 DATA    NAME  / 4HTRNS,4HP   /
 
 IF (nscrh /= 8) CALL conmsg (iparm,2,0)
 iat(1) = 0
 iat(2) = 0
 incr1  = 1
 ii     = 1
 IF (itypat == 0) itypat = itypa
 iotyp  = MIN0(itypat,itypa)
 iotypa = iotyp
 iotyp1 = iotyp
 IF (iforma == 4 .OR. iforma == 5) GO TO 50
!               LOWER            UPPER TRIANG. MATRICES
 
 j = MOD(ksys94,10000)/1000
 IF (j == 1) GO TO 50
 CALL sswtch (41,j)
 IF (j == 1) GO TO 50
 
!     NASTRAN MAINTENANCE WORK IS DONE ON VAX
 
 IF (mach /= 5 .OR. iforma < 3 .OR. iforma == 6) GO TO 40
!               VAX      NOT SQUARE, RECTANG., AND SYMM.
 CALL fname (namea,a)
 WRITE  (4,30) a,iforma
 30 FORMAT (40X,'MATRIX ',2A4,', FORM =',i2,' ===>TRNSPS')
 40 ncolat = 0
 nscrth = 1
 CALL trnsps (core,core)
 GO TO 500
 
 50 iparm1 = namea
 nscrth = nscrh
 im1    = 1
 ncalat = ncolat
 ncolat = 0
 ij1    = 0
 last   = 1
 ntype  = iotypa
 IF (ntype == 3) ntype = 2
 lcore  = lcare
 ibuf1  = lcore - sysbuf
 ibuf   = ibuf1 - sysbuf
 lcore  = ibuf  - 1
 IF (lcore > 0) THEN
   GO TO    70
 ELSE
   GO TO   440
 END IF
 
!     COMMENT FROM G.CHAN/UNISYS    1/91
!     ABOUT THE SQUARE OR RECTANGULAR MATRIX TRANSPOSE BY THE VAX -
!     DATA, 1.0**-10 OR SMALLER, ON THE TRANSPOSED MATRIX MAY DIFFER
!     FROM THE ORIGINAL VALUES. CAN NOT EXPLAIN WHY.
!     THE NORMAL DATA, 1.0**+5 OR LARGER, ARE ALL OK.
!     (NO CHECK ON THE OTHER MACHINES)
 
 70 nrowo  = MIN0(nrowat,ncola)
 nbrut  = lcore/(nrowo*ntype)
 IF (nbrut == 0) GO TO 440
 nrem   = nbrut
 IF (nbrut > ncalat) GO TO 380
 k = AMAX1(FLOAT(nrowat)*SQRT(FLOAT(ntype)/FLOAT(lcore)),1.0)
 80 nrow2 = nbrut*k
 nrow  = MIN0(nscrth*nrow2,ncalat)
 km    = (ncalat+nrow-1)/nrow
 icol  = nbrut*ntype
 IF (lcore < nrow*ntype+(nscrth-1)*sysbuf) GO TO 440
 
!     THERE ARE NROW2 ROWS IN EACH SUBMATRIX
!     WE GENERATE NROW ROWS PER PASS OF FULL MATRIX
!     THERE WILL BE KM SUCH PASSES
 
 ioloop = 1
 90 IF (ij1 == 0) THEN
   GO TO   100
 ELSE
   GO TO   210
 END IF
 100 IF (ioloop == km) GO TO 390
 nrow1 = nrow*ioloop
 110 is1   = nrow1 - nrow + 1
 IF (ioloop /= 1) GO TO 120
 iparm1= namea
 CALL OPEN (*410,namea,core(ibuf1),0)
 120 CALL fwdrec (*420,namea)
 nl = nrow*ntype
 
!     OPEN SCRATCHES
 
 j = ibuf
 DO  i = 1,nscrth
   iparm1 = scrth(i)
   CALL OPEN (*410,scrth(i),core(j),1)
   j = j - sysbuf
   DO  iii = 1,7
     trb1(iii,i) = 0
   END DO
 END DO
 DO  iloop = 1,nrowo
   CALL unpack (*180,namea,core)
   150 ik = 1
   jj = nrow2
   incr = 1
   DO  i = 1,nscrth
     CALL pack (core(ik),scrth(i),trb1(1,i))
     ik = ik + nrow2*ntype
   END DO
   
!     END LOOP ON BUILDING 1 COL OF SUB MATRICES
   
   CYCLE
   
   180 DO  i = 1,nl
     core(i) = 0.0
   END DO
   GO TO 150
 END DO
 CALL REWIND (namea)
 
!     END LOOP ON BUILDING NSCRATH SUB MATRICES
 
 DO  i = 1,nscrth
   CALL CLOSE (scrth(i),1)
 END DO
 loop350:  210 DO  j = 1,nscrth
   IF (ij1 == 0) THEN
     GO TO   220
   ELSE
     GO TO   230
   END IF
   220 IF (ioloop /= km .OR. j /= nscrth) GO TO 230
   last = 0
   230 DO  m = 1,k
     iparm1 = scrth(j)
     CALL OPEN (*410,scrth(j),core(ibuf),0)
     IF (last == 1 .OR. ncalat-ncolat >= nrem) GO TO 240
     nbrut = ncalat - ncolat
     icol  = nbrut*ntype
     is1   = (m-1)*nrem + 1
     nrow1 = is1 + nbrut
     GO TO 270
     240 IF (ij1 == 0) THEN
       GO TO   260
     END IF
     250 CALL fwdrec (*420,scrth(j))
     260 is1   = (m-1)*nbrut + 1
     nrow1 = nbrut*m
     270 l = 1
     DO  i = 1,nrowo
       CALL unpack (*280,scrth(j),core(l))
       GO TO 300
       280 DO  nl = 1,icol
         m2 = nl + l - 1
         core(m2) = 0.0
       END DO
       300 l = l + icol
     END DO
     CALL CLOSE (scrth(j),1)
     iparm1 = nameat
     CALL OPEN (*410,nameat,core(ibuf),im1)
     IF (im1 == 3) GO TO 320
     CALL fname (nameat,a(1))
     CALL WRITE (nameat,a(1),2,1)
     im1  = 3
     320 incr = nbrut
     jj = nrowo
     DO  l = 1,nbrut
       m2 = ntype*(l-1) + 1
       CALL pack (core(m2),nameat,nameat)
     END DO
     CALL CLOSE (nameat,2)
     
!     END LOOP ON SUBMATRIX
     
     IF (ncolat >= ncalat) CYCLE loop350
   END DO
   
!     END LOOP ON EACH SCRATCH
   
 END DO loop350
 
!     END LOOP ON EACH PASS THROUGH LARGE MATRIX
 
 ioloop = ioloop + 1
 IF (ioloop <= km) GO TO 90
 iparm1 = nameat
 CALL OPEN  (*410,nameat,core(ibuf),3)
 CALL CLOSE (nameat,1)
 CALL CLOSE (namea, 1)
 GO TO 500
 
!     ONE PASS ONLY
 
 380 nscrth  = 1
 scrth(1)= namea
 nbrut   = ncalat
 k   = 1
 ij1 = 1
 iotyp = itypa
 GO TO 80
 390 iover = ncalat - (km-1)*nrow
 nbrut = MIN0(nbrut,iover)
 icol  = nbrut*ntype
 nrow  = iover
 nrow2 = MIN0(nbrut*k,nrow)
 k     = (nrow2+nbrut-1)/nbrut
 nscrth= MIN0((iover+k*nbrut-1)/(k*nbrut),nscrth)
 IF (nscrth == 0) nscrth =  1
 nrow1 = ncalat
 GO TO 110
 
!     ERROR MESSAGES
 
 410 n1 = -1
 GO TO 450
 420 n1 = -2
 GO TO 450
 440 n1 = -8
 450 CALL mesage (n1,iparm1,NAME)
 
!     ONE FINAL CHECK BEFORE RETURN
 
 500 IF (iforma == 3 .OR. iforma == 7) GO TO 520
 IF (ncolat == nrowa .AND. nrowat == ncola) GO TO 520
 CALL fname (namea,a)
 WRITE (otpe,510) swm,a,iforma,ncola,nrowa,iforat,ncolat,nrowat
 510 FORMAT (a27,' FORM TRNSP. TRANSPOSED MATRIX APPEARS IN ERROR',  &
     /5X,'ORIGINAL ',2A4, ' - FORM =',i3,',  (',i6,' X',i6,')',  &
     /5X,'TRNASPOSED MATRIX - FORM =',i3,',  (',i6,' X',i6,')')
 520 IF (nscrh /= 8) CALL conmsg (iparm,2,0)
 RETURN
END SUBROUTINE trnsp
