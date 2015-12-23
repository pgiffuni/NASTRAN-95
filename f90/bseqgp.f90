SUBROUTINE bseqgp (norig,ild,jump)
     
 
 INTEGER, INTENT(IN)                      :: norig(2)
 INTEGER, INTENT(IN)                      :: ild(1)
 INTEGER, INTENT(IN OUT)                  :: jump
 EXTERNAL        orf
 INTEGER :: geom1,    geom2,    seqgp(3), eof(3),   sub(2),  &
     two,      orf,      obw,      op,       rd,  &
     rdrew,    wrt,      wrtrew,   rew,      grid(8), z
 DIMENSION  isys(100)
 COMMON /banda / ibuf1,    dum2a(2), nopch,    dum1a,    method,  &
     icrit,    ngpts,    nspts
 COMMON /bandb / nbit,     kore,     dum1b,    ngrd
 COMMON /bandd / obw,      nbw,      op,       np,       ncm,  &
     nzero,    nel,      neq,      neqr
 COMMON /bands / nn,       mm,       dum2(2),  ngrid,    dum3(3),  &
     mindeg,   nedge
 COMMON /bandw / maxw0,    rms0,     maxw1,    rms1,     i77, brms0,    brms1
 COMMON /two   / two(1)
 COMMON /system/ ibuf,     nout
 COMMON /names / rd,       rdrew,    wrt,      wrtrew,   rew, norew
 COMMON /geomx / geom1,    geom2
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ibuf,isys(1)),     (nlpp,isys(9)),  &
     (lpch,isys(91)),    (iecho,isys(19))
 DATA            sub           ,  eof    ,  seqgp          /  &
     4HSSEQ, 4HGP  ,  3*2147483647,  5301, 53, 4    /
 
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     NORIG(I) = ORIGINAL GRID POINT CORRESPONDING TO BANDIT INTERNAL
!                LABLE I
!     ILD(I)   = NEW RESEQUENCED LABEL CORRESPONDING TO BANDIT INTERNAL
!                LABLE I
!     NN       = NUMBER OF GRID POINTS
!     NGRD     .LT.0, INSUFF. WORKING CORE, OR SCRATCH ARRAY FOR BANDIT
 
 j77 = 0
 IF (nn <= 0 .OR. ngrd < 0) GO TO 145
 
!     PRINT BANDIT SUMMARY.
 
 IF (nlpp <= 48 .AND. method == 0) CALL page1
 WRITE  (nout,10)
 10 FORMAT (//53X,'*** BANDIT SUMMARY ***',/, /72X,'BEFORE',5X,'AFTER')
 
 WRITE  (nout,20) obw,nbw,op,np,maxw0,maxw1
 20 FORMAT (40X,'BANDWIDTH (B)',15X,2I10, /40X,'PROFILE (P)', 17X,2I10,  &
           /40X,'MAXIMUM WAVEFRONT (C-MAX)',3X,2I10)
 
 ann = FLOAT(nn)
 av1 = FLOAT(op)/ann
 av2 = FLOAT(np)/ann
 WRITE  (nout,30) av1,av2,rms0,rms1,brms0,brms1,ngpts
 30 FORMAT (40X,'AVERAGE WAVEFRONT (C-AVG)',3X,2F10.3,  &
           /40X,'RMS WAVEFRONT (C-RMS)',7X,2F10.3,  &
           /40X,'RMS BANDWIDTH (B-RMS)',7X,2F10.3,  &
           /40X,'NUMBER OF GRID POINTS (N)',15X,i8)
 
 IF (nspts > 0) WRITE (nout,35) nspts
 35 FORMAT (40X,'NUMBER OF SCALAR POINTS',17X,i8)
 
 WRITE  (nout,40) nel,neqr,neq
 40 FORMAT (40X,'NUMBER OF ELEMENTS (NON-RIGID)',10X,I8,   &
           /40X,'NUMBER OF RIGID ELEMENTS PROCESSED*',5X,I8,&
           /40X,'NUMBER OF MPC  EQUATIONS PROCESSED*',5X,I8)
 
 WRITE  (nout,50) ncm,mm,mindeg
 50 FORMAT (40X,'NUMBER OF COMPONENTS',20X,I8, &
           /40X,'MAXIMUM NODAL DEGREE',20X,I8, &
           /40X,'MINIMUM NODAL DEGREE',20X,I8)
 
 nonz = 2*nedge + nn
 an   = nn*nn
 den  = FLOAT(nonz)*100./an
 WRITE  (nout,60) nedge,den,nzero,kore
 60 FORMAT (40X,'NUMBER OF UNIQUE EDGES',18X,I8,  &
           /40X,'MATRIX DENSITY, PERCENT', 16X,F9.3, &
           /40X,'NUMBER OF POINTS OF ZERO DEGREE',9X,I8, &
           /40X,'BANDIT OPEN CORE',24X,I8)
 
 IF (icrit == 1) WRITE (nout,61)
 IF (icrit == 2) WRITE (nout,62)
 IF (icrit == 3) WRITE (nout,63)
 IF (icrit == 4) WRITE (nout,64)
 61 FORMAT (40X,'CRITERION*',25X,'RMS WAVEFRONT')
 62 FORMAT (40X,'CRITERION*',29X,'BANDWIDTH')
 63 FORMAT (40X,'CRITERION*',31X,'PROFILE')
 64 FORMAT (40X,'CRITERION*',25X,'MAX WAVEFRONT')
 
 IF (method == -1) WRITE (nout,66)
 IF (method == +1) WRITE (nout,67)
 IF (method ==  0) WRITE (nout,68)
 66 FORMAT (40X,'METHOD USED*',34X,'CM')
 67 FORMAT (40X,'METHOD USED*',33X,'GPS')
 68 FORMAT (40X,'METHOD USED*',26X,'CM AND GPS')
 
 IF (jump == 0) GO TO 90
 WRITE  (nout,75)
 75 FORMAT (/31X,'(* THESE DEFAULT OPTIONS CAN BE OVERRIDDEN BY THE',  &
     ' NASTRAN CARD)')
 WRITE  (nout,80)
 80 FORMAT (//31X,'BANDIT FINDS GRID POINT RE-SEQUENCING NOT ', 'NECESSARY')
 GO TO 142
 
!     GENERATE SEQGP ARRAY AND OUTPUT SEQGP CARDS
 
 90 j = 0
 DO  i = 1,nn
   z(j+1) = norig(i)
   z(j+2) = ild(i)
   j = j + 2
 END DO
 CALL sort (0,0,2,1,z(1),j)
 
!     CHECK AGAINST ORIGINAL GRID POINT DATA, AND SEE ANY UNUSED GRIDS
!     (SUCH AS THE THIRD GRID ON CBAR CARD). IF THEY EXIST, BRING THEM
!     IN, AND RE-SORT TABLE.  (GEOM1 IS READY HERE, SEE BGRID)
 
 CALL OPEN (*160,geom1,z(ibuf1),rd)
 nnx = nn
 IF (nn == ngrid) GO TO 106
 CALL READ (*104,*104,geom1,grid,3,0,k)
 102 CALL READ (*104,*104,geom1,grid,8,0,k)
 CALL bisloc (*103,grid(1),z,2,nnx,k)
 GO TO 102
 103 nn = nn + 1
 z(j+1) = grid(1)
 z(j+2) = nn
 j = j + 2
 GO TO 102
 
!     DO THE SAME CHECK IF SCALAR POINTS ARE PRESENT
 
 104 IF (nspts == 0) GO TO 1045
 nonz = j + 2*nspts + 2
 CALL preloc (*1045,z(nonz),geom2)
 grid(1) = 5551
 grid(2) = 49
 CALL locate (*1044,z(nonz),grid,k)
 1042 CALL READ (*1044,*1044,geom2,i,1,0,k)
 CALL bisloc (*1043,i,z,2,nnx,k)
 GO TO 1042
 1043 nn = nn + 1
 z(j+1) = i
 z(j+2) = nn
 j = j + 2
 GO TO 1042
 1044 CALL CLOSE (geom2,rew)
 1045 i = nn - nnx
 IF (i > 0) WRITE (nout,105) i
 105 FORMAT (40X,'NO. OF NON-ACTIVE GRID POINTS',11X,I8)
 106 i = (j+7)/8
 WRITE  (nout,107) i
 107 FORMAT (40X,'NO. OF SEQGP CARDS GENERATED',12X,I8)
 WRITE  (nout,75)
 IF (nopch == +9) GO TO 147
 IF (nnx   /= nn) CALL sort (0,0,2,1,z(1),j)
 IF (iecho == -1) GO TO 125
 CALL page1
 WRITE  (nout,110)
 110 FORMAT (//35X,'S Y S T E M  G E N E R A T E D  S E Q G P  C A R D S',/)
 WRITE  (nout,120) (z(i),i=1,j)
 120 FORMAT (25X,'SEQGP   ',8I8)
 121 FORMAT (    'SEQGP   ',8I8)
 125 IF (nopch <= 0) GO TO 130
 WRITE (lpch,121) (z(i),i=1,j)
 127 j77 = -2
 GO TO 141
 
!     BEEF UP INTERNAL GRID NOS. BY 1000 AS REQUIRED BY NASTRAN
 
 130 DO  i = 2,j,2
   z(i) = z(i)*1000
 END DO
 
!     REWIND AND SKIP FORWARDS TO THE END OF GEOM1 FILE.
!     OVERWRITE THE OLD SEQGP RECORD IF NECESSARY.
!     (WARNING - IF SEQGP IS NOT THE VERY LAST ITEM IN GEOM1 FILE, THE
!      FOLLOWING LOGIC OF INSERTING SEQGP CARDS NEEDS MODIFICATION -
!      BECAUSE GEOM1 IS IN ALPHA-NUMERIC SORTED ORDER).
 
 CALL REWIND (geom1)
 CALL skpfil (geom1,+1)
 CALL skpfil (geom1,-1)
 CALL bckrec (geom1)
 CALL READ (*150,*150,geom1,norig(1),3,1,i)
 IF (norig(1) == seqgp(1) .AND. norig(2) == seqgp(2)) CALL bckrec (geom1)
 CALL CLOSE (geom1,norew)
 
!     ADD SEQGP CARDS TO THE END OF GEOM1 FILE
!     SET GEOM1 TRAILER, AND CLEAR /SYSTEM/ 76TH WORD
 
 CALL OPEN  (*160,geom1,z(ibuf1),wrt)
 CALL WRITE (geom1,seqgp(1),3,0)
 CALL WRITE (geom1,z(1),j,1)
 CALL WRITE (geom1,eof(1),3,1)
 
 z(1) = geom1
 CALL rdtrl (z(1))
 i = (seqgp(2)+31)/16
 j = seqgp(2)-i*16 + 48
 z(i) = orf(z(i),two(j))
 CALL wrttrl (z(1))
 141 CALL CLOSE (geom1,rew)
 142 DO  i = 1,kore
   z(i) = 0
 END DO
 145 isys(i77) = j77
 IF (ngrd < 0) RETURN
 CALL page2 (-2)
 WRITE  (nout,146)
 146 FORMAT ('0',9X,'**NO ERRORS FOUND - EXECUTE NASTRAN PROGRAM**')
 RETURN
 
!     SPECIAL PUNCH OPTION (BANDTPCH=+9)
!     TO PUNCH OUT EXTERNAL GRIDS IN RE-SEQUENCED INTERNAL ORDER
 
 147 CALL sort (0,0,2,2,z(1),j)
 WRITE  (nout,148) (z(i),i=1,j,2)
 148 FORMAT ('1',35X,'LIST OF EXTERNAL GRID POINTS IN INTERNAL RE-SEQUENCED ORDER', &
           /4X,31(4H----),/,(/5X,15I8))
 WRITE  (lpch,149) (z(i),i=1,j,2)
 149 FORMAT (10I7)
 GO TO 127
 
!     FILE ERROR
 
 150 k = -2
 GO TO 170
 160 k = -1
 170 CALL mesage (k,geom1,sub)
 RETURN
END SUBROUTINE bseqgp
