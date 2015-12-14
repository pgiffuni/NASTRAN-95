SUBROUTINE frrd1f (mhh,bhh,khh,frl,frqset,nload,nfreq,ph,uhv)
     
!     ROUTINE  SOLVES DIRECTLY FOR UNCOUPLED MODAL FORMULATION
 
 
 INTEGER, INTENT(IN)                      :: mhh
 INTEGER, INTENT(IN)                      :: bhh
 INTEGER, INTENT(IN)                      :: khh
 INTEGER, INTENT(IN OUT)                  :: frl
 INTEGER, INTENT(IN OUT)                  :: frqset
 INTEGER, INTENT(IN)                      :: nload
 INTEGER, INTENT(IN)                      :: nfreq
 INTEGER, INTENT(IN OUT)                  :: ph
 INTEGER, INTENT(IN OUT)                  :: uhv
 INTEGER :: sysbuf,FILE,mcb(7),NAME(2)
 COMMON /system/ sysbuf
 COMMON /zblpkx/ b(4),jj
 COMMON /zntpkx/ a(4),ii,ieol,IEOR
 COMMON /zzzzzz/ core(1)
 DATA    NAME  / 4HFRRD,4H1F  /
 
 
 ibuf1 = korsz(core) - sysbuf + 1
 
!     PICK UP FREQUENCY LIST
 
 CALL gopen  (frl,core(ibuf1),0)
 CALL skprec (frl,frqset-1)
 IF (ibuf1-1 < nfreq) GO TO 170
 CALL fread (frl,core,nfreq,1)
 CALL CLOSE (frl,1)
 
!     BRING IN  MODAL MATRICES
 
 imhh   = nfreq
 mcb(1) = mhh
 CALL rdtrl (mcb)
 lhset  = mcb(2)
 IF (ibuf1-1 < nfreq+3*lhset) GO TO 170
 ibhh   = imhh + lhset
 ikhh   = ibhh + lhset
 
!     BRING IN MHH
 
 matnam = mhh
 ASSIGN 30 TO iret
 ipnt   = imhh
 GO TO 110
 
!     BRING  IN  BHH
 
 30 matnam = bhh
 ASSIGN 40 TO iret
 ipnt   = ibhh
 GO TO 110
 
!     BRING IN KHH
 
 40 matnam =  khh
 ASSIGN 50 TO iret
 ipnt   = ikhh
 GO TO 110
 
!     READY LOADS
 
 50 CALL gopen (ph,core(ibuf1),0)
 
!     READY SOLUTIONS
 
 ibuf2 = ibuf1 - sysbuf
 CALL gopen  (uhv,core(ibuf2),1)
 CALL makmcb (mcb,uhv,lhset,2,3)
 
!     COMPUTE  SOLUTIONS
 
 DO  i = 1,nload
   DO   j = 1,nfreq
     
!     PICK  UP  FREQ
     
     w  = core(j)
     w2 = -w*w
     CALL bldpk (3,3,uhv,0,0)
     CALL intpk (*80,ph,0,3,0)
     60 IF (ieol == 0) THEN
       GO TO    70
     ELSE
       GO TO    80
     END IF
     70 CALL zntpki
     
!     COMPUTE  REAL AND COMPLEX PARTS OF DENOMINATOR
     
     ik   = ikhh + ii
     ib   = ibhh + ii
     im   = imhh + ii
     rdem = w2*core(im) + core(ik)
     cdem = core(ib)*w
!IBMD DEM  = RDEM*RDEM + CDEM*CDEM
!IBMR IF (DEM .NE. 0.0) GO TO 71
     IF (rdem /= 0.0 .OR. cdem /= 0.0) GO TO 71
     CALL mesage (5,j,NAME)
     b(1) = 0.0
     b(2) = 0.0
     GO TO 72
     71 CONTINUE
     
!     COMPUTE REAL AND COMPLEX PHI-S
     
!IBMD B(1) = (A(1)*RDEM + A(2)*CDEM)/DEM
!IBMD B(2) = (A(2)*RDEM - A(1)*CDEM)/DEM
!IBMNB
     IF (rdem == 0.0) GO TO 715
     ratio = cdem/rdem
     factr = 1.0 / (rdem + ratio*cdem)
     b(1) = (a(1) + a(2)*ratio) * factr
     b(2) = (a(2) - a(1)*ratio) * factr
     GO TO 72
     715 ratio = rdem/cdem
     factr = 1.0 / (ratio*rdem + cdem)
     b(1) = (a(1)*ratio + a(2)) * factr
     b(2) = (a(2)*ratio - a(1)) * factr
!IBMNE
     72 jj   = ii
     CALL zblpki
     GO TO 60
     
!     END  COLUMN
     
     80 CALL bldpkn (uhv,0,mcb)
   END DO
 END DO
 CALL CLOSE  (uhv,1)
 CALL CLOSE  (ph,1)
 CALL wrttrl (mcb)
 RETURN
 
!     INTERNAL SUBROUTINE TO BRING IN  H MATRICES
 
 110 FILE = matnam
 CALL OPEN (*132,matnam,core(ibuf1),0)
 CALL skprec (matnam,1)
 DO  i = 1,lhset
   ipnt = ipnt + 1
   CALL intpk (*120,matnam,0,1,0)
   CALL zntpki
   IF (ii /= i .OR. ieol /= 1) GO TO 180
   core(ipnt) = a(1)
   CYCLE
   
!     NULL COLUMN
   
   120 core(ipnt) = 0.0
 END DO
 CALL CLOSE (matnam,1)
 131 GO TO iret, (30,40,50)
 
!     ZERO CORE FOR PURGED MATRIX
 
 132 DO  i = 1,lhset
   ipnt = ipnt + 1
   core(ipnt) = 0.0
 END DO
 GO TO 131
 
!     ERROR MESAGES
 
 150 CALL mesage (ip1,FILE,NAME)
 170 ip1 = -8
 GO TO  150
 180 ip1 = -7
 GO TO  150
END SUBROUTINE frrd1f
