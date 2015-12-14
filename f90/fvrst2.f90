SUBROUTINE fvrst2
     
!    1. ENTRY POINT - FVRST2
 
!    2. PURPOSE -  THIS MODULE IS USED DURING A FORCED VIBRATION
!                  RESPONSE ANALYSIS OF ROTATING CYCLIC STRUCTURES
!                  TO GENERATE TABLE DATA BLOCKS FRL AND FOL AND TO
!                  GENERATE MATRIX DATA BLOCKS REORDER1 AND REORDER2.
!                  FVRSTR2 ALSO COMPUTES PARAMETERS LMAX, NTSTEPS,
!                  FLMAX, NORO1 AND NORO2.
 
!    3. DMAP CALLING SEQUENCE -
 
!         FVRSTR2  TOL,,,,,,, / FRL,FOL,REORDER1,REORDER2,,,, /
!                  V,Y,NSEGS/ V,Y,CYCIO/ V,Y,LMAX=-1/ V,N,FKMAX/
!                  V,N,FLMAX/ V,N,NTSTEPS/ V,N,NORO1/ V,N,NORO2  $
 
!    4. INPUT DATA BLOCKS -
 
!         TOL    - TIME OUTPUT LIST.
 
!         NOTE   - (1) TOL MUST BE PRESENT.
 
!    5. OUTPUT DATA BLOCKS -
 
!         FRL      - FREQUENCY RESPONSE LIST.
!         FOL      - FREQUENCY OUTPUT LIST.
!         REORDER1 - LOAD REORDERING MATRIX FO TIME-DEPENDENT PROBLEMS.
!         REORDER2 - LOAD REORDERING MATRIX FO TIME-DEPENDENT PROBLEMS.
 
!         NOTE     - (1) FRL AND FOL CANNOT BE PURGED.
!                    (2) REORDER1 AND REORDER2 SHOULD NOT BE PURGED.
 
!    6. PARAMETERS -
 
!        (A) NSEGS   - INPUT-INTEGER-NO DEFAULT.  THE NUMBER OF
!                      IDENTICAL SEGMENTS IN THE STRUCTURAL MODEL.
!        (B) CYCIO   - INPUT-INTEGER-NO DEFAULT.  THE INTEGER VALUE
!                      OF THIS PARAMETER SPECIFIES THE FORM OF THE INPUT
!                      AND OUTPUT DATA FOR CYCLIC STRUCTURES. A VALUE
!                      OF +1 IS USED TO SPECIFY PHYSICAL SEGMENT REPRE-
!                      SENTATION AND A VALUE OF -1 FOR CYCLIC TRANSFOR-
!                      MATION REPRESENTATION.
!        (C) LMAX    - INPUT/OUTPUT-INTEGER.  THE INTEGER VALUE OF THIS
!                      PARAMETER SPECIFIES THE MAXIMUM TIME HARMONIC
!                      INDEX FOR CYCLIC STRUCTURES. THE DEFAULT VALUE
!                      IS NTSTEPS/2, WHERE NTSTEPS IS THE NUMBER OF
!                      TIME STEPS DEFINED BELOW.
!        (D) FKMAX   - INPUT-INTEGER-NO DEFAULT.  FUNCTION OF KMAX.
!        (E) FLMAX   - OUTPUT-INTEGER-NO DEFAULT.  FUNCTION OF LMAX.
!        (F) NTSTEPS - OUTPUT-INTEGER-NO DEFAULT.  THE NUMBER OF
!                      TIME STEPS FOUND IN DATA BLOCK TOL.
!        (G) NORO1   - OUTPUT-INTEGER-NO DEFAULT.  NORO1 =-1 IF DATA
!                      BLOCK REORDER1 IS NOT GENERATED.
!        (H) NORO2   - OUTPUT-INTEGER-NO DEFAULT.  NORO2 =-1 IF DATA
!                      BLOCK REORDER2 IS NOT GENERATED.
 
!    7. METHOD -
 
!         DATA BLOCK TOL IS READ AND THE LIST OF SOLUTION TIMES IS
!         STORED. SET NTSTEPS TO THE NUMBER OF SOLUTION TIMES READ.
!         IF NECESSARY COMPUTE THE DEFAULT VALUE OF LMAX AND THEN
!         COMPUTE FLMAX.
!         GENERATE TABLE DATA BLOCKS FOL AND FRL.
!         GENERATE MATRIX DATA BLOCKS REORDER1 AND REORDER2 AND
!         PARAMETERS NORO1 AND NORO2.
 
!    8. SUBROUTINES - FVRST2 CALLS SUBROUTINE FVRS2A AND OTHER
!                     STANDARD NASTRAN UTILITY ROUTINES.
 
!    9. DESIGN REQUIREMENTS -
 
!         (1) OPEN CORE IS DEFINED AT /ZZFVR2/.
!         (2) NO SCRATCH FILES ARE USED.
!         (3) FVRST2 RESIDES IN LINKNS07.
!         (4) OPEN CORE FOR ONE BUFFER+1 IS REQUIRED.
 
!   10. DIAGNOSTIC MESSAGES -
 
!         THE FOLLOWING MESSAGES MAY BE ISSUED - 3001,3002,3003,3008.
 
 
 INTEGER :: modnam(2),FILE,fnam(2),trl(7),dum(2),sysbuf,tol,  &
     frl,fol,reord1,reord2,cycio,fkmax,flmax
 DOUBLE PRECISION :: period,freq,fact,dpi,dtwopi,dradeg,ddegra,d4pisq
 COMMON /BLANK /  nsegs,cycio,lmax,fkmax,flmax,ntstps,noro1,noro2
 COMMON /zzzzzz/  z(1)
 COMMON /system/  sysbuf,nout
 COMMON /condad/  dpi,dtwopi,dradeg,ddegra,d4pisq
 DATA    modnam/  4HFVRS,4HTR2  /
 DATA       tol,  frl, fol, reord1, reord2 / 101,  201, 202, 203,    204    /
 
 
!     DETERMINE LENGTH OF OPEN CORE AND ALLOCATE BUFFERS.
 
 nz    = korsz(z)
 ibuf1 = nz - sysbuf
 nz    = ibuf1 - 1
 IF (nz <= 0) GO TO 9908
 
!     READ DATA BLOCK TOL (TIME OUTPUT LIST).
!     LIST OF OUTPUT TIME VALUES ARE STORED IN TOL HEADER.
 
 FILE = tol
 itol = 1
 CALL fname (FILE,fnam)
 CALL OPEN (*9901,FILE,z(ibuf1),0)
 CALL fread (FILE,dum,2,0)
 CALL READ (*9902,*10,FILE,z(itol),nz,1,ntimes)
 
!     INSUFFICIENT CORE TO HOLD ALL TIMES.
 
 GO TO 9908
 
 10 CALL CLOSE (FILE,1)
 
 nz   = nz - ntimes
 next = ntimes + 1
 IF (nz <= 0) GO TO 9908
 
!     DEFINE PARAMETER NTSTEPS.
 
!     IF (CYCIO .EQ. -1) NTSTEPS = (NTIMES*FKMAX)/FKMAX
!     IF (CYCIO .EQ. +1) NTSTEPS = (NTIMES*NSEGS)/NSEGS
 
 ntstps = ntimes
 
!     SET DEFAULT VALUE OF PARAMETER LMAX.
 
 IF (lmax < 0) lmax = ntstps/2
 
!     DEFINE PARAMETER FLMAX
 
 kk = (ntstps/2)*2
 IF (kk /= ntstps) GO TO 20
 
!     NTSTPS IS EVEN.
 
 IF (lmax /= ntstps/2) GO TO 20
 flmax = ntstps
 GO TO 30
 
!     NTSTPS IS ODD.
 
 20 flmax = 2*lmax + 1
 
 30 CONTINUE
 
!     GENERATE DATA BLOCKS FRL AND FOL BY CONVERTING TOL TIMES
!     TO THE FREQUENCY DOMAIN.
 
 nfreq= flmax
 ifol = next
 next = ifol + nfreq
 nz   = nz   - nfreq
 IF (nz <= 0) GO TO 9908
 
!     GENERATE FREQUENCY LIST FROM TOL TIME LIST.
 
 z(ifol) = 0.0
 IF (nfreq <= 1) GO TO 60
 
 period = DBLE(z(itol+1)) + DBLE(z(itol+ntimes-1))
 freq   = 1.0D0/period
 fact   = 1.0D0
 
 ifreq1 = ifol + 1
 ifreq2 = ifol + nfreq - 1
 
 DO  ifreq = ifreq1,ifreq2,2
   z(ifreq)   = fact*freq
   z(ifreq+1) = z(ifreq)
   fact       = fact + 1.0D0
 END DO
 
 kk = (nfreq/2)*2
 IF (kk /= nfreq) GO TO 60
 z(ifreq2) = fact*freq
 
 60 CONTINUE
 
!     OUTPUT FOL TABLE (FREQUENCY OUTPUT RESPONSE LIST).
 
 FILE = fol
 CALL fname (FILE,fnam)
 CALL OPEN (*9901,FILE,z(ibuf1),1)
 CALL WRITE (FILE,fnam,2,0)
 CALL WRITE (FILE,z(ifol),nfreq,1)
 CALL CLOSE (FILE,1)
 
 trl(1) = FILE
 trl(2) = nfreq
 trl(3) = 1
 trl(4) = 0
 trl(5) = 0
 trl(6) = 0
 trl(7) = 0
 CALL wrttrl (trl)
 
!     GENERATE DATA BLOCK FRL FROM FOL (W = F*2*PI).
!     USE SAME CORE WHERE FOL IS STORED.
 
 DO  ifreq = ifreq1,ifreq2
   z(ifreq) = z(ifreq)*dtwopi
 END DO
 
!     OUTPUT FRL TABLE (FREQUENCY RESPONSE LIST).
 
 FILE = frl
 CALL fname (FILE,fnam)
 CALL OPEN (*9901,FILE,z(ibuf1),1)
 CALL WRITE (FILE,fnam,2,0)
 CALL WRITE (FILE,1,1,1)
 CALL WRITE (FILE,z(ifol),nfreq,1)
 CALL CLOSE (FILE,1)
 
 trl(1) = FILE
 trl(2) = 1
 trl(3) = 0
 trl(4) = 0
 trl(5) = 0
 trl(6) = 0
 trl(7) = 0
 CALL wrttrl (trl)
 
!     GENERATE MATRIX DATA BLOCKS REORDER1 AND REORDER2 USED FOR
!     REORDERING COLUMNS OF A MATRIX BY POST-MULTIPLYING THE MATRIX
!     WHOSE COLUMNS ARE TO BE REORDERED.
 
 k1 = ntstps
 k3 = flmax
 IF (cycio == -1) k2 = fkmax
 IF (cycio == +1) k2 = nsegs
 
!     GENERATE MATRIX REORDER1
 
 CALL fvrs2a (reord1,k1,k2,noro1,z(ibuf1))
 
!     GENERATE MATRIX REORDER2
 
 CALL fvrs2a (reord2,k2,k3,noro2,z(ibuf1))
 
 RETURN
 
!     ERROR PROCESSING
 
!     DATA SET NOT IN FIST
 
 9901 ip1 = -1
 GO TO 9999
 
!     E-O-F ENCOUNTERED
 
 9902 ip1 = -2
 GO TO 9999
 
!     E-O-L ENCOUNTERED
 
 9908 ip1 = -8
 GO TO 9999
 9999 CALL mesage (ip1,FILE,modnam)
 CALL mesage (-37,0,modnam)
 
 RETURN
END SUBROUTINE fvrst2
