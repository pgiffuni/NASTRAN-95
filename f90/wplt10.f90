SUBROUTINE wplt10 (a,opt)
     
!     TO WRITE PLOTTER COMMANDS FOR NASTRAN GENERAL PURPOSE PLOTTER
!     REF - NASTRAN PROGRAMMER'S MANUAL P.3.4-111
 
!     REVISED  9/1990 BY G.CHAN/UNISYS
!     SEE SGINO FOR IMPLEMENTATION OF PLT1 FILE
 
!     INPUT -
!        OPT = 0 IF ARRAY A IS A PLOT COMMAND.
!        OPT = 1 IF CURRENT SERIES OF PLOT COMMANDS IS TO BE TERMINATED
 
!     OUTPUT -
!       A(1) = PLOT MODE DIGIT
!       A(2) = CONTROL DIGIT
!       A(3) = X1 = X-COORDINATE
!       A(4) = Y1 = Y-COORDINATE
!       A(5) = X2 = X-COORDINATE
!       A(6) = Y2 = Y-COORDINATE
 
!     A PLT2 FILE PLOTTER COMMAND IS OF THE FOLLOWING FORMAT
 
!         MC1111122222333334444400000000
!            WHERE M = MODE            1 BYTE
!                  C = CONTROL         1 BYTE
!                  1 = DIGIT OF X1     5 BYTES
!                  2 = ..... .. Y1     5 BYTES
!                  3 = ..... .. X2     5 BYTES
!                  4 = ..... .. Y2     5 BYTES
!                  0 = ZERO            8 BYTES
!                              ---------------
!                              TOTAL  30 BYTES
 
!     SEE SGINO FOR PLT1 FILE PLOTTER COMMAND FORMAT
 
!     /PLTDAT/
!     EDGE = SIZE OF THE BORDERS (X,Y) IN PLOTTER UNITS,    REAL - INPUT
!     PLOT = GINO FILE NAME OF THE PLOT TAPE TO BE WRITTEN,  BCD - INPUT
!     MAXCHR = PLOT TAPE BUFFER SIZE (NUMBER OF CHARACTERS), INT - INPUT
!              (AN INTEGER MULTIPLE OF THE NUMBER OF CHARACTERS
!              PER WORD ON THE COMPUTER ON WHICH THE PLOT TAPE IS
!              BEING READ)
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(IN)                      :: a(6)
 INTEGER, INTENT(IN OUT)                  :: opt
 INTEGER :: c(30),ten(5),zero(30)

 REAL :: edge
 
 COMMON /pltdat/  skpplt(8),edge(12),skpa(10),plot,maxchr
 
 EQUIVALENCE      (c1,c(1))
 
 DATA    nchr  ,  ten,pzero / 0, 10000, 1000, 100, 10, 1, +0 /,  &
         plt2  ,  nc,zero,c / 4HPLT2,   30,   30*0,    30*0  /
 
 
 IF (plot == plt2) GO TO 100
 
!     PLT1 FILE - NON BYTE PACKING LOGIC
!     A FORMAT OF (5(2I3,4I5)) IS COMPOSED IN SGINO
!     =============================================
 
 nc = 6
 IF (opt /= 0) GO TO 40
 
!     SET UP THE MODE AND CONTROL CHARACTERS IN THE COMMAND.
 
 c1   = a(1)
 c(2) = a(2)
 
 i3   = IFIX(edge(1) + .1)
 i4   = IFIX(edge(2) + .1)
 c(3) = a(3) + i3
 c(4) = a(4) + i4
 c(5) = a(5)
 c(6) = a(6)
 IF (c1 == 4 .OR. c1 == 14) GO TO 20
 c(5) = a(5) + i3
 c(6) = a(6) + i4
 20 CALL swrite (plot,c,nc,0)
 GO TO 200
 
!     TERMINATE A SET OF PLOT COMMANDS
!     SEND A RECORD OF ALL ZERO-S TO SWRITE
 
 40 CALL swrite (plot,zero,nc,0)
 CALL swrite (plot,0,0,1)
 GO TO 200
 
!     PLT2 FILE - WITH BYTE PACKING LOGIC
!     A FORMAT OF (10(180A4)) IS COMPOSED IN SGINO
!     ============================================
 
 100 IF (opt /= 0) GO TO 140
 
!     SET UP THE MODE + CONTROL CHARACTERS IN THE COMMAND.
 
 c1   = a(1)
 c(2) = a(2)
 
!     SEPARATE THE DECIMAL DIGITS OF THE X + Y COORDINATES.
 
 DO  j = 1,4
   i = 1
   IF (j == 2 .OR. j == 4) i = 2
   n = a(j+2)
   IF (j < 3 .OR. (c1 /= 4 .AND. c1 /= 14)) n = n + IFIX(edge(i)+.1)
   k = 5*(j-1)
   DO  i = 1,5
     m = n/ten(i)
     
!   . M MAY BE A -0 (UNIVAC), SET IT TO +0 FOR SURE
     
     IF (m == 0) m = pzero
     c(k+3) = m
     k = k + 1
     n = n - m*ten(i)
   END DO
 END DO
 
 CALL swrite (plot,c,nc,0)
 nchr = nchr + nc
 IF (nchr == maxchr) nchr = 0
 GO TO 200
 
!     TERMINATE A SET OF PLOT COMMANDS (FILL THE RECORD WITH ZERO-S).
 
 140 IF (nchr == 0) GO TO 160
 150 CALL swrite (plot,zero,nc,0)
 nchr = nchr + nc
 IF (nchr /= maxchr) GO TO 150
 nchr = 0
 160 CALL swrite (plot,0,0,1)
 
 200 RETURN
END SUBROUTINE wplt10
