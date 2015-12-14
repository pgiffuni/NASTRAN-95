SUBROUTINE sgino
     
!     REVISED  9/90 BY G.CHAN/UNISYS. TO REACTIVATE PLT1 FILE
 
!     THE HIGH POINTS OF PLT1 FILE ARE
!       NEW 130 COLUMN FORMAT RECORD
!       MACHINE PORTABLE FILE
!       NO DATA RECONSTRUCTION REQUIRED WHEN PLT1 IS USED BY AN EXTERNAL
!          TRANSLATOR PROGRAM
!                                          PLT2 FILE       PLT1 FILE
!     ---------------------------------  -------------  --------------
!     FILE TYPE, SEQUENTIAL FORMATTED    NO CARRIGE CTRL CARRIAGE CTRL
!     RECORD  TYPE                        ASSCII/BINARY*    ASCII
!     RECORD  LENGTH                       3000 BYTES     130 COLUMNS
!     FORTRAN FORMAT                      (10(180A4))*   (5(2I3,4I5))
!     PLOT COMMANDS PER PHYSICAL RECORD        100            5
!     DATA TYPE PER COMMAND (TOTAL)         30 BYTES       26 DECIMALS
!          COMMAND, P   (SEE USER'S MANUAL   1 BYTE         3 DIGITS
!          CONTROL, C    PAGE 4.4-2)         1 BYTE         3 DITITS
!          FIRST  VALUE, R                   5 BYTES        5 DIGITS
!          SECOND VALUE, S                   5 BYTES        5 DIGITS
!          THIRD  VALUE, T                   5 BYTES        5 DIGITS
!          FOURTH VALUE, U                   5 BYTES        5 DIGITS
!          FILLER (ALL ZEROS)                8 BYTES          NONE
!     DATA BYTE PACKING                        YES            NO
!     FILE - EDITED, PRINTED, SCREEN VIEWING   NO             YES
!     PORTABLE FILE AMONG MACHINES             NO             YES
!     FORTRAN UNIT NUMBER                      13             12
!     DISC STORAGE REQUIREMENT                  -           25% LESS
!     IF MAGNETIC TAPE - TRACK AND PARITY     9,ODD          9,ODD
!     (* 1. ASCII RECORD, BUT DATA STORED IN BINARY BYTES.
!           (IN EARLY NASTRAN PLOT TAPE DESIGN, A BYTE HAD 6
!           BITS. BUT IT IS NO LONGER TRUE. NOW, A BYTE CAN
!           BE 6, 8 OR 9 BITS, DEPENDING ON THE MACHINE)
!        2. SINCE THE RECORD LENGTH IS 3000 BYTES, A FORMAT
!           OF (750A4) IS SUFFICIENT)
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: OPEN,nopack
 INTEGER :: buf(1),lbuf(1),a(1),NAME(2),FORMAT(3),formtx(3)
 CHARACTER (LEN=7) :: fortn,NONE
 COMMON /system/ idum1,nout,skpsys(36),nbpc,nbpw,ncpw
!WKBNB
 CHARACTER (LEN=80) :: dsnames
 COMMON / dsname / dsnames(80)
!WKBNE
 DATA    OPEN  /.false. /, NAME  / 4H sgi, 2HNO /
 DATA    plt1,plt2,pltx /  4HPLT1, 4HPLT2,    0 /
 DATA    FORMAT/ 4H(10( ,  4H180A, 4H4))        /,  &
     formtx/ 4H(5(2 ,  4HI3,4, 4HI5))       /
 DATA    fortn , NONE   / 'FORTRAN', 'NONE   '  /
 DATA    shift , nbits  /  0, 0                 /
 
 GO TO 250
 
 
 ENTRY sopen (*,pltape,buf,ibfsz)
!     ================================
 
!     PLT2 FILE -
!     IBFSZ (FIRST WORD OF /XXPARM/), IS THE PLOT FILE BUFFER SIZE. IT
!     IS SET EQUAL TO PDATA(12,1)/NCPW IN PLTSET. PDATA(12,1) IS
!     INITIALIZED  IN PLOTBD VIA DATA(12,1) WHICH IS EQUIVALENT TO
!     PBFSIZ(1,1). COMPLICATED ISN'T IT?
 
!     (PBFSIZ(1,1)=3000, NCPW=4, IBFSZ AND BFSZ ARE THEREFORE =750 EACH
!     EACH PHYSICAL RECORD HOLDS 100 (=3000/30) PLOT COMMANDS)
 
!     NOTE - BOTH PLT2 AND PLT1 ARE SEQUENTIAL FILES, NOT DIRECT ACCESS
!     FILES. THE RECORD LENGTH, IF USED, IS BASED ON NO. CHARACTERS PER
!     WORD
 
 ptape = 10
 IF (pltape /= plt1 .AND. pltape /= plt2) RETURN 1
 pltx   = pltape
 nopack = pltx == plt1
 IF (nopack) GO TO 10
 
!     PLT2 -
 
 bfsz   = ibfsz
 irecsz = ncpw*bfsz
 noff   = locfx(buf(1)) - locfx(lbuf(1))
 GO TO 20
 
!     PLT1 -
 
!WKBR 10   PTAPE  = 12
 10   CONTINUE
 nopack = .true.
 bfsz   = 30
 irecsz = (bfsz/6)*(2*3 + 4*5)
 NONE   = fortn
 FORMAT(1) = formtx(1)
 FORMAT(2) = formtx(2)
 FORMAT(3) = formtx(3)
 
!     NOFF CAN BE SET TO ZERO IF LBUF IS LOCALLY DIMENSIONED TO 30 WORDS
!     AND OPEN CORE IS NOT USED
 
 noff   = locfx(buf(1)) - locfx(lbuf(1))
 
!     OPEN STATEMENT ADDED TO SET OUTPUT RECORDSIZE GREATER THAN DEFAULT
!     (COMMENTS FORM G.C./UNISYS 1989 -
!     RECORDSIZE IS NOT ALLOWED FOR SEQUENTIAL FILE IN SOME COMPILERS,
!     (e.g. DEC/ULTRIX(RISC), AND BLOCKSIZE AND ORGANINZATION ARE NOT
!     DEFINED. RECORDTYPE='FIXED' IS ALSO NOT ALLOWED FOR SEQUENTIAL
!     FORMATTED FILE.
!     FOR UNICOS, RECL IS NOT ALLOWED IF ASSCESS=SEQUENTIAL)
 
!     FOR MACHINES THAT DO NOT HAVE 'APPEN' FEATURE
 
 20   IF (OPEN) GO TO 80
!     MA = 'A'
!     IF (NONE .EQ. 'NONE') MA = 'M'
!     IF (MACH .EQ IBM) CALL FILEDEF (PTAPE,RECFM,FB(MA))
 OPEN (UNIT   = ptape,
!WKBI  &
 FILE = dsnames(10), STATUS = 'OLD',  &
     FORM   = 'FORMATTED', ACCESS = 'SEQUENTIAL',  &
     IOSTAT = j
!HP  5      ,CARRIAGECONTROL = NONE
!HP  6      ,RECL  = IRECSZ
!            RECL IS NEEDED BY VAX, AND POSSIBLY OTHER MACHINES)  &
 )
 IF (j /= 0) GO TO 60
 30   READ   (ptape,40,END=50) j
 40   FORMAT (a1)
 GO TO 30
 50   BACKSPACE ptape
 GO TO 80
 
 60   OPEN (UNIT   = ptape,
!WKBI  &
 FILE = dsnames(10), STATUS = 'NEW',  &
     FORM   = 'FORMATTED', ACCESS = 'SEQUENTIAL',  &
     IOSTAT = j
!HP  5      ,CARRIAGECONTROL = NONE
!HP  6      ,RECL  = IRECSZ
!            RECL IS NEEDED BY VAX, AND POSSIBLY OTHER MACHINES)  &
 )
 IF (j == 0) GO TO 80
 WRITE  (nout,70) pltx,ptape
 70   FORMAT ('0*** SYSTEM FATAL ERROR. SGINO CAN NOT OPEN ',a4,  &
     ' FILE, FORTRAN UNIT',i5)
 CALL mesage (-61,0,0)
 
 80   OPEN  = .true.
 nb    = 1
 IF (nopack) GO TO 210
 ASSIGN 100 TO tra
 word  = 0
 nbits = nbpw - nbpc
 shift = nbits
 GO TO 250
 
 
 ENTRY swrite (pltape,a,n,eorx)
!     ==============================
 
!     SWRITE IS CALLED ONLY BY WPLT10
 
 IF (pltape /= pltx) GO TO 180
 eor = eorx
 nw  = 1
 90  IF (nopack) GO TO 120
 
!     ORIGINAL BYTE PACKING LOGIC
 
 100  IF (nw > n) GO TO 110
!UNIX IF (A(NW) .NE. 0) WORD =  OR(ISHFT(A(NW),SHIFT),WORD)
 IF (a(nw) /= 0) word = ior(ishft(a(nw),shift),word)
 nw  = nw + 1
 IF (shift == 0) GO TO 105
 shift = shift - nbpc
 GO TO 100
 105  lbuf(nb+noff) = word
 IF (nb == bfsz) GO TO 200
 word  = 0
 nb    = nb + 1
 shift = nbits
 GO TO 100
 
 110  IF (eor == 0) GO TO 250
 eor = 0
 IF (shift /= nbits) GO TO 115
 nb  = nb - 1
 IF (nb > 0) THEN
   GO TO   200
 ELSE
   GO TO   190
 END IF
 
 115  lbuf(nb+noff) = word
 GO TO 200
 
!     NON BYTE PACKING LOGIC
 
 120  IF (nw > n) GO TO 125
 lbuf(nb+noff) = a(nw)
 nw = nw + 1
 nb = nb + 1
 IF (nb <= bfsz) GO TO 120
 nb = nb - 1
 GO TO 200
 
 125  IF (eor == 0) GO TO 250
 eor = 0
 IF (nb >= bfsz) GO TO 135
 DO  j = nb,bfsz
   lbuf(j+noff) = 0
 END DO
 135  nb = bfsz
 GO TO 200
 
 
 ENTRY sclose (pltape)
!     =====================
 
 eof = 0
 
 150  IF (pltape /= pltx) GO TO 180
 IF (.NOT.nopack .AND. shift /= nbits) GO TO 155
 nb = nb - 1
 IF (nb > 0) THEN
   GO TO   160
 ELSE
   GO TO   170
 END IF
 155  lbuf(nb+noff) = word
 160  ASSIGN 165 TO tra
 GO TO 200
 165  ASSIGN 100 TO tra
 IF (nopack) ASSIGN 120 TO tra
 170  IF (eof == 0) GO TO 175
 ENDFILE ptape
 GO TO 190
 175  pltx = 0
 GO TO 190
 
 180  WRITE  (nout,185) pltx,pltape
 185  FORMAT ('0*** SYSTEM FATAL ERROR FROM SGINO. ',a4,' FILE OR ',a4,  &
     ' FILE GOT LOST')
 CALL errtrc (NAME)
 
 
 ENTRY seof (pltape)
!     ===================
 
 eof = 1
 GO TO 150
 
 190  nb = 1
 GO TO 250
 
 200  WRITE (ptape,FORMAT) (lbuf(noff+j),j=1,nb)
 nb = 1
 word  = 0
 shift = nbits
 GO TO tra, (100,120,165)
 
 210  ASSIGN 120 TO tra
 
 250  RETURN
END SUBROUTINE sgino
