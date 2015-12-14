SUBROUTINE sswtch (nbit,l)
     
!     PURPOSE OF THIS ROUTINE IS TO SET L = 1 IF SENSE SWITCH BIT IS
!     ON, OTHERWISE L = 0.
 
!     SENSE SWITCH DEFINITION
!      1 = DUMP CORE WHEN SUBROUTINE DUMP OR PDUMP(NO ARGUMENTS) IS
!          CALLED
!      2 = DUMP FIAT TABLE AFTER ALLOCATION
!      3 = DUMP DATA POOL DICTIONARY AFTER ALLOCATION
!      4 = DUMP OSCAR FILE AT END OF XGPI
!      5 = CONSOLE MESSAGE DESIRED (BEGIN)
!      6 = CONSOLE MESSAGE DESIRED (END)
!      7 = EIGENVALUE EXTRACTION DIAGNOSTICS
!          (DETERMINANT AND INVERSE POWER)
!      8 = TRACES NPTP ON 1108
!      9 = TURNS ON PRINTER PLOTTER FOR ANY XYPLOT REQUESTS
!     10 = USES ALTERNATE ALGORITHUM FOR NON LINEAR LOADS SEE SPR 153
!     11 = ACTIVE ROW AND COLUMN TIME PRINTS
!     12 = CONPLEX EIGENVALUE EXTRACTION DIAGNOSTICS
!          (INVERSE POWER)
!     28 = PUNCHES OUT LINK SPECIFICATION TABLE - DECK XBSBD
!     29 = PROCESS LINK SPECIFICATION UPDATE DECK
!     30 = PUNCHES OUT ALTERS TO XSEM-S FOR SWITCHES 1-15
!     31 = PRINT LINK SPECIFICATION TABLE
 
 
 INTEGER, INTENT(IN)                      :: nbit
 INTEGER, INTENT(OUT)                     :: l
 EXTERNAL        lshift,rshift,andf
 INTEGER :: switch,andf,rshift,renter
 COMMON /system/ xsys(78),switch(3)
 COMMON /xlink / lxlink,maxlnk
 COMMON /sem   / dummy(3),ns(1)
 DATA    renter/ 4HREEN /
 
 l = 0
 IF (iret ==  1) RETURN
 IF (nbit > 31) GO TO 10
 IF (andf(lshift(1,nbit-1),switch(1)) /= 0) l = 1
 RETURN
 
 10 nbit2 = nbit - 31
 IF (andf(lshift(1,nbit2-1),switch(2)) /= 0) l = 1
 RETURN
 
 
 ENTRY pressw (nbit,l)
!     =====================
 
!     PRESSW IS CALLED ONLY BY BGNSYS AND XCSA TO SETUP DIAGNOSTIC BITS
!     FOR A PARTICULAR LINK.
!     BITS  0 THRU 47 ARE USED FOR 48 DIAGNOSTICS
!     BITS 49 THRU 63 ARE RESERVED FOR 15 LINK NOS.
!     NBIT HERE (INPUT) CONTAINS BCD WORD NSXX WHERE XX IS LINK NO.
 
 iret = 0
 IF (nbit == renter) RETURN
 IF (switch(3)+switch(2) == 0) iret = 1
 i = 32 - maxlnk
 IF (rshift(switch(2),i) == 0) GO TO 40
 DO  i = 1,maxlnk
   IF (nbit == ns(i)) GO TO 30
 END DO
 i = 0
 IF (nbit == ns (16)) i = 5
 IF (nbit == ns (17)) i = 8
 IF (nbit == ns (18)) i = 13
 IF (nbit == ns (19)) i = 6
 IF (nbit == ns (20)) i = 2
 IF (nbit == ns (21)) i = 9
 IF (nbit == ns (22)) i = 11
 IF (nbit == ns (23)) i = 15
 IF (i /= 0) GO TO 30
 GO TO 40
 30 nbit2 = i + 31 - maxlnk
 IF (andf(lshift(1,nbit2),switch(2)) == 0) iret = 1
 40 IF (iret == 0) switch(1) = switch(3)
 RETURN
END SUBROUTINE sswtch
