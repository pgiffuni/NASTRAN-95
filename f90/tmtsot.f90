SUBROUTINE tmtsot
     
!     THIS SUBROUTINE PRINTS THE CONTENTS OF COMMON /NTIME/
 
 COMMON /ntime / nitems, tgino , tbldpk, tintpk, tpack ,  &
     tunpak, tgetst, tputst, ttlrsp, ttlrdp, ttlcsp, ttlcdp,  &
     tllrsp, tllrdp, tllcsp, tllcdp, tgetsb  &
     , rgino , rbldpk  , rintpk , rpack  , runpak, rgetst  , rputst
 COMMON /system/ isysbf, nout  , dummy(74)     , isy77
 
 WRITE (nout,2000) nitems
 WRITE (nout,2010) tgino
 WRITE (nout,2020) tbldpk
 WRITE (nout,2030) tintpk
 WRITE (nout,2040) tpack
 WRITE (nout,2050) tunpak
 WRITE (nout,2060) tgetst
 WRITE (nout,2070) tputst
 WRITE (nout,2080) ttlrsp
 WRITE (nout,2090) ttlrdp
 WRITE (nout,2100) ttlcsp
 WRITE (nout,2110) ttlcdp
 WRITE (nout,2120) tllrsp
 WRITE (nout,2130) tllrdp
 WRITE (nout,2140) tllcsp
 WRITE (nout,2150) tllcdp
 WRITE (nout,2160) tgetsb
 WRITE (nout,2210) rgino
 WRITE (nout,2220) rbldpk
 WRITE (nout,2230) rintpk
 WRITE (nout,2240) rpack
 WRITE (nout,2250) runpak
 WRITE (nout,2260) rgetst
 WRITE (nout,2270) rputst
 IF (isy77 /= -3) WRITE (nout,2200)
 RETURN
 2000 FORMAT (1H1,23X,  &
     ' DIAG 35 OUTPUT OF TIMING CONSTANTS IN COMMON /NTIME/'/ 24X,  &
     ' ----------------------------------------------------'//  &
     ' NUMBER OF TIMING CONSTANTS IN COMMON /NTIME/   ',  &
     '                         --- ', i11                   / )
 2010 FORMAT (' READ + WRITE + BACKWARD READ                   ',  &
     ' (AVERAGE PER WORD     ) --- ', e11.4, ' MICROSECONDS'/ )
 2020 FORMAT (' BLDPK  - PACK   SUCCESSIVE ELEMENTS OF A COLUMN',  &
     ' (AVERAGE PER WORD     ) --- ', e11.4, ' MICROSECONDS'/ )
 2030 FORMAT (' INTPK  - UNPACK SUCCESSIVE ELEMENTS OF A COLUMN',  &
     ' (AVERAGE PER WORD     ) --- ', e11.4, ' MICROSECONDS'/ )
 2040 FORMAT (' PACK   - PACK   AN ENTIRE COLUMN               ',  &
     ' (AVERAGE PER WORD     ) --- ', e11.4, ' MICROSECONDS'/ )
 2050 FORMAT (' UNPACK - UNPACK AN ENTIRE COLUMN               ',  &
     ' (AVERAGE PER WORD     ) --- ', e11.4, ' MICROSECONDS'/ )
 2060 FORMAT (' GETSTR - FORWARD READ  A STRING OF DATA        ',  &
     ' (AVERAGE PER WORD     ) --- ', e11.4, ' MICROSECONDS'/ )
 2070 FORMAT (' PUTSTR - WRITE A STRING OF DATA                ',  &
     ' (AVERAGE PER WORD     ) --- ', e11.4, ' MICROSECONDS'/ )
 2080 FORMAT (' TIGHT-LOOP MULTIPLY - REAL    SINGLE PRECISION ',  &
     ' (AVERAGE PER OPERATION) --- ', e11.4, ' MICROSECONDS'/ )
 2090 FORMAT (' TIGHT-LOOP MULTIPLY - REAL    DOUBLE PRECISION ',  &
     ' (AVERAGE PER OPERATION) --- ', e11.4, ' MICROSECONDS'/ )
 2100 FORMAT (' TIGHT-LOOP MULTIPLY - COMPLEX SINGLE PRECISION ',  &
     ' (AVERAGE PER OPERATION) --- ', e11.4, ' MICROSECONDS'/ )
 2110 FORMAT (' TIGHT-LOOP MULTIPLY - COMPLEX DOUBLE PRECISION ',  &
     ' (AVERAGE PER OPERATION) --- ', e11.4, ' MICROSECONDS'/ )
 2120 FORMAT (' LOOSE-LOOP MULTIPLY - REAL    SINGLE PRECISION ',  &
     ' (AVERAGE PER OPERATION) --- ', e11.4, ' MICROSECONDS'/ )
 2130 FORMAT (' LOOSE-LOOP MULTIPLY - REAL    DOUBLE PRECISION ',  &
     ' (AVERAGE PER OPERATION) --- ', e11.4, ' MICROSECONDS'/ )
 2140 FORMAT (' LOOSE-LOOP MULTIPLY - COMPLEX SINGLE PRECISION ',  &
     ' (AVERAGE PER OPERATION) --- ', e11.4, ' MICROSECONDS'/ )
 2150 FORMAT (' LOOSE-LOOP MULTIPLY - COMPLEX DOUBLE PRECISION ',  &
     ' (AVERAGE PER OPERATION) --- ', e11.4, ' MICROSECONDS'/ )
 2160 FORMAT (' GETSTB - BACKWARD READ  A STRING OF DATA       ',  &
     ' (AVERAGE PER WORD     ) --- ', e11.4, ' MICROSECONDS'/ )
 2200 FORMAT ('0*** NASTRAN INFORMATION MESSAGE, TO INCORPORATE THESE ',  &
     'TIMING CONSTANTS INTO NASTRAN PERMANENTLY', /5X,  &
     'RE-RUN JOB WITH ''NASTRAN BULKDATA=-3'' FOR MORE ', 'INSTRUCTIONS',/)
 2210 FORMAT (' READ + WRITE + BACKWARD READ                   ',  &
     ' (AVERAGE PER RECORD   ) --- ', e11.4, ' MICROSECONDS'/ )
 2220 FORMAT (' BLDPK  - PACK   SUCCESSIVE ELEMENTS OF A COLUMN',  &
     ' (AVERAGE PER RECORD   ) --- ', e11.4, ' MICROSECONDS'/ )
 2230 FORMAT (' INTPK    UNPACK SUCCESSIVE ELEMENTS OF A COLUMN',  &
     ' (AVERAGE PER RECORD   ) --- ', e11.4, ' MICROSECONDS'/ )
 2240 FORMAT (' PACK   - PACK   AN ENTIRE COLUMN               ',  &
     ' (AVERAGE PER RECORD   ) --- ', e11.4, ' MICROSECONDS'/ )
 2250 FORMAT (' UNPACK - UNPACK AN ENTIRE COLUMN               ',  &
     ' (AVERAGE PER RECORD   ) --- ', e11.4, ' MICROSECONDS'/ )
 2260 FORMAT (' GETSTR - READ  A STRING OF DATA                ',  &
     ' (AVERAGE PER RECORD   ) --- ', e11.4, ' MICROSECONDS'/ )
 2270 FORMAT (' PUTSTR - WRITE A STRING OF DATA                ',  &
     ' (AVERAGE PER RECORD   ) --- ', e11.4, ' MICROSECONDS'/ )
END SUBROUTINE tmtsot
