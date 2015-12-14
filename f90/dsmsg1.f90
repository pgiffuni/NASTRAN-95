SUBROUTINE dsmsg1 ( BLOCK )
     
 INTEGER, INTENT(IN OUT)                  :: BLOCK(15)
 COMMON / zblpkx / a1(4), irow1
 COMMON / zntpkx / a2(4), irow2, ieol2, ieor2
 COMMON / packx  / itin3, itout3, irow3, nrow3, incr3
 COMMON / unpakx / itout4,irow4, nrow4, incr4
 COMMON / system / NONE,  iwr
 
 WRITE( iwr, 9000 )
 WRITE( iwr, 9010 )
 WRITE( iwr, 9015 ) BLOCK
 WRITE( iwr, 9020 )
 WRITE( iwr, 9015 ) a1, irow1
 WRITE( iwr, 9030 )
 WRITE( iwr, 9015 ) a2, irow2, ieol2, ieor2
 WRITE( iwr, 9040 )
 WRITE( iwr, 9015 ) itin3, itout3, irow3, nrow3, incr3
 WRITE( iwr, 9050 )
 WRITE( iwr, 9015 ) itout4, irow4, nrow4, incr4
 9000    FORMAT(' *** ERROR OCCURRED IN PAKUNPK I/O SUBSYSTEM ***')
 9010    FORMAT(' CONTENTS OF THE STRING CONTROL BLOCK')
 9015    FORMAT(10(5(1X,z8),/))
 9020    FORMAT(' CONTENTS OF /ZBLPKX/')
 9030    FORMAT(' CONTENTS OF /ZNTPKX/')
 9040    FORMAT(' CONTENTS OF /PACKX/ ')
 9050    FORMAT(' CONTENTS OF /UNPAKX/')
 RETURN
 
END SUBROUTINE dsmsg1
