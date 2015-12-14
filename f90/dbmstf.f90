SUBROUTINE dbmstf
     INCLUDE   'DSIOF.COM'
 COMMON / system / isysbf, iwr
 COMMON / logout / lout
 iblksz = isysbf - 4
 IF ( maxblk /= 0 )  perc1 = maxblk*1.0 / maxalc
 iperc1 = perc1 * 100.
 imemnu = ( maxalc - maxblk ) * lenalc
 WRITE( lout, 901 ) lenopc, idblen, maxblk, maxalc, iperc1, maxdsk  &
     , iblksz, numopn, numcls, numwri, numrea
 IF ( idbdir /= 0 ) WRITE( lout, 902 ) imemnu
 901   FORMAT(1H1  &
     ,5X,'STATISTICS ON IN-MEMORY DATA BASE AND DISK I/O USAGE',/  &
     ,/,8X,' LENGTH (IN WORDS) OF OPEN CORE ALLOCATED          ',i16  &
     ,/,8X,' LENGTH (IN WORDS) OF IN-MEMORY DATA BASE ALLOCATED',i16  &
     ,/,8X,' NUMBER OF BLOCKS USED IN THE IN-MEMORY DATA BASE  ',i16  &
     ,/,8X,' NUMBER OF BLOCKS ALLOCATED FOR THE IN-MEMORY DATA ',i16  &
     ,/,8X,' PERCENTAGE OF IN-MEMORY DATA USED                 ',i16 ,'%'  &
     ,/,8X,' TOTAL BLOCKS WRITTEN TO DISK                      ',i16  &
     ,/,8X,' BLOCK SIZE (IN WORDS)                             ',i16  &
     ,/,8X,' NUMBER OF OPENS TO DISK FILES                     ',i16  &
     ,/,8X,' NUMBER OF CLOSES TO DISK FILES                    ',i16  &
     ,/,8X,' NUMBER OF WRITES TO DISK FILES                    ',i16  &
     ,/,8X,' NUMBER OF READS FROM DISK FILES                   ',i16)
 902   FORMAT( 8X,' MEMORY (IN WORDS) NOT USED BY IN-MEM. DATA BASE   ',i16  &
     )
 RETURN
END SUBROUTINE dbmstf
