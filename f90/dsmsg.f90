SUBROUTINE dsmsg ( iflag )
     INCLUDE 'DSIOF.COM'
 
 INTEGER, INTENT(IN)                      :: iflag
 COMMON / dsio   / mssoft
 INCLUDE 'GINOX.COM'
 COMMON / logout / lout
 INCLUDE 'XNSTRN.COM'
 COMMON / ddiosv / iflpos(2,maxpri)
 COMMON / sofcom / nfiles,filnam(10),filsiz(10)
 COMMON / sys    / blksiz,dirsiz,supsiz,avblks,hiblk
 COMMON / system / isystm(175)
 
 INTEGER :: xname(2), BLANK
 INTEGER :: filnam,filsiz,blksiz,dirsiz,supsiz,avblks,hiblk
 INTEGER :: gino(52)
 
 EQUIVALENCE ( IEOR, gino(1) )
 EQUIVALENCE (isystm(  2), iwr   ), (isystm(151), nllog ),  &
     (isystm(152), loglin)
 
 DATA BLANK / 1H /
 DATA iname /4HDSMS/
 
 CALL fname ( NAME, xname )
 IF ( xname( 1 ) /= 0 ) GO TO 4
 xname( 1 ) = BLANK
 xname( 2 ) = BLANK
 4 CONTINUE
 IF ( IABS(iflag) == 777 ) GO TO 7770
 IF ( iflag /= 1 .AND. iflag /= 2 .AND. iflag /= 8 ) WRITE( iwr, 5 ) iflag
 5 FORMAT(' I/O SUBSYSTEM ERROR NUMBER',i10)
 IF ( iflag > 100 ) GO TO 1000
 SELECT CASE ( iflag )
   CASE (    1)
     GO TO  10
   CASE (    2)
     GO TO  20
   CASE (    3)
     GO TO  30
   CASE (    4)
     GO TO  40
   CASE (    5)
     GO TO  50
   CASE (    6)
     GO TO  60
   CASE (    7)
     GO TO  70
   CASE (    8)
     GO TO  80
   CASE (    9)
     GO TO  90
 END SELECT
 10 CONTINUE
 WRITE ( lout, 15 ) 'OPEN ',xname, iocode
 loglin = loglin + 1
 GO TO 90000
 20 CONTINUE
 WRITE ( lout, 15 ) 'CLOSE ', xname, iocode
 loglin = loglin + 1
 GO TO 90000
 30 WRITE ( iwr, 35 ) xname, ifilex
 GO TO 99999
 40 WRITE ( iwr, 45 ) xname, ifilex
 GO TO 7000
 50 WRITE ( iwr, 55 ) xname, ifilex
 GO TO 99999
 60 WRITE ( iwr, 65 ) xname, ifilex
 GO TO 99999
 70 WRITE ( iwr, 75 ) xname, ifilex
 GO TO 7000
 80 CONTINUE
 WRITE ( lout, 85 ) xname, ifilex, idsn
 GO TO 90000
 90 WRITE ( iwr, 95 ) nblock
 GO TO 99999
 100 CONTINUE
 GO TO 7000
 1000 itemp = iflag - 100
 GO TO ( 1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080  &
     ,1090, 1100, 1110, 1120, 1130, 1140, 1150, 1160  &
     ,1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240 ), itemp
!1010 WRITE ( IWR, 1015 ) IOERR
 1010 GO TO 7000
 1020 WRITE ( iwr, 1025 ) xname, ifilex
 GO TO 7000
 1030 WRITE ( iwr, 1035 ) xname, ifilex
 GO TO 7000
 1040 WRITE ( iwr, 1045 ) xname, ifilex
 GO TO 7000
 1050 WRITE ( iwr, 1055 ) xname, ifilex
 GO TO 7000
 1060 WRITE ( iwr, 1065 ) xname, ifilex
 GO TO 7000
 1070 WRITE ( iwr, 1075 ) xname, ifilex
 GO TO 7000
 1080 WRITE ( iwr, 1085 ) xname, ifilex
 GO TO 7000
 1090 WRITE ( iwr, 1095 ) xname, ifilex
 GO TO 7000
 1100 WRITE ( iwr, 1105 ) xname, ifilex
 GO TO 7000
 1110 WRITE ( iwr, 1115 ) xname, ifilex
 GO TO 7000
 1120 WRITE ( iwr, 1125 ) xname, ifilex
 GO TO 7000
 1130 WRITE ( iwr, 1135 ) xname, ifilex
 GO TO 7000
 1140 WRITE ( iwr, 1145 ) xname, ifilex
 GO TO 7000
 1150 WRITE ( iwr, 1155 ) xname, ifilex
 GO TO 7000
 1160 WRITE ( iwr, 1165 ) xname, ifilex
 GO TO 7000
 1170 WRITE ( iwr, 1175 ) xname, ifilex
 GO TO 7000
 1180 WRITE ( iwr, 1185 ) xname, ifilex
 GO TO 7000
 1190 WRITE ( iwr, 1195 ) xname, ifilex
 GO TO 7000
 1200 WRITE ( iwr, 1205 ) xname, ifilex
 GO TO 7000
 1210 WRITE ( iwr, 1215 ) xname, ifilex
 GO TO 7000
 1220 WRITE ( iwr, 1225 ) xname, ifilex
 GO TO 7000
 1230 CONTINUE
 1240 CONTINUE
 7000 CONTINUE
 7770  CONTINUE
 WRITE ( iwr, 91000 ) ioerr, NAME, xname, ifilex
 WRITE ( iwr, 92000 )
 DO  i = 1, maxfcb
   CALL dshxdd ( i, mdsfcb( 1, i ), 3 )
 END DO
 WRITE( iwr, 92001 )
 DO  i = 1, 80
!WKBR NCL93007 11/94
!      WRITE ( IWR, 92003 ) I, ( FCB(K,I),K=1,15)
   WRITE ( iwr, 92003 ) i, ( fcb(k,i),k=1,17)
!WKBR NCL93007 11/04
!92003 FORMAT(I3,'-',I3,I7,4I5,I12,I2,4I7,1X,2A4,I4)
   92003 FORMAT(i3,'-',i3,i7,4I5,i12,i2,4I7,1X,2A4,i4,2I8)
 END DO
 WRITE ( iwr, 92002)idbbas, idbfre, idbdir, indbas, indclr, indcbp  &
     ,                  nblock, lenalc, iocode, ifilex, NAME,   maxalc  &
     ,                  maxblk, maxdsk, idblen, idbadr, ibasbf, inddir  &
     ,                  numopn, numcls, numwri, numrea, lenopc
 92002 FORMAT(/,' CONTENTS OF / DBM / FOLLOW:'  &
     ,/,' IDBBAS =',i8,' IDBFRE =',i8,' IDBDIR =',i8,' INDBAS =',i8  &
     ,/,' INDCLR =',i8,' INDCBP =',i8,' NBLOCK =',i8,' LENALC =',i8  &
     ,/,' IOCODE =',i8,' IFILEX =',i8,' NAME   =',i8,' MAXALC =',i8  &
     ,/,' MAXBLK =',i8,' MAXDSK =',i8,' IDBLEN =',i8,' IDBADR =',i8  &
     ,/,' IBASBF =',i8,' INDDIR =',i8,' NUMOPN =',i8,' NUMCLS =',i8  &
     ,/,' NUMWRI =',i8,' NUMREA =',i8,' LENOPC =',i8)
 IF ( iflag >= 118 .AND. iflag <= 120 ) GO TO 946
 WRITE ( iwr, 95020 )
 CALL dshxdp ( nfiles, 16 )
 WRITE ( iwr, 95030 )
 CALL dshxdp ( blksiz, 1 )
 WRITE ( iwr, 96000 )
 CALL dshxdp ( IEOR, 59 )
 WRITE ( iwr, 97000 )
 DO  i = 1, maxpri
   WRITE ( iwr, 97001 ) i, iflpos(1,i), iflpos(2,i)
   97001 FORMAT(i5,2I10)
 END DO
 loop = (nbuff+lendsp) / 8 + 4
 INDEX = indbas
 WRITE ( iwr, 94510 )
 DO  i = 1, loop
   iii = (i-1) * 8 + 1
   CALL dshxdd ( iii, ibase( INDEX ), 8 )
   INDEX = INDEX + 8
 END DO
 CALL dbmdia
 946 IF ( iflag /= 777 ) GO TO 99999
!     CALL TRBK( IWR )
 RETURN
 91000 FORMAT(' I/O ERROR #',i6,' ON FILE ',z8,' NAME=',2A4,' UNIT=',i4)
 92000 FORMAT(//' CONTENTS OF MDSFCB' )
 92001 FORMAT(//' CONTENTS OF FCB' )
 94510 FORMAT(//' CONTENTS OF I/O BUFFER' )
 95020 FORMAT(//' CONTENTS OF SOFCOM ')
 95030 FORMAT(//' CONTENTS OF SYS ')
 96000 FORMAT(//' CONTENTS OF /DSIO/')
 97000 FORMAT(//' CONTENTS OF /DDIOSV/')
 99999 CALL mesage (-61, 0, 0)
 90000 CONTINUE
 RETURN
 15 FORMAT( 40X, a6, 2A4, 2X, i2 )
 35 FORMAT( ' BUFFER CONFLICTS WITH EXISTING BUFFERS',  &
     ' ON FILE ',2A4, ' LOGICAL UNIT', i4 )
 45 FORMAT(' ATTEMPT TO READ FILE OPENED FOR WRITE',  &
     ' FILE=',2A4,' UNIT=',i4 )
 55 FORMAT(' FILE IS ALREADY OPENED-FILE ',2A4, ' UNIT=', i4 )
 65 FORMAT(' ATTEMPT TO WRITE LESS THAN ONE WORD',  &
     ' ON FILE ',2A4,' UNIT= ',i4 )
 75 FORMAT(' ATTEMPT TO WRITE ON FILE OPENED FOR READ ',  &
     '-FILE=',2A4,' UNIT =',i4)
85 FORMAT(//,' ****** GINO SUBSYSTEM WILL EXTEND FILE ',2A4,  &
    ' ON UNIT',i4,' TO UNIT',i4,' ******' )
95 FORMAT(//,' INSUFFICIENT SPACE ALLOCATION ON FILE NPTP',  &
    ' -, NUMBER OF BLOCKS WRITTEN WERE ',i10)
1015 FORMAT(' ERROR DURING I/O REQUEST - ERROR FLAG=',z8)
1025 FORMAT(' INCORRECT BLOCK NUMBER ENCOUNTERED',  &
    ' ON FILE ',2A4,' UNIT=',i4)
1035 FORMAT(' EXPECTED RH, SB, EF, OR EB CONTROL WORD',  &
    ' ON FILE ',2A4,' UNIT=',i4)
1045 FORMAT(' EXPECTED RT CONTROL WORD ON FILE ',2A4, ' UNIT=',i4)
1055 FORMAT(' EXPECTED RH OR EF CONTROL WORD ON FILE ',2A4, ' UNIT=',i4)
1065 FORMAT(' EXPECTED RH, EB OR SB CONTROL WORD ON FILE ',2A4, ' UNIT=',i4)
1075 FORMAT(' REFERENCE IS MADE TO FILE ',2A4, ' THAT IS NOT OPENED-UNIT=',i4)
1085 FORMAT(' INSUFFICIENT SPACE FOR I/O CONTROL WORDS ON FILE '  &
    ,2A4,' UNIT=',i4)
1095 FORMAT(' TOO MANY TERMS WRITTEN TO STRING ON FILE ',2A4, ' UNIT=',i4)
1105 FORMAT(' EXPECTED A SB OR EB CONTROL WORD ON FILE ',2A4, ' UNIT=',i4)
1115 FORMAT(' EXPECTED A CH CONTROL WORD ON FILE ',2A4, ' UNIT=',i4)
1125 FORMAT(' EXPECTED A SE, SD, CT, OR SH CONTROL WORD ON FILE ',  &
    2A4,' UNIT=',i4)
1135 FORMAT(' ERROR  - CLR.GT. LCW  ON FILE ',2A4, ' UNIT=',i4)
1145 FORMAT(' EXPECTED A RT CONTROL WORD ON FILE ',2A4, ' UNIT=',i4)
1155 FORMAT(' EXPECTED A CH CONTROL WORD ON FILE ',2A4,' UNIT=',i4)
1165 FORMAT(' EXPECTED A CH,ST,SH,SD,RT, OR SE CONTROL WORD ON FILE ',  &
    2A4,' UNIT=',i4)
1175 FORMAT(' EXPECTED A ST CONTROL WORD ON FILE ',2A4, ' UNIT=',i4)
1185 FORMAT(' TYPIN OR TYPOUT FOR MATRIX PACK IS OUT OF RANGE ON',  &
    ' FILE ',2A4,' UNIT=',i4)
1195 FORMAT(' NON-ASCENDING ROW NUMBER GIVEN', ' ON FILE ',2A4, ' UNIT=',i10)
1205 FORMAT(' FILE NAME DOES NOT MATCH STRING CONTROL BLOCK FOR ',  &
    'FILE ',2A4,' UNIT=',i4)
1215 FORMAT(' INVALID UNIT NUMBER IN MDSFCB FOR FILE ',2A4,' UNIT=',i4)
1225 FORMAT(' INSUFFICIENT NUMBER OF FILES AVAILABLE FOR FILE ',  &
    2A4,' UNIT=',i4)
END SUBROUTINE dsmsg
