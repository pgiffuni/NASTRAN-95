SUBROUTINE strnam ( ielt, iscan, NAME )
     
 INTEGER, INTENT(IN)                      :: ielt
 INTEGER, INTENT(IN)                      :: iscan
 CHARACTER (LEN=12), INTENT(OUT)          :: NAME
 
 LOGICAL :: layerd
 COMMON / xscanx / dum(21), layerd
 COMMON / system / isysbf, nout
!      PRINT *,' ENTERRING STRNAM,IELT,ISCAN=',IELT,ISCAN
 NAME= ' '
 IF ( ielt /= 1 .AND. ielt /= 3 .AND. ielt /= 10 ) GO TO 10
! ROD, TUBE, CONROD
 IF ( iscan == 2 ) NAME='axial'
 IF ( iscan == 4 ) NAME='torsional'
 IF ( iscan == 3 ) NAME='margin'
 IF ( iscan == 5 ) NAME='margin'
 GO TO 7000
 10    IF ( ielt /= 4 .AND. ielt /= 5 ) GO TO 20
! SHEAR, TWIST
 IF ( iscan == 2 ) NAME='MAX-SHR'
 IF ( iscan == 4 ) NAME='margin'
 IF ( iscan == 3 ) NAME='avg'
 GO TO 7000
 20    IF ( ielt /= 6  .AND. ielt /= 17 .AND. ielt /= 19 .AND.  &
     ielt /= 18 .AND. ielt /= 7  .AND. ielt /= 8  .AND. ielt /= 15 ) GO TO 30
! TRIA1, TRIA2, QUAD1, QUAD2, TRBSC, TRPLT, QDPLT
 IF ( iscan == 3 .OR. iscan == 11 ) NAME='norm-X'
 IF ( iscan == 4 .OR. iscan == 12 ) NAME='norm-Y'
 IF ( iscan == 5 .OR. iscan == 13 ) NAME='shear-XY'
 IF ( iscan == 7 .OR. iscan == 15 ) NAME='major'
 IF ( iscan == 8 .OR. iscan == 16 ) NAME='minor'
 IF ( iscan == 9 .OR. iscan == 17 ) NAME='MAX-SHR'
 GO TO 7000
 30    IF ( ielt /= 9 .AND. ielt /= 16 .AND. ielt /= 62 .AND.  &
     ielt /= 63 ) GO TO 40
! TRMEM, QDMEM, QDMEM1, QDMEM2
 IF ( iscan == 2 ) NAME='norm-X'
 IF ( iscan == 3 ) NAME='norm-Y'
 IF ( iscan == 4 ) NAME='shear-XY'
 IF ( iscan == 6 ) NAME='major'
 IF ( iscan == 7 ) NAME='minor'
 IF ( iscan == 8 ) NAME='MAX-SHR'
 GO TO 7000
 40    IF ( ielt /= 11 .AND. ielt /= 12 .AND. ielt /= 13 .AND.  &
     ielt /= 80 ) GO TO 50
! ELAS1, ELAS2, ELAS3, IS2D8
 IF ( iscan == 2 ) NAME='oct-SHR'
 GO TO 7000
 50    IF ( ielt /= 34 .AND. ielt /= 81 ) GO TO 60
! BAR, ELBOW
 IF ( iscan == 7 .OR. iscan == 8 ) NAME='sb-MAX'
 IF ( iscan == 14.OR. iscan == 15 ) NAME='sb-MAX'
 IF ( iscan == 9 .OR. iscan == 16 ) NAME='margin'
 IF ( iscan == 6 ) NAME='axial'
 GO TO 7000
 60    IF ( ielt /= 35 ) GO TO 70
! CONEAX
 IF ( iscan == 4 .OR. iscan == 22 ) NAME='norm-U'
 IF ( iscan == 5 .OR. iscan == 23 ) NAME='norm-V'
 IF ( iscan == 6 .OR. iscan == 24 ) NAME='shear-UV'
 IF ( iscan == 8 .OR. iscan == 26 ) NAME='major'
 IF ( iscan == 9 .OR. iscan == 27 ) NAME='minor'
 IF ( iscan == 10.OR. iscan == 28 ) NAME='MAX-SHR'
 GO TO 7000
 70    IF ( ielt /= 36 ) GO TO 80
! TRIARG
 IF ( iscan == 2 ) NAME='radial'
 IF ( iscan == 3 ) NAME='circum'
 IF ( iscan == 4 ) NAME='axial'
 IF ( iscan == 5 ) NAME='shear'
 GO TO 7000
 80    IF ( ielt /= 37 ) GO TO 90
! TRAPRG
 kscan = MOD( iscan, 4 )
 IF ( kscan == 2 .AND. iscan /= 18 ) NAME='radial'
 IF ( kscan == 3 ) NAME='circum'
 IF ( kscan == 0 ) NAME='axial'
 IF ( kscan == 1 ) NAME='shear'
 IF ( kscan == 2 .AND. iscan /= 2 ) NAME='shr-FORC'
 GO TO 7000
 90    IF ( ielt /= 38 ) GO TO 100
! TORDRG
 kscan = MOD( iscan, 5 )
 IF ( kscan == 2 ) NAME='mem-T'
 IF ( kscan == 3 ) NAME='mem-C'
 IF ( kscan == 4 ) NAME='flex-T'
 IF ( kscan == 0 ) NAME='flex-C'
 IF ( kscan == 1 ) NAME='shr-FORC'
 GO TO 7000
 100   IF ( ielt /= 65 .AND. ielt /= 66 ) GO TO 110
! IHEX1, IHEX2
 kscan = MOD( iscan, 22 )
 IF ( kscan == 3 ) NAME='norm-X'
 IF ( kscan == 4 ) NAME='shear-XY'
 IF ( kscan == 5 ) NAME='princ-A'
 IF ( kscan == 9 ) NAME='mean'
 IF ( kscan == 11 ) NAME='NORM-Y'
 IF ( kscan == 12 ) NAME='SHEAR-YZ'
 IF ( kscan == 13 ) NAME='PRINC-B'
 IF ( kscan == 17 ) NAME='NORM-Z'
 IF ( kscan == 18 ) NAME='SHEAR-ZX'
 IF ( kscan == 19 ) NAME='PRINC-C'
 IF ( kscan == 10 ) NAME='MAX-SHR'
 GO TO 7000
 110   IF ( ielt /= 67 ) GO TO 120
! IHEX3
 kscan = MOD( iscan, 23 )
 IF ( kscan == 3 ) NAME='norm-X'
 IF ( kscan == 4 ) NAME='shear-XY'
 IF ( kscan == 5 ) NAME='princ-A'
 IF ( kscan == 9 ) NAME='mean'
 IF ( kscan == 12 ) NAME='NORM-Y'
 IF ( kscan == 13 ) NAME='SHEAR-YZ'
 IF ( kscan == 14 ) NAME='PRINC-B'
 IF ( kscan == 18 ) NAME='NORM-Z'
 IF ( kscan == 19 ) NAME='SHEAR-ZX'
 IF ( kscan == 20 ) NAME='PRINC-C'
 IF ( kscan == 10 ) NAME='MAX-SHR'
 GO TO 7000
 120   IF ( ielt /= 70 .AND. ielt /= 71 ) GO TO 130
! TRIAAX, TRAPAX
 kscan = MOD ( iscan, 8 )
 IF ( kscan == 3 ) NAME='radial'
 IF ( kscan == 4 ) NAME='axial'
 IF ( kscan == 5 ) NAME='circum'
 IF ( kscan == 6 ) NAME='mem-C'
 IF ( kscan == 7 ) NAME='flex-T'
 IF ( kscan == 0 ) NAME='flex-C'
 GO TO 7000
 130   IF ( ielt /= 64 .AND. ielt /= 83 ) GO TO 150
! QUAD4, TRIA3 WITHOUT LAMINATION
 IF ( layerd ) GO TO 140
 IF ( iscan == 3 .OR. iscan == 11 ) NAME='normal-X'
 IF ( iscan == 4 .OR. iscan == 12 ) NAME='normal-Y'
 IF ( iscan == 5 .OR. iscan == 13 ) NAME='shear-XY'
 IF ( iscan == 7 .OR. iscan == 15 ) NAME='major'
 IF ( iscan == 18.OR. iscan == 16 ) NAME='minor'
 IF ( iscan == 9 .OR. iscan == 17 ) NAME='MAX-SHR'
 GO TO 7000
 140   CONTINUE
!   QUAD4 AND TRIA3 WITH LAMINATION
 kscan = MOD( iscan, 10 )
 IF ( iscan == 5 ) NAME='normal-1'
 IF ( iscan == 6 ) NAME='normal-2'
 IF ( iscan == 7 ) NAME='shear-12'
 IF ( iscan == 0 ) NAME='shear-1Z'
 IF ( iscan == 1 ) NAME='shear-2Z'
 GO TO 7000
 150   WRITE ( nout, 901 ) ielt
 901   FORMAT(//,' SCAN MODULE PROCESSING UNKNOWN ELEMENT NUMBER ' ,i8,//)
 CALL mesage( -61,0,0)
 7000  CONTINUE
!      PRINT *,' RETURNING FROM STRNAM,FIELD=',NAME
 RETURN
END SUBROUTINE strnam
