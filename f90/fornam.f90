SUBROUTINE fornam ( ielt, iscan, NAME )
     
 INTEGER, INTENT(IN)                      :: ielt
 INTEGER, INTENT(IN)                      :: iscan
 CHARACTER (LEN=12), INTENT(OUT)          :: NAME
 COMMON / system / isysbf, nout
 
 
 NAME= ' '
 IF ( ielt /= 1 .AND. ielt /= 3 .AND. ielt /= 10 ) GO TO 10
! ROD, TUBE, CONROD
 IF ( iscan == 2 ) NAME='axial'
 IF ( iscan == 4 ) NAME='torque'
 GO TO 7000
 10    IF ( ielt /= 4 .AND. ielt /= 5 ) GO TO 20
! SHEAR, TWIST
 IF ( iscan == 2 ) NAME='force-1'
 IF ( iscan == 3 ) NAME='force-2'
 GO TO 7000
 20    IF ( ielt /= 6  .AND. ielt /= 17 .AND. ielt /= 19 .AND.  &
     ielt /= 18 .AND. ielt /= 7  .AND. ielt /= 8  .AND. ielt /= 15 ) GO TO 30
! TRIA1, TRIA2, QUAD1, QUAD2, TRBSC, TRPLT, QDPLT
 IF ( iscan == 2 ) NAME='moment-X'
 IF ( iscan == 3 ) NAME='moment-Y'
 IF ( iscan == 5 ) NAME='shear-X'
 IF ( iscan == 6 ) NAME='shear-Y'
 IF ( iscan == 4 ) NAME='twist'
 GO TO 7000
 30    IF ( ielt /= 9 .AND. ielt /= 16 .AND. ielt /= 62 .AND.  &
     ielt /= 63 ) GO TO 40
! TRMEM, QDMEM, QDMEM1, QDMEM2
 IF ( iscan == 3 .OR. iscan == 4 ) NAME='force-12'
 IF ( iscan == 5 .OR. iscan == 6 ) NAME='force-23'
 IF ( iscan == 7 .OR. iscan == 8 ) NAME='force-34'
 IF ( iscan == 2 .OR. iscan == 9 ) NAME='force-41'
 IF ( iscan == 10) NAME='kick ON1'
 IF ( iscan == 12) NAME='kick ON2'
 IF ( iscan == 14) NAME='kick ON3'
 IF ( iscan == 16) NAME='kick ON4'
 IF ( iscan == 11) NAME='shear-XY'
 IF ( iscan == 13) NAME='shear-YZ'
 IF ( iscan == 15) NAME='shear-ZX'
 IF ( iscan == 17) NAME='shear'
 GO TO 7000
 40    IF ( ielt /= 11 .AND. ielt /= 12 .AND. ielt /= 13 .AND.  &
     ielt /= 80 ) GO TO 50
! ELAS1, ELAS2, ELAS3, IS2D8
 IF ( iscan == 2 ) NAME='circum'
 IF ( iscan == 4 .AND. iscan == 9 ) NAME='force-1'
 IF ( iscan == 3 .AND. iscan == 6 ) NAME='force-2'
 IF ( iscan == 5 .AND. iscan == 8 ) NAME='force-3'
 IF ( iscan == 2 .AND. iscan == 7 ) NAME='force-4'
 GO TO 7000
 50    IF ( ielt /= 34 .AND. ielt /= 81 ) GO TO 60
! BAR, ELBOW
 IF ( iscan == 5 .OR. iscan == 6 ) NAME='shear'
 IF ( iscan == 2 .OR. iscan == 3 ) NAME='moment-A'
 IF ( iscan == 4 .OR. iscan == 5 ) NAME='moment-B'
 IF ( iscan == 8 ) NAME='axial'
 IF ( iscan == 9 ) NAME='torque'
 GO TO 7000
 60    IF ( ielt /= 35 ) GO TO 70
! CONEAX
 IF ( iscan == 3 ) NAME='moment-U'
 IF ( iscan == 4 ) NAME='moment-V'
 IF ( iscan == 6 ) NAME='shear-XY'
 IF ( iscan == 7 ) NAME='shear-YZ'
 GO TO 7000
 70    IF ( ielt /= 36 ) GO TO 80
! TRIARG
 kscan = MOD ( iscan, 3 )
 IF ( kscan == 2 ) NAME='radial'
 IF ( kscan == 3 ) NAME='circum'
 IF ( kscan == 1 ) NAME='axial'
 GO TO 7000
 80    IF ( ielt /= 37 ) GO TO 90
! TRAPRG
 kscan = MOD( iscan, 3 )
 IF ( kscan == 2 ) NAME='radial'
 IF ( kscan == 3 ) NAME='circum'
 IF ( kscan == 1 ) NAME='axial'
 GO TO 7000
 90    IF ( ielt /= 38 ) GO TO 120
! TORDRG
 kscan = MOD( iscan, 6 )
 IF ( kscan == 2 ) NAME='radial'
 IF ( kscan == 3 ) NAME='circum'
 IF ( kscan == 4 ) NAME='axial'
 IF ( kscan == 5 ) NAME='moment'
 IF ( kscan == 1 ) NAME='curv'
 GO TO 7000
 120   IF ( ielt /= 70 .AND. ielt /= 71 ) GO TO 130
! TRIAAX, TRAPAX
 kscan = MOD ( iscan, 4 )
 IF ( kscan == 3 ) NAME='radial'
 IF ( kscan == 0 ) NAME='circum'
 IF ( kscan == 1 ) NAME='axial'
 GO TO 7000
 130   IF ( ielt /= 64 .AND. ielt /= 83 ) GO TO 150
! QUAD4, TRIA3
 IF ( iscan == 2 .OR. iscan == 3 ) NAME='fx+FY'
 IF ( iscan == 4                   ) NAME='fxy'
 IF ( iscan == 5 .OR. iscan == 6 ) NAME='mx+MY'
 IF ( iscan == 7                   ) NAME='mxy'
 IF ( iscan == 8 .OR. iscan == 9 ) NAME='vx+VY'
 GO TO 7000
 150   WRITE ( nout, 901 ) ielt
 901   FORMAT(//' SCAN MODULE PROCESSING UNKNOWN ELEMENT NUMBER ' ,i8,//)
 CALL mesage( -61,0,0)
 7000  RETURN
END SUBROUTINE fornam
