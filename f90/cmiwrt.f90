SUBROUTINE cmiwrt (icode,name1,name2,loc,nw,a,iz)
     
!     THIS SUBROUTINE WRITES FORMATTED SOF ITEMS.
!     ICODE = 1 FOR EQSS    ICODE = 2 FOR BGSS    ICODE = 3 FOR CSTM
!     ICODE = 4 FOR PLTS    ICODE = 5 FOR LODS    ICODE = 7 FOR LOAP
!     NAME1 IS PSEUDOSTRUCTURE NAME, NAME2 IS COMPONENT NAME
 
 EXTERNAL        andf
 
 INTEGER, INTENT(IN OUT)                  :: icode
 INTEGER, INTENT(IN)                      :: name1(2)
 INTEGER, INTENT(IN OUT)                  :: name2(2)
 INTEGER, INTENT(IN)                      :: loc
 INTEGER, INTENT(IN)                      :: nw
 INTEGER, INTENT(OUT)                     :: iz(1)
 INTEGER :: outt,andf

 REAL, INTENT(IN OUT)                     :: a(1)

 DIMENSION  ibits(32),ipl(6), ih1(96),ih2(96),ih3(96),ih4(96),ih5(96),ih6(96)

 COMMON /system/ xxx,outt,junk1(6),nlpp,junk2(2),nline
 COMMON /output/ ititl(96),ihead(96)

 DATA ih1 / 9*4H    ,4H EQS,4HS IT,4HEM F,4HOR S,4HUBST,4HRUCT,     &
            4HURE ,2*4H    ,4H COM,4HPONE,4HNT  ,11*4H    ,4HGRID,  &
            4H POI,4HNT  ,4H INT,4HERNA,4HL   ,4H  CO,4HMPON,	    &
            4HENT ,2*4H    ,4H GRI,4HD PO,4HINT ,4H  IN,4HTERN,	    &
            4HAL  ,4H   C,4HOMPO,4HNENT,2*4H    ,4H  GR,4HID P,	    &
            4HOINT,4H   I,4HNTER,4HNAL ,4H    ,4HCOMP,4HONEN,	    &
            4HT   ,4H    ,4HID  ,4H    ,4H POI,4HNT I,4HD   ,	    &
            4H    ,4H DOF,4*4H    ,4H ID ,4H    ,4H  PO,4HINT ,	    &
            4HID  ,4H    ,4H  DO,4HF   ,3*4H    ,4H  ID,4H    ,	    &
            4H   P,4HOINT,4H ID ,4H    ,4H   D,4HOF  ,4H     /	    
 DATA ih2 / 11*4H    ,4HBGSS,4H ITE,4HM FO,4HR SU,4HBSTR,4HUCTU,    &
            4HRE  ,21*4H    ,4HINTE,4HRNAL,4H    ,4H CST,4HM ID,    &
            4*4H    ,4H  C ,4HO O ,4HR D ,4HI N ,4HA T ,4HE S ,	    &
            17*4H    ,4HPOIN,4HT ID,4H    ,4H   N,4HO.  ,3*4H    ,  &
            4HX1  ,3*4H    ,4HX2  ,3*4H    ,4HX3  ,8*4H     /	    
 DATA ih3 / 12*4H    ,4HCSTM,4H ITE,4HM FO,4HR SU,4HBSTR,	    &
            4HUCTU,4HRE  ,13*4H    ,2*4H    ,4H CST,4HM   ,4HTYPE,  &
            2*4H    ,4HC O ,4HO R ,4HD I ,4HN A ,4HT E ,4HS   ,	    &
            4HO F ,4H  O ,4HR I ,4HG I ,4HN   ,3*4H    ,4H   T,	    &
            4H R A,4H N S,4H F O,4H R M,4H A T,4H I O,4H N  ,	    &
            5*4H    ,4H  ID,5*4H    ,4HX1  ,3*4H    ,4HX2  ,	    &
            3*4H    ,4HX3  ,6*4H    ,4H   M,4H A T,4H R I,4H X  ,   &
            5*4H      /						    
 DATA ih4 / 12*4H    ,4HPLTS,4H ITE,4HM FO,4HR SU,4HBSTR,	    &
            4HUCTU,4HRE  ,13*4H    ,2*4H    ,4HCOMP,4HONEN,4HT   ,  &
            4H    ,4H C O,4H O R,4H D I,4H N A,4H T E,4H S  ,	    &
            4H O F,4H   O,4HR I ,4HG I ,4HN   ,3*4H    ,4H   T,	    &
            4H R A,4H N S,4H F O,4H R M,4H A T,4H I O,4H N  ,	    &
            6*4H    ,4H  NA,4HME  ,3*4H    ,4H X1 ,3*4H    ,	    &
            4H X2 ,3*4H    ,4H X3 ,6*4H    ,4H   M,4H A T,4H R I,   &
            4H X  ,6*4H      /					    
 DATA ih5 / 12*4H    ,4HLODS,4H ITE,4HM FO,4HR SU,4HBSTR,	    &
            4HUCTU,4HRE  ,18*4H    ,4H COM,4HPONE,4HNT  ,4H  NU,    &
            4HMBER,4H OF ,21*4H    ,5*4H    ,4H   N,4HAME ,	    &
            4H    ,4H  LO,4HAD S,4HETS ,4H  L ,4HO A ,4HD   ,	    &
            4HS E ,4HT   ,4HI D ,4HE N ,4HT I ,4HF I ,4HC A ,	    &
            4HT I ,4HO N ,4H  N ,4HU M ,4HB E ,4HR S ,5*4H      /   
 DATA ih6 / 9*4H         ,4HEQSS,4H ITE,4HM - ,4HSCAL,4HAR I,	    &
            4HNDEX,4H LIS,4HT FO,4HR SU,4HBSTR,4HUCTU,4HRE  ,	    &
            11*4H        ,4H INT,4HERNA,4HL   ,4H INT,4HERNA,	    &
            4HL   ,4H  CO,4HMPON,4HENT ,2*4H         ,4H  IN,	    &
            4HTERN,4HAL  ,4H  IN,4HTERN,4HAL  ,4H   C,4HOMPO,	    &
            4HNENT,2*4H         ,4H   I,4HNTER,4HNAL ,4H   I,	    &
            4HNTER,4HNAL ,4H    ,4HCOMP,4HONEN,4HT   ,4H POI,	    &
            4HNT I,4HD   ,4H  SI,4HL ID,2*4H         ,4H DOF,	    &
            3*4H         ,4H  PO,4HINT ,4HID  ,4H   S,4HIL I,	    &
            4HD   ,4H    ,4H  DO,4HF   ,2*4H         ,4H   P,	    &
            4HOINT,4H ID ,4H    ,4HSIL ,4HID  ,4H    ,4H   D,	    &
            4HOF  ,4H     /
 DATA loap/ 4HLOAP/
 
 ist  = loc
 ifin = loc + nw - 1
 SELECT CASE ( icode )
   CASE (    1)
     GO TO 1
   CASE (    2)
     GO TO 2
   CASE (    3)
     GO TO 3
   CASE (    4)
     GO TO 4
   CASE (    5)
     GO TO 5
   CASE (    6)
     GO TO 6
   CASE (    7)
     GO TO 5
   CASE (    8)
     GO TO 8
 END SELECT
 
!     EQSS ITEM
 
 1 DO  i = 1,96
   ihead(i) = ih1(i)
 END DO
 
!     INSERT NAMES INTO HEADING
 
 ihead(17) = name1(1)
 ihead(18) = name1(2)
 ihead(22) = name2(1)
 ihead(23) = name2(2)
 CALL page
 IF (nw /= 0) GO TO 140
 WRITE (outt,1009)
 GO TO 700
 
 140 DO  i = ist,ifin,9
   nline = nline + 1
   IF (nline <= nlpp) GO TO 150
   CALL page
   nline = nline + 1
   150 CONTINUE
   icomp = andf(iz(i+2),63)
   CALL bitpat (icomp,ibits(1))
   i2 = 3
   IF (i+5 > ifin) GO TO 151
   icomp = andf(iz(i+5),63)
   CALL bitpat (icomp,ibits(4))
   i2 = 6
   IF (i+8 > ifin) GO TO 151
   icomp = andf(iz(i+8),63)
   CALL bitpat (icomp,ibits(7))
   i2 = 9
   151 CONTINUE
   WRITE (outt,1000) (iz(i+j-1),iz(i+j),ibits(j),ibits(j+1),j=1,i2,3)
 END DO
 GO TO 700
 
!     EQSS - SCALER INDEX LIST
 
 8 DO  i = 1,96
   ihead(i) = ih6(i)
 END DO
 ihead(22) = name1(1)
 ihead(23) = name1(2)
 CALL page
 
 ip = 0
 DO  i = ist,ifin,6
   nline = nline + 1
   IF (nline <= nlpp) GO TO 601
   CALL page
   nline = nline + 1
   601 CONTINUE
   kcode = iz(i+1)
   CALL bitpat (kcode,ibits(1))
   i2 = 2
   ipl(1) = ip + 1
   IF (i+3 > ifin) GO TO 602
   kcode = iz(i+3)
   CALL bitpat (kcode,ibits(3))
   i2 = 4
   ipl(3) = ip + 2
   IF (i+5 > ifin) GO TO 602
   kcode = iz(i+5)
   CALL bitpat (kcode,ibits(5))
   i2 = 6
   ipl(5) = ip + 3
   602 CONTINUE
   WRITE (outt,1000) (ipl(j),iz(i+j-1),ibits(j),ibits(j+1),j=1,i2,2)
   ip = ip + 3
 END DO
 GO TO 700
 
!     BGSS ITEM
 
 2 DO  i = 1,96
   ihead(i)  = ih2(i)
 END DO
 ihead(20) = name1(1)
 ihead(21) = name1(2)
 CALL page
 j = 0
 DO  i = ist,ifin,4
   j = j + 1
   nline = nline + 1
   IF (nline <= nlpp) GO TO 250
   CALL page
   nline = nline + 1
   250 CONTINUE
   WRITE (outt,1001) j,iz(i),a(i+1),a(i+2),a(i+3)
 END DO
 GO TO 700
 
!     CSTM ITEM
 
 3 DO  i = 1,96
   ihead(i)  = ih3(i)
 END DO
 ihead(20) = name1(1)
 ihead(21) = name1(2)
 CALL page
 DO  i = ist,ifin,14
   nline = nline + 4
   IF (nline <= nlpp) GO TO 350
   CALL page
   nline = nline + 4
   350 CONTINUE
   i1 = i + 2
   i2 = i + 13
   WRITE (outt,1002) iz(i),iz(i+1),(a(kk),kk= i1,i2)
 END DO
 GO TO 700
 
 4 DO  i = 1,96
   
!     PLTS ITEM
   
   ihead(i)  = ih4(i)
 END DO
 ihead(20) = name1(1)
 ihead(21) = name1(2)
 CALL page
 DO  i = ist,ifin,14
   nline = nline + 4
   IF (nline <= nlpp) GO TO 450
   CALL page
   nline = nline + 4
   450 CONTINUE
   i1 = i + 2
   i2 = i + 13
   WRITE (outt,1004) iz(i),iz(i+1),(a(j),j=i1,i2)
 END DO
 GO TO 700
 
!     LODS AND LOAP ITEMS
 
 5 DO  i = 1,96
   ihead(i)  = ih5(i)
 END DO
 ihead(20) = name1(1)
 ihead(21) = name1(2)
 IF (icode == 7) ihead(13) = loap
 CALL page
 6 IF (nw == 0 .OR. nw == 1) GO TO 520
 nl = nw/5 + 3
 nline = nline + nl
 IF (nline <= nlpp) GO TO 550
 CALL page
 nline = nline + nl
 550 CONTINUE
 ist1 = ist + 1
 WRITE (outt,1006) name2(1),name2(2),iz(ist),(iz(j),j=ist1,ifin)
 GO TO 700
 
 520 nline = nline + 2
 IF (nline <= nlpp) GO TO 560
 CALL page
 nline = nline + 2
 560 CONTINUE
 WRITE (outt,1008) name2(1),name2(2)
 700 RETURN
 
 1000 FORMAT (6X,i8,4X,i8,6X,a4,a2,2(13X,i8,4X,i8,6X,a4,a2))
 1001 FORMAT (33X,i8,4X,i8,3X,3(3X,e13.6))
 1002 FORMAT (/10X,i8,3X,i4,3X,3(3X,e13.6),4X,3(3X,e13.6),  &
              /80X,3(3X,e13.6), /80X,3(3X,e13.6))
 1004 FORMAT (/14X,2A4,3X,3(3X,e13.6),4X,3(3X,e13.6)  &
              /77X,3(3X,e13.6), /77X,3(3X,e13.6))
 1006 FORMAT (/26X,2A4,3X,i8,5X,6(2X,i8)/(50X,2X,i8,2X,i8,2X,i8,2X,i8,  &
              2X,i8,2X,i8,/))
 1008 FORMAT (/26X,2A4,17X,32HNO load sets for this component. )
 1009 FORMAT (/30X,64HALL degrees of freedom for this component have been &
                  &reduced out. )
     
END SUBROUTINE cmiwrt
