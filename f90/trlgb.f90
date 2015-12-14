SUBROUTINE trlgb (usetd,ap,gmd,god,phidh,as,ad,ah,iflag1,scr1,  &
        scr2,scr3,scr4)
     
!     THE PURPOSE OF THIS ROUTINE IS TO REDUCE THE SCALE FACTOR MATRIX
!     AP TO  A TRANS FORMATION MATRIX  AS, AD, AH
 
!     INPUTS (5)
!         USETD
!         AP     SCALE MATRIX --P SIZE
!         GMD    M- SET TRASNFORMATION MATRIX
!         GOD    0- SET TRASNFORMATION MATRIX
!         PHIDH  H- SET TRASNFORMATION MATRIX
 
!     OUTPUTS(3)
!         AS     SCALE MATRIX --S SET
!         AD     SCALE MATRIX --D SET
!         AH     SCALE MATRIX --H SET
 
!     NOTE  IFLAG1 WILL BE SET  TO -1  IF  AP = AD (N0 M,S,O)
 
 
 
 INTEGER, INTENT(IN)                      :: usetd
 INTEGER, INTENT(IN)                      :: ap
 INTEGER, INTENT(IN OUT)                  :: gmd
 INTEGER, INTENT(IN OUT)                  :: god
 INTEGER, INTENT(IN)                      :: phidh
 INTEGER, INTENT(IN OUT)                  :: as
 INTEGER, INTENT(IN OUT)                  :: ad
 INTEGER, INTENT(IN OUT)                  :: ah
 INTEGER, INTENT(OUT)                     :: iflag1
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER, INTENT(IN)                      :: scr4
 EXTERNAL        andf
 INTEGER :: mcb(7), uset1,anbar,am,an,af,adbar,ao,andf,multi,single,  &
     omit,SIGN,trnsp,prec, um,us,uo
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,ue,up,une,ufe, ud
 COMMON /zzzzzz/ iz(1)
 COMMON /system/ iskip(54),iprec
 COMMON /patx  / nz,n1,n2,n3,uset1
 COMMON /two   / two1(32)
 
 
 anbar = scr2
 am    = scr3
 an    = scr4
 af    = scr2
 adbar = scr3
 ao    = scr4
 
!     SET FLAGS FOR PRESCENCE OF SETS
 
 mcb(1) = usetd
 CALL rdtrl (mcb)
 uset1  = usetd
 multi  = andf(mcb(5),two1(um))
 single = andf(mcb(5),two1(us))
 omit   = andf(mcb(5),two1(uo))
 modal  = 0
 mcb(1) = phidh
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) modal = 1
 nz     = korsz(iz)
 SIGN   = 1
 trnsp  = 1
 prec   = iprec
 
!     REMOVE EACH CONSTRAINT
 
 IF (multi == 0) GO TO 10
 IF (single == 0 .AND. omit == 0) an = ad
 CALL calcv (scr1,up,une,um,iz)
 CALL ssg2a (ap,anbar,am,scr1)
 CALL ssg2b (gmd,am,anbar,an,trnsp,prec,SIGN,scr1)
 GO TO 20
 
!     NO MULTI-POINT CONSTRAINTS
 
 10 an = ap
 
!     REMOVE SINGLES
 
 20 IF (single == 0) GO TO 30
 IF (omit  ==  0) af = ad
 CALL calcv (scr1,une,ufe,us,iz)
 CALL ssg2a (an,af,as,scr1)
 GO TO 40
 
!     NO SINGLES
 
 30 af = an
 40 IF (omit == 0) GO TO 50
 
!     REMOVE OMITS
 
 CALL calcv (scr1,ufe,ud,uo,iz)
 IF (af == ao) ao = scr2
 CALL ssg2a (af,adbar,ao,scr1)
 CALL ssg2b (god,ao,adbar,ad,trnsp,prec,SIGN,scr1)
 GO TO 60
 
!     NO OMITS
 
 50 ad = af
 
!     REMOVE TO H SET
 
 60 IF (modal /= 0) GO TO 70
 CALL ssg2b (phidh,ad,0,ah,trnsp,prec,SIGN,scr1)
 70 iflag1 = multi + single + omit
 IF (iflag1 == 0) iflag1 = -1
 RETURN
END SUBROUTINE trlgb
