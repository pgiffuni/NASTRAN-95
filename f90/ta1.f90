SUBROUTINE ta1
     
!     TA1 CONTROLS THE EXECUTION OF THE TABLE ASSEMBLER.
 
!     DMAP CALL IS
 
!     TA1   ECT,EPT,BGPDT,SIL,GPTT,CSTM,MPT,EQEXIN/EST,GEI,GPECT,
!           ECPT,GPCT,MPTX,PCOMPS,EPTX/V,N,LUSET/V,N,NOSIMP=-1/
!           V,N,NOSUP=-1,1,2/V,N,NOGENEL=-1/V,N,GENEL/V,N,COMPS=1 $
 
 
!     EITHER THE GPECT OR BOTH GPECT AND ECPT, GPCT MAY BE GENERATED.
!     IF NOSUP .EQ. 1, GENERATE GPECT. IF NOSUP .EQ. 2 , GENERATE ALL.
!     IF NOSUP .LT. 0, GENERATE NONE.
 
!   1. TA1 EXECUTES TA1A WHICH BUILDS THE ELEMENT SUMMARY TABLE (EST)
!   2. TA1 EXECUTES TA1B WHICH BUILDS THE ELEMENT CONNECTION AND
!      PROPERTIES TABLE (ECPT) AND THE GRID POINT CONNECTION TABLE(GPCT)
!   3. IF GENERAL ELEMENTS ARE PRESENT, TA1 EXECUTES TA1C WHICH BUILDS
!      THE GENERAL ELEMENT INPUT (GEI).
!   4. IF LAMINATED COMPOSITE ELEMENTS ARE PRESENT, TA1 EXECUTES
!      TA1CPS/D WHICH -
!      (1) CREATES PCOMPS DATA, WHICH INCLUDES THE ECHOING OF
!          INTRINSIC LAYER PROPERTIES, AND
!      (2) CALCULATES OVERALL MATERIAL PROPERTIES.
 
 
 EXTERNAL        andf
 INTEGER :: genl  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,  &
     ecpt  ,gpct  ,scr1  ,scr2  ,two   ,est   ,  &
     andf  ,scr3  ,scr4  ,gei   ,cstm  ,gpect ,  &
     pcomps,eptx  ,comps ,eqexin,genel(2)
 DIMENSION       mcb(7)
 COMMON /BLANK / luset ,nosimp,nosup ,nogenl,genl  ,comps
 COMMON /system/ isystm(54)   ,iprec
 COMMON /ta1com/ nsil  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,cstm  ,  &
     mpt   ,est   ,gei   ,gpect ,ecpt  ,gpct  ,mptx  ,  &
     pcomps,eptx  ,scr1  ,scr2  ,scr3  ,scr4  ,eqexin
 COMMON /two   / two(32)
 DATA    genel / 4301 , 43 /
 
!     INITIALIZE
 
 CALL delset
 ect    = 101
 ept    = 102
 bgpdt  = 103
 sil    = 104
 gptt   = 105
 cstm   = 106
 mpt    = 107
 eqexin = 108
 
 est    = 201
 gei    = 202
 gpect  = 203
 ecpt   = 204
 gpct   = 205
 mptx   = 206
 pcomps = 207
 eptx   = 208
 
 scr1   = 301
 scr2   = 302
 scr3   = 303
 scr4   = 304
 
!     TEST FOR PRESENCE OF GENERAL ELEMENTS
 
 nogenl = -1
 mcb(1) = ect
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) GO TO 100
 j = (genel(2)-1)/16
 k =  genel(2)-16*j
 IF (andf(mcb(j+2),two(k+16)) /= 0) nogenl = 1
 
!     EXECUTE TA1A FOR ALL PROBLEMS
 
 100 CALL ta1a
 
!     EXECUTE TA1CPD/S TO BUILD PCOMPS DATA
 
 IF (nosup == 0) GO TO 300
 IF (comps /= -1) GO TO 200
 IF (iprec == 1) CALL ta1cps
 IF (iprec == 2) CALL ta1cpd
 200 IF (nosup == 1) GO TO 400
 
!     EXECUTE TA1B IF SIMPLE ELEMENTS ARE PRESENT
 
 300 IF (nosimp > 0) CALL ta1b
 IF (nosup == 0) GO TO 500
 
!     CALL TA1H TO GENERATE GPECT
 
 400 IF (nosimp > 0) CALL ta1h
 
!     EXECUTE TA1C IF GENERAL ELEMENTS ARE PRESENT
 
 500 IF (nogenl > 0) CALL ta1c
 genl = -nogenl
 
 RETURN
END SUBROUTINE ta1
