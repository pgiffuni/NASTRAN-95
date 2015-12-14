SUBROUTINE eandm (itype,ido,nextz,lcore,nbdys,all,nelout)
     
!     COMPUTES ADDITIONAL LOAD IN ZIEKIEWICZ PAPER DUE TO SPECIFIED
!     MAGNETIC FIELD OR CURRENT LOOP
 
!     ITYPE = 20  SPCFLD
!     ITYPE = 21  CEMLOOP
!     ITYPE = 22  GEMLOOP
!     ITYPE = 23  MDIPOLE
!     ITYPE = 24  REMFLUX
!     IDO   = NUMBER OF CARDS OF PRESENT TYPE
!     NEXTZ = NEXT AVAILABLE POINTER INTO OPEN CORE
!     LAST AVAILABLE POINTER INTO OPEN CORE
!     *** ALL CEMLOOP, SPCFLD, GEMLOOP, AND MDIPOLE CARDS WERE COMBINED
!     INTO ONE SPCFLD-TYPE CARD WITH 3*NROWSP WORDS-HCX, HCY, HCZ AT
!     EACH POINT AND IS INDICATED BY ITYPE =-20. THESE 3*NROWSP WORDS
!     ARE WRITTEN TO HCFLDS FOR LATER USE. THE OTHER CARDS ARE STILL ON
!     SLT FOR USE IN THE NUMERICAL INTEGRATION.
 
 
 INTEGER, INTENT(IN)                      :: itype
 INTEGER, INTENT(IN)                      :: ido
 INTEGER, INTENT(IN)                      :: nextz
 INTEGER, INTENT(IN)                      :: lcore
 INTEGER, INTENT(IN OUT)                  :: nbdys
 REAL, INTENT(IN OUT)                     :: all
 INTEGER, INTENT(IN OUT)                  :: nelout
 LOGICAL :: done
 INTEGER :: FILE,buf1,sysbuf,est,slt,eltype,estwds,outpt,scr6,  &
     hcflds,mcb(7),remfls,mcb1(7),mcb2(7)
 DIMENSION       iz(1),nam(2),necpt(1),NAME(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / nrowsp
 COMMON /system/ ksystm(64)
 COMMON /emecpt/ ecpt(200)
 COMMON /packx / ita,itb,ii,jj,incur
 COMMON /zblpkx/ a(4),irow
 COMMON /zzzzzz/ z(1)
 COMMON /gpta1 / nelems,last,incr,NE(1)
 EQUIVALENCE     (ksystm(1),sysbuf),(ksystm(2),outpt),(z(1),iz(1)),  &
     (ksystm(56),ithrml),(ecpt(1),necpt(1))
 DATA    nam   / 4HEAND,4HM              /
 DATA    est   , slt, hcflds,remfls,scr6 / 105   , 205, 304,   305,   306  /
 DATA    mcb   / 304, 0, 0, 2, 1, 0, 0   /
 DATA    mcb1  / 305, 0, 0, 2, 1, 0, 0   /
 DATA    done  / .false.                 /
 
!     CHECK IF THERMAL FORMULATION
 
 IF (ithrml == 0) RETURN
 
!     READ A CARD TYPE FROM SLT. TYPE=-20 IS THE COMBINATION HC FOR ALL
!     CARD TYPES EXCEPT REMFLUX AND SIGNIFIES END OF A SUBCASE. READ AND
!     PACK IT. SAME FOR TYPE = 24(REMFLUX). FOR TYPES 20-24, COMPUTE
!     LOAD = INTEGRAL(GRAD(NI)*MU*HC)*D(VOL). WE WILL USE NUMERICAL
!     INTEGRATION FOR ITYPE=24-REMFLUX-ONLY ONE CARD GIVING FLUX IN
!     EACH ELEMENT COMPUTE INTEGRAL(GRAD NI*BR)*D(VOL)
 
 buf1  = lcore - sysbuf + 1
 icore = buf1  - 1
 IF (nextz > buf1) GO TO 420
 
 IF (itype /= -20 .AND. itype /= 24) GO TO 40
 IF (ido /= 1) GO TO 300
 
!     END OF SUBCASE-WRAP UP SCR6 AND CALL HCCOM TO COMBINE CENTROID
!     RESULTS IF NOT REMFLUX, CALL HCCOM NOW.IF REMFLUX, WAIT UNTIL
!     KCOUNT IS SET.
 
 CALL CLOSE (scr6,1)
 IF (itype == -20) CALL hccom (itype,lcore,icore,nextz,kcount)
 jj = nrowsp
 
!     ITYPE=-20 OR +24--END OF SUBCASE. IF +24, WRITE ZEROS TO HCFLDS
!     AND HCCENS AND REMFLUX VECTOR TO REMFLS. THEN CONTINUE ON TO
!     COMPUTE LOADS. IF ITYPE=-20, WRITE ZRROS TO REMFLS, GRID POINT
!     HC VALUES TO HCFLDS AND CENTROIDAL VALUES TO HCCENS (ALREADY DONE
!     IN HCCOM). FOR ITYPE=-20, NO FURTHER PROCESSING IS DONE SINCE
!     LOADS HAVE ALREADY BEEN COMPUTED.
 
 ita = 1
 itb = 1
 ii  = 1
 jj  = 3*nrowsp
 incur  = 1
 mcb(3) = jj
 mcb2(1)= est
 CALL rdtrl (mcb2)
 nel = mcb2(2)
 jj1 = 3*nel
 mcb1(3)= jj1
 
!     READ IN THE ONE SPCFLD OR REMFLUX-TYPE CARD
 
 nwords = 3*nrowsp
 IF (itype /= 24) GO TO 10
 nwords = 3*nel
 jj  = nwords
 jj1 = 3*nrowsp
 10 istart = nextz
 IF (nextz+nwords-1 > icore) GO TO 420
 CALL fread (slt,z(nextz),nwords,0)
 
!     CREATE A ZERO VECTOR FOR EITHER REMFLS OR HCFLDS(WHICHEVER IS NOT
!     USED IN THIS SET ID-REMEMBER THAT SPCFLD AND REMFLUX CANNOT HAVE
!     THE SAME SET ID
 
!     PACK THE 3*NROWSP HC FIELD OUT TO BE USED LATER BY EMFLD. HCFLDS
!     WILL CONTAIN ONE COLUMN PER CASE CONTROL SIMPLE SELECTION
!     (SIMPLE LOADS ON LOAD CARDS ARE INCLUDED). COMBIN WILL COMBINE
!     FOR LOAD BULK DATA CARDS AND PUT LOADS IN ORDER OF SELECTION ONTO
!     HCFL (SAME HOLDS FOR 3*NEL WORDS OF REMFLS)
 
 IF (itype == 24) GO TO 20
 CALL pack (z(nextz),hcflds,mcb)
 CALL wrttrl (mcb)
 jj = jj1
 CALL bldpk  (1,1,remfls,0,0)
 CALL bldpkn (remfls,0,mcb1)
 CALL wrttrl (mcb1)
 GO TO 30
 20 CALL pack(z (nextz),remfls,mcb1)
 CALL wrttrl (mcb1)
 jj = jj1
 CALL bldpk  (1,1,hcflds,0,0)
 CALL bldpkn (hcflds,0,mcb)
 CALL wrttrl (mcb)
 
!     RETURN JJ TO VALUE EXPECTED IN EXTERN
 
 30 jj = nrowsp
 IF (itype == -20) RETURN
 
!     GET INFO FROM EST
 
 40 FILE = est
 CALL gopen (est,z(buf1),0)
 ncount = 0
 IF (.NOT.done) kcount = 0
 
!     READ IN ALL CARDS OF THIS TYPE FOR THIS SUBCASE. NO NEED TO READ
!     IN THE ONE REMFLUX CARD SINCE IT WAS DONE ABOVE.
 
 ijk = itype - 19
 SELECT CASE ( ijk )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 60
   CASE (    3)
     GO TO 70
   CASE (    4)
     GO TO 80
   CASE (    5)
     GO TO 100
 END SELECT
 50 iwords = 3*nrowsp
 IF (ido /= 1) GO TO 300
 GO TO 90
 60 iwords = 12
 GO TO 90
 70 iwords = 48
 GO TO 90
 80 iwords = 9
 90 nwords = iwords*ido
 IF (nextz+nwords-1 > icore) GO TO 420
 CALL fread (slt,z(nextz),nwords,0)
 istart = nextz
 
 100 CALL READ (*260,*410,est,eltype,1,0,iflag)
 idx = (eltype-1)*incr
 estwds = NE(idx+12)
 ngrids = NE(idx+10)
 NAME(1)= NE(idx+1)
 NAME(2)= NE(idx+2)
 
 120 CALL READ (*400,*100,est,ecpt,estwds,0,iflag)
 ncount = ncount + 1
 IF (done) GO TO 130
 IF (eltype < 65) kcount = kcount + 3
 IF (eltype == 65) kcount = kcount + 27
 IF (eltype == 66 .OR. eltype == 67) kcount = kcount + 63
 IF (eltype == 80) kcount = kcount + 27
 
 130 IF (eltype > 80) GO TO 230
 GO TO (200,230,200,230,230,210,230,230,210,200,  &
     230,230,230,230,230,210,210,210,210,230,  &
     230,230,230,230,230,230,230,230,230,230,  &
     230,230,230,200,230,210,210,230,220,220,  &
     220,220,230,230,230,230,230,230,230,230,  &
     230,230,230,230,230,230,230,230,230,230,  &
     230,230,230,230,220,220,220,230,230,230,  &
     230,230,230,230,230,230,230,230,230,210), eltype
 
 200 CALL em1d (eltype,istart,itype,ncount,ido,iwords,nbdys,all,nelout)
 GO TO 120
 210 CALL em2d (eltype,istart,itype,ncount,ido,iwords,nbdys,all,nelout)
 GO TO 120
 220 CALL em3d (eltype,istart,itype,ncount,ido,iwords,nbdys,all,nelout)
 GO TO 120
 
 230 WRITE  (outpt,240) ufm,NAME
 240 FORMAT (a23,', ELEMENT TYPE ',2A4,' WAS USED IN AN E AND M ',  &
     'PROBLEM. NOT A LEGAL TYPE')
 250 CALL mesage (-61,0,0)
 
!     DONE
 
 260 CALL CLOSE (est,1)
 IF (itype == 24) GO TO 270
 CALL WRITE (scr6,0,0,1)
 GO TO 280
 270 CALL hccom (itype,lcore,icore,nextz,kcount)
 jj = nrowsp
 280 done =.true.
 RETURN
 
!     FATAL ERROR MESSAGES
 
 300 WRITE  (outpt,310) ufm,nam
 310 FORMAT (a23,', LOGIC ERROR IN SUBROUTINE ',2A4,  &
     '. ONLY ONE SPCFLD OR REMFLUX SHOULD NOW EXIST')
 GO TO 250
 
 400 n = -2
 GO TO 430
 410 n = -3
 GO TO 430
 420 n = -8
 FILE = 0
 430 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE eandm
