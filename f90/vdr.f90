SUBROUTINE vdr
     
!     VDR IS THE CONTROL PROGRAM FOR THE VECTOR DATA RECOVERY MODULE
 
!                                                          OPHID
!                                  PHID                    OUDVC1
!                                  UDVF CLAMA              OUDV1
!             CASECC  EQDYN  USETD UDVT PPF          PHLD  OPHIH  OPNL1
!     VDR     CASEXX,HEQDYN,HUSETD,PHIH,TOL   ,XYCBD,PNLH /OUHVC1,HOPNL1
!                                  UHVT HTOL         HPNLD OUHV1
!                                  HUDVT                   HOUVD1
 
!                   TRANRESP     DIRECT
!              /C,N,FREQRESP/C,N,MODAL /V,N,SORT2/V,N,OUTPUT/V,N,SDR2
!                   CEIGN
 
!              /V,N,FMODE  $      PROGRAMMER'S MANUAL PP. 4.60-1 TRHU -7
 
 
 INTEGER :: pnl   ,outfle,opnl1 ,app   ,trn   ,vdrreq,sort2 ,  &
     output,sdr2  ,sscell,buf   ,casecc
 DIMENSION       nam(2),buf(50)      ,masks(6)     ,mcb(7),cei(2),  &
     frq(2),trn(2),modal(2)     ,DIRECT(2)
 COMMON /vdrcom/ vdrcom,idisp ,ivel  ,iacc  ,ispcf ,iloads,istr  ,  &
     ielf  ,iadisp,iavel ,iaacc ,ipnl  ,ittl  ,ilsym ,  &
     ifrout,idload,casecc,eqdyn ,usetd ,infile,oeigs ,  &
     pp    ,xycdb ,pnl   ,outfle,opnl1 ,scr1  ,scr2  ,  &
     buf1  ,buf2  ,buf3  ,nam   ,buf   ,masks ,cei   ,  &
     frq   ,trn   ,DIRECT,xset0 ,vdrreq,modal
 COMMON /BLANK / app(2),FORM(2),sort2,output,sdr2  ,imode
 COMMON /system/ dumi(68),sscell
 
!     EXECUTE THE PHASES OF VDR.
 
 DO  i = 1,50
   buf(i) = 0
 END DO
 casecc = 101
 output = -1
 sort2  = -1
 CALL vdra
 IF (sscell /= 0) sdr2 = 1
 IF (vdrreq == 0) RETURN
 mcb(1) = infile
 CALL rdtrl (mcb)
 IF (mcb(1) /= infile) GO TO 20
 CALL vdrb (infile,outfle,iadisp)
 20 IF (app(1) /= trn(1)) RETURN
 mcb(1) = pnl
 CALL rdtrl (mcb)
 IF (mcb(1) /= pnl) RETURN
 CALL vdrb (pnl,opnl1,ipnl)
 RETURN
END SUBROUTINE vdr
