SUBROUTINE frlg
     
!     FREQUENCE RESPONSE LOAD GENERATOR
 
!     INPUTS  - CASEXX,USETD,DLT,FRL,GMD,GOD,DIT,PHIDH
 
!     OUTPUTS - PPF,PSF,PDF,FOL,PHF
 
!     4 SCRATHCES
 
 
 EXTERNAL        andf
 INTEGER :: casexx,usetd,dlt,frl,gmd,god,dit,phidh,ppf,psf,  &
     pdf,fol,phf,scr1,scr2,scr3,scr4,mcb(7),andf, single,omit,ifreq(2),itran(2)
 COMMON /two   / itwo(32)
 COMMON /BLANK / modal(2),notrd,iapp(2)
 COMMON /bitpos/ ium,iuo,skp(6),ius
 DATA    casexx, usetd,dlt,frl,gmd,god,dit,phidh      /  &
     101   , 102  ,103,104,105,106,107,108        /
 DATA    ppf   , psf,pdf,fol,phf, scr1,scr2,scr3,scr4 /  &
     201   , 202,203,204,205, 301 ,302 ,303 ,304  /
 DATA    moda  / 4HMODA  /
 DATA    ifreq , itran   /4HFREQ,1H ,4HTRAN,1H /
 
!     DETERMINE USET DATA
 
 mcb(1) = usetd
 CALL rdtrl (mcb)
 lusetd = mcb(2)
 multi  =-1
 IF (andf(mcb(5),itwo(ium)) /= 0) multi =  1
 single =-1
 IF (andf(mcb(5),itwo(ius)) /= 0) single =  1
 omit =-1
 IF (andf(mcb(5),itwo(iuo)) /= 0) omit =  1
 
 iapp(1) = ifreq(1)
 iapp(2) = ifreq(2)
 
!     BUILD LOADS ON P SET
 
!     ORDER IS ALL FREQUENCIES FOR GIVEN LOAD TOGETHER
 
 CALL frlga (dlt,frl,casexx,dit,ppf,lusetd,nfreq,nload,frqset,fol, notrd)
 IF (notrd == -1) GO TO 10
 iapp(1) = itran(1)
 iapp(2) = itran(2)
 10  CONTINUE
 
!     REDUCE LOADS TO D OR H SET
 
 IF (multi < 0 .AND. single < 0 .AND. omit < 0 .AND. modal(1) /= moda) RETURN
 CALL frlgb (ppf,usetd,gmd,god,multi,single,omit,modal,phidh,pdf,  &
     psf,phf,scr1,scr2,scr3,scr4)
 RETURN
END SUBROUTINE frlg
