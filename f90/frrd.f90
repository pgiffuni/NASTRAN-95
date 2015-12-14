SUBROUTINE frrd
     
!     FREQUENCY AND RANDOM RESPONSE MODULE
 
!     INPUTS CASECC,USETD,DLT,FRL,GMD,GOD,KDD,BDD,MDD,PHIDH,DIT
 
!     OUTPUTS UDV,PS,PD,PP
 
!     8 SCRATCHES
 
 INTEGER :: single,omit,casecc,usetd,dlt,frl,gmd,god,bdd,  &
     phidh,dit,scr1,scr2,scr3,scr4,scr5,scr6,scr7,scr8,  &
     udv,ps,pd,pp,pdd,fol,NAME(2),mcb(7)
 COMMON /BLANK / app(2),modal(2),lusetd,multi,single,omit, noncup,frqset
 COMMON /frrdst/ ovf(150),icnt,ifrst,itl(3),idit,ifrd,k4dd
 COMMON /cdcmpx/ dum32(32),ib,ibbar
 DATA    casecc, usetd,dlt,frl,gmd,god,kdd,bdd,mdd,phidh,dit /  &
     101   , 102,  103,104,105,106,107,108,109,110,  111 /
 DATA    udv   , ps, pd, pp ,pdd,fol / 201   , 202,203,204,203,205 /
 DATA    scr1  , scr2,scr3,scr4,scr5,scr6,scr7,scr8 /  &
     301   , 302, 303, 304, 305, 306, 307, 308  /
 DATA    moda  / 4HMODA /,      NAME /4HFRRD,4H     /
 
 pdd  = 203
 scr6 = 306
 ib   = 0
 
!     BUILD LOADS ON P SET ORDER IS ALL FREQ. FOR LOAD TOGETHER
!     FRRD1A IS AN ENTRY POINT IN FRLGA
 
 CALL frrd1a (dlt,frl,casecc,dit,pp,lusetd,nfreq,nload,frqset,fol, notrd)
 IF (multi < 0 .AND. single < 0 .AND. omit < 0 .AND.  &
     modal(1) /= moda) GO TO 60
 
!     REDUCE LOADS TO D OR H SET
!     FRRD1B IS AN ENTRY POINT IN FRLGB
 
 CALL frrd1b (pp,usetd,gmd,god,multi,single,omit,modal(1),phidh,pd,  &
     ps,scr5,scr1,scr2,scr3,scr4)
 
!     SCR5 HAS PH IF MODAL FORMULATION
 
 IF (modal(1) == moda) pdd = scr5
 
!     SOLVE PROBLEM FOR EACH FREQUENCY
 
 IF (noncup < 0 .AND. modal(1) == moda) GO TO 50
 10 IF (nfreq == 1  .OR.  nload == 1) scr6 = udv
 DO  i = 1,nfreq
   CALL klock (itime1)
   
!     FORM AND DECOMPOSE MATRICES
!     IF MATRIX IS SINGULAR, IGOOD IS SET TO 1 IN FRRD1C. ZERO OTHERWISE
   
   CALL frrd1c (frl,frqset,mdd,bdd,kdd,i,scr1,scr2,scr3,scr4,scr8, scr7,igood)
   
!     ULL IS ON SCR1 -- LLL IS IN SCR2
   
!     SOLVE FOR PD LOADS STACK ON SCR6
   
   CALL frrd1d (pdd,scr1,scr2,scr3,scr4,scr6,i,nload,igood,nfreq)
   CALL klock  (itime2)
   CALL tmtogo (itleft)
   IF (2*(itime2-itime1) > itleft .AND. i /= nfreq) GO TO 70
 END DO
 
 i = nfreq
 30 IF (nfreq == 1 .OR. nload == 1) GO TO 40
 
!     RESORT SOLUTION VECTORS INTO SAME ORDER AS LOADS
!     FRRD1E IS AN ENTRY POINT IN FRRD1D
 
 CALL frrd1e (scr6,udv,nload,i)
 40 RETURN
 
!     UNCOUPLED MODAL
 
 50 CALL frrd1f (mdd,bdd,kdd,frl,frqset,nload,nfreq,pdd,udv)
 GO TO 40
 60 pdd = pp
 GO TO 10
 
!     INSUFFICIENT TIME TO COMPLETE ANOTHER LOOP
 
 70 CALL mesage (45,nfreq-i,NAME)
 mcb(1) = scr6
 CALL rdtrl (mcb(1))
 ndone  = mcb(2)
 mcb(1) = pp
 CALL rdtrl (mcb(1))
 mcb(2) = ndone
 CALL wrttrl (mcb(1))
 IF (single < 0) GO TO 80
 mcb(1) = ps
 CALL rdtrl (mcb(1))
 mcb(2) = ndone
 CALL wrttrl (mcb(1))
 80 mcb(1) = pd
 CALL rdtrl( mcb(1))
 mcb(2) = ndone
 CALL wrttrl (mcb(1))
 GO TO 30
END SUBROUTINE frrd
