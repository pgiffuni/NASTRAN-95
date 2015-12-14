SUBROUTINE frlgb (pp,usetd,gmd,god,multi,single,omit,modal,phidh,  &
        pd,ps,ph,scr1,scr2,scr3,scr4)
     
!     THIS ROUTINE REDUCES LOADS FROM P SET TO D SET
 
!     ENTRY POINT - FRRD1B
!                   ======
 
 
 INTEGER, INTENT(IN)                      :: pp
 INTEGER, INTENT(IN)                      :: usetd
 INTEGER, INTENT(IN OUT)                  :: gmd
 INTEGER, INTENT(IN OUT)                  :: god
 INTEGER, INTENT(IN OUT)                  :: multi
 INTEGER, INTENT(IN OUT)                  :: single
 INTEGER, INTENT(IN)                      :: omit
 INTEGER, INTENT(IN OUT)                  :: modal
 INTEGER, INTENT(IN OUT)                  :: phidh
 INTEGER, INTENT(IN OUT)                  :: pd
 INTEGER, INTENT(IN OUT)                  :: ps
 INTEGER, INTENT(IN)                      :: ph
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER, INTENT(IN)                      :: scr4
 INTEGER :: po,  uset,pn,pnbar,pm,pf,pdbar
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,ue,up,une,ufe, ud
 COMMON /patx  / nz,n1,n2,n3,uset
 COMMON /zzzzzz/ core(1)
 DATA    moda  / 4HMODA /
 
 GO TO 5
 
 
 ENTRY frrd1b (pp,usetd,gmd,god,multi,single,omit,modal,phidh,  &
     pd,ps,ph,scr1,scr2,scr3,scr4)
!     =============================================================
 
!     SET UP INITIAL VALUES
 
 5 nz    = korsz(core)
 uset  = usetd
 pnbar = scr2
 pm    = scr3
 pn    = scr4
 pf    = scr2
 pdbar = scr3
 po    = ph
 
!     REMOVE EACH TYPE OF CONSTRAINT
 
 IF (multi < 0) GO TO 10
 
!     REMOVE MULTIPOINT CONSTRAINTS
 
 IF (single < 0 .AND. omit < 0) pn = pd
 CALL calcv (scr1,up,une,um,core(1))
 CALL ssg2a (pp,pnbar,pm,scr1)
 CALL ssg2b (gmd,pm,pnbar,pn,1,1,1,scr1)
 GO TO 20
 
!     NO M-S
 
 10 pn = pp
 20 IF (single < 0) GO TO 30
 
!     REMOVE SINGLE POINT CONSTRAINTS
 
 IF (omit < 0) pf = pd
 CALL calcv (scr1,une,ufe,us,core(1))
 CALL ssg2a (pn,pf,ps,scr1)
 GO TO 40
 
!     NO SINGLE POINT CONSTRAINTS
 
 30 pf = pn
 40 IF (omit < 0) GO TO 50
 
!     REMOVE OMITS
 
 CALL calcv (scr1,ufe,ud,uo,core(1))
 CALL ssg2a (pf,pdbar,po,scr1)
 CALL ssg2b (god,po,pdbar,pd,1,1,1,scr1)
 GO TO 60
 50 pd = pf
 60 IF (modal /= moda) GO TO 70
 
!     TRANSFORM TO MODAL COORDINATES
 
 CALL ssg2b (phidh,pd,0,ph,1,1,1,scr1)
 70 RETURN
END SUBROUTINE frlgb
