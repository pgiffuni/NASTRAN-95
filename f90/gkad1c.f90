SUBROUTINE gkad1c (xmd,xod,xcr1,xcr2,xcr3,xcr4,xcr5,xcr6,xsetd)
     
!     GKAD1C SETS UP TO REDUCE STRUCTURAL MODAL
 
 
 INTEGER, INTENT(IN)                      :: xmd
 INTEGER, INTENT(IN)                      :: xod
 INTEGER, INTENT(IN)                      :: xcr1
 INTEGER, INTENT(IN)                      :: xcr2
 INTEGER, INTENT(IN)                      :: xcr3
 INTEGER, INTENT(IN)                      :: xcr4
 INTEGER, INTENT(IN)                      :: xcr5
 INTEGER, INTENT(IN)                      :: xcr6
 INTEGER, INTENT(IN)                      :: xsetd
 INTEGER :: gmd,god,scr1,scr2,scr3,scr4,scr5,scr6,usetd,  &
     omit,single,check,NAME(2)
!NV  3                MCB(7),T
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,ue,up,une,ufe, ud
 COMMON /BLANK / TYPE(2),app(2),modal(2),g,w3,w4,ik2pp,  &
     im2pp,ib2pp,multi,single,omit,noue
 DATA    NAME  / 4HGKAD,4H1C      /
 
 gmd   = xmd
 god   = xod
 scr1  = xcr1
 scr2  = xcr2
 scr3  = xcr3
 scr4  = xcr4
 scr5  = xcr5
 scr6  = xcr6
 usetd = xsetd
 check = 123456789
 RETURN
 
 
 ENTRY gkad1d (k2pp,k2dd)
!     ========================
 
 IF (check /= 123456789) CALL mesage (-37,0,NAME)
 
!     NAVY'S FIX (MARKED BY CNV) TO FORCE K2NN BE SYMMETRIC IF K2PP IS
!     SYMMETRIC. A PARAMETER OF -6 IS PASSED TO SSG2B TO FLAG THE FORM
!     OF THE MATRIX TO BE SYMMETRIC.
!     ALSO, IN SSG2B, ABOUT LINE 55, ADD FOLLOWING 2 LINES
!           IF (T1 .EQ. -6) T = 1
!           IF (T1 .EQ. -6) FILED(4) = SYMM
 
!     (THE FIX IS NOT ADOPTED HERE. A MORE GENERAL FIX IS ADDED IN SSG2B
!     WHICH SHOULD TAKE CARE OF THE PROBLEM HERE   G.C/UNISYS 3/93)
 
!NV   MCB(1) = K2PP
!NV   CALL RDTRL (MCB)
!NV   T = 1
!NV   IF (MCB(4) .EQ. 6) T = -6
 
 k2ff = k2dd
 IF (multi < 0) GO TO 20
 IF (omit < 0 .AND. single < 0) GO TO 10
 k2nn = scr4
 IF (single < 0) k2nn = k2dd
 GO TO 30
 10 k2nn = k2dd
 GO TO 30
 20 k2nn = k2pp
 30 IF (single >= 0) GO TO 40
 k2ff = k2nn
 40 IF (multi < 0) GO TO 50
 
!     MULTI POINT CONSTRAINTS
 
 CALL upart (usetd,scr1,up,une,um)
 CALL mpart (k2pp,scr2,scr3,scr5,scr4)
 CALL ssg2b (scr4,gmd,scr3,scr1,0,2,1,scr6)
 CALL ssg2b (scr5,gmd,scr2,scr3,0,2,1,scr6)
 
!NV   CALL SSG2B (GMD,SCR1,SCR3,K2NN,T,2,1,SCR6)
 CALL ssg2b (gmd,scr1,scr3,k2nn,1,2,1,scr6)
 
 50 IF (single < 0) GO TO 60
 CALL upart (usetd,scr1,une,ufe,us)
 CALL mpart (k2nn,k2ff,0,0,0)
 60 IF (omit < 0) GO TO 70
 CALL upart (usetd,scr1,ufe,ud,uo)
 CALL mpart (k2ff,scr2,scr3,scr5,scr4)
 CALL ssg2b (scr4,god,scr3,scr1,0,2,1,scr6)
 CALL ssg2b (scr5,god,scr2,scr3,0,2,1,scr6)
 
!NV   CALL SSG2B (GOD,SCR1,SCR3,K2DD,T,2,1,SCR6)
 CALL ssg2b (god,scr1,scr3,k2dd,1,2,1,scr6)
 
 70 RETURN
END SUBROUTINE gkad1c
