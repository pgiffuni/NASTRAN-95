SUBROUTINE smp2
!*****
! THIS MODULE, WHICH IS CALLED ONLY FOR DIFFERENTIAL STIFFNESS, PARTI-
! TIONS KDFF AND THEN COMPUTES KDAA AS FOLLOWS ....
!                       I
!               -D      I     D
!               K       I    K
!   D            AA     I     AO
!  K   =     ----------------------
!   FF                  I
!                D  T   I     D
!              (K  )    I    K
!                AO     I     OO
!                       I
 
!   D     -D         D                 D         T             D
!  K   =  K     +   K    X  G    +   (K    X  G )    +   G  X K    X  G
!   AA     AA        AO      O         AO      O          O    OO      O
 
!*****
 
! DMAP CALL ...
 
!     SMP2     USET,GO,KDFF/KDAA/
 
 INTEGER :: uset               ,GO  &
     ,                  scr1               ,scr2  &
     ,                  uf                 ,ua  &
     ,                  uo                 ,mcb(7)
 
 
 
 COMMON /BLANK/ icom
 
! INPUT FILES
 
 DATA     uset,GO,kdff /101,102,103/
 
! OUTPUT FILE
 
 DATA     kdaa /201/
 
! SCRATCH FILES
 
 DATA     scr1,scr2,kdaab,kdao,kdoo /301,302,303,304,305/
 
! USET BIT POSITIONS
 
 DATA     uf,ua,uo/26,25,30/
 
!  TEST FOR PRESENCE OF KDFF
 
 mcb(1)=kdff
 CALL rdtrl(mcb)
 IF(mcb(1) < 0) RETURN
 
! PARTITION KFF
 
 CALL upart (uset,scr1,uf,ua,uo)
 CALL mpart (kdff,kdaab,kdao,o,kdoo)
 
! COMPUTE KDAA
 
 CALL elim(kdaab,kdao,kdoo,GO,kdaa,scr1,scr2,306)
 RETURN
END SUBROUTINE smp2
