SUBROUTINE gkam1b(usetd,scr1,scr2,phidh,phidh1,modes,core,  &
        lhset,noue,scr3)
     
 
 INTEGER, INTENT(IN)                      :: usetd
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: phidh
 INTEGER, INTENT(IN)                      :: phidh1
 INTEGER, INTENT(IN)                      :: modes
 INTEGER, INTENT(IN OUT)                  :: core(1)
 INTEGER, INTENT(OUT)                     :: lhset
 INTEGER, INTENT(IN)                      :: noue
 INTEGER, INTENT(IN OUT)                  :: scr3
 INTEGER :: uset
 INTEGER :: mcb(7),sysbuf
 
 
 
 COMMON  /system/ sysbuf
 COMMON  /patx/lc,n1,n2,n3,uset
 COMMON  /bitpos/um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,ue,up,une,ufe, ud
 COMMON  /zblpkx/a(4),ii
 COMMON  /parmeg/ ia(7),ia11(7),ia12(7),ia21(7),ia22(7),lcore,irule
 
! ----------------------------------------------------------------------
 
 lhset =lhset  + noue
 
!     BUILD  MERGE  VECTOR
 
 uset = usetd
 lc  = korsz(core(modes))
 lcore = lc
 CALL  calcv(scr1,ud,ua,ue,core(modes))
 
!     BUILD  EXE  IDENTY   MATRIX
 
 nz = lc-sysbuf
 CALL gopen(scr2,core(nz+1),1)
 CALL makmcb(mcb,scr2,noue,6,1)
 a(1) = 1.0
 DO  i=1,noue
   CALL bldpk(1,1,scr2,0,0)
   ii = i
   CALL zblpki
   CALL bldpkn(scr2,0,mcb)
 END DO
 CALL  CLOSE(scr2, 1)
 CALL wrttrl ( mcb )
 
!     SET  UP  FOR  MERGE
 
 irule =0
 ia22(1)= scr2
 CALL rdtrl(ia22)
 ia(1) = phidh
 ia(2)= lhset
 ia(3) =  n1+n2+n3
 ia(4) =  2
 ia(5) = 1
 ia21 (1) = 0
 ia12 (1) = 0
 ia11 (1) = phidh1
 CALL makmcb(core(modes),scr3,lhset,2,1)
 CALL rdtrl (ia11)
 
!     BUILD  VECTOR IN CORE
 
 CALL gopen(scr3,core(nz+1),1)
 CALL bldpk(1,1,scr3,0,0)
 ii = modes -1
 DO  i=1,noue
   ii = ii+1
   CALL zblpki
 END DO
 CALL bldpkn( scr3, 0, core(modes) )
 CALL CLOSE(scr3,1)
 CALL wrttrl(core(modes))
 CALL merge(scr3,scr1,core(modes))
 CALL wrttrl( ia )
 RETURN
END SUBROUTINE gkam1b
