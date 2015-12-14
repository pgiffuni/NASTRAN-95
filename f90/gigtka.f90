SUBROUTINE gigtka(multi,single,omit)
     
 
 LOGICAL, INTENT(IN OUT)                  :: multi
 LOGICAL, INTENT(IN OUT)                  :: single
 LOGICAL, INTENT(IN OUT)                  :: omit
 
 INTEGER :: core,uset1,gm,GO,gka,gkg,gknb,gkm,scr1,gkab,gkf, gks,gko,useta,gkn
 INTEGER :: um,uo,ur,usg,usb,ul,ua,uf,us,un,ug
 
 COMMON /patx/  lc,n,no,ny,uset1,ibc
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug
 COMMON /zzzzzz/ core(1)
 COMMON/gicom/ spline,useta,cstm,bgpt,sila,eqaero,gm,GO,gka,  &
     ksize,gsize,scr1,gkg,gknb,gkm,gkab
 
!-----------------------------------------------------------------------
 
 lc = korsz(core)
 gkf = gknb
 gks = gkm
 gko = gks
 uset1 = useta
 
!     REDUCE TO N SET IF MULTI POINT CONSTRAINTS
 
 gkn = gkg
 IF(.NOT.multi) GO TO 20
 IF(.NOT.single.AND..NOT.omit) gkn = gka
 CALL calcv(scr1,ug,un,um,core)
 CALL ssg2a(gkg,gknb,gkm,scr1)
 CALL ssg2b(gm,gkm,gknb,gkn,1,1,1,scr1)
 
!     PARTITION INTO F SET IF SINGLE POINT CONSTRAINTS
 
 20 IF(.NOT.single) GO TO 30
 IF(.NOT.omit) gkf = gka
 CALL calcv(scr1,un,uf,us,core)
 CALL ssg2a(gkn,gkf,  0,scr1)
 GO TO 40
 
!     REDUCE TO A SET IF OMITS
 
 30 gkf = gkn
 40 IF(.NOT.omit) GO TO 50
 CALL calcv(scr1,uf,ua,uo,core)
 CALL ssg2a(gkf,gkab,gko,scr1)
 CALL ssg2b(GO,gko,gkab,gka,1,1,1,scr1)
 50 RETURN
END SUBROUTINE gigtka
