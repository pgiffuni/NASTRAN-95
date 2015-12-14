SUBROUTINE ssg2
     
!     MULTI  = 0  IMPLIES NO MULTI-POINT CONSTRAINTS PN = PG
 
!     SINGLE = 0  IMPLIES NO SINGLE POINT CONSTRAINTS PF = PN
 
!     OMIT   = 0  IMPLIES NO OMITTED POINTS PA = PF
 
!     REACT  = 0  IMPLIES NO FREE BODY PROBLEM PL = PA
 
 
 EXTERNAL        andf
 INTEGER :: uset,gm,kfs,GO,pnbar,pg,pm,po,pa,single,omit,  &
     pvect,ps,d,pl,pr,qr,react,um,un,ug,us,uf,uo,ua,  &
     ul,ur,pf,pabar,pn,pfbar,andf,ys,ia(7),uset1,sr4
 DIMENSION       core(1)
 COMMON /system/ dum54(54),iprec
 COMMON /patx  / lc,n,no,n4,uset1,ibc
 COMMON /BLANK / single
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug
 COMMON /two   / two1(32)
 COMMON /zzzzzz/ core
 
 DATA    uset  , gm    ,kfs   ,GO    ,pnbar ,pm    ,po    ,pn    /  &
     101   , 102   ,104   ,105   ,302   ,303   ,202   ,204   /
 DATA    pfbar , pf    ,pabar ,pa    ,ps    ,d     ,pl    ,pr    /  &
     302   , 204   ,302   ,204   ,203   ,106   ,204   ,302   /
 DATA    qr    , pvect ,ys    ,pg    ,sr4   /  &
     201   , 301   ,103   ,107   ,304   /
 
 
 pnbar = 302
 pn    = 204
 pr    = 302
 pf    = 204
 pa    = 204
 lc    = korsz(core)
 
!     DECIDE IF MULTI,SINGLE,OMIT,REACT ARE 1 OR ZERO
 
 ia(1) = uset
 uset1 = uset
 CALL rdtrl (ia)
 multi = andf(ia(5),two1(um))
 single= andf(ia(5),two1(us))
 omit  = andf(ia(5),two1(uo))
 react = andf(ia(5),two1(ur))
 IF (react <= 0) GO TO 10
 IF (.NOT.(multi > 0 .AND. single == 0 .AND. omit == 0)) GO TO 5
 pnbar = 204
 pn    = 302
 pr    = 303
 GO TO 20
 5 CONTINUE
 pf    = 201
 pa    = 303
 10 IF (multi == 0) THEN
   GO TO    30
 END IF
 
 20 CALL calcv (pvect,ug,un,um,core(1))
 CALL ssg2a (pg,pnbar,pm,pvect)
 CALL ssg2b (gm,pm,pnbar,pn,1,iprec,1,sr4)
 GO TO 40
 
 30 pn = pg
 40 IF (single == 0.0) THEN
   GO TO    70
 END IF
 50 CALL calcv (pvect,un,uf,us,core(1))
 CALL ssg2a (pn,pfbar,ps,pvect)
 CALL ssg2b (kfs,ys,pfbar,pf,0,iprec,0,sr4)
 GO TO 80
 70 pf = pn
 80 IF (omit == 0.0) THEN
   GO TO   100
 END IF
 
 90 CALL calcv (pvect,uf,ua,uo,core(1))
 CALL ssg2a (pf,pabar,po,pvect)
 CALL ssg2b (GO,po,pabar,pa,1,iprec,1,sr4)
 GO TO 110
 100 pa = pf
 110 IF (react == 0.0) THEN
   GO TO   130
 END IF
 
 120 CALL calcv (pvect,ua,ul,ur,core(1))
 CALL ssg2a (pa,pl,pr,pvect)
 CALL ssg2b (d,pl,pr,qr,1,iprec,-1,sr4)
 130 RETURN
END SUBROUTINE ssg2
