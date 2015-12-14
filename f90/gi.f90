SUBROUTINE gi
     
 EXTERNAL        andf
 LOGICAL :: multi,single,omit
 INTEGER :: andf,two1,ia(7)
 INTEGER :: um,uo,ur,usg,usb,ul,ua,uf,us,un,ug
 INTEGER :: spline,useta,cstm,bagpdt,sila,ecta,gm,GO,scr1,  &
     scr2,scr3,scr4,scr5,ksize,gsize,gtka
 COMMON /gicom / spline,useta,cstm,bagpdt,sila,ecta,gm,GO,gtka,  &
     ksize,gsize,scr1,scr2,scr3,scr4,scr5
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug
 COMMON /two   / two1(32)
 COMMON /BLANK / nk,ng
 DATA    single/ .true./, multi /.true./, omit /.true./
 DATA       ia / 7*0  /
 
 spline = 101
 useta  = 102
 cstm   = 103
 bagpdt = 104
 sila   = 105
 ecta   = 106
 gm     = 107
 GO     = 108
 gtka   = 201
 ksize  = nk
 gsize  = ng
 IF (gsize > 0) GO TO 5
 ia(1)  = sila
 CALL rdtrl (ia)
 gsize  = ia(3)
 5 CONTINUE
 scr1   = 301
 scr2   = 302
 scr3   = 303
 scr4   = 304
 scr5   = 305
 CALL giggks
 ia(1)  = useta
 CALL rdtrl (ia)
 IF (andf(ia(5),two1(um)) == 0) multi  = .false.
 IF (andf(ia(5),two1(us)) == 0) single = .false.
 IF (andf(ia(5),two1(uo)) == 0) omit   = .false.
 IF (multi .OR. single .OR. omit) GO TO 10
 scr2 = gtka
 10 CALL gigtkg
 CALL gipsst
 IF (multi .OR. single .OR. omit) GO TO 20
 GO TO 30
 20 CALL gigtka (multi,single,omit)
 30 RETURN
END SUBROUTINE gi
