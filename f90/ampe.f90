SUBROUTINE ampe (phidh,gtka,gkh,scr1,scr2,useta)
     
    !     THE PURPOSE OF THIS ROUTINE IS TO COMPUTE GKH
 
 
    INTEGER, INTENT(IN)                      :: phidh
    INTEGER, INTENT(IN OUT)                  :: gtka
    INTEGER, INTENT(IN OUT)                  :: gkh
    INTEGER, INTENT(IN OUT)                  :: scr1
    INTEGER, INTENT(IN)                      :: scr2
    INTEGER, INTENT(IN)                      :: useta
    INTEGER :: phiah
    COMMON /BLANK / noue
    COMMON /patx  / lc,ns0,ns1,ns2,iuset
    COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,ue,up,une,ufe,  &
        ud,ups,usa,uk,upa
    COMMON /zzzzzz/ z(1)
    COMMON /system/ sysbuf,nout,skip(52),iprec
    COMMON /patx  / xxx,nrow1,nrow2
 
    phiah = phidh
 
    !     DETERMINE IF PHIDH MUST BE MODIFIED
 
    IF (noue == -1) GO TO 20
 
    !     BUILD PARTITIONING VECTORS
 
    iuset = useta
    lc = korsz(z)
    CALL calcv (scr1,ud,ua,ue,z)
 
    !     PERFORM PARTITION
 
    nrow1 = ns0
    nrow2 = ns1
    phiah = scr2
    CALL ssg2a (phidh,phiah,0,scr1)
 
!     COMPUTE GKH
 
20 CONTINUE
   CALL ssg2b (gtka,phiah,0,gkh,1,iprec,1,scr1)

   RETURN
END SUBROUTINE ampe
