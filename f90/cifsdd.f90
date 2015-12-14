SUBROUTINE cifsdd
     
    !     THIS SUBROUTINE INITIALIZES THE CIFS1P, CIFS2P, CIFS3P,
    !     CIFS4P, AND CIFS5P COMMON BLOCKS.
 
    INTEGER :: b1, bardf2, bardf5, bardf6, bardf7, bardf8
    INTEGER :: g1
 
    LOGICAL :: slot
    LOGICAL :: fphys, fphys1, dmiflg, fthru, fphys2
    LOGICAL :: grdmsg, lh, idfreq
    LOGICAL :: lflsym, ffphys
 
    COMMON /cifs1p/ b1, bardf2, bardf5, bardf6, bardf7, bardf8, km1, slot, idrdl
    COMMON /cifs2p/ fphys, fphys1, km2, dmiflg, ibcds, fthru, fphys2
    COMMON /cifs3p/ grdmsg, la1, l7, km3, l0, g1, lh,  &
                    igdst2, igdst6, igdst7, igdst8, iddsf,  &
                    idfreq, idrad, nvar, ids, jms, kms, lplf
    COMMON /cifs4p/ j(20), km4, lflsym, ffphys
    COMMON /cifs5p/ km5, ic, ip, icont, iaero, ipopt
 
    DATA icc/1HC/, ipp/1HP/
 
    b1 = 1
    bardf2 = 0
    bardf5 = 0
    bardf6 = 0
    bardf7 = 0
    bardf8 = 0
    km1 = 0
    slot = .false.
    idrdl = 0
 
    fphys  = .true.
    fphys1 = .true.
    km2 = 0
    dmiflg = .false.
    ibcds = 0
    fthru  = .false.
    fphys2 = .true.
 
    grdmsg = .false.
    la1 = 0
    l7  = 0
    km3 = 0
    l0 = 1
    g1 = 1
    lh     = .true.
    igdst2 = 0
    igdst6 = 0
    igdst7 = 0
    igdst8 = 0
    iddsf = 0
    idfreq = .true.
    idrad = 0
    nvar  = 0
    ids = 0
    jms = 0
    kms = 0
    lplf = 0
 
    DO  i=3,20
        j(i)  = 0
    END DO
    j(1) = 20
    j(2) = 2
    km4 = 0
    lflsym = .false.
    ffphys = .true.
 
    km5 = 0
    ic = icc
    ip = ipp
    icont = 0
    iaero = 0
    ipopt = 0
 
    RETURN
END SUBROUTINE cifsdd
