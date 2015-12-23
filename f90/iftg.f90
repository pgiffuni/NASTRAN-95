SUBROUTINE iftg(tha,rp,cp)

    CALL ifte2(tha,r,c)
    CALL ifte4(tha,r1,c1)
    rp = 2.*r - r1
    cp = 2.*c - c1

    RETURN
END SUBROUTINE iftg
