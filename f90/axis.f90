SUBROUTINE axis (xa,ya,xb,yb,penx,opt)
     
    !     (XA,YA) = STARTING POINT OF THE AXIS.
    !     (XB,YB) = TERMINAL POINT OF THE AXIS.
    !     PENX    = PEN NUMBER OR LINE DENSITY (DEPENDS ON PLOTTER).
    !     OPT     = -1 TO INITIATE  THE LINE MODE.
    !             = +1 TO TERMINATE THE LINE MODE.
    !             =  0 TO DRAW A LINE.
 
 
    REAL, INTENT(IN OUT)                     :: xa
    REAL, INTENT(IN OUT)                     :: ya
    REAL, INTENT(IN OUT)                     :: xb
    REAL, INTENT(IN OUT)                     :: yb
    INTEGER, INTENT(IN OUT)                  :: penx
    INTEGER, INTENT(IN OUT)                  :: opt
    INTEGER :: pen, ploter
    COMMON /pltdat/ model,ploter,skpplt(18),skpa(6),npens
 
    IF (opt /= 0) GO TO 110
    pen = MAX0(penx,1)
    pen = pen - npens*((pen-1)/npens)
 
110 CALL axis10 (xa,ya,xb,yb,pen,opt)
    RETURN
END SUBROUTINE axis
