SUBROUTINE fndplt (ploter, model, pmodel)
     
!     PLOTER = PLOTTER INDEX.
!     MODEL = MODEL INDEX.
!     PMODEL = PLOTTER MODEL ID.
 
!...  DATA FOR PLOTTER + MODEL RECOGNITION.
 
 
 INTEGER, INTENT(OUT)                     :: ploter
 INTEGER, INTENT(OUT)                     :: model
 INTEGER, INTENT(IN OUT)                  :: pmodel(2)
 INTEGER :: pltter(2,6), pltmdl(2,6)
 
 DATA pltmdl / &
!       NASTRAN GENERAL PURPOSE PLOTTER  &
 1HM,1, 1HT,1, 1HD,1, 1HM,0, 1HT,0, 1HD,0/
 DATA pltter / 1,-1,  2,-2,  2,-3,  &
     1,+1,  2,+2,  2,+3/
 
!     FIND THE MODEL ID.
 
 n = -1
 n1 = pmodel(2)
 DO  i = 1, 6
   IF (pmodel(1) /= pltmdl(1,i))  CYCLE
   IF (n <= 0)  n=i
   IF (n1 == pltmdl(2,i)) n = i
 END DO
 
!     SETUP THE PLOTTER + MODEL INDICES.
 
 i2 = pmodel(2)
 IF (n < 0) i2 = 0
 n = IABS (n)
 DO  i = 1,2
   IF (pltmdl(i,n) /= 0)  pmodel(i)=pltmdl(i,n)
 END DO
 ploter = pltter(1,n)
 model  = pltter(2,n)
 
 RETURN
END SUBROUTINE fndplt
