   for %%f in (*.o) do (
    copy %%~nf.f90 ..\%%~nf.f90
    del %%~nf.f90
)