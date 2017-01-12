@echo off
echo Use this batch file if you do not have a schema validator but
echo do have an XSLT processor. Please ensure that the path to the
echo XSLT Processor is set correctly in ecb-run-nocheck.bat
echo If you have a schema validator as well as an XSLT processor 
echo try demo.bat
echo If you do not have an XSLT processor try demo-raw.bat 
echo.
pause
:menu
echo Choose the sample file to be run
echo  1. American 1
echo  2. American 2
echo  3. Convertible 1
echo  4. Convertible 2
echo.
echo  Q. Quit
choice /c:1234q>nul
if errorlevel 1 set file=American-1.xml
if errorlevel 2 set file=American-2.xml
if errorlevel 3 set file=Convertible-1.xml
if errorlevel 4 set file=Convertible-2.xml
if errorlevel 5 goto exit
command /e:1024 /c ecb-run-nocheck samples/%file%
pause
goto menu
:exit