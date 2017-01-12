@echo off
echo Use this batch file if you have both a schema validator and
echo an XSLT processor. Please ensure that the paths to these
echo are set correctly in ecb-run.bat
echo If you have only a XSLT processor try demo-nocheck.bat
echo If you do not have that either try demo-raw.bat 
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
command /e:1024 /c ecb-run samples/%file%
pause
goto menu
:exit