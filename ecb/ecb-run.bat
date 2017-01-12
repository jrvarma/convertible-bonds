@echo off

rem Usage and Help
if %1/==/ goto usage
if %1/==--help/ goto usage
goto main
:usage
echo Usage %0 [--help] infile [outfile]
goto end
:main

rem set path to schema validation program
rem I use xsv from the W3 Consortium
rem http://www.w3.org/XML/Schema or
rem http://www.ltg.ed.ac.uk/software/xml/ or
rem ftp://ftp.cogsci.ed.ac.uk/pub/XSV/XSV12.EXE
set validate="C:\Program Files\xsv\xsv"

rem set path to xslt processor program
rem I use sabcmd from the Ginger Alliance:
rem http://www.gingerall.com
set xslate=sabcmd

echo Checking data file (%1) against schema
%validate% -o xsv.out %1 ecbdata.xsd
%xslate% xsv-test.xsl xsv.out | find /c "error" >nul
if not errorlevel 1 goto error

rem there are no errors, so we run ecb
del xsv.out
echo Data file is OK. Running ecb
ecb-run-nocheck %1 %2
goto end

:error
rem there are errors so we display error messages from the validator
%xslate% xsv-test.xsl xsv.out | more
del xsv.out
echo Please correct errors and run again

:end

