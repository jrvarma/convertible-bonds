@echo off

rem Usage and Help
if %1/==/ goto usage
if %1/==--help/ goto usage
goto main
:usage
echo Usage %0 [--help] infile [outfile]
goto end
:main

rem set path to xslt processor program
rem I use sabcmd from the Ginger Alliance:
rem http://www.gingerall.com
set xslate=sabcmd

if %2/==/ goto screen
ecb %1 | %xslate% xmlout.xsl >%2
goto end
:screen
ecb %1 | %xslate% xmlout.xsl | more
:end