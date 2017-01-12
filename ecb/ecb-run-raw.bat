@echo off

rem Usage and Help
if %1/==/ goto usage
if %1/==--help/ goto usage
goto main
:usage
echo Usage %0 [--help] infile [outfile]
goto end
:main

if %2/==/ goto screen
ecb %1 >%2
goto end
:screen
ecb %1 | more
:end