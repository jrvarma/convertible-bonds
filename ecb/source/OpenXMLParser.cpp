/*
    OpenXMLFile.cpp 
    Additional functions defined by J. R. Varma for AdvXMLParser
    Part of ecb, a program for Valuation of Convertible Bonds
    Copyright (C) 2001  Prof. Jayanth R. Varma, jrvarma@iimahd.ernet.in,
    Indian Institute of Management, Ahmedabad 380 015, INDIA

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (see file COPYING); if not, write to the 
    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
    Boston, MA  02111-1307  USA
*/

// ===================================================================
// AdvXMLParser
// -------------------------------------------------------------------
// Version:     1.1.4
// Date:        November 19, 2000
// OS:          Windows 2000 SP1
// Compiler:    Microsoft Visual C++ 6.0 SP4
// STL:         STLport 4.0
// -------------------------------------------------------------------
// Sebastien Andrivet grants Licensee a non-exclusive, non-transferable, 
// royalty-free license to use AdvXMLParser and its documentation (the
// 'Software') without fee.
// 
// By downloading, using, or copying the Software or any portion thereof, 
// Licensee agrees to abide by the intellectual property laws and all other 
// applicable laws, and to all of the terms and conditions of this Agreement.
// 
// Licensee shall maintain the following copyright and permission notices on 
// the Software sources and its documentation unchanged :
// 
// Copyright © 1999,2000 Sebastien Andrivet
// 
// THE SOFTWARE IS PROVIDED ''AS IS'', WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL SEBASTIEN ANDRIVET BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// 
// Permission to use or copy this software for any purpose is hereby granted
// without fee, provided the above notices are retained on all copies. 
// Permission to modify the code and to distribute modified code is granted,
// provided the above notices are retained, and a notice that the code was 
// modified is included with the above copyright notice. 
// 
// The Licensee may distribute binaries compiled with the Software (whether 
// original or modified) without any royalties or restrictions.
// 
// The Licensee may distribute original or modified the Software sources,
// provided that the conditions indicated in the above permission notice are met.
// 
// Except as contained in this notice, the name of Sebastien Andrivet
// shall not be used in advertising or otherwise to promote the sale, 
// use or other dealings in this Software without prior written 
// authorization from Sebastien Andrivet.
// ===================================================================

#include <stdio.h>
using namespace std;

#include "OpenXMLParser.h"
using namespace AdvXMLParser;

// This is a modification of a similar function defined in
// the sample files of AdvXMLParser

char* OpenXmlFile(FILE* pFile, long& nSize)
// ---------------------------------------------
// Load a XML file a return it in a buffer
//
// pFile:  	[in]  input XML file
// nSize:       [out] size of the data returned
// Return:      pointer to XML file content
//
// Note:        The caller must destroy the
//              returned buffer with delete[]
// ---------------------------------------------
{
    // Compute the size of the file
    fseek(pFile, 0, SEEK_END);
    nSize = ftell(pFile);
    // Put the file pointer at the beginning
    fseek(pFile, 0, SEEK_SET);

    // Allocate a buffer big enough
    char* pBuffer = new char[nSize + 1];
    // Put the XML file data in the buffer
    fread(pBuffer, nSize, 1, pFile);
    // Put a 0 char at the end
    pBuffer[nSize] = 0;
    // Close the file
    fclose(pFile), pFile = 0;

    // Return the buffer
    return(pBuffer);
}

int count_children(const Element &elem)
{
  ConstIterator<Element> itBeg = elem.Begin();
  ConstIterator<Element> itEnd = elem.End();
  int n = 0;
  for(ConstIterator<Element> it = itBeg; it < itEnd; ++it)
    n++;
  return n;
}

void operator<<(bool& bValue, const ConstValue& value)
{
    bValue = atoi(value.Get().c_str());
}






