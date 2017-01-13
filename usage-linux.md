USAGE
=====

1. Run make in the folder `ecb/source` to create the executable file `ecb.exe`. Rename this file to `ecb` in accordance with unix conventions. Alternatively, you may want to modify the Makefile to name the executable as `ecb` instead of `ecb.exe`. 

2. (Optional) Install `xmllint` from `libxml2`

3. (Optional) Install `xsltproc` from `libxslt`

4. For the following steps `cd` to the `ecb` folder.

5. (Optional) Validate input file (say `convertible-1.xml`) against the schema `ecbdata.xsd` using `xmllint`:

    `xmllint --noout --schema ecbdata.xsd convertible-1.xml`

6. Run `ecb` on input file (say `convertible-1.xml`):

    `./ecb convertible-1.xml`
    
7. (Optional) For more human readable output, transform the `XML` output of `ecb` using the stylesheet `xmlout.xsl` with `xsltproc`:

    `./ecb convertible-1.xml | xsltproc xmlout.xsl -`
    

COPYRIGHT
---------

The program `ecb` is copyrighted and distributed under GNU GPL. Copyright (C) 2001 Prof. Jayanth R. Varma, jrvarma@iima.ac.in, Indian Institute of Management, Ahmedabad 380 015, INDIA

The program ecb uses the software `AdvXMLParser` Copyright 1999,2000 Sebastien Andrivet which is covered by a separate license (see the `License` file in the `AdvXMLParser` folder)

