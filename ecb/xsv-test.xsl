<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
 xmlns="http://www.w3.org/TR/REC-html40" 
 xmlns:xsv="http://www.w3.org/2000/05/xsv" version="1.0">
<xsl:output method="text" indent="no" omit-xml-declaration="yes"/>

<xsl:variable name="_blanks">
a
</xsl:variable>
<xsl:variable name="newline" select="substring($_blanks,1,1)"/>

<xsl:template match="//xsv:invalid">
    <xsl:text>Schema error in line </xsl:text>
    <xsl:value-of select="@line"/>
    <xsl:text> of </xsl:text>
    <xsl:value-of select="@resource"/>
    <xsl:value-of select="$newline"/>
    <xsl:value-of select="string(./text())"/>
    <xsl:value-of select="$newline"/>
</xsl:template>

<xsl:template match="//xsv:XMLMessages">
    <xsl:text>Low level XML error</xsl:text>
    <xsl:value-of select="string()"/>
</xsl:template>

<xsl:template match="xsv:xsv">
    <xsl:apply-templates select="//xsv:XMLMessages"/>
    <xsl:apply-templates select="//xsv:invalid"/>
</xsl:template>

</xsl:stylesheet>

