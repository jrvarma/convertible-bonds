<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 version ="1.0">
<xsl:output method="text" indent="no" omit-xml-declaration="yes"/>

<!-- 
     variables used for indentation and for new lines 
     do not change indentation of the line defining _blanks
-->
<xsl:variable name="_blanks">
"                                                       "
</xsl:variable>
<xsl:variable name="blanks" select="substring($_blanks,3,53)"/>
<xsl:variable name="newline" select="substring($_blanks,1,1)"/>


<!-- 
     main template
-->
<xsl:template match="Output//*">
    <xsl:variable name="level" select="count(ancestor::*) - 1"/>
    <xsl:variable name="indent" select="substring($blanks,1,3*$level)"/>
    <xsl:variable name="indent2" select="substring($blanks,1,3*$level+3)"/>
    <xsl:value-of select="$indent"/>
    <xsl:value-of select="name()"/>
    <xsl:value-of select="$newline"/>
    <xsl:for-each select="@*">
	<xsl:value-of select="$indent2"/>
	<xsl:value-of select="name()"/>
	<xsl:text> = </xsl:text>
	<xsl:value-of select="string()"/>
	<xsl:value-of select="$newline"/>
    </xsl:for-each>
    <xsl:apply-templates select="*"/>
</xsl:template>

<!-- 
     template to handle root node
-->
<xsl:template match="Output">
    <xsl:text>All interest rates below are on the basis of </xsl:text>
    <xsl:value-of select="@Output_Interest_Rate_Compounding"/>
    <xsl:text> compounding.</xsl:text>
    <xsl:value-of select="$newline"/>
    <xsl:apply-templates select="*"/>
</xsl:template>

<!-- 
     template to handle STAGE node of full lattice
     all attribute-value pairs are written out in a single line
-->
<xsl:template match="//Stage">
    <xsl:variable name="level" select="count(ancestor::*) - 1"/>
    <xsl:variable name="indent" select="substring($blanks,1,3*$level)"/>
    <xsl:value-of select="$indent"/>
    <xsl:for-each select="@*">
	<xsl:value-of select="name()"/>
	<xsl:text> = </xsl:text>
	<xsl:value-of select="string()"/>
	<xsl:text>.  </xsl:text>
    </xsl:for-each>
    <xsl:value-of select="$newline"/>
    <xsl:apply-templates select="*"/>
</xsl:template>

<!-- 
     template to handle following nodes - NODE, CALL, PUT
     all attribute-value pairs are written out in a single line inside square brackets []
-->
<xsl:template match="//Node|//Call|//Put">
    <xsl:variable name="level" select="count(ancestor::*) - 1"/>
    <xsl:variable name="indent" select="substring($blanks,1,3*$level)"/>
    <xsl:value-of select="$indent"/>
    <xsl:text>[</xsl:text>
    <xsl:for-each select="@*">
	<xsl:value-of select="name()"/>
	<xsl:text> = </xsl:text>
	<xsl:value-of select="string()"/>
	<xsl:if test="not(position()= last())">
	    <xsl:text>  </xsl:text>
	</xsl:if>
    </xsl:for-each>
    <xsl:text>]</xsl:text>
    <xsl:value-of select="$newline"/>
    <xsl:apply-templates select="*"/>
</xsl:template>

</xsl:stylesheet>

