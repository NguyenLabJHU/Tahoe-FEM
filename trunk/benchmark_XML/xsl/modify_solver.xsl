<?xml version="1.0"?>
<!-- $Id: modify_solver.xsl,v 1.1 2004-07-12 16:25:33 paklein Exp $ -->
<!-- This style sheet tranforms tahoe XML input file to another tahoe XML
     input file, modifying the solver element by converting matrix_type attribute
     an element with the name <[value of the matrix_type attribute]_matrix> -->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output method="xml" indent="yes"/> 

<!-- identity transformation -->
<xsl:template match="/ | @* | node()">
	<xsl:copy>
		<xsl:apply-templates select="@* | node()" />
	</xsl:copy>
</xsl:template>

<!-- modify solver node -->
<xsl:template match="nonlinear_solver|linear_solver|PCG_solver|nonlinear_solver_LS">
	<xsl:copy>

		<!-- process attributes -->
		<xsl:apply-templates select="@*" />

		<!-- convert matrix_type to an element or profile_matrix if missing -->
		<xsl:choose>
			<xsl:when test="string-length(@matrix_type) > 0">
				<xsl:element name="{@matrix_type}_matrix"/>
			</xsl:when>
			<xsl:otherwise>
				<xsl:element name="profile_matrix"/>
			</xsl:otherwise>
		</xsl:choose>

		<!-- process elements -->
		<xsl:apply-templates select="node()" />
	</xsl:copy>
</xsl:template>

<!-- don't copy matrix_type attribute -->
<xsl:template match="@matrix_type"/>

</xsl:stylesheet> 
