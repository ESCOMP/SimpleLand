<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="namelist_definition">
  <head>
    <title>SLIM Namelist Definition</title>
  </head>
  <body>
<p>
</p>
<hr/>
<p>
</p>
    <h1>Definition of SLIM namelist variables</h1>
    <p>We list all of the relevant namelist variables for SLIM cases. This includes
    SLIM Namelist items.</p>
<hr/>
    <h2>Definition of SLIM namelist variables</h2>
    <p>Note, these all would go into the user_nl_slim file</p>
    <p>Included in the table are the following pieces of information:</p>
    <ul>
    <li>Variable name.</li>
    <li>Variable type (<code>char</code>, <code>integer</code>,
    <code>real</code>, or <code>logical</code>).  The type
    <code>char</code> has the length appended
    following an asterisk, e.g., <code>char*256</code>.  Variables that are
    arrays have their dimension specifier appended inside parentheses.  For
    example <code>char*1(6)</code> denotes a array of six
    <code>char*1</code> values.
    </li>
    <li>Variable description (includes information on defaults).</li>
    <li>Valid values (if restricted).</li>
    </ul>

    <table border="1" cellpadding="10">
    <caption>SLIM Namelist Physics Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[category='slim_physics']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>SLIM Namelist Datasets</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description
      </th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[category='datasets']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>SLIM Namelist History output settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[category='history']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>SLIM Namelist Performance Tuning</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[category='slim_performance']"/>
    </table>

<p>
<hr/>

</p>
</body>
</xsl:template>

<xsl:template match="entry">
  <tr>
    <td rowspan="2"><font color="#ff0000"><xsl:value-of select="@id"/></font></td>
    <td rowspan="2"><xsl:value-of select="type"/></td>
    <td><xsl:value-of select="desc"/><xsl:apply-templates/></td>
  </tr>
  <tr>
    <td colspan="1"><xsl:if test="string-length(valid_values)>0"><b>Valid Values: </b>
         <xsl:value-of select="valid_values"/></xsl:if></td>
  </tr>
</xsl:template>

</xsl:stylesheet>
