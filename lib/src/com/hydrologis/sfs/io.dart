part of dart_sfs;

/// An enumeration of possible Well-Known-Text or Well-Known-Binary ordinates.
/// <p>
/// Intended to be used as an {@code List<Ordinate>}, optimized create methods have been provided for {@link #createXY()}, {@link #createXYM()}, {@link #createXYZ()} and {@link #createXYZM()}.
enum Ordinate {
  /// X-ordinate
  X,

  /// Y-ordinate
  Y,

  /// Z-ordinate
  Z,

  /// Measure-ordinate
  M
}

List<Ordinate> XY = List.from([Ordinate.X, Ordinate.Y]);
List<Ordinate> XYZ = List.from([Ordinate.X, Ordinate.Y, Ordinate.Z]);
List<Ordinate> XYM = List.from([Ordinate.X, Ordinate.Y, Ordinate.M]);
List<Ordinate> XYZM = List.from([Ordinate.X, Ordinate.Y, Ordinate.Z, Ordinate.M]);

/// A filter implementation to test if a coordinate sequence actually has
/// meaningful values for an ordinate bit-pattern
class CheckOrdinatesFilter implements CoordinateSequenceFilter {
  List<Ordinate> checkOrdinateFlags;
  List<Ordinate> outputOrdinates;

  /// Creates an instance of this class
  /// @param checkOrdinateFlags the index for the ordinates to test.
  CheckOrdinatesFilter(List<Ordinate> checkOrdinateFlags) {
    this.outputOrdinates = List.from([Ordinate.X, Ordinate.Y]);
    this.checkOrdinateFlags = checkOrdinateFlags;
  }

  /// @see org.locationtech.jts.geom.CoordinateSequenceFilter#isGeometryChanged */
  void filter(CoordinateSequence seq, int i) {
    if (checkOrdinateFlags.contains(Ordinate.Z) && !outputOrdinates.contains(Ordinate.Z)) {
      if (!seq.getZ(i).isNaN) outputOrdinates.add(Ordinate.Z);
    }

    if (checkOrdinateFlags.contains(Ordinate.M) && !outputOrdinates.contains(Ordinate.M)) {
      if (!seq.getM(i).isNaN) outputOrdinates.add(Ordinate.M);
    }
  }

  /// @see org.locationtech.jts.geom.CoordinateSequenceFilter#isGeometryChanged */
  bool isGeometryChanged() {
    return false;
  }

  /// @see org.locationtech.jts.geom.CoordinateSequenceFilter#isDone */
  bool isDone() {
    return CollectionsUtils.areEqual(outputOrdinates, checkOrdinateFlags);
  }

  /// Gets the evaluated ordinate bit-pattern
  ///
  /// @return A bit-pattern of ordinates with valid values masked by {@link #checkOrdinateFlags}.
  List<Ordinate> getOutputOrdinates() {
    return outputOrdinates;
  }
}

/// Writes the Well-Known Text representation of a {@link Geometry}.
/// The Well-Known Text format is defined in the
/// OGC <a href="http://www.opengis.org/techno/specs.htm">
/// <i>Simple Features Specification for SQL</i></a>.
/// See {@link WKTReader} for a formal specification of the format syntax.
/// <p>
/// The <code>WKTStringBuffer</code> outputs coordinates rounded to the precision
/// model. Only the maximum number of decimal places
/// necessary to represent the ordinates to the required precision will be
/// output.
/// <p>
/// The SFS WKT spec does not define a special tag for {@link LinearRing}s.
/// Under the spec, rings are output as <code>LINESTRING</code>s.
/// In order to allow precisely specifying constructed geometries,
/// JTS also supports a non-standard <code>LINEARRING</code> tag which is used
/// to output LinearRings.
///
/// @version 1.7
/// @see WKTReader
class WKTWriter {
  /// Generates the WKT for a <tt>POINT</tt>
  /// specified by a {@link Coordinate}.
  ///
  /// @param p0 the point coordinate
  ///
  /// @return the WKT
  static String toPoint(Coordinate p0) {
    return "POINT ( ${p0.x}  ${p0.y} )";
  }

  /// Generates the WKT for a <tt>LINESTRING</tt>
  /// specified by a {@link CoordinateSequence}.
  ///
  /// @param seq the sequence to write
  ///
  /// @return the WKT string
  static String toLineStringFromSequence(CoordinateSequence seq) {
    StringBuffer buf = StringBuffer();
    buf.write("LINESTRING ");
    if (seq.size() == 0)
      buf.write(" EMPTY");
    else {
      buf.write("(");
      for (int i = 0; i < seq.size(); i++) {
        if (i > 0) {
          buf.write(", ");
        }
        buf.write(seq.getX(i));
        buf.write(" ");
        buf.write(seq.getY(i));
      }
      buf.write(")");
    }
    return buf.toString();
  }

  /// Generates the WKT for a <tt>LINESTRING</tt>
  /// specified by a {@link CoordinateSequence}.
  ///
  /// @param coord the sequence to write
  ///
  /// @return the WKT string
  static String toLineString(List<Coordinate> coord) {
    StringBuffer buf = new StringBuffer();
    buf.write("LINESTRING ");
    if (coord.length == 0)
      buf.write(" EMPTY");
    else {
      buf.write("(");
      for (int i = 0; i < coord.length; i++) {
        if (i > 0) buf.write(", ");
        buf.write(coord[i].x);
        buf.write(" ");
        buf.write(coord[i].y);
      }
      buf.write(")");
    }
    return buf.toString();
  }

  /// Generates the WKT for a <tt>LINESTRING</tt>
  /// specified by two {@link Coordinate}s.
  ///
  /// @param p0 the first coordinate
  /// @param p1 the second coordinate
  ///
  /// @return the WKT
  static String toLineStringFromCoords(Coordinate p0, Coordinate p1) {
    return "LINESTRING ( ${p0.x} ${p0.y}, ${p1.x} ${p1.y} )";
  }

  static final int INDENT = 2;
  static final int OUTPUT_DIMENSION = 2;

  ///  Creates the <code>NumberFormat</code> used to write <code>double</code>s
  ///  with a sufficient number of decimal places.
  ///
  ///@param  precisionModel  the <code>PrecisionModel</code> used to determine
  ///      the number of decimal places to write.
  ///@return                 a <code>NumberFormat</code> that write <code>double</code>
  ///      s without scientific notation.
  static NumberFormat createFormatter(PrecisionModel precisionModel) {
    // the default number of decimal places is 16, which is sufficient
    // to accommodate the maximum precision of a double.
    int decimalPlaces = precisionModel.getMaximumSignificantDigits();
    // specify decimal separator explicitly to avoid problems in other locales
//    NumberFormatSymbols symbols = new NumberFormatSymbols();
//    symbols.setDecimalSeparator('.');
    String fmtString = "0" + (decimalPlaces > 0 ? "." : "") + stringOfChar('#', decimalPlaces);
    return NumberFormat(fmtString);
  }

  ///  Returns a <code>String</code> of repeated characters.
  ///
  ///@param  ch     the character to repeat
  ///@param  count  the number of times to repeat the character
  ///@return        a <code>String</code> of characters
  static String stringOfChar(String ch, int count) {
    StringBuffer buf = new StringBuffer();
    for (int i = 0; i < count; i++) {
      buf.write(ch);
    }
    return buf.toString();
  }

  List<Ordinate> outputOrdinates;
  int outputDimension;
  PrecisionModel precisionModel = null;
  bool isFormatted = false;
  int coordsPerLine = -1;
  String indentTabStr;

  /// Creates a new WKTStringBuffer with default settings
  WKTWriter() : this.withDimension(OUTPUT_DIMENSION);

  /// Creates a writer that writes {@link Geometry}s with
  /// the given output dimension (2 to 4).
  /// The output follows the following rules:
  /// <ul>
  ///   <li>If the specified <b>output dimension is 3</b> and the <b>z is measure flag
  ///   is set to true</b>, the Z value of coordinates will be written if it is present
  /// (i.e. if it is not <code>Double.NaN</code>)</li>
  ///   <li>If the specified <b>output dimension is 3</b> and the <b>z is measure flag
  ///   is set to false</b>, the Measure value of coordinates will be written if it is present
  /// (i.e. if it is not <code>Double.NaN</code>)</li>
  ///   <li>If the specified <b>output dimension is 4</b>, the Z value of coordinates will
  ///   be written even if it is not present when the Measure value is present.The Measrue
  ///   value of coordinates will be written if it is present
  /// (i.e. if it is not <code>Double.NaN</code>)</li>
  /// </ul>
  ///
  /// @param outputDimension the coordinate dimension to output (2 to 4)
  WKTWriter.withDimension(int outputDimension) {
    setTab(INDENT);
    this.outputDimension = outputDimension;

    if (outputDimension < 2 || outputDimension > 4) throw ArgumentError("Invalid output dimension (must be 2 to 4)");

    this.outputOrdinates = List.from([Ordinate.X, Ordinate.Y]);
    if (outputDimension > 2) {
      outputOrdinates.add(Ordinate.Z);
    }
    if (outputDimension > 3) {
      outputOrdinates.add(Ordinate.M);
    }
  }

  /// Sets whether the output will be formatted.
  ///
  /// @param isFormatted true if the output is to be formatted
  void setFormatted(bool isFormatted) {
    this.isFormatted = isFormatted;
  }

  /// Sets the maximum number of coordinates per line
  /// written in formatted output.
  /// If the provided coordinate number is &lt;= 0,
  /// coordinates will be written all on one line.
  ///
  /// @param coordsPerLine the number of coordinates per line to output.
  void setMaxCoordinatesPerLine(int coordsPerLine) {
    this.coordsPerLine = coordsPerLine;
  }

  /// Sets the tab size to use for indenting.
  ///
  /// @param size the number of spaces to use as the tab string
  /// @throws IllegalArgumentException if the size is non-positive
  void setTab(int size) {
    if (size <= 0) throw ArgumentError("Tab count must be positive");
    this.indentTabStr = stringOfChar(' ', size);
  }

  /// Sets the {@link Ordinate} that are to be written. Possible members are:
  /// <ul>
  /// <li>{@link Ordinate#X}</li>
  /// <li>{@link Ordinate#Y}</li>
  /// <li>{@link Ordinate#Z}</li>
  /// <li>{@link Ordinate#M}</li>
  /// </ul>
  /// Values of {@link Ordinate#X} and {@link Ordinate#Y} are always assumed and not
  /// particularly checked for.
  ///
  /// @param outputOrdinates A set of {@link Ordinate} values
  void setOutputOrdinates(List<Ordinate> outputOrdinates) {
    this.outputOrdinates.remove(Ordinate.Z);
    this.outputOrdinates.remove(Ordinate.M);

    if (this.outputDimension == 3) {
      if (outputOrdinates.contains(Ordinate.Z))
        this.outputOrdinates.add(Ordinate.Z);
      else if (outputOrdinates.contains(Ordinate.M)) this.outputOrdinates.add(Ordinate.M);
    }
    if (this.outputDimension == 4) {
      if (outputOrdinates.contains(Ordinate.Z)) this.outputOrdinates.add(Ordinate.Z);
      if (outputOrdinates.contains(Ordinate.M)) this.outputOrdinates.add(Ordinate.M);
    }
  }

  /// Gets a bit-pattern defining which ordinates should be
  /// @return an ordinate bit-pattern
  /// @see #setOutputOrdinates(List)
  List<Ordinate> getOutputOrdinates() {
    return this.outputOrdinates;
  }

  /// Sets a {@link PrecisionModel} that should be used on the ordinates written.
  /// <p>If none/{@code null} is assigned, the precision model of the {@link Geometry#getFactory()}
  /// is used.</p>
  /// <p>Note: The precision model is applied to all ordinate values, not just x and y.</p>
  /// @param precisionModel
  ///    the flag indicating if {@link Coordinate#z}/{} is actually a measure value.
  void setPrecisionModel(PrecisionModel precisionModel) {
    this.precisionModel = precisionModel;
  }

  ///  Converts a <code>Geometry</code> to its Well-known Text representation.
  ///
  ///@param  geometry  a <code>Geometry</code> to process
  ///@return           a &lt;Geometry Tagged Text&gt; string (see the OpenGIS Simple
  ///      Features Specification)
  String write(Geometry geometry) {
    StringBuffer sw = new StringBuffer();

    // determine the precision model
    PrecisionModel pm = this.precisionModel;
    if (pm == null) pm = geometry.getFactory().getPrecisionModel();

    try {
      writeFormatted4Args(geometry, false, sw, pm);
    } catch (ex) {
      print(ex);
    }
    return sw.toString();
  }

  ///  Converts a <code>Geometry</code> to its Well-known Text representation.
  ///
  ///@param  geometry  a <code>Geometry</code> to process
  void writeWithStringBuffer(Geometry geometry, StringBuffer writer) {
    // determine the precision model
    PrecisionModel pm = this.precisionModel;
    if (pm == null) pm = geometry.getFactory().getPrecisionModel();

    // write the geometry
    writeFormatted4Args(geometry, isFormatted, writer, pm);
  }

  ///  Same as <code>write</code>, but with newlines and spaces to make the
  ///  well-known text more readable.
  ///
  ///@param  geometry  a <code>Geometry</code> to process
  ///@return           a &lt;Geometry Tagged Text&gt; string (see the OpenGIS Simple
  ///      Features Specification), with newlines and spaces
  String writeFormatted(Geometry geometry) {
    StringBuffer sw = new StringBuffer();
    try {
      writeFormatted4Args(geometry, true, sw, precisionModel);
    } catch (ex) {
      print(ex);
    }
    return sw.toString();
  }

  ///  Same as <code>write</code>, but with newlines and spaces to make the
  ///  well-known text more readable.
  ///
  ///@param  geometry  a <code>Geometry</code> to process
  void writeFormattedWithStringBuffer(Geometry geometry, StringBuffer writer) {
    writeFormatted4Args(geometry, true, writer, precisionModel);
  }

  ///  Converts a <code>Geometry</code> to its Well-known Text representation.
  ///
  ///@param  geometry  a <code>Geometry</code> to process
  void writeFormatted4Args(Geometry geometry, bool useFormatting, StringBuffer writer, PrecisionModel precisionModel) {
    // ensure we have a precision model
    if (precisionModel == null) {
      precisionModel = geometry.getPrecisionModel();
    }

    // create the formatter
    NumberFormat formatter = createFormatter(precisionModel);

    // append the WKT
    appendGeometryTaggedText(geometry, useFormatting, writer, formatter);
  }

  ///  Converts a <code>Geometry</code> to &lt;Geometry Tagged Text&gt; format,
  ///  then appends it to the writer.
  ///
  /// @param  geometry           the <code>Geometry</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  writer             the output writer to append to
  /// @param  formatter       the <code>NumberFormatter</code> to use to convert
  ///      from a precise coordinate to an external coordinate
  void appendGeometryTaggedText(Geometry geometry, bool useFormatting, StringBuffer writer, NumberFormat formatter) {
    // evaluate the ordinates actually present in the geometry
    CheckOrdinatesFilter cof = new CheckOrdinatesFilter(this.outputOrdinates);
    geometry.applyCSF(cof);

    // Append the WKT
    appendGeometryTaggedText6Args(geometry, cof.getOutputOrdinates(), useFormatting, 0, writer, formatter);
  }

  ///  Converts a <code>Geometry</code> to &lt;Geometry Tagged Text&gt; format,
  ///  then appends it to the writer.
  ///
  /// @param  geometry           the <code>Geometry</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  level              the indentation level
  /// @param  writer             the output writer to append to
  /// @param  formatter       the <code>NumberFormatter</code> to use to convert
  ///      from a precise coordinate to an external coordinate
  void appendGeometryTaggedText6Args(
      Geometry geometry, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    indent(useFormatting, level, writer);

    if (geometry is Point) {
      appendPointTaggedText(geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is LinearRing) {
      appendLinearRingTaggedText(geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is LineString) {
      appendLineStringTaggedText(geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is Polygon) {
      appendPolygonTaggedText(geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is MultiPoint) {
      appendMultiPointTaggedText(geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is MultiLineString) {
      appendMultiLineStringTaggedText(geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is MultiPolygon) {
      appendMultiPolygonTaggedText(geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is GeometryCollection) {
      appendGeometryCollectionTaggedText(geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else {
      throw ArgumentError("Unsupported Geometry implementation:  ${geometry.runtimeType}");
    }
  }

  ///  Converts a <code>Coordinate</code> to &lt;Point Tagged Text&gt; format,
  ///  then appends it to the writer.
  ///
  /// @param  point           the <code>Point</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  level              the indentation level
  /// @param  writer             the output writer to append to
  /// @param  formatter          the formatter to use when writing numbers
  void appendPointTaggedText(Point point, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    writer.write("POINT ");
    appendOrdinateText(outputOrdinates, writer);
    appendSequenceText(point.getCoordinateSequence(), outputOrdinates, useFormatting, level, false, writer, formatter);
  }

  ///  Converts a <code>LineString</code> to &lt;LineString Tagged Text&gt;
  ///  format, then appends it to the writer.
  ///
  /// @param  lineString  the <code>LineString</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  level              the indentation level
  /// @param  writer             the output writer to append to
  /// @param  formatter       the <code>NumberFormatter</code> to use to convert
  ///      from a precise coordinate to an external coordinate
  void appendLineStringTaggedText(
      LineString lineString, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    writer.write("LINESTRING ");
    appendOrdinateText(outputOrdinates, writer);
    appendSequenceText(lineString.getCoordinateSequence(), outputOrdinates, useFormatting, level, false, writer, formatter);
  }

  ///  Converts a <code>LinearRing</code> to &lt;LinearRing Tagged Text&gt;
  ///  format, then appends it to the writer.
  ///
  /// @param  linearRing  the <code>LinearRing</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  level              the indentation level
  /// @param  writer             the output writer to append to
  /// @param  formatter       the <code>NumberFormatter</code> to use to convert
  ///      from a precise coordinate to an external coordinate
  void appendLinearRingTaggedText(
      LinearRing linearRing, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    writer.write("LINEARRING ");
    appendOrdinateText(outputOrdinates, writer);
    appendSequenceText(linearRing.getCoordinateSequence(), outputOrdinates, useFormatting, level, false, writer, formatter);
  }

  ///  Converts a <code>Polygon</code> to &lt;Polygon Tagged Text&gt; format,
  ///  then appends it to the writer.
  ///
  /// @param  polygon  the <code>Polygon</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  level              the indentation level
  /// @param  writer             the output writer to append to
  /// @param  formatter       the <code>NumberFormatter</code> to use to convert
  ///      from a precise coordinate to an external coordinate
  void appendPolygonTaggedText(Polygon polygon, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    writer.write("POLYGON ");
    appendOrdinateText(outputOrdinates, writer);
    appendPolygonText(polygon, outputOrdinates, useFormatting, level, false, writer, formatter);
  }

  ///  Converts a <code>MultiPoint</code> to &lt;MultiPoint Tagged Text&gt;
  ///  format, then appends it to the writer.
  ///
  /// @param  multipoint  the <code>MultiPoint</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  level              the indentation level
  /// @param  writer             the output writer to append to
  /// @param  formatter       the <code>NumberFormatter</code> to use to convert
  ///      from a precise coordinate to an external coordinate
  void appendMultiPointTaggedText(
      MultiPoint multipoint, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    writer.write("MULTIPOINT ");
    appendOrdinateText(outputOrdinates, writer);
    appendMultiPointText(multipoint, outputOrdinates, useFormatting, level, writer, formatter);
  }

  ///  Converts a <code>MultiLineString</code> to &lt;MultiLineString Tagged
  ///  Text&gt; format, then appends it to the writer.
  ///
  /// @param  multiLineString  the <code>MultiLineString</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  level              the indentation level
  /// @param  writer             the output writer to append to
  /// @param  formatter       the <code>NumberFormatter</code> to use to convert
  ///      from a precise coordinate to an external coordinate
  void appendMultiLineStringTaggedText(
      MultiLineString multiLineString, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    writer.write("MULTILINESTRING ");
    appendOrdinateText(outputOrdinates, writer);
    appendMultiLineStringText(multiLineString, outputOrdinates, useFormatting, level, /*false, */ writer, formatter);
  }

  ///  Converts a <code>MultiPolygon</code> to &lt;MultiPolygon Tagged Text&gt;
  ///  format, then appends it to the writer.
  ///
  /// @param  multiPolygon  the <code>MultiPolygon</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  level              the indentation level
  /// @param  writer             the output writer to append to
  /// @param  formatter       the <code>NumberFormatter</code> to use to convert
  ///      from a precise coordinate to an external coordinate
  void appendMultiPolygonTaggedText(
      MultiPolygon multiPolygon, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    writer.write("MULTIPOLYGON ");
    appendOrdinateText(outputOrdinates, writer);
    appendMultiPolygonText(multiPolygon, outputOrdinates, useFormatting, level, writer, formatter);
  }

  ///  Converts a <code>GeometryCollection</code> to &lt;GeometryCollection
  ///  Tagged Text&gt; format, then appends it to the writer.
  ///
  /// @param  geometryCollection  the <code>GeometryCollection</code> to process
  /// @param  useFormatting      flag indicating that the output should be formatted
  /// @param  level              the indentation level
  /// @param  writer             the output writer to append to
  /// @param  formatter       the <code>NumberFormatter</code> to use to convert
  ///      from a precise coordinate to an external coordinate
  void appendGeometryCollectionTaggedText(
      GeometryCollection geometryCollection, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    writer.write("GEOMETRYCOLLECTION ");
    appendOrdinateText(outputOrdinates, writer);
    appendGeometryCollectionText(geometryCollection, outputOrdinates, useFormatting, level, writer, formatter);
  }

  /// Appends the i'th coordinate from the sequence to the writer
  /// <p>If the {@code seq} has coordinates that are {@link double.NAN}, these are not written, even though
  /// {@link #outputDimension} suggests this.
  ///
  /// @param  seq        the <code>CoordinateSequence</code> to process
  /// @param  i          the index of the coordinate to write
  /// @param  writer     the output writer to append to
  /// @param  formatter  the formatter to use for writing ordinate values
  static void appendCoordinate(CoordinateSequence seq, List<Ordinate> outputOrdinates, int i, StringBuffer writer, NumberFormat formatter) {
    writer.write(writeNumber(seq.getX(i), formatter) + " " + writeNumber(seq.getY(i), formatter));

    if (outputOrdinates.contains(Ordinate.Z)) {
      double z = seq.getZ(i);
      if (!z.isNaN) {
        writer.write(" ");
        writer.write(writeNumber(seq.getZ(i), formatter));
      } else {
        writer.write(" NaN");
      }
    }

    if (outputOrdinates.contains(Ordinate.M)) {
      writer.write(" ");
      writer.write(writeNumber(seq.getM(i), formatter));
    }
  }

  ///  Converts a <code>double</code> to a <code>String</code>, not in scientific
  ///  notation.
  ///
  ///@param  d  the <code>double</code> to convert
  ///@return    the <code>double</code> as a <code>String</code>, not in
  ///      scientific notation
  static String writeNumber(double d, NumberFormat formatter) {
    return formatter.format(d);
  }

  /// Appends additional ordinate information. This function may
  /// <ul>
  ///   <li>append 'Z' if in {@code outputOrdinates} the
  ///   {@link Ordinate#Z} value is included
  ///   </li>
  ///   <li>append 'M' if in {@code outputOrdinates} the
  ///   {@link Ordinate#M} value is included
  ///   </li>
  ///   <li> append 'ZM' if in {@code outputOrdinates} the
  ///   {@link Ordinate#Z} and
  ///   {@link Ordinate#M} values are included/li>
  /// </ul>
  ///
  /// @param outputOrdinates  a bit-pattern of ordinates to write.
  /// @param writer         the output writer to append to.
  /// @throws IOException   if an error occurs while using the writer.
  void appendOrdinateText(List<Ordinate> outputOrdinates, StringBuffer writer) {
    if (outputOrdinates.contains(Ordinate.Z)) {
      writer.write('Z');
    }
    if (outputOrdinates.contains(Ordinate.M)) {
      writer.write('M');
    }
  }

  ///  Appends all members of a <code>CoordinateSequence</code> to the stream. Each {@code Coordinate} is separated from
  ///  another using a colon, the ordinates of a {@code Coordinate} are separated by a space.
  ///
  /// @param  seq             the <code>CoordinateSequence</code> to process
  /// @param  useFormatting   flag indicating that
  /// @param  level           the indentation level
  /// @param  indentFirst     flag indicating that the first {@code Coordinate} of the sequence should be indented for
  ///                         better visibility
  /// @param  writer          the output writer to append to
  /// @param  formatter       the formatter to use for writing ordinate values.
  void appendSequenceText(
      CoordinateSequence seq, List<Ordinate> outputOrdinates, bool useFormatting, int level, bool indentFirst, StringBuffer writer, NumberFormat formatter) {
    if (seq.size() == 0) {
      writer.write("EMPTY");
    } else {
      if (indentFirst) indent(useFormatting, level, writer);
      writer.write("(");
      for (int i = 0; i < seq.size(); i++) {
        if (i > 0) {
          writer.write(", ");
          if (coordsPerLine > 0 && i % coordsPerLine == 0) {
            indent(useFormatting, level + 1, writer);
          }
        }
        appendCoordinate(seq, outputOrdinates, i, writer, formatter);
      }
      writer.write(")");
    }
  }

  ///  Converts a <code>Polygon</code> to &lt;Polygon Text&gt; format, then
  ///  appends it to the writer.
  ///
  /// @param  polygon         the <code>Polygon</code> to process
  /// @param  useFormatting   flag indicating that
  /// @param  level           the indentation level
  /// @param  indentFirst     flag indicating that the first {@code Coordinate} of the sequence should be indented for
  ///                         better visibility
  /// @param  writer          the output writer to append to
  /// @param  formatter       the formatter to use for writing ordinate values.
  void appendPolygonText(
      Polygon polygon, List<Ordinate> outputOrdinates, bool useFormatting, int level, bool indentFirst, StringBuffer writer, NumberFormat formatter) {
    if (polygon.isEmpty) {
      writer.write("EMPTY");
    } else {
      if (indentFirst) indent(useFormatting, level, writer);
      writer.write("(");
      appendSequenceText(polygon._exterior.getCoordinateSequence(), outputOrdinates, useFormatting, level, false, writer, formatter);
      for (int i = 0; i < polygon._interiors.length; i++) {
        writer.write(", ");
        appendSequenceText(polygon._interiors[i].getCoordinateSequence(), outputOrdinates, useFormatting, level + 1, true, writer, formatter);
      }
      writer.write(")");
    }
  }

  ///  Converts a <code>MultiPoint</code> to &lt;MultiPoint Text&gt; format, then
  ///  appends it to the writer.
  ///
  /// @param  multiPoint      the <code>MultiPoint</code> to process
  /// @param  useFormatting   flag indicating that
  /// @param  level           the indentation level
  /// @param  writer          the output writer to append to
  /// @param  formatter       the formatter to use for writing ordinate values.
  void appendMultiPointText(MultiPoint multiPoint, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    if (multiPoint.isEmpty) {
      writer.write("EMPTY");
    } else {
      writer.write("(");
      for (int i = 0; i < multiPoint.getNumGeometries(); i++) {
        if (i > 0) {
          writer.write(", ");
          indentCoords(useFormatting, i, level + 1, writer);
        }
        appendSequenceText(multiPoint.getGeometryN(i).getCoordinateSequence(), outputOrdinates, useFormatting, level, false, writer, formatter);
      }
      writer.write(")");
    }
  }

  ///  Converts a <code>MultiLineString</code> to &lt;MultiLineString Text&gt;
  ///  format, then appends it to the writer.
  ///
  /// @param  multiLineString  the <code>MultiLineString</code> to process
  /// @param  useFormatting    flag indicating that
  /// @param  level            the indentation level
  /// //@param  indentFirst      flag indicating that the first {@code Coordinate} of the sequence should be indented for
  /// //                         better visibility
  /// @param  writer           the output writer to append to
  /// @param  formatter        the formatter to use for writing ordinate values.
  void appendMultiLineStringText(
      MultiLineString multiLineString,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      /*bool indentFirst, */ StringBuffer writer,
      NumberFormat formatter) {
    if (multiLineString.isEmpty) {
      writer.write("EMPTY");
    } else {
      int level2 = level;
      bool doIndent = false;
      writer.write("(");
      for (int i = 0; i < multiLineString.getNumGeometries(); i++) {
        if (i > 0) {
          writer.write(", ");
          level2 = level + 1;
          doIndent = true;
        }
        appendSequenceText(multiLineString.getGeometryN(i).getCoordinateSequence(), outputOrdinates, useFormatting, level2, doIndent, writer, formatter);
      }
      writer.write(")");
    }
  }

  ///  Converts a <code>MultiPolygon</code> to &lt;MultiPolygon Text&gt; format,
  ///  then appends it to the writer.
  ///
  /// @param  multiPolygon  the <code>MultiPolygon</code> to process
  /// @param  useFormatting   flag indicating that
  /// @param  level           the indentation level
  /// @param  writer          the output writer to append to
  /// @param  formatter       the formatter to use for writing ordinate values.
  void appendMultiPolygonText(
      MultiPolygon multiPolygon, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    if (multiPolygon.isEmpty) {
      writer.write("EMPTY");
    } else {
      int level2 = level;
      bool doIndent = false;
      writer.write("(");
      for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
        if (i > 0) {
          writer.write(", ");
          level2 = level + 1;
          doIndent = true;
        }
        appendPolygonText(multiPolygon.getGeometryN(i), outputOrdinates, useFormatting, level2, doIndent, writer, formatter);
      }
      writer.write(")");
    }
  }

  ///  Converts a <code>GeometryCollection</code> to &lt;GeometryCollectionText&gt;
  ///  format, then appends it to the writer.
  ///
  /// @param  geometryCollection  the <code>GeometryCollection</code> to process
  /// @param  useFormatting   flag indicating that
  /// @param  level           the indentation level
  /// @param  writer          the output writer to append to
  /// @param  formatter       the formatter to use for writing ordinate values.
  void appendGeometryCollectionText(
      GeometryCollection geometryCollection, List<Ordinate> outputOrdinates, bool useFormatting, int level, StringBuffer writer, NumberFormat formatter) {
    if (geometryCollection.isEmpty) {
      writer.write("EMPTY");
    } else {
      int level2 = level;
      writer.write("(");
      for (int i = 0; i < geometryCollection.getNumGeometries(); i++) {
        if (i > 0) {
          writer.write(", ");
          level2 = level + 1;
        }
        appendGeometryTaggedText6Args(geometryCollection.getGeometryN(i), outputOrdinates, useFormatting, level2, writer, formatter);
      }
      writer.write(")");
    }
  }

  void indentCoords(bool useFormatting, int coordIndex, int level, StringBuffer writer) {
    if (coordsPerLine <= 0 || coordIndex % coordsPerLine != 0) return;
    indent(useFormatting, level, writer);
  }

  void indent(bool useFormatting, int level, StringBuffer writer) {
    if (!useFormatting || level <= 0) return;
    writer.write("\n");
    for (int i = 0; i < level; i++) {
      writer.write(indentTabStr);
    }
  }
}

/// Converts a geometry in Well-Known Text format to a {@link Geometry}.
/// <p>
/// <code>WKTReader</code> supports
/// extracting <code>Geometry</code> objects from either {@link Reader}s or
///  {@link String}s. This allows it to function as a parser to read <code>Geometry</code>
///  objects from text blocks embedded in other data formats (e.g. XML). <P>
/// <p>
///  A <code>WKTReader</code> is parameterized by a <code>GeometryFactory</code>,
///  to allow it to create <code>Geometry</code> objects of the appropriate
///  implementation. In particular, the <code>GeometryFactory</code>
///  determines the <code>PrecisionModel</code> and <code>SRID</code> that is
///  used. <P>
///
///  The <code>WKTReader</code> converts all input numbers to the precise
///  internal representation.
///
/// <h3>Notes:</h3>
/// <ul>
/// <li>Keywords are case-insensitive.
/// <li>The reader supports non-standard "LINEARRING" tags.
/// <li>The reader uses <tt>Double.parseDouble</tt> to perform the conversion of ASCII
/// numbers to floating point.  This means it supports the Java
/// syntax for floating point literals (including scientific notation).
/// </ul>
///
/// <h3>Syntax</h3>
/// The following syntax specification describes the version of Well-Known Text
/// supported by JTS.
/// (The specification uses a syntax language similar to that used in
/// the C and Java language specifications.)
/// <p>
/// As of version 1.15, JTS can read (but not write) WKT Strings including Z, M or ZM
/// in the name of the geometry type (ex. POINT Z, LINESTRINGZM).
/// Note that it only makes the reader more flexible, but JTS could already read
/// 3D coordinates from WKT String and still can't read 4D coordinates.
///
/// <blockquote><pre>
/// <i>WKTGeometry:</i> one of<i>
///
///       WKTPoint  WKTLineString  WKTLinearRing  WKTPolygon
///       WKTMultiPoint  WKTMultiLineString  WKTMultiPolygon
///       WKTGeometryCollection</i>
///
/// <i>WKTPoint:</i> <b>POINT</b><i>[Dimension]</i> <b>( </b><i>Coordinate</i> <b>)</b>
///
/// <i>WKTLineString:</i> <b>LINESTRING</b><i>[Dimension]</i> <i>CoordinateSequence</i>
///
/// <i>WKTLinearRing:</i> <b>LINEARRING</b><i>[Dimension]</i> <i>CoordinateSequence</i>
///
/// <i>WKTPolygon:</i> <b>POLYGON</b><i>[Dimension]</i> <i>CoordinateSequenceList</i>
///
/// <i>WKTMultiPoint:</i> <b>MULTIPOINT</b><i>[Dimension]</i> <i>CoordinateSingletonList</i>
///
/// <i>WKTMultiLineString:</i> <b>MULTILINESTRING</b><i>[Dimension]</i> <i>CoordinateSequenceList</i>
///
/// <i>WKTMultiPolygon:</i>
///         <b>MULTIPOLYGON</b><i>[Dimension]</i> <b>(</b> <i>CoordinateSequenceList {</i> , <i>CoordinateSequenceList }</i> <b>)</b>
///
/// <i>WKTGeometryCollection: </i>
///         <b>GEOMETRYCOLLECTION</b><i>[Dimension]</i> <b> (</b> <i>WKTGeometry {</i> , <i>WKTGeometry }</i> <b>)</b>
///
/// <i>CoordinateSingletonList:</i>
///         <b>(</b> <i>CoordinateSingleton {</i> <b>,</b> <i>CoordinateSingleton }</i> <b>)</b>
///         | <b>EMPTY</b>
///
/// <i>CoordinateSingleton:</i>
///         <b>(</b> <i>Coordinate</i> <b>)</b>
///         | <b>EMPTY</b>
///
/// <i>CoordinateSequenceList:</i>
///         <b>(</b> <i>CoordinateSequence {</i> <b>,</b> <i>CoordinateSequence }</i> <b>)</b>
///         | <b>EMPTY</b>
///
/// <i>CoordinateSequence:</i>
///         <b>(</b> <i>Coordinate {</i> , <i>Coordinate }</i> <b>)</b>
///         | <b>EMPTY</b>
///
/// <i>Coordinate:
///         Number Number Number<sub>opt</sub></i>
///
/// <i>Number:</i> A Java-style floating-point number (including <tt>NaN</tt>, with arbitrary case)
///
/// <i>Dimension:</i>
///         <b>Z</b>|<b> Z</b>|<b>M</b>|<b> M</b>|<b>ZM</b>|<b> ZM</b>
///
/// </pre></blockquote>
///
///
///@version 1.7
/// @see WKTStringBuffer
class WKTReader {
  static final String EMPTY = "EMPTY";
  static final String COMMA = ",";
  static final String L_PAREN = "(";
  static final String R_PAREN = ")";
  static final String NAN_SYMBOL = "NaN";

  GeometryFactory geometryFactory;
  CoordinateSequenceFactory csFactory;
  static CoordinateSequenceFactory csFactoryXYZM = CoordinateArraySequenceFactory();
  PrecisionModel precisionModel;

  /// Flag indicating that the old notation of coordinates in JTS
  /// is supported.
  static final bool ALLOW_OLD_JTS_COORDINATE_SYNTAX = true;
  bool isAllowOldJtsCoordinateSyntax = ALLOW_OLD_JTS_COORDINATE_SYNTAX;

  /// Flag indicating that the old notation of MultiPoint coordinates in JTS
  /// is supported.
  static final bool ALLOW_OLD_JTS_MULTIPOINT_SYNTAX = true;
  bool isAllowOldJtsMultipointSyntax = ALLOW_OLD_JTS_MULTIPOINT_SYNTAX;

  /// Creates a reader that creates objects using the default {@link GeometryFactory}.
  WKTReader() : this.withFactory(GeometryFactory.defaultPrecision());

  ///  Creates a reader that creates objects using the given
  ///  {@link GeometryFactory}.
  ///
  ///@param  geometryFactory  the factory used to create <code>Geometry</code>s.
  WKTReader.withFactory(GeometryFactory geometryFactory) {
    this.geometryFactory = geometryFactory;
    this.csFactory = geometryFactory.getCoordinateSequenceFactory();
    this.precisionModel = geometryFactory.getPrecisionModel();
  }

  /// Sets a flag indicating, that coordinates may have 3 ordinate values even though no Z or M ordinate indicator
  /// is present. The default value is {@link #ALLOW_OLD_JTS_COORDINATE_SYNTAX}.
  ///
  /// @param value a bool value
  void setIsOldJtsCoordinateSyntaxAllowed(bool value) {
    isAllowOldJtsCoordinateSyntax = value;
  }

  /// Sets a flag indicating, that point coordinates in a MultiPoint geometry must not be enclosed in paren.
  /// The default value is {@link #ALLOW_OLD_JTS_MULTIPOINT_SYNTAX}
  /// @param value a bool value
  void setIsOldJtsMultiPointSyntaxAllowed(bool value) {
    isAllowOldJtsMultipointSyntax = value;
  }

  /// Reads a Well-Known Text representation of a {@link Geometry}
  /// from a {@link String}.
  ///
  /// @param wellKnownText
  ///            one or more &lt;Geometry Tagged Text&gt; strings (see the OpenGIS
  ///            Simple Features Specification) separated by whitespace
  /// @return a <code>Geometry</code> specified by <code>wellKnownText</code>
  /// @throws ParseException
  ///             if a parsing problem occurs
  Geometry read(String wellKnownText) {
    return null; // TODO
  }
}
