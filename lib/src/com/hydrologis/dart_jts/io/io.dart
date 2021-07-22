part of dart_jts;

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

final List<Ordinate> OrdinateSet_XY = List.from([Ordinate.X, Ordinate.Y]);
final List<Ordinate> OrdinateSet_XYZ =
    List.from([Ordinate.X, Ordinate.Y, Ordinate.Z]);
final List<Ordinate> OrdinateSet_XYM =
    List.from([Ordinate.X, Ordinate.Y, Ordinate.M]);
final List<Ordinate> OrdinateSet_XYZM =
    List.from([Ordinate.X, Ordinate.Y, Ordinate.Z, Ordinate.M]);

/// A filter implementation to test if a coordinate sequence actually has
/// meaningful values for an ordinate bit-pattern
class CheckOrdinatesFilter implements CoordinateSequenceFilter {
  late List<Ordinate> checkOrdinateFlags;
  late List<Ordinate> outputOrdinates;

  /// Creates an instance of this class
  /// @param checkOrdinateFlags the index for the ordinates to test.
  CheckOrdinatesFilter(List<Ordinate> checkOrdinateFlags) {
    this.outputOrdinates = List.from([Ordinate.X, Ordinate.Y]);
    this.checkOrdinateFlags = checkOrdinateFlags;
  }

  /// @see org.locationtech.jts.geom.CoordinateSequenceFilter#isGeometryChanged */
  void filter(CoordinateSequence seq, int i) {
    if (checkOrdinateFlags.contains(Ordinate.Z) &&
        !outputOrdinates.contains(Ordinate.Z)) {
      if (!seq.getZ(i).isNaN) outputOrdinates.add(Ordinate.Z);
    }

    if (checkOrdinateFlags.contains(Ordinate.M) &&
        !outputOrdinates.contains(Ordinate.M)) {
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
    StringBuffer buf = StringBuffer();
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
//    NumberFormatSymbols symbols = NumberFormatSymbols();
//    symbols.setDecimalSeparator('.');
    String fmtString =
        "0" + (decimalPlaces > 0 ? "." : "") + stringOfChar('#', decimalPlaces);
    return NumberFormat(fmtString);
  }

  ///  Returns a <code>String</code> of repeated characters.
  ///
  ///@param  ch     the character to repeat
  ///@param  count  the number of times to repeat the character
  ///@return        a <code>String</code> of characters
  static String stringOfChar(String ch, int count) {
    StringBuffer buf = StringBuffer();
    for (int i = 0; i < count; i++) {
      buf.write(ch);
    }
    return buf.toString();
  }

  late List<Ordinate> outputOrdinates;
  int outputDimension = 0;
  PrecisionModel? precisionModel = null;
  bool isFormatted = false;
  int coordsPerLine = -1;
  String? indentTabStr;

  /// Creates a WKTStringBuffer with default settings
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

    if (outputDimension < 2 || outputDimension > 4)
      throw ArgumentError("Invalid output dimension (must be 2 to 4)");

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
      else if (outputOrdinates.contains(Ordinate.M))
        this.outputOrdinates.add(Ordinate.M);
    }
    if (this.outputDimension == 4) {
      if (outputOrdinates.contains(Ordinate.Z))
        this.outputOrdinates.add(Ordinate.Z);
      if (outputOrdinates.contains(Ordinate.M))
        this.outputOrdinates.add(Ordinate.M);
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
    StringBuffer sw = StringBuffer();

    // determine the precision model
    PrecisionModel? pm = this.precisionModel;
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
    PrecisionModel? pm = this.precisionModel;
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
    StringBuffer sw = StringBuffer();
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
  void writeFormatted4Args(Geometry geometry, bool useFormatting,
      StringBuffer writer, PrecisionModel? precisionModel) {
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
  void appendGeometryTaggedText(Geometry geometry, bool useFormatting,
      StringBuffer writer, NumberFormat formatter) {
    // evaluate the ordinates actually present in the geometry
    CheckOrdinatesFilter cof = CheckOrdinatesFilter(this.outputOrdinates);
    geometry.applyCSF(cof);

    // Append the WKT
    appendGeometryTaggedText6Args(geometry, cof.getOutputOrdinates(),
        useFormatting, 0, writer, formatter);
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
      Geometry geometry,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    indent(useFormatting, level, writer);

    if (geometry is Point) {
      appendPointTaggedText(
          geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is LinearRing) {
      appendLinearRingTaggedText(
          geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is LineString) {
      appendLineStringTaggedText(
          geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is Polygon) {
      appendPolygonTaggedText(
          geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is MultiPoint) {
      appendMultiPointTaggedText(
          geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is MultiLineString) {
      appendMultiLineStringTaggedText(
          geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is MultiPolygon) {
      appendMultiPolygonTaggedText(
          geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else if (geometry is GeometryCollection) {
      appendGeometryCollectionTaggedText(
          geometry, outputOrdinates, useFormatting, level, writer, formatter);
    } else {
      throw ArgumentError(
          "Unsupported Geometry implementation:  ${geometry.runtimeType}");
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
  void appendPointTaggedText(
      Point point,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    writer.write("POINT ");
    appendOrdinateText(outputOrdinates, writer);
    appendSequenceText(point.getCoordinateSequence(), outputOrdinates,
        useFormatting, level, false, writer, formatter);
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
      LineString lineString,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    writer.write("LINESTRING ");
    appendOrdinateText(outputOrdinates, writer);
    appendSequenceText(lineString.getCoordinateSequence(), outputOrdinates,
        useFormatting, level, false, writer, formatter);
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
      LinearRing linearRing,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    writer.write("LINEARRING ");
    appendOrdinateText(outputOrdinates, writer);
    appendSequenceText(linearRing.getCoordinateSequence(), outputOrdinates,
        useFormatting, level, false, writer, formatter);
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
  void appendPolygonTaggedText(
      Polygon polygon,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    writer.write("POLYGON ");
    appendOrdinateText(outputOrdinates, writer);
    appendPolygonText(polygon, outputOrdinates, useFormatting, level, false,
        writer, formatter);
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
      MultiPoint multipoint,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    writer.write("MULTIPOINT ");
    appendOrdinateText(outputOrdinates, writer);
    appendMultiPointText(
        multipoint, outputOrdinates, useFormatting, level, writer, formatter);
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
      MultiLineString multiLineString,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    writer.write("MULTILINESTRING ");
    appendOrdinateText(outputOrdinates, writer);
    appendMultiLineStringText(multiLineString, outputOrdinates, useFormatting,
        level, /*false, */ writer, formatter);
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
      MultiPolygon multiPolygon,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    writer.write("MULTIPOLYGON ");
    appendOrdinateText(outputOrdinates, writer);
    appendMultiPolygonText(
        multiPolygon, outputOrdinates, useFormatting, level, writer, formatter);
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
      GeometryCollection geometryCollection,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    writer.write("GEOMETRYCOLLECTION ");
    appendOrdinateText(outputOrdinates, writer);
    appendGeometryCollectionText(geometryCollection, outputOrdinates,
        useFormatting, level, writer, formatter);
  }

  /// Appends the i'th coordinate from the sequence to the writer
  /// <p>If the {@code seq} has coordinates that are {@link double.NAN}, these are not written, even though
  /// {@link #outputDimension} suggests this.
  ///
  /// @param  seq        the <code>CoordinateSequence</code> to process
  /// @param  i          the index of the coordinate to write
  /// @param  writer     the output writer to append to
  /// @param  formatter  the formatter to use for writing ordinate values
  static void appendCoordinate(
      CoordinateSequence seq,
      List<Ordinate> outputOrdinates,
      int i,
      StringBuffer writer,
      NumberFormat formatter) {
    writer.write(writeNumber(seq.getX(i), formatter) +
        " " +
        writeNumber(seq.getY(i), formatter));

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
  /// @   if an error occurs while using the writer.
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
      CoordinateSequence seq,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      bool indentFirst,
      StringBuffer writer,
      NumberFormat formatter) {
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
      Polygon polygon,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      bool indentFirst,
      StringBuffer writer,
      NumberFormat formatter) {
    if (polygon.isEmpty()) {
      writer.write("EMPTY");
    } else {
      if (indentFirst) indent(useFormatting, level, writer);
      writer.write("(");
      appendSequenceText(polygon.getExteriorRing().getCoordinateSequence(),
          outputOrdinates, useFormatting, level, false, writer, formatter);
      for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
        writer.write(", ");
        appendSequenceText(polygon.getInteriorRingN(i).getCoordinateSequence(),
            outputOrdinates, useFormatting, level + 1, true, writer, formatter);
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
  void appendMultiPointText(
      MultiPoint multiPoint,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    if (multiPoint.isEmpty()) {
      writer.write("EMPTY");
    } else {
      writer.write("(");
      for (int i = 0; i < multiPoint.getNumGeometries(); i++) {
        if (i > 0) {
          writer.write(", ");
          indentCoords(useFormatting, i, level + 1, writer);
        }
        appendSequenceText(
            (multiPoint.getGeometryN(i) as Point).getCoordinateSequence(),
            outputOrdinates,
            useFormatting,
            level,
            false,
            writer,
            formatter);
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
    if (multiLineString.isEmpty()) {
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
        appendSequenceText(
            (multiLineString.getGeometryN(i) as LineString)
                .getCoordinateSequence(),
            outputOrdinates,
            useFormatting,
            level2,
            doIndent,
            writer,
            formatter);
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
      MultiPolygon multiPolygon,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    if (multiPolygon.isEmpty()) {
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
        appendPolygonText(
            multiPolygon.getGeometryN(i) as Polygon,
            outputOrdinates,
            useFormatting,
            level2,
            doIndent,
            writer,
            formatter);
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
      GeometryCollection geometryCollection,
      List<Ordinate> outputOrdinates,
      bool useFormatting,
      int level,
      StringBuffer writer,
      NumberFormat formatter) {
    if (geometryCollection.isEmpty()) {
      writer.write("EMPTY");
    } else {
      int level2 = level;
      writer.write("(");
      for (int i = 0; i < geometryCollection.getNumGeometries(); i++) {
        if (i > 0) {
          writer.write(", ");
          level2 = level + 1;
        }
        appendGeometryTaggedText6Args(geometryCollection.getGeometryN(i),
            outputOrdinates, useFormatting, level2, writer, formatter);
      }
      writer.write(")");
    }
  }

  void indentCoords(
      bool useFormatting, int coordIndex, int level, StringBuffer writer) {
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

class WKTTokenType {
  static const QUOTED_NAME = WKTTokenType._(0, "QUOTED_NAME");
  static const KEYWORD = WKTTokenType._(1, "KEYWORD");
  static const NUMERIC_LITERAL = WKTTokenType._(2, "NUMERIC_LITERAL");
  static const WHITESPACE = WKTTokenType._(3, "WHITESPACE");
  static const LPAREN = WKTTokenType._(4, "LPAREN");
  static const RPAREN = WKTTokenType._(5, "RPAREN");
  static const LBRACKET = WKTTokenType._(6, "LBRACKET");
  static const RBRACKET = WKTTokenType._(7, "RBRACKET");
  static const ERROR = WKTTokenType._(8, "ERROR");
  static const EOS = WKTTokenType._(9, "EOS");
  static const COMMA = WKTTokenType._(10, "COMMA");
  static const NOTHING = WKTTokenType._(-1, "NOTHING");

  final _type;
  final _displayName;

  const WKTTokenType._(this._type, this._displayName);

  String toString() => _displayName;
}

class WKTToken {
  final type;
  final value;

  const WKTToken(this.type, this.value);

  String toString() => "{${type}: $value}";

  bool get isKeyword => type == WKTTokenType.KEYWORD;

  bool get isLParen => type == WKTTokenType.LPAREN;

  bool get isRParen => type == WKTTokenType.RPAREN;

  bool get isComma => type == WKTTokenType.COMMA;

  bool get isNumber => type == WKTTokenType.NUMERIC_LITERAL;

  bool get isEOS => type == WKTTokenType.EOS;

  ensureKeyword([kw = null]) {
    if (!isKeyword) {
      throw ArgumentError("expected a keyword, got '$value'");
    }
    if (kw != null && value.toLowerCase() != kw.toLowerCase()) {
      throw ArgumentError("unexpected keyword: expected '$kw', got '$value'");
    }
  }

  ensureLParen() {
    if (!isLParen) throw ArgumentError("expected left paren '(', got '$value'");
  }

  ensureRParen() {
    if (!isRParen)
      throw ArgumentError("expected right paren ')', got '$value'");
  }

  ensureNumber() {
    if (!isNumber)
      throw ArgumentError("expected numeric literal, got '$value'");
  }

  ensureNotEOS() {
    if (isEOS) throw ArgumentError("unexpected end of input");
  }

  bool matchKeyword(kw) => isKeyword && value.toLowerCase() == kw.toLowerCase();
}

class WKTTokenizer {
  static final CC_0 = "0".codeUnitAt(0);
  static final CC_9 = "9".codeUnitAt(0);
  static final CC_a = "a".codeUnitAt(0);
  static final CC_z = "z".codeUnitAt(0);
  static final CC_A = "A".codeUnitAt(0);
  static final CC_Z = "Z".codeUnitAt(0);
  static final CC_SPACE = " ".codeUnitAt(0);
  static final CC_CR = "\n".codeUnitAt(0);
  static final CC_LF = "\r".codeUnitAt(0);
  static final CC_LPAREN = "(".codeUnitAt(0);
  static final CC_RPAREN = ")".codeUnitAt(0);
  static final CC_DQUOTE = "\"".codeUnitAt(0);
  static final CC_QUOTE = "'".codeUnitAt(0);
  static final CC_MINUS = "-".codeUnitAt(0);
  static final CC_PLUS = "+".codeUnitAt(0);
  static final CC_UNDERSCORE = "_".codeUnitAt(0);
  static final CC_PERIOD = ".".codeUnitAt(0);
  static final CC_COMMA = ",".codeUnitAt(0);
  static final CC_DOT = ".".codeUnitAt(0);
  static final CC_E = "E".codeUnitAt(0);
  static final CC_e = "e".codeUnitAt(0);

  static bool isDigit(int c) => c >= CC_0 && c <= CC_9;

  static bool isLowerCaseLetter(int c) => c >= CC_a && c <= CC_z;

  static bool isUpperCaseLetter(int c) => c >= CC_A && c <= CC_Z;

  static bool isLetter(int c) => isLowerCaseLetter(c) || isUpperCaseLetter(c);

  static bool isWS(int c) => c == CC_SPACE || c == CC_CR || c == CC_LF;

  static bool isParen(int c) => c == CC_LPAREN || c == CC_RPAREN;

  static bool isDQuote(int c) => c == CC_DQUOTE;

  static bool isSpecial(int c) =>
      c == CC_RPAREN ||
      c == CC_LPAREN ||
      c == CC_MINUS ||
      c == CC_UNDERSCORE ||
      c == CC_PERIOD ||
      c == CC_QUOTE ||
      c == CC_SPACE;

  static bool isComma(int c) => c == CC_COMMA;

  static bool isSign(int c) => c == CC_MINUS || c == CC_PLUS;

  static bool isDot(int c) => c == CC_DOT;

  static bool isE(int c) => c == CC_E || c == CC_e;

  static bool isDelimiter(int c) =>
      isWS(c) || isParen(c) || isComma(c) || isDQuote(c) || isDot(c);

  static bool isLParen(int c) => c == CC_LPAREN;

  static bool isRParen(int c) => c == CC_RPAREN;

  static errorToken(value) => WKTToken(WKTTokenType.ERROR, value);

  static numericLiteralToken(value) =>
      WKTToken(WKTTokenType.NUMERIC_LITERAL, value);

  static eosToken() => WKTToken(WKTTokenType.EOS, null);

  static wsToken(value) => WKTToken(WKTTokenType.WHITESPACE, value);

  static keywordToken(value) => WKTToken(WKTTokenType.KEYWORD, value);

  static commaToken(value) => WKTToken(WKTTokenType.COMMA, value);

  static quotedNameToken(value) => WKTToken(WKTTokenType.QUOTED_NAME, value);

  late Iterable<int> source;
  int head = 0;
  bool eos = false;
  bool skipWhitespace = true;

  WKTToken? _lastReadToken;
  bool pushedBack = false;

  WKTTokenizer(source, {this.skipWhitespace: true}) {
    if (source is Iterable<int>) {
      this.source = source;
    } else if (source is String) {
      this.source = source.codeUnits;
    }
    try {
      this.source.elementAt(0);
    } on RangeError catch (e) {
      eos = true;
    }
  }

  advance() {
    head++;
    try {
      source.elementAt(head);
      return true;
    } on RangeError catch (e) {
      eos = true;
      return false;
    }
  }

  pushBack() {
    if (_lastReadToken != null) pushedBack = true;
  }

  int get cur => source.elementAt(head);

  currentToken() => String.fromCharCodes(source.take(head));

  consumeToken() {
    var token = currentToken();
    source = source.skip(head);
    head = 0;
    return token;
  }

  consumeWhiteSpace() {
    assert(isWS(cur));
    while (advance()) {
      if (!isWS(cur)) break;
    }
    var token = consumeToken();
    return wsToken(token);
  }

  consumeComma() {
    assert(isComma(cur));
    advance();
    var token = consumeToken();
    return commaToken(token);
  }

  consumeKeyword() {
    assert(isLetter(cur));
    while (advance()) {
      var c = cur;
      if (isLetter(c)) continue;
      if (isWS(c)) break;
      if (isParen(c)) break;
      if (isComma(c)) break;
      return errorToken(
          "unexpected character in keyword at <${currentToken()}>");
    }
    var token = consumeToken();
    return keywordToken(token);
  }

  consumeQuotedName() {
    assert(isDQuote(cur));
    while (advance()) {
      var c = cur;
      if (isWS(c) || isLetter(c) || isDigit(c) || isSpecial(c)) continue;
      if (isDQuote(c)) {
        advance();
        var token = consumeToken();
        return quotedNameToken(token);
      }
      return errorToken(
          "unexpected character in quoted name at  <${currentToken()}>");
    }
    if (eos) {
      return errorToken("quoted name not terminated at <${currentToken()}>");
    }
  }

  scanUnsignedInteger() {
    assert(isDigit(cur));
    while (advance()) {
      if (!isDigit(cur)) break;
    }
  }

  scanSignedInteger() {
    assert(isDigit(cur) || isSign(cur));
    if (isSign(cur)) {
      if (!advance()) throw "unexpected end of stream";
      if (!isDigit(cur)) throw "unexpected character after sign";
    }
    scanUnsignedInteger();
  }

  scanExactNumericLiteral() {
    if (isDot(cur)) {
      if (!advance()) return;
      scanUnsignedInteger();
    } else if (isDigit(cur)) {
      scanUnsignedInteger();
      if (eos) return;
      if (isDot(cur)) {
        if (!advance()) return;
        scanUnsignedInteger();
      }
    }
  }

  scanUnsignedNumericLiteral() {
    scanExactNumericLiteral();
    if (eos) return;
    if (isE(cur)) {
      if (!advance()) throw "unexpected end of stream";
      if (!isDigit(cur) && !isSign(cur)) throw "unexpected character after 'E'";
      scanSignedInteger();
    }
  }

  scanSignedNumericLiteral() {
    if (isSign(cur)) {
      if (!advance()) throw "unexpected end of stream";
      if (!(isDot(cur) || isDigit(cur)))
        throw "unexpected character after sign";
    }
    scanUnsignedNumericLiteral();
  }

  consumeSignedNumericLiteral() {
    try {
      scanSignedNumericLiteral();
      var token = consumeToken();
      return numericLiteralToken(token);
    } catch (e, st) {
      print(st);
      return errorToken(e.toString());
    }
  }

  consumeLParen() {
    assert(isLParen(cur));
    advance();
    var token = consumeToken();
    return WKTToken(WKTTokenType.LPAREN, "(");
  }

  consumeRParen() {
    assert(isRParen(cur));
    advance();
    var token = consumeToken();
    return WKTToken(WKTTokenType.RPAREN, ")");
  }

  WKTToken? next() {
    if (pushedBack) {
      pushedBack = false;
      return _lastReadToken;
    }
    if (eos) _lastReadToken = eosToken();
    while (!eos) {
      var c = cur;
      if (isWS(c)) {
        var token = consumeWhiteSpace();
        if (!skipWhitespace) {
          _lastReadToken = token;
          break;
        }
      } else if (isLetter(c)) {
        _lastReadToken = consumeKeyword();
        break;
      } else if (isDQuote(c)) {
        _lastReadToken = consumeQuotedName();
        break;
      } else if (isComma(c)) {
        _lastReadToken = consumeComma();
        break;
      } else if (isLParen(c)) {
        _lastReadToken = consumeLParen();
        break;
      } else if (isRParen(c)) {
        _lastReadToken = consumeRParen();
        break;
      } else if (isDigit(c) || isDot(c) || isSign(c)) {
        _lastReadToken = consumeSignedNumericLiteral();
        break;
      } else {
        _lastReadToken =
            errorToken("unexpected character at <${currentToken()}>");
      }
    }
    return _lastReadToken;
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
/// @see WKTWriter
class WKTReader {
  static final String EMPTY = "EMPTY";
  static final String COMMA = ",";
  static final String L_PAREN = "(";
  static final String R_PAREN = ")";
  static final String NAN_SYMBOL = "NaN";

  GeometryFactory geometryFactory;
  late CoordinateSequenceFactory csFactory;
  static CoordinateSequenceFactory csFactoryXYZM =
      CoordinateArraySequenceFactory();
  late PrecisionModel precisionModel;

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
  WKTReader.withFactory(this.geometryFactory) {
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
  /// @
  Geometry? read(String wellKnownText) {
    return readGeometryTaggedText(WKTTokenizer(wellKnownText));
//  StringReader reader = StringReader(wellKnownText);
//  try {
//  return read(reader);
//  }
//  finally {
//  reader.close();
//  }
  }

  ///**
// * Reads a Well-Known Text representation of a {@link Geometry}
// * from a {@link Reader}.
// *
// *@param  reader           a Reader which will return a &lt;Geometry Tagged Text&gt;
// *      string (see the OpenGIS Simple Features Specification)
// *@return                  a <code>Geometry</code> read from <code>reader</code>
// *@throws  ParseException  if a parsing problem occurs
// */
// Geometry read(Reader reader)  {
//WKTTokenizer tokenizer = createTokenizer(reader);
//try {
//return readGeometryTaggedText(tokenizer);
//}
//catch (IOException e) {
//throw ParseException(e.toString());
//}
//}

  ///**
// * Utility function to create the tokenizer
// * @param reader a reader
// *
// * @return a WKT Tokenizer.
// */
// static WKTTokenizer createTokenizer(Reader reader) {
//  WKTTokenizer tokenizer = WKTTokenizer(reader);
//  // set tokenizer to NOT parse numbers
//  tokenizer.resetSyntax();
//  tokenizer.wordChars('a', 'z');
//  tokenizer.wordChars('A', 'Z');
//  tokenizer.wordChars(128 + 32, 255);
//  tokenizer.wordChars('0', '9');
//  tokenizer.wordChars('-', '-');
//  tokenizer.wordChars('+', '+');
//  tokenizer.wordChars('.', '.');
//  tokenizer.whitespaceChars(0, ' ');
//  tokenizer.commentChar('#');
//
//  return tokenizer;
//}

  /// Reads a <code>Coordinate</Code> from a stream using the given {@link WKTTokenizer}.
  /// <p>
  ///   All ordinate values are read, but -depending on the {@link CoordinateSequenceFactory} of the
  ///   underlying {@link GeometryFactory}- not necessarily all can be handled. Those are silently dropped.
  /// </p>
  /// @param tokenizer the tokenizer to use
  /// @param ordinateFlags a bit-mask defining the ordinates to read.
  /// @param tryParen a value indicating if a starting {@link #L_PAREN} should be probed.
  /// @return a {@link CoordinateSequence} of length 1 containing the read ordinate values
  ///
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  CoordinateSequence getCoordinate(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags, bool tryParen) {
    bool opened = false;
    if (tryParen && isOpenerNext(tokenizer)) {
      tokenizer.next();
      opened = true;
    }

// create a sequence for one coordinate
    int offsetM = ordinateFlags.contains(Ordinate.Z) ? 1 : 0;
    CoordinateSequence sequence = csFactory.createSizeDimMeas(1,
        toDimension(ordinateFlags), ordinateFlags.contains(Ordinate.M) ? 1 : 0);
    sequence.setOrdinate(0, CoordinateSequence.X,
        precisionModel.makePrecise(getNextNumber(tokenizer)));
    sequence.setOrdinate(0, CoordinateSequence.Y,
        precisionModel.makePrecise(getNextNumber(tokenizer)));

// additionally read other vertices
    if (ordinateFlags.contains(Ordinate.Z)) {
      sequence.setOrdinate(0, CoordinateSequence.Z, getNextNumber(tokenizer));
    }
    if (ordinateFlags.contains(Ordinate.M)) {
      sequence.setOrdinate(
          0, CoordinateSequence.Z + offsetM, getNextNumber(tokenizer));
    }

    if (ordinateFlags.length == 2 &&
        this.isAllowOldJtsCoordinateSyntax &&
        isNumberNext(tokenizer)) {
      sequence.setOrdinate(0, CoordinateSequence.Z, getNextNumber(tokenizer));
    }

// read close token if it was opened here
    if (opened) {
      getNextCloser(tokenizer);
    }

    return sequence;
  }

  /// Reads a <code>Coordinate</Code> from a stream using the given {@link WKTTokenizer}.
  /// <p>
  ///   All ordinate values are read, but -depending on the {@link CoordinateSequenceFactory} of the
  ///   underlying {@link GeometryFactory}- not necessarily all can be handled. Those are silently dropped.
  /// </p>
  /// <p>
  ///
  /// </p>
  /// @param tokenizer the tokenizer to use
  /// @param ordinateFlags a bit-mask defining the ordinates to read.
  /// @return a {@link CoordinateSequence} of length 1 containing the read ordinate values
  ///
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  CoordinateSequence getCoordinateSequence(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags) {
    return getCoordinateSequenceTryParen(tokenizer, ordinateFlags, false);
  }

  /// Reads a <code>CoordinateSequence</Code> from a stream using the given {@link WKTTokenizer}.
  /// <p>
  ///   All ordinate values are read, but -depending on the {@link CoordinateSequenceFactory} of the
  ///   underlying {@link GeometryFactory}- not necessarily all can be handled. Those are silently dropped.
  /// </p>
  /// @param tokenizer the tokenizer to use
  /// @param ordinateFlags a bit-mask defining the ordinates to read.
  /// @param tryParen a value indicating if a starting {@link #L_PAREN} should be probed for each coordinate.
  /// @return a {@link CoordinateSequence} of length 1 containing the read ordinate values
  ///
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  CoordinateSequence getCoordinateSequenceTryParen(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags, bool tryParen) {
    if (getNextEmptyOrOpener(tokenizer) == EMPTY)
      return this.csFactory.createSizeDimMeas(0, toDimension(ordinateFlags),
          ordinateFlags.contains(Ordinate.M) ? 1 : 0);

    List<CoordinateSequence> coordinates = [];
    do {
      coordinates.add(getCoordinate(tokenizer, ordinateFlags, tryParen));
    } while (getNextCloserOrComma(tokenizer) == COMMA);

    return mergeSequences(coordinates, ordinateFlags);
  }

  /// Computes the required dimension based on the given ordinate values.
  /// It is assumed that {@link Ordinate#X} and {@link Ordinate#Y} are included.
  ///
  /// @param ordinateFlags the ordinate bit-mask
  /// @return the number of dimensions required to store ordinates for the given bit-mask.
  int toDimension(List<Ordinate> ordinateFlags) {
    int dimension = 2;
    if (ordinateFlags.contains(Ordinate.Z)) dimension++;
    if (ordinateFlags.contains(Ordinate.M)) dimension++;

    if (dimension == 2 && this.isAllowOldJtsCoordinateSyntax) dimension++;

    return dimension;
  }

  /// Merges an array of one-coordinate-{@link CoordinateSequence}s into one {@link CoordinateSequence}.
  ///
  /// @param sequences an array of coordinate sequences. Each sequence contains <b>exactly one</b> coordinate.
  /// @param ordinateFlags a bit-mask of required ordinates.
  /// @return a coordinate sequence containing all coordinate
  CoordinateSequence mergeSequences(
      List<CoordinateSequence> sequences, List<Ordinate> ordinateFlags) {
    // if the sequences array is empty or null create an empty sequence
    if (sequences == null || sequences.isEmpty)
      return csFactory.createSizeDim(0, toDimension(ordinateFlags));

    if (sequences.length == 1) return sequences[0];

    List<Ordinate> mergeOrdinates;
    if (this.isAllowOldJtsCoordinateSyntax && ordinateFlags.length == 2) {
      mergeOrdinates = []
        ..addAll(ordinateFlags); // TODO check if this clone is enough
      for (int i = 0; i < sequences.length; i++) {
        if ((sequences[i]).hasZ()) {
          mergeOrdinates.add(Ordinate.Z);
          break;
        }
      }
    } else
      mergeOrdinates = ordinateFlags;

    // create and fill the result sequence // TODO check
    CoordinateSequence sequence = this.csFactory.createSizeDimMeas(
        sequences.length,
        toDimension(mergeOrdinates),
        mergeOrdinates.contains(Ordinate.M) ? 1 : 0);

    int offsetM =
        CoordinateSequence.Z + (mergeOrdinates.contains(Ordinate.Z) ? 1 : 0);
    for (int i = 0; i < sequences.length; i++) {
      CoordinateSequence item = sequences[i];
      sequence.setOrdinate(
          i, CoordinateSequence.X, item.getOrdinate(0, CoordinateSequence.X));
      sequence.setOrdinate(
          i, CoordinateSequence.Y, item.getOrdinate(0, CoordinateSequence.Y));
      if (mergeOrdinates.contains(Ordinate.Z)) {
        sequence.setOrdinate(
            i, CoordinateSequence.Z, item.getOrdinate(0, CoordinateSequence.Z));
      }
      if (mergeOrdinates.contains(Ordinate.M)) {
        sequence.setOrdinate(i, offsetM, item.getOrdinate(0, offsetM));
      }
    }

    // return it
    return sequence;
  }

  /// Returns the next array of <code>Coordinate</code>s in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next element returned by the stream should be L_PAREN (the
  ///      beginning of "(x1 y1, x2 y2, ..., xn yn)") or EMPTY.
  ///@return                  the next array of <code>Coordinate</code>s in the
  ///      stream, or an empty array if EMPTY is the next element returned by
  ///      the stream.
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  ///
  ///@deprecated in favor of functions returning {@link CoordinateSequence}s
  List<Coordinate> getCoordinates(WKTTokenizer tokenizer) {
    String nextToken = getNextEmptyOrOpener(tokenizer);
    if (nextToken == EMPTY) {
      return <Coordinate>[];
    }
    List<Coordinate> coordinates = [];
    coordinates.add(getPreciseCoordinate(tokenizer));
    nextToken = getNextCloserOrComma(tokenizer);
    while (nextToken == COMMA) {
      coordinates.add(getPreciseCoordinate(tokenizer));
      nextToken = getNextCloserOrComma(tokenizer);
    }
    return coordinates;
  }

  /// Returns the next array of <code>Coordinate</code>s in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next element returned by the stream should be a number.
  ///@return                  the next array of <code>Coordinate</code>s in the
  ///      stream.
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  ///
  ///@deprecated in favor of functions returning {@link CoordinateSequence}s
  List<Coordinate> getCoordinatesNoLeftParen(WKTTokenizer tokenizer) {
    List<Coordinate> coordinates = [];
    coordinates.add(getPreciseCoordinate(tokenizer));
    String nextToken = getNextCloserOrComma(tokenizer);
    while (nextToken == COMMA) {
      coordinates.add(getPreciseCoordinate(tokenizer));
      nextToken = getNextCloserOrComma(tokenizer);
    }
    return coordinates;
  }

  /// Returns the next precise <code>Coordinate</code> in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next element returned by the stream should be a number.
  ///@return                  the next array of <code>Coordinate</code>s in the
  ///      stream.
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  ///
  ///@deprecated in favor of functions returning {@link CoordinateSequence}s
  Coordinate getPreciseCoordinate(WKTTokenizer tokenizer) {
    Coordinate coord = Coordinate.empty2D();
    coord.x = getNextNumber(tokenizer);
    coord.y = getNextNumber(tokenizer);
    if (isNumberNext(tokenizer)) {
      coord.z = getNextNumber(tokenizer);
    }
    if (isNumberNext(tokenizer)) {
      getNextNumber(tokenizer); // ignore M value
    }
    precisionModel.makeCoordinatePrecise(coord);
    return coord;
  }

  /// Tests if the next token in the stream is a number
  ///
  /// @param tokenizer the tokenizer
  /// @return {@code true} if the next token is a number, otherwise {@code false}
  /// @throws  IOException     if an I/O error occurs
  static bool isNumberNext(WKTTokenizer tokenizer) {
    WKTToken token = tokenizer.next()!;

    tokenizer.pushBack();
    return token.type == WKTTokenType.KEYWORD; // TT_WORD;
  }

  /// Tests if the next token in the stream is a left opener ({@link #L_PAREN})
  ///
  /// @param tokenizer the tokenizer
  /// @return {@code true} if the next token is a {@link #L_PAREN}, otherwise {@code false}
  /// @throws  IOException     if an I/O error occurs
  static bool isOpenerNext(WKTTokenizer tokenizer) {
    WKTToken token = tokenizer.next()!;
    tokenizer.pushBack();
    return token.type == WKTTokenType.LPAREN;
  }

  /// Parses the next number in the stream.
  /// Numbers with exponents are handled.
  /// <tt>NaN</tt> values are handled correctly, and
  /// the case of the "NaN" symbol is not significant.
  ///
  /// @param  tokenizer        tokenizer over a stream of text in Well-known Text
  /// @return                  the next number in the stream
  /// @throws  ParseException  if the next token is not a valid number
  /// @throws  IOException     if an I/O error occurs
  double getNextNumber(WKTTokenizer tokenizer) {
    WKTToken token = tokenizer.next()!;
    switch (token.type) {
      case WKTTokenType.KEYWORD:
        {
          if (StringUtils.equalsIgnoreCase(token.value, NAN_SYMBOL)) {
            return double.nan;
          } else {
            try {
              return double.parse(token.value);
            } catch (ex) {
              throw ArgumentError("Invalid number: ${token.value}");
            }
          }
          break;
        }
      case WKTTokenType.NUMERIC_LITERAL:
        {
          try {
            return double.parse(token.value);
          } catch (ex) {
            throw ArgumentError("Invalid number: ${token.value}");
          }
        }
    }
    throw ArgumentError(
        "Expected number token but found ${token.type}: ${token.value}");
  }

  ///  Returns the next EMPTY or L_PAREN in the stream as uppercase text.
  ///
  ///@return                  the next EMPTY or L_PAREN in the stream as uppercase
  ///      text.
  ///@throws  ParseException  if the next token is not EMPTY or L_PAREN
  ///@throws  IOException     if an I/O error occurs
  /// @param  tokenizer        tokenizer over a stream of text in Well-known Text
  static String getNextEmptyOrOpener(WKTTokenizer tokenizer) {
    String nextWord = getNextWord(tokenizer);
    if (StringUtils.equalsIgnoreCase(nextWord, "Z")) {
//z = true;
      nextWord = getNextWord(tokenizer);
    } else if (StringUtils.equalsIgnoreCase(nextWord, "M")) {
//m = true;
      nextWord = getNextWord(tokenizer);
    } else if (StringUtils.equalsIgnoreCase(nextWord, "ZM")) {
//z = true;
//m = true;
      nextWord = getNextWord(tokenizer);
    }
    if (nextWord == EMPTY || nextWord == L_PAREN) {
      return nextWord;
    }
    throw ArgumentError(
        "Expected $EMPTY or $L_PAREN token but found $nextWord");
  }

  ///  Returns the next ordinate flag information in the stream as uppercase text.
  ///  This can be Z, M or ZM.
  ///
  ///@return                  the next EMPTY or L_PAREN in the stream as uppercase
  ///      text.
  ///@throws  ParseException  if the next token is not EMPTY or L_PAREN
  ///@throws  IOException     if an I/O error occurs
  /// @param  tokenizer        tokenizer over a stream of text in Well-known Text
  static List<Ordinate> getNextOrdinateFlags(WKTTokenizer tokenizer) {
    List<Ordinate> result = List.from([Ordinate.X, Ordinate.Y]);

    String nextWord = lookAheadWord(tokenizer).toUpperCase();
    if (StringUtils.equalsIgnoreCase(nextWord, "Z")) {
      tokenizer.next();
      result.add(Ordinate.Z);
    } else if (StringUtils.equalsIgnoreCase(nextWord, "M")) {
      tokenizer.next();
      result.add(Ordinate.M);
    } else if (StringUtils.equalsIgnoreCase(nextWord, "ZM")) {
      tokenizer.next();
      result.add(Ordinate.Z);
      result.add(Ordinate.M);
    }
    return result;
  }

  ///  Returns the next word in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next token must be a word.
  ///@return                  the next word in the stream as uppercase text
  ///@throws  ParseException  if the next token is not a word
  ///@throws  IOException     if an I/O error occurs
  static String lookAheadWord(WKTTokenizer tokenizer) {
    String nextWord = getNextWord(tokenizer);
    tokenizer.pushBack();
    return nextWord;
  }

  ///  Returns the next {@link #R_PAREN} or {@link #COMMA} in the stream.
  ///
  ///@return                  the next R_PAREN or COMMA in the stream
  ///@throws  ParseException  if the next token is not R_PAREN or COMMA
  ///@throws  IOException     if an I/O error occurs
  /// @param  tokenizer        tokenizer over a stream of text in Well-known Text
  static String getNextCloserOrComma(WKTTokenizer tokenizer) {
    String nextWord = getNextWord(tokenizer);
    if (nextWord == COMMA || nextWord == R_PAREN) {
      return nextWord;
    }
    throw ArgumentError(
        "Expected $COMMA or $R_PAREN token but found $nextWord");
  }

  ///  Returns the next {@link #R_PAREN} in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next token must be R_PAREN.
  ///@return                  the next R_PAREN in the stream
  ///@throws  ParseException  if the next token is not R_PAREN
  ///@throws  IOException     if an I/O error occurs
  String getNextCloser(WKTTokenizer tokenizer) {
    String nextWord = getNextWord(tokenizer);
    if (nextWord == R_PAREN) {
      return nextWord;
    }
    throw ArgumentError("Expected $R_PAREN token but found $nextWord");
  }

  ///  Returns the next word in the stream.
  ///
  ///@return                  the next word in the stream as uppercase text
  ///@throws  ParseException  if the next token is not a word
  ///@throws  IOException     if an I/O error occurs
  /// @param  tokenizer        tokenizer over a stream of text in Well-known Text
  static String getNextWord(WKTTokenizer tokenizer) {
    WKTToken token = tokenizer.next()!;
    switch (token.type) {
      case WKTTokenType.KEYWORD:
        String word = token.value;
        if (StringUtils.equalsIgnoreCase(word, EMPTY)) {
          return EMPTY;
        }
        return word;

      case WKTTokenType.LPAREN:
        return L_PAREN; //'('
      case WKTTokenType.RPAREN:
        return R_PAREN; //')'
      case WKTTokenType.COMMA:
        return COMMA; //','
    }
    throw ArgumentError(
        "Expected word token but found ${token.type}: ${token.value}");
  }

  /// Gets a description of the current token type
  /// @param tokenizer the tokenizer
  /// @return a description of the current token
  static String tokenString(WKTToken token) {
    switch (token.type) {
      case WKTTokenType.NUMERIC_LITERAL:
        return "<NUMBER>";
//      case WKTTokenType.TT_EOL: // TODO check
//        return "End-of-Line";
      case WKTTokenType.EOS:
        return "End-of-Stream";
      case WKTTokenType.KEYWORD:
        return "'${token.value}'";
    }
    return "'${token.type}'";
  }

  ///  Creates a <code>Geometry</code> using the next token in the stream.
  ///
  ///@return                  a <code>Geometry</code> specified by the next token
  ///      in the stream
  ///@throws  ParseException  if the coordinates used to create a <code>Polygon</code>
  ///      shell and holes do not form closed linestrings, or if an unexpected
  ///      token was encountered
  ///@throws  IOException     if an I/O error occurs
  /// @param  tokenizer        tokenizer over a stream of text in Well-known Text
  Geometry? readGeometryTaggedText(WKTTokenizer tokenizer) {
    String type;

    List<Ordinate> ordinateFlags = List.from([Ordinate.X, Ordinate.Y]);
    try {
      type = getNextWord(tokenizer).toUpperCase();
      if (type.endsWith("ZM")) {
        ordinateFlags.add(Ordinate.Z);
        ordinateFlags.add(Ordinate.M);
      } else if (type.endsWith("Z")) {
        ordinateFlags.add(Ordinate.Z);
      } else if (type.endsWith("M")) {
        ordinateFlags.add(Ordinate.M);
      }
    } catch (e) {
      return null;
    }

    return readGeometryTaggedTextWithOpts(tokenizer, type, ordinateFlags);
  }

  Geometry readGeometryTaggedTextWithOpts(
      WKTTokenizer tokenizer, String type, List<Ordinate> ordinateFlags) {
    if (ordinateFlags.length == 2) {
      ordinateFlags = getNextOrdinateFlags(tokenizer);
    }

    // if we can create a sequence with the required dimension everything is ok, otherwise
    // we need to take a different coordinate sequence factory.
    // It would be good to not have to try/catch this but if the CoordinateSequenceFactory
    // exposed a value indicating which min/max dimension it can handle or even an
    // ordinate bit-flag.
    try {
      csFactory.createSizeDimMeas(0, toDimension(ordinateFlags),
          ordinateFlags.contains(Ordinate.M) ? 1 : 0);
    } catch (e) {
      geometryFactory = GeometryFactory(geometryFactory.getPrecisionModel(),
          geometryFactory.getSRID(), csFactoryXYZM);
    }

    if (type.startsWith("POINT")) {
      return readPointText(tokenizer, ordinateFlags);
    } else if (type.startsWith("LINESTRING")) {
      return readLineStringText(tokenizer, ordinateFlags);
    } else if (type.startsWith("LINEARRING")) {
      return readLinearRingText(tokenizer, ordinateFlags);
    } else if (type.startsWith("POLYGON")) {
      return readPolygonText(tokenizer, ordinateFlags);
    } else if (type.startsWith("MULTIPOINT")) {
      return readMultiPointText(tokenizer, ordinateFlags);
    } else if (type.startsWith("MULTILINESTRING")) {
      return readMultiLineStringText(tokenizer, ordinateFlags);
    } else if (type.startsWith("MULTIPOLYGON")) {
      return readMultiPolygonText(tokenizer, ordinateFlags);
    } else if (type.startsWith("GEOMETRYCOLLECTION")) {
      return readGeometryCollectionText(tokenizer, ordinateFlags);
    }
    var currentToken = tokenizer.currentToken();
    throw ArgumentError(
        "Unknown geometry type: ${currentToken.type}  -> ${currentToken.value}");
  }

  ///  Creates a <code>Point</code> using the next token in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next tokens must form a &lt;Point Text&gt;.
  ///@return                  a <code>Point</code> specified by the next token in
  ///      the stream
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  Point readPointText(WKTTokenizer tokenizer, List<Ordinate> ordinateFlags) {
    Point point = geometryFactory
        .createPointSeq(getCoordinateSequence(tokenizer, ordinateFlags));
    return point;
  }

  ///  Creates a <code>LineString</code> using the next token in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next tokens must form a &lt;LineString Text&gt;.
  ///@return                  a <code>LineString</code> specified by the next
  ///      token in the stream
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  LineString readLineStringText(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags) {
    return geometryFactory
        .createLineStringSeq(getCoordinateSequence(tokenizer, ordinateFlags));
  }

  ///  Creates a <code>LinearRing</code> using the next token in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next tokens must form a &lt;LineString Text&gt;.
  ///@return                  a <code>LinearRing</code> specified by the next
  ///      token in the stream
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if the coordinates used to create the <code>LinearRing</code>
  ///      do not form a closed linestring, or if an unexpected token was
  ///      encountered
  LinearRing readLinearRingText(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags) {
    return geometryFactory
        .createLinearRingSeq(getCoordinateSequence(tokenizer, ordinateFlags));
  }

  ///  Creates a <code>MultiPoint</code> using the next tokens in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next tokens must form a &lt;MultiPoint Text&gt;.
  ///@return                  a <code>MultiPoint</code> specified by the next
  ///      token in the stream
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  MultiPoint readMultiPointText(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags) {
    return geometryFactory.createMultiPointSeq(getCoordinateSequenceTryParen(
        tokenizer, ordinateFlags, this.isAllowOldJtsMultipointSyntax));
  }

  ///  Creates a <code>Polygon</code> using the next token in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next tokens must form a &lt;Polygon Text&gt;.
  ///@return                  a <code>Polygon</code> specified by the next token
  ///      in the stream
  ///@throws  ParseException  if the coordinates used to create the <code>Polygon</code>
  ///      shell and holes do not form closed linestrings, or if an unexpected
  ///      token was encountered.
  ///@throws  IOException     if an I/O error occurs
  Polygon readPolygonText(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags) {
    String nextToken = getNextEmptyOrOpener(tokenizer);
    if (nextToken == EMPTY) {
      return geometryFactory.createPolygonEmpty();
    }
    List<LinearRing> holes = [];
    LinearRing shell = readLinearRingText(tokenizer, ordinateFlags);
    nextToken = getNextCloserOrComma(tokenizer);
    while (nextToken == COMMA) {
      LinearRing hole = readLinearRingText(tokenizer, ordinateFlags);
      holes.add(hole);
      nextToken = getNextCloserOrComma(tokenizer);
    }
    // List<LinearRing> array = [];//..length = (holes.length);
    return geometryFactory.createPolygon(shell, holes);
  }

  ///  Creates a <code>MultiLineString</code> using the next token in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next tokens must form a &lt;MultiLineString Text&gt;.
  ///@return                  a <code>MultiLineString</code> specified by the
  ///      next token in the stream
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  MultiLineString readMultiLineStringText(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags) {
    String nextToken = getNextEmptyOrOpener(tokenizer);
    if (nextToken == EMPTY) {
      return geometryFactory.createMultiLineStringEmpty();
    }

    List<LineString> lineStrings = [];
    do {
      LineString lineString = readLineStringText(tokenizer, ordinateFlags);
      lineStrings.add(lineString);
      nextToken = getNextCloserOrComma(tokenizer);
    } while (nextToken == COMMA);

    return geometryFactory.createMultiLineString(lineStrings);
  }

  ///  Creates a <code>MultiPolygon</code> using the next token in the stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next tokens must form a &lt;MultiPolygon Text&gt;.
  ///@return                  a <code>MultiPolygon</code> specified by the next
  ///      token in the stream, or if if the coordinates used to create the
  ///      <code>Polygon</code> shells and holes do not form closed linestrings.
  ///@throws  IOException     if an I/O error occurs
  ///@throws  ParseException  if an unexpected token was encountered
  MultiPolygon readMultiPolygonText(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags) {
    String nextToken = getNextEmptyOrOpener(tokenizer);
    if (nextToken == EMPTY) {
      return geometryFactory.createMultiPolygonEmpty();
    }
    List<Polygon> polygons = [];
    do {
      Polygon polygon = readPolygonText(tokenizer, ordinateFlags);
      polygons.add(polygon);
      nextToken = getNextCloserOrComma(tokenizer);
    } while (nextToken == COMMA);
    return geometryFactory.createMultiPolygon(polygons);
  }

  ///  Creates a <code>GeometryCollection</code> using the next token in the
  ///  stream.
  ///
  ///@param  tokenizer        tokenizer over a stream of text in Well-known Text
  ///      format. The next tokens must form a &lt;GeometryCollection Text&gt;.
  ///@return                  a <code>GeometryCollection</code> specified by the
  ///      next token in the stream
  ///@throws  ParseException  if the coordinates used to create a <code>Polygon</code>
  ///      shell and holes do not form closed linestrings, or if an unexpected
  ///      token was encountered
  ///@throws  IOException     if an I/O error occurs
  GeometryCollection readGeometryCollectionText(
      WKTTokenizer tokenizer, List<Ordinate> ordinateFlags) {
    String nextToken = getNextEmptyOrOpener(tokenizer);
    if (nextToken == EMPTY) {
      return geometryFactory.createGeometryCollectionEmpty();
    }
    List<Geometry> geometries = [];
    do {
      Geometry geometry = readGeometryTaggedText(tokenizer)!;
      geometries.add(geometry);
      nextToken = getNextCloserOrComma(tokenizer);
    } while (nextToken == COMMA);

    return geometryFactory.createGeometryCollection(geometries);
  }
}
