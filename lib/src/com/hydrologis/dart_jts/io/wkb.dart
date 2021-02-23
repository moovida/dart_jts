part of dart_jts;

/**
 * Constant values used by the WKB format
 */
class WKBConstants {
  static const wkbXDR = 0;
  static const wkbNDR = 1;

  static const wkbPoint = 1;
  static const wkbLineString = 2;
  static const wkbPolygon = 3;
  static const wkbMultiPoint = 4;
  static const wkbMultiLineString = 5;
  static const wkbMultiPolygon = 6;
  static const wkbGeometryCollection = 7;
}

/**
 * Allows reading a stream of Java primitive datatypes from an underlying
 * {@link InStream},
 * with the representation being in either common byte ordering.
 */
class ByteOrderDataInStream {
  Endian _byteOrder = Endian.big;
  Uint8List _bytesList;

  int _readOffset = 0;

  ByteOrderDataInStream(List<int> bytesList) {
    if (bytesList is Uint8List) {
      this._bytesList = bytesList;
    } else {
      this._bytesList = Uint8List.fromList(bytesList);
    }
  }

  int get readOffset => _readOffset;

//  /**
//   * Allows a single ByteOrderDataInStream to be reused
//   * on multiple InStreams.
//   *
//   * @param stream
//   */
//  void setInStream(Uint8List stream) {
//    this._stream = stream;
//  }

  void setOrder(Endian byteOrder) {
    this._byteOrder = byteOrder;
  }

  /**
   * Reads a byte value
   *
   * @return the byte read
   */
  int readByte() {
    int b = _bytesList[_readOffset++];
    return b;
  }

  int readInt() {
    var sublist = _bytesList.sublist(_readOffset, _readOffset + 4);
    var bdata = ByteData.view(sublist.buffer);
    var int32 = bdata.getInt32(0, _byteOrder);

    _readOffset += 4;
    return int32;
  }

  int readLong() {
    var sublist = _bytesList.sublist(_readOffset, _readOffset + 8);
    var bdata = ByteData.view(sublist.buffer);
    var int64 = bdata.getInt64(0, _byteOrder);
    _readOffset += 8;
    return int64;
  }

  double readDouble() {
    var sublist = _bytesList.sublist(_readOffset, _readOffset + 8);
    var bdata = ByteData.view(sublist.buffer);
    var float64 = bdata.getFloat64(0, _byteOrder);
    _readOffset += 8;
    return float64;
  }
}

/**
 * Reads a {@link Geometry}from a byte stream in Well-Known Binary format.
 * Supports use of an {@link InStream}, which allows easy use
 * with arbitrary byte stream sources.
 * <p>
 * This class reads the format describe in {@link WKBWriter}.
 * It also partially handles
 * the <b>Extended WKB</b> format used by PostGIS,
 * by parsing and storing SRID values.
 * The reader repairs structurally-invalid input
 * (specifically, LineStrings and LinearRings which contain
 * too few points have vertices added,
 * and non-closed rings are closed).
 * <p>
 * This class is designed to support reuse of a single instance to read multiple
 * geometries. This class is not thread-safe; each thread should create its own
 * instance.
 * <p>
 * As of version 1.15, the reader can read geometries following OGC 06-103r4
 * speification used by Spatialite/Geopackage.
 * <p>
 * The difference between PostGIS EWKB format and the new OGC specification is
 * that Z and M coordinates are detected with a bit mask on the higher byte in
 * the former case (0x80 for Z and 0x40 for M) while new OGC specification use
 * specif int ranges for 2D gemetries, Z geometries (2D code+1000), M geometries
 * (2D code+2000) and ZM geometries (2D code+3000).
 * <p>
 * Note that the {@link WKBWriter} is not changed and still write PostGIS WKB
 * geometries
 * @see WKBWriter for a formal format specification
 */
class WKBReader {
  /**
   * Converts a hexadecimal string to a byte array.
   * The hexadecimal digit symbols are case-insensitive.
   *
   * @param hex a string containing hex digits
   * @return an array of bytes with the value of the hex string
   */
  static List<int> hexToBytes(String hex) {
    List<int> bytes = HEX.decode(hex);
    return bytes;
//
//    int byteLen = hex.length ~/ 2;
//    List<int> bytes = List(byteLen);
//
//    for (int i = 0; i < byteLen; i++) {
//      int i2 = 2 * i;
//      if (i2 + 1 > hex.length)
//        throw ArgumentError("Hex string has odd length");
//
//      int nib1 = hexToInt(hex.charAt(i2));
//      int nib0 = hexToInt(hex.charAt(i2 + 1));
//      int b = (nib1 << 4) + nib0;
//      bytes[i] = b;
//    }
//    return bytes;
  }

//  static int hexToInt(String hex) {
//    int nib = Character.digit(hex, 16);
//    if (nib < 0)
//      throw new IllegalArgumentException("Invalid hex digit: '" + hex + "'");
//    return nib;
//  }

  static final String INVALID_GEOM_TYPE_MSG =
      "Invalid geometry type encountered in ";

  GeometryFactory factory;
  CoordinateSequenceFactory csFactory;
  PrecisionModel precisionModel;

  // default dimension - will be set on read
  int inputDimension = 2;
  bool hasSRID = false;
  int SRID = 0;

  /**
   * true if structurally invalid input should be reported rather than repaired.
   * At some point this could be made client-controllable.
   */
  bool isStrict = false;
  List<double> ordValues;

  ByteOrderDataInStream dis;

  WKBReader() : this.withFactory(new GeometryFactory.defaultPrecision());

  WKBReader.withFactory(GeometryFactory geometryFactory) {
    this.factory = geometryFactory;
    precisionModel = factory.getPrecisionModel();
    csFactory = factory.getCoordinateSequenceFactory();
  }

//  /**
//   * Reads a single {@link Geometry} in WKB format from a byte array.
//   *
//   * @param bytes the byte array to read from
//   * @return the geometry read
//   * @throws ParseException if the WKB is ill-formed
//   */
//  Geometry read(List<int> bytes) {
//    // possibly reuse the ByteArrayInStream?
//    // don't throw IOExceptions, since we are not doing any I/O
//    try {
//      return read(new ByteArrayInStream(bytes));
//    } catch (ex) {
//      throw IOException("Unexpected IOException caught: ${ex.getMessage()}");
//    }
//  }

  /**
   * Reads a {@link Geometry} in binary WKB format from an {@link InStream}.
   *
   * @param is the stream to read from
   * @return the Geometry read
   * @throws IOException if the underlying stream creates an error
   * @throws ParseException if the WKB is ill-formed
   */
  Geometry read(List<int> ins) {
    var bytesList = ins;
    if (!(ins is Uint8List)) {
      bytesList = Uint8List.fromList(ins);
    }
    dis = ByteOrderDataInStream(bytesList);
    Geometry g = readGeometry();
    return g;
  }

  Geometry readGeometry() {
// determine byte order
    int byteOrderWKB = dis.readByte();

// always set byte order, since it may change from geometry to geometry
    if (byteOrderWKB == WKBConstants.wkbNDR) {
      dis.setOrder(Endian.little);
    } else if (byteOrderWKB == WKBConstants.wkbXDR) {
      dis.setOrder(Endian.big);
    } else if (isStrict) {
      throw new ParseException(
          "Unknown geometry byte order (not NDR or XDR): $byteOrderWKB");
    }
//if not strict and not XDR or NDR, then we just use the dis default set at the
//start of the geometry (if a multi-geometry).  This  allows WBKReader to work
//with Spatialite native BLOB WKB, as well as other WKB variants that might just
//specify endian-ness at the start of the multigeometry.

    int typeInt = dis.readInt();
// Adds %1000 to make it compatible with OGC 06-103r4
    int geometryType = (typeInt & 0xffff) % 1000;

// handle 3D and 4D WKB geometries
// geometries with Z coordinates have the 0x80 flag (postgis EWKB)
// or are in the 1000 range (Z) or in the 3000 range (ZM) of geometry type (OGC 06-103r4)
    bool hasZ = ((typeInt & 0x80000000) != 0 ||
        (typeInt & 0xffff) / 1000 == 1 ||
        (typeInt & 0xffff) / 1000 == 3);
// geometries with M coordinates have the 0x40 flag (postgis EWKB)
// or are in the 1000 range (M) or in the 3000 range (ZM) of geometry type (OGC 06-103r4)
    bool hasM = ((typeInt & 0x40000000) != 0 ||
        (typeInt & 0xffff) / 1000 == 2 ||
        (typeInt & 0xffff) / 1000 == 3);
//System.out.println(typeInt + " - " + geometryType + " - hasZ:" + hasZ);
    inputDimension = 2 + (hasZ ? 1 : 0) + (hasM ? 1 : 0);

// determine if SRIDs are present
    hasSRID = (typeInt & 0x20000000) != 0;
    int SRID = 0;
    if (hasSRID) {
      SRID = dis.readInt();
    }

// only allocate ordValues buffer if necessary
    if (ordValues == null || ordValues.length < inputDimension)
      ordValues = List(inputDimension);

    Geometry geom = null;
    switch (geometryType) {
      case WKBConstants.wkbPoint:
        geom = readPoint();
        break;
      case WKBConstants.wkbLineString:
        geom = readLineString();
        break;
      case WKBConstants.wkbPolygon:
        geom = readPolygon();
        break;
      case WKBConstants.wkbMultiPoint:
        geom = readMultiPoint();
        break;
      case WKBConstants.wkbMultiLineString:
        geom = readMultiLineString();
        break;
      case WKBConstants.wkbMultiPolygon:
        geom = readMultiPolygon();
        break;
      case WKBConstants.wkbGeometryCollection:
        geom = readGeometryCollection();
        break;
      default:
        throw new ParseException("Unknown WKB type $geometryType");
    }
    setSRID(geom, SRID);
    return geom;
  }

  /**
   * Sets the SRID, if it was specified in the WKB
   *
   * @param g the geometry to update
   * @return the geometry with an updated SRID value, if required
   */
  Geometry setSRID(Geometry g, int SRID) {
    if (SRID != 0) g.setSRID(SRID);
    return g;
  }

  Point readPoint() {
    CoordinateSequence pts = readCoordinateSequence(1);
    return factory.createPointSeq(pts);
  }

  LineString readLineString() {
    int size = dis.readInt();
    CoordinateSequence pts = readCoordinateSequenceLineString(size);
    return factory.createLineStringSeq(pts);
  }

  LinearRing readLinearRing() {
    int size = dis.readInt();
    CoordinateSequence pts = readCoordinateSequenceRing(size);
    return factory.createLinearRingSeq(pts);
  }

  Polygon readPolygon() {
    int numRings = dis.readInt();
    List<LinearRing> holes = null;
    if (numRings > 1) holes = List(numRings - 1);

    LinearRing shell = readLinearRing();
    for (int i = 0; i < numRings - 1; i++) {
      holes[i] = readLinearRing();
    }
    return factory.createPolygon(shell, holes);
  }

  MultiPoint readMultiPoint() {
    int numGeom = dis.readInt();
    List<Point> geoms = List(numGeom);
    for (int i = 0; i < numGeom; i++) {
      Geometry g = readGeometry();
      if (!(g is Point))
        throw new ParseException(INVALID_GEOM_TYPE_MSG + "MultiPoint");
      geoms[i] = g;
    }
    return factory.createMultiPoint(geoms);
  }

  MultiLineString readMultiLineString() {
    int numGeom = dis.readInt();
    List<LineString> geoms = List(numGeom);
    for (int i = 0; i < numGeom; i++) {
      Geometry g = readGeometry();
      if (!(g is LineString))
        throw new ParseException(INVALID_GEOM_TYPE_MSG + "MultiLineString");
      geoms[i] = g;
    }
    return factory.createMultiLineString(geoms);
  }

  MultiPolygon readMultiPolygon() {
    int numGeom = dis.readInt();
    List<Polygon> geoms = List(numGeom);

    for (int i = 0; i < numGeom; i++) {
      Geometry g = readGeometry();
      if (!(g is Polygon))
        throw new ParseException(INVALID_GEOM_TYPE_MSG + "MultiPolygon");
      geoms[i] = g;
    }
    return factory.createMultiPolygon(geoms);
  }

  GeometryCollection readGeometryCollection() {
    int numGeom = dis.readInt();
    List<Geometry> geoms = List(numGeom);
    for (int i = 0; i < numGeom; i++) {
      geoms[i] = readGeometry();
    }
    return factory.createGeometryCollection(geoms);
  }

  CoordinateSequence readCoordinateSequence(int size) {
    CoordinateSequence seq = csFactory.createSizeDim(size, inputDimension);
    int targetDim = seq.getDimension();
    if (targetDim > inputDimension) targetDim = inputDimension;
    for (int i = 0; i < size; i++) {
      readCoordinate();
      for (int j = 0; j < targetDim; j++) {
        seq.setOrdinate(i, j, ordValues[j]);
      }
    }
    return seq;
  }

  CoordinateSequence readCoordinateSequenceLineString(int size) {
    CoordinateSequence seq = readCoordinateSequence(size);
    if (isStrict) return seq;
    if (seq.size() == 0 || seq.size() >= 2) return seq;
    return CoordinateSequences.extend(csFactory, seq, 2);
  }

  CoordinateSequence readCoordinateSequenceRing(int size) {
    CoordinateSequence seq = readCoordinateSequence(size);
    if (isStrict) return seq;
    if (CoordinateSequences.isRing(seq)) return seq;
    return CoordinateSequences.ensureValidRing(csFactory, seq);
  }

  /**
   * Reads a coordinate value with the specified dimensionality.
   * Makes the X and Y ordinates precise according to the precision model
   * in use.
   */
  void readCoordinate() {
    for (int i = 0; i < inputDimension; i++) {
      if (i <= 1) {
        ordValues[i] = precisionModel.makePrecise(dis.readDouble());
      } else {
        ordValues[i] = dis.readDouble();
      }
    }
  }
}

/**
 * Writes a {@link Geometry} into Well-Known Binary format.
 * Supports use of an {@link OutStream}, which allows easy use
 * with arbitrary byte stream sinks.
 * <p>
 * The WKB format is specified in the
 * OGC <A HREF="http://www.opengis.org/techno/specs.htm"><i>Simple Features for SQL</i></a>
 * specification.
 * This implementation also supports the <b>Extended WKB</b>
 * standard. Extended WKB allows writing 3-dimensional coordinates
 * and including the geometry SRID value.
 * The presence of 3D coordinates is signified
 * by setting the high bit of the <tt>wkbType</tt> word.
 * The presence of an SRID is signified
 * by setting the third bit of the <tt>wkbType</tt> word.
 * EWKB format is upward compatible with the original SFS WKB format.
 * <p>
 * Empty Points cannot be represented in WKB; an
 * {@link IllegalArgumentException} will be thrown if one is
 * written.
 * <p>
 * The WKB specification does not support representing {@link LinearRing}s;
 * they will be written as {@link LineString}s.
 * <p>
 * This class is designed to support reuse of a single instance to read multiple
 * geometries. This class is not thread-safe; each thread should create its own
 * instance.
 *
 * <h3>Syntax</h3>
 * The following syntax specification describes the version of Well-Known Binary
 * supported by JTS.
 * <p>
 * <i>The specification uses a syntax language similar to that used in
 * the C language.  Bitfields are specified from hi-order to lo-order bits.</i>
 * <p>
 * <blockquote><pre>
 *
 * <b>byte</b> = 1 byte
 * <b>uint32</b> = 32 bit unsigned integer (4 bytes)
 * <b>double</b> = double precision number (8 bytes)
 *
 * abstract Point { }
 *
 * Point2D extends Point {
 * 	<b>double</b> x;
 * 	<b>double</b> y;
 * }
 *
 * Point3D extends Point {
 * 	<b>double</b> x;
 * 	<b>double</b> y;
 * 	<b>double</b> z;
 * }
 *
 * LinearRing {
 * 	<b>uint32</b> numPoints;
 * 	Point points[numPoints];
 * }
 *
 * enum wkbGeometryType {
 * 	wkbPoint = 1,
 * 	wkbLineString = 2,
 * 	wkbPolygon = 3,
 * 	wkbMultiPoint = 4,
 * 	wkbMultiLineString = 5,
 * 	wkbMultiPolygon = 6,
 * 	wkbGeometryCollection = 7
 * }
 *
 * enum byteOrder {
 * 	wkbXDR = 0,	// Big Endian
 * 	wkbNDR = 1 	// Little Endian
 * }
 *
 * WKBType {
 * 	<b>uint32</b> wkbGeometryType : 8; // values from enum wkbGeometryType
 * }
 *
 * EWKBType {
 * 	<b>uint32</b> is3D : 1; 	// 0 = 2D, 1 = 3D
 * 	<b>uint32</b> noData1 : 1;
 * 	<b>uint32</b> hasSRID : 1;  	// 0, no, 1 = yes
 * 	<b>uint32</b> noData2 : 21;
 * 	<b>uint32</b> wkbGeometryType : 8; // values from enum wkbGeometryType
 * }
 *
 * abstract WKBGeometry {
 * 	<b>byte</b> byteOrder;		// values from enum byteOrder
 * 	EWKBType wkbType
 * 	[ <b>uint32</b> srid; ] 	// only if hasSRID = yes
 * }
 *
 * WKBPoint extends WKBGeometry {
 * 	Point point;
 * }
 *
 * WKBLineString extends WKBGeometry {
 * 	<b>uint32</b> numCoords;
 * 	Point points[numCoords];
 * }
 *
 * WKBPolygon extends WKBGeometry {
 * 	<b>uint32</b> numRings;
 * 	LinearRing rings[numRings];
 * }
 *
 * WKBMultiPoint extends WKBGeometry {
 * 	<b>uint32</b> numElems;
 * 	WKBPoint elems[numElems];
 * }
 *
 * WKBMultiLineString extends WKBGeometry {
 * 	<b>uint32</b> numElems;
 * 	WKBLineString elems[numElems];
 * }
 *
 * wkbMultiPolygon extends WKBGeometry {
 * 	<b>uint32</b> numElems;
 * 	WKBPolygon elems[numElems];
 * }
 *
 * WKBGeometryCollection extends WKBGeometry {
 * 	<b>uint32</b> numElems;
 * 	WKBGeometry elems[numElems];
 * }
 *
 * </pre></blockquote>
 * @see WKBReader
 */
class WKBWriter {
  /**
   * Converts a byte array to a hexadecimal string.
   *
   * @param bytes a byte array
   * @return a string of hexadecimal digits
   */
  static String toHex(List<int> bytes) {
    if (!(bytes is Uint8List)) {
      bytes = Uint8List.fromList(bytes);
    }
    return HEX.encode(bytes);
//    StringBuffer buf = new StringBuffer();
//    for (int i = 0; i < bytes.length; i++) {
//      int b = bytes[i];
//      buf.write(toHexDigit((b >> 4) & 0x0F));
//      buf.write(toHexDigit(b & 0x0F));
//    }
//    return buf.toString();
//  }
//
//
//  static String toHexDigit(int n) {
//    if (n < 0 || n > 15) throw new ArgumentError("Nibble value out of range: $n");
//    if (n <= 9) return '0$n';
//    return 'A${n - 10}';
  }

  int outputDimension = 2;
  Endian byteOrder;
  bool includeSRID = false;
  List<int> byteArrayOutStream;

  // holds output data values
  List<int> buf = List.from([0, 0, 0, 0, 0, 0, 0, 0]);

  /**
   * Creates a writer that writes {@link Geometry}s with
   * output dimension = 2 and Endian.big byte order
   */
  WKBWriter() : this.withDimOrder(2, Endian.big);

  /**
   * Creates a writer that writes {@link Geometry}s with
   * the given dimension (2 or 3) for output coordinates
   * and {@link ByteOrderValues#Endian.big} byte order.
   * If the input geometry has a small coordinate dimension,
   * coordinates will be padded with {@link Coordinate#NULL_ORDINATE}.
   *
   * @param outputDimension the coordinate dimension to output (2 or 3)
   */
  WKBWriter.withDim(int outputDimension)
      : this.withDimOrder(outputDimension, Endian.big);

  /**
   * Creates a writer that writes {@link Geometry}s with
   * the given dimension (2 or 3) for output coordinates
   * and {@link ByteOrderValues#Endian.big} byte order. This constructor also
   * takes a flag to control whether srid information will be
   * written.
   * If the input geometry has a smaller coordinate dimension,
   * coordinates will be padded with {@link Coordinate#NULL_ORDINATE}.
   *
   * @param outputDimension the coordinate dimension to output (2 or 3)
   * @param includeSRID indicates whether SRID should be written
   */
  WKBWriter.withDimSrid(int outputDimension, bool includeSRID)
      : this.withDimOrderSrid(outputDimension, Endian.big, includeSRID);

  /**
   * Creates a writer that writes {@link Geometry}s with
   * the given dimension (2 or 3) for output coordinates
   * and byte order
   * If the input geometry has a small coordinate dimension,
   * coordinates will be padded with {@link Coordinate#NULL_ORDINATE}.
   *
   * @param outputDimension the coordinate dimension to output (2 or 3)
   * @param byteOrder the byte ordering to use
   */
  WKBWriter.withDimOrder(int outputDimension, Endian byteOrder)
      : this.withDimOrderSrid(outputDimension, byteOrder, false);

  /**
   * Creates a writer that writes {@link Geometry}s with
   * the given dimension (2 or 3) for output coordinates
   * and byte order. This constructor also takes a flag to
   * control whether srid information will be written.
   * If the input geometry has a small coordinate dimension,
   * coordinates will be padded with {@link Coordinate#NULL_ORDINATE}.
   *
   * @param outputDimension the coordinate dimension to output (2 or 3)
   * @param byteOrder the byte ordering to use
   * @param includeSRID indicates whether SRID should be written
   */
  WKBWriter.withDimOrderSrid(
      int outputDimension, Endian byteOrder, bool includeSRID) {
    this.outputDimension = outputDimension;
    this.byteOrder = byteOrder;
    this.includeSRID = includeSRID;

    if (outputDimension < 2 || outputDimension > 3)
      throw ArgumentError("Output dimension must be 2 or 3");
  }

  /**
   * Writes a {@link Geometry} into a byte array.
   *
   * @param geom the geometry to write
   * @return the byte array containing the WKB
   */
  List<int> write(Geometry geom) {
    try {
      byteArrayOutStream = [];
      writeStream(geom, byteArrayOutStream);
    } catch (ex) {
      throw RuntimeException("Unexpected IO exception: $ex");
    }
    return byteArrayOutStream;
  }

  /**
   * Writes a {@link Geometry} to an {@link OutStream}.
   *
   * @param geom the geometry to write
   * @param os the out stream to write to
   * @throws IOException if an I/O error occurs
   */
  void writeStream(Geometry geom, List<int> os) {
    if (geom is Point)
      writePoint(geom, os);
    // LinearRings will be written as LineStrings
    else if (geom is LineString)
      writeLineString(geom, os);
    else if (geom is Polygon)
      writePolygon(geom, os);
    else if (geom is MultiPoint)
      writeGeometryCollection(WKBConstants.wkbMultiPoint, geom, os);
    else if (geom is MultiLineString)
      writeGeometryCollection(WKBConstants.wkbMultiLineString, geom, os);
    else if (geom is MultiPolygon)
      writeGeometryCollection(WKBConstants.wkbMultiPolygon, geom, os);
    else if (geom is GeometryCollection)
      writeGeometryCollection(WKBConstants.wkbGeometryCollection, geom, os);
    else {
      Assert.shouldNeverReachHere("Unknown Geometry type");
    }
  }

  void writePoint(Point pt, List<int> os) {
    if (pt.getCoordinateSequence().size() == 0)
      throw ArgumentError("Empty Points cannot be represented in WKB");
    writeByteOrder(os);
    writeGeometryType(WKBConstants.wkbPoint, pt, os);
    writeCoordinateSequence(pt.getCoordinateSequence(), false, os);
  }

  void writeLineString(LineString line, List<int> os) {
    writeByteOrder(os);
    writeGeometryType(WKBConstants.wkbLineString, line, os);
    writeCoordinateSequence(line.getCoordinateSequence(), true, os);
  }

  void writePolygon(Polygon poly, List<int> os) {
    writeByteOrder(os);
    writeGeometryType(WKBConstants.wkbPolygon, poly, os);
    writeInt(poly.getNumInteriorRing() + 1, os);
    writeCoordinateSequence(
        poly.getExteriorRing().getCoordinateSequence(), true, os);
    for (int i = 0; i < poly.getNumInteriorRing(); i++) {
      writeCoordinateSequence(
          poly.getInteriorRingN(i).getCoordinateSequence(), true, os);
    }
  }

  void writeGeometryCollection(
      int geometryType, GeometryCollection gc, List<int> os) {
    writeByteOrder(os);
    writeGeometryType(geometryType, gc, os);
    writeInt(gc.getNumGeometries(), os);
    for (int i = 0; i < gc.getNumGeometries(); i++) {
      writeStream(gc.getGeometryN(i), os);
    }
  }

  void writeByteOrder(List<int> os) {
    if (byteOrder == Endian.little)
      os.add(WKBConstants.wkbNDR);
    else
      os.add(WKBConstants.wkbXDR);
  }

  void writeGeometryType(int geometryType, Geometry g, List<int> os) {
    int flag3D = (outputDimension == 3) ? 0x80000000 : 0;
    int typeInt = geometryType | flag3D;
    typeInt |= includeSRID ? 0x20000000 : 0;
    writeInt(typeInt, os);
    if (includeSRID) {
      writeInt(g.getSRID(), os);
    }
  }

  void writeInt(int intValue, List<int> os) {
    Byteutils.putInt32(intValue, buf, byteOrder);
    os.addAll(buf.sublist(0, 4));
  }

  void writeCoordinateSequence(
      CoordinateSequence seq, bool writeSize, List<int> os) {
    if (writeSize) writeInt(seq.size(), os);

    for (int i = 0; i < seq.size(); i++) {
      writeCoordinate(seq, i, os);
    }
  }

  void writeCoordinate(CoordinateSequence seq, int index, List<int> os) {
    Byteutils.putFloat64(seq.getX(index), buf, byteOrder);
    os.addAll(buf);
    Byteutils.putFloat64(seq.getY(index), buf, byteOrder);
    os.addAll(buf);

// only write 3rd dim if caller has requested it for this writer
    if (outputDimension >= 3) {
// if 3rd dim is requested, only write it if the CoordinateSequence provides it
      double ordVal = Coordinate.NULL_ORDINATE;
      if (seq.getDimension() >= 3) ordVal = seq.getOrdinate(index, 2);
      Byteutils.putFloat64(ordVal, buf, byteOrder);
      os.addAll(buf);
    }
  }
}

///**
// * Methods to read and write primitive Java datatypes from/to byte
// * sequences, allowing the byte order to be specified
// * <p>
// * Similar to the standard Java <code>ByteBuffer</code> class.
// */
//class ByteOrderValues {
//  static int getInt(List<int> buf, Endian byteOrder) {
//    if (byteOrder == Endian.big) {
//      return ((buf[0] & 0xff) << 24) | ((buf[1] & 0xff) << 16) | ((buf[2] & 0xff) << 8) | ((buf[3] & 0xff));
//    } else {
//      // LITTLE_ENDIAN
//      return ((buf[3] & 0xff) << 24) | ((buf[2] & 0xff) << 16) | ((buf[1] & 0xff) << 8) | ((buf[0] & 0xff));
//    }
//  }
//
//  static void putInt(int intValue, List<int> buf, Endian byteOrder) {
//    if (byteOrder == Endian.big) {
//      buf[0] = (intValue >> 24);
//      buf[1] = (intValue >> 16);
//      buf[2] = (intValue >> 8);
//      buf[3] = intValue;
//    } else {
//      // LITTLE_ENDIAN
//      buf[0] = intValue;
//      buf[1] = (intValue >> 8);
//      buf[2] = (intValue >> 16);
//      buf[3] = (intValue >> 24);
//    }
//  }
//
//  static int getLong(List<int> buf, Endian byteOrder) {
//    if (byteOrder == Endian.big) {
//      return (buf[0] & 0xff) << 56 |
//          (buf[1] & 0xff) << 48 |
//          (buf[2] & 0xff) << 40 |
//          (buf[3] & 0xff) << 32 |
//          (buf[4] & 0xff) << 24 |
//          (buf[5] & 0xff) << 16 |
//          (buf[6] & 0xff) << 8 |
//          (buf[7] & 0xff);
//    } else {
//      // LITTLE_ENDIAN
//      return (buf[7] & 0xff) << 56 |
//          (buf[6] & 0xff) << 48 |
//          (buf[5] & 0xff) << 40 |
//          (buf[4] & 0xff) << 32 |
//          (buf[3] & 0xff) << 24 |
//          (buf[2] & 0xff) << 16 |
//          (buf[1] & 0xff) << 8 |
//          (buf[0] & 0xff);
//    }
//  }
//
//  static void putLong(int longValue, List<int> buf, Endian byteOrder) {
//    if (byteOrder == Endian.big) {
//      buf[0] = (longValue >> 56);
//      buf[1] = (longValue >> 48);
//      buf[2] = (longValue >> 40);
//      buf[3] = (longValue >> 32);
//      buf[4] = (longValue >> 24);
//      buf[5] = (longValue >> 16);
//      buf[6] = (longValue >> 8);
//      buf[7] = longValue;
//    } else {
//      // LITTLE_ENDIAN
//      buf[0] = longValue;
//      buf[1] = (longValue >> 8);
//      buf[2] = (longValue >> 16);
//      buf[3] = (longValue >> 24);
//      buf[4] = (longValue >> 32);
//      buf[5] = (longValue >> 40);
//      buf[6] = (longValue >> 48);
//      buf[7] = (longValue >> 56);
//    }
//  }
//
//  static double getDouble(List<int> buf, Endian byteOrder) {
//    var bdata = new ByteData.view(Uint8List.fromList(buf).buffer);
//    return bdata.getFloat64(0, byteOrder);
////    int longVal = getLong(buf, byteOrder);
////    return Double.longBitsToDouble(longVal);
//  }
//
//  static void putDouble(double doubleValue, List<int> buf, Endian byteOrder) {
//    Float64List float64List = Float64List.fromList([doubleValue]);
//    Uint8List bytes = new Uint8List.view(float64List.buffer);
//    buf.setAll(0, bytes);
////    long longVal = Double.doubleToLongBits(doubleValue);
////    putLong(longVal, buf, byteOrder);
//  }
//}
