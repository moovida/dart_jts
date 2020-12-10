part of dart_jts;

/// This file contains code that has not been ported
/// from the JTS project.
/// Geometry types used by the utility.
///
/// @author Andrea Antonello (www.hydrologis.com)
class EGeometryType {
//  static const NONE = const EGeometryType._(0, 0);
  static const POINT = const EGeometryType._(Point, MultiPoint, "Point");
  static const MULTIPOINT =
      const EGeometryType._(MultiPoint, MultiPoint, "MultiPoint");
  static const LINESTRING =
      const EGeometryType._(LineString, MultiLineString, "LineString");
  static const MULTILINESTRING = const EGeometryType._(
      MultiLineString, MultiLineString, "MultiLineString");
  static const POLYGON =
      const EGeometryType._(Polygon, MultiPolygon, "Polygon");
  static const MULTIPOLYGON =
      const EGeometryType._(MultiPolygon, MultiPolygon, "MultiPolygon");
  static const GEOMETRYCOLLECTION = const EGeometryType._(
      GeometryCollection, GeometryCollection, "GeometryCollection");
  static const GEOMETRY = const EGeometryType._(Geometry, Geometry, "GEOMETRY");
  static const UNKNOWN = const EGeometryType._(null, null, "Unknown");

  final clazz;
  final multiClazz;
  final String typeName;

  const EGeometryType._(this.clazz, this.multiClazz, this.typeName);

  dynamic getClazz() {
    return clazz;
  }

  dynamic getMultiClazz() {
    return multiClazz;
  }

  bool isMulti() {
    switch (this) {
      case MULTILINESTRING:
      case MULTIPOINT:
      case MULTIPOLYGON:
        return true;
      default:
        return false;
    }
  }

  bool isPoint() {
    switch (this) {
      case MULTIPOINT:
      case POINT:
        return true;
      default:
        return false;
    }
  }

  bool isLine() {
    switch (this) {
      case MULTILINESTRING:
      case LINESTRING:
        return true;
      default:
        return false;
    }
  }

  bool isPolygon() {
    switch (this) {
      case MULTIPOLYGON:
      case POLYGON:
        return true;
      default:
        return false;
    }
  }

  bool isCompatibleWith(EGeometryType geometryType) {
    switch (geometryType) {
      case LINESTRING:
        return this == LINESTRING;
      case MULTILINESTRING:
        return this == LINESTRING || this == MULTILINESTRING;
      case POINT:
        return this == POINT;
      case MULTIPOINT:
        return this == POINT || this == MULTIPOINT;
      case POLYGON:
        return this == POLYGON;
      case MULTIPOLYGON:
        return this == POLYGON || this == MULTIPOLYGON;
      default:
        return false;
    }
  }

  /// Returns the [EGeometryType] for a given [geometry].
  static EGeometryType forGeometry(Geometry geometry) {
    if (geometry is LineString) {
      return EGeometryType.LINESTRING;
    } else if (geometry is MultiLineString) {
      return EGeometryType.MULTILINESTRING;
    } else if (geometry is Point) {
      return EGeometryType.POINT;
    } else if (geometry is MultiPoint) {
      return EGeometryType.MULTIPOINT;
    } else if (geometry is Polygon) {
      return EGeometryType.POLYGON;
    } else if (geometry is MultiPolygon) {
      return EGeometryType.MULTIPOLYGON;
    } else if (geometry is GeometryCollection) {
      return EGeometryType.GEOMETRYCOLLECTION;
    } else {
      return EGeometryType.GEOMETRY;
    }
  }

  /// Returns the [EGeometryType] for a given [wktName].
  static EGeometryType forWktName(String wktName) {
    if (StringUtils.equalsIgnoreCase(wktName, POINT.getTypeName())) {
      return POINT;
    } else if (StringUtils.equalsIgnoreCase(
        wktName, MULTIPOINT.getTypeName())) {
      return MULTIPOINT;
    } else if (StringUtils.equalsIgnoreCase(
        wktName, LINESTRING.getTypeName())) {
      return LINESTRING;
    } else if (StringUtils.equalsIgnoreCase(
        wktName, MULTILINESTRING.getTypeName())) {
      return MULTILINESTRING;
    } else if (StringUtils.equalsIgnoreCase(wktName, POLYGON.getTypeName())) {
      return POLYGON;
    } else if (StringUtils.equalsIgnoreCase(
        wktName, MULTIPOLYGON.getTypeName())) {
      return MULTIPOLYGON;
    } else if (StringUtils.equalsIgnoreCase(
        wktName, GEOMETRYCOLLECTION.getTypeName())) {
      return GEOMETRYCOLLECTION;
    } else if (StringUtils.equalsIgnoreCase(wktName, GEOMETRY.getTypeName())) {
      return GEOMETRY;
    }
    return UNKNOWN;
  }

  static EGeometryType forTypeName(String typeName) {
    return forWktName(typeName);
  }

  /// Checks if the given [geometry] is a [LineString] (or [MultiLineString]) geometry.
  static bool isGeomLine(Geometry geometry) {
    if (geometry is LineString || geometry is MultiLineString) {
      return true;
    }
    return false;
  }

  /// Checks if the given [geometry] is a [Polygon] (or [MultiPolygon]) geometry.
  static bool isGeomPolygon(Geometry geometry) {
    if (geometry is Polygon || geometry is MultiPolygon) {
      return true;
    }
    return false;
  }

  /// Checks if the given [geometry] is a [Point] (or [MultiPoint]) geometry.
  static bool isGeomPoint(Geometry geometry) {
    if (geometry is Point || geometry is MultiPoint) {
      return true;
    }
    return false;
  }

  int getTypeCode([Geometry geometry]) {
    bool hasZ = false;
    bool hasM = false;
    if (geometry != null) {
      var coord = geometry.getCoordinate();
      hasZ = coord.z != null && !coord.z.isNaN;
      hasM = coord.getM() != null && !coord.getM().isNaN;
    }

    if (hasZ && !hasM) {
      switch (this) {
        case GEOMETRY:
          return 1000;
        case POINT:
          return 1001;
        case LINESTRING:
          return 1002;
        case POLYGON:
          return 1003;
        case MULTIPOINT:
          return 1004;
        case MULTILINESTRING:
          return 1005;
        case MULTIPOLYGON:
          return 1006;
        case GEOMETRYCOLLECTION:
          return 1007;
        default:
          return 1000;
      }
    } else if (!hasZ && hasM) {
      switch (this) {
        case GEOMETRY:
          return 2000;
        case POINT:
          return 2001;
        case LINESTRING:
          return 2002;
        case POLYGON:
          return 2003;
        case MULTIPOINT:
          return 2004;
        case MULTILINESTRING:
          return 2005;
        case MULTIPOLYGON:
          return 2006;
        case GEOMETRYCOLLECTION:
          return 2007;
        default:
          return 2000;
      }
    } else if (hasZ && hasM) {
      switch (this) {
        case GEOMETRY:
          return 3000;
        case POINT:
          return 3001;
        case LINESTRING:
          return 3002;
        case POLYGON:
          return 3003;
        case MULTIPOINT:
          return 3004;
        case MULTILINESTRING:
          return 3005;
        case MULTIPOLYGON:
          return 3006;
        case GEOMETRYCOLLECTION:
          return 3007;
        default:
          return 3000;
      }
    } else {
      switch (this) {
        case GEOMETRY:
          return 0;
        case POINT:
          return 1;
        case LINESTRING:
          return 2;
        case POLYGON:
          return 3;
        case MULTIPOINT:
          return 4;
        case MULTILINESTRING:
          return 5;
        case MULTIPOLYGON:
          return 6;
        case GEOMETRYCOLLECTION:
          return 7;
        default:
          return 0;
      }
    }
  }

  /// Returns the [EGeometryType] for a spatialite geometries types [value].
  static EGeometryType fromGeometryTypeCode(int value) {
    switch (value) {
      case 0:
        return GEOMETRY;
      case 1:
        return POINT;
      case 2:
        return LINESTRING;
      case 3:
        return POLYGON;
      case 4:
        return MULTIPOINT;
      case 5:
        return MULTILINESTRING;
      case 6:
        return MULTIPOLYGON;
      case 7:
        return GEOMETRYCOLLECTION;
      /*
         * XYZ
         */
      case 1000:
        return GEOMETRY;
      case 1001:
        return POINT;
      case 1002:
        return LINESTRING;
      case 1003:
        return POLYGON;
      case 1004:
        return MULTIPOINT;
      case 1005:
        return MULTILINESTRING;
      case 1006:
        return MULTIPOLYGON;
      case 1007:
        return GEOMETRYCOLLECTION;
      /*
         * XYM
         */
      case 2000:
        return GEOMETRY;
      case 2001:
        return POINT;
      case 2002:
        return LINESTRING;
      case 2003:
        return POLYGON;
      case 2004:
        return MULTIPOINT;
      case 2005:
        return MULTILINESTRING;
      case 2006:
        return MULTIPOLYGON;
      case 2007:
        return GEOMETRYCOLLECTION;
      /*
         * XYZM
         */
      case 3000:
        return GEOMETRY;
      case 3001:
        return POINT;
      case 3002:
        return LINESTRING;
      case 3003:
        return POLYGON;
      case 3004:
        return MULTIPOINT;
      case 3005:
        return MULTILINESTRING;
      case 3006:
        return MULTIPOLYGON;
      case 3007:
        return GEOMETRYCOLLECTION;
      default:
        break;
    }
    return UNKNOWN;
  }

  String getTypeName() {
    return typeName;
  }
}
