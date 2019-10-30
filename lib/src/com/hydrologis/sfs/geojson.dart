part of dart_sfs;

/// Parses GeoJSON.
///
/// Returns a [Geometry], a [Feature], or a [FeatureCollection].
///
//TODO: more checks, throw FormatException on error
//TODO: provide an optional factory object, in particular for Features
//  and FeatureCollection
parseGeoJson(String geoJson) {
  var value = json.jsonDecode(geoJson);
  assert(value is Map);

  Point pos(coord) => Point(coord[0], coord[1]);

  List<Point> poslist(l) {
    var map = l.map(pos);
    List<dynamic> list2 = map.toList();
    List<Point> list = [];
    list2.forEach((dyn) {
      list.add(dyn as Point);
    });
    return list;
  }

  Point deserializePoint(Map gj) => pos(gj["coordinates"]);

  MultiPoint deserializeMultiPoint(Map gj) =>
      MultiPoint(poslist(gj["coordinates"]));

  LineString deserializeLineString(Map gj) =>
      LineString(poslist(gj["coordinates"]));

  MultiLineString deserializeMultiLineString(Map gj) => MultiLineString(
      gj["coordinates"].map((ls) => LineString(poslist(ls))).toList());

  Polygon polygonFromCoordinates(coords) {
    var rings = coords
        .map((l) => poslist(l))
        .map((poslist) => LineString(poslist))
        .toList();
    return Polygon(rings.first, rings.skip(1));
  }

  Polygon deserializePolygon(Map gj) =>
      polygonFromCoordinates(gj["coordinates"]);

  MultiPolygon deserializeMultipolygon(Map gj) => MultiPolygon(gj["coordinates"]
      .map((coords) => polygonFromCoordinates(coords))
      .toList());

  var deserialize;

  GeometryCollection deserializeGeometryCollection(Map gj) {
   var list = gj["geometries"].map((o) => deserialize(o)).toList();
   List<Geometry> geoms = [];
   list.forEach((dyn) {
     geoms.add(dyn as Geometry);
   });
   return GeometryCollection(geoms);
  }

  Feature deserializeFeature(Map gj) {
    var geometry = deserialize(gj["geometry"]);
    var properties = gj["properties"];
    return Feature(geometry, properties);
  }

  FeatureCollection deserializeFeatureCollection(Map gj) {
    List<dynamic> featuresDyn =
        gj["features"].map((f) => deserializeFeature(f)).toList();
    List<Feature> features = [];
    featuresDyn.forEach((dyn) {
      features.add(dyn as Feature);
    });
    return FeatureCollection(features);
  }

  deserialize = (Map gj) {
    switch (gj["type"]) {
      case "Point":
        return deserializePoint(gj);
      case "MultiPoint":
        return deserializeMultiPoint(gj);
      case "LineString":
        return deserializeLineString(gj);
      case "MultiLineString":
        return deserializeMultiLineString(gj);
      case "Polygon":
        return deserializePolygon(gj);
      case "MultiPolygon":
        return deserializeMultipolygon(gj);
      case "GeometryCollection":
        return deserializeGeometryCollection(gj);
      case "Feature":
        return deserializeFeature(gj);
      case "FeatureCollection":
        return deserializeFeatureCollection(gj);
      default:
        throw FormatException("unknown GeoJson object type '${gj['type']}");
    }
  };

  return deserialize(value);
}
