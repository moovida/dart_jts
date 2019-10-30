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

  List<Coordinate> poslist(l) => l.map(pos).toList();

  deserializePoint(Map gj) => pos(gj["coordinates"]);

  deserializeMultiPoint(Map gj) =>
      MultiPoint(poslist(gj["coordinates"]).map((dp) => Point(dp.x, dp.y)));

  deserializeLineString(Map gj) =>
      LineString(poslist(gj["coordinates"]).map((dp) => Point(dp.x, dp.y)));

  deserializeMultiLineString(Map gj) => MultiLineString(gj["coordinates"]
      .map((ls) => LineString(poslist(ls).map((dp) => Point(dp.x, dp.y))))
      .toList());

  polygonFromCoordinates(coords) {
    var rings = coords
        .map((l) => poslist(l))
        .map((poslist) => LineString(poslist))
        .toList();
    return new Polygon(rings.first, rings.skip(1));
  }

  deserializePolygon(Map gj) => polygonFromCoordinates(gj["coordinates"]);

  deserializeMultipolygon(Map gj) => MultiPolygon(gj["coordinates"]
      .map((coords) => polygonFromCoordinates(coords))
      .toList());

  var deserialize;

  deserializeGeometryCollection(Map gj) =>
      GeometryCollection(gj["geometries"].map((o) => deserialize(o)).toList());

  deserializeFeature(Map gj) {
    var geometry = deserialize(gj["geometry"]);
    var properties = gj["properties"];
    return Feature(geometry, properties);
  }

  deserializeFeatureCollection(Map gj) {
    var features = gj["features"].map((f) => deserializeFeature(f)).toList();
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
