part of dart_sfs;

/// A [Feature] isn't a standard class in the SFS.
///
/// It represents an arbirtray object with a set of non-spatial
/// properties and one dedicated spatial property called [geometry].
///
/// It is is mainly used to
/// deserialize GeoJSON objects with `type: "Feature"`.
abstract class Feature {
  /// Creates a feature with a specific [geometry] and an (optional) set
  /// of properties.
  ///
  /// [geometry] must not be null, otherwise throws an [ArgumentError].
  ///
  /// [properties] are optional, and if present, they can be null or empty.
  factory Feature(Geometry geometry, [Map<String, dynamic> properties = null]) {
    return _FeatureImpl(geometry, properties);
  }

  Geometry get geometry;

  Map<String, dynamic> get properties;
}

class _FeatureImpl implements Feature {
  Map<String, dynamic> _properties;
  Geometry _geometry;

  _FeatureImpl(this._geometry, [Map<String, dynamic> properties]) {
    _require(_geometry != null);
    if (properties == null || properties.isEmpty) return;
    _properties = Map<String, dynamic>.from(properties);
  }

  Map<String, dynamic> get properties => _properties == null ? {} : _properties;

  Geometry get geometry => _geometry;
}

/// A [FeatureCollection] isn't a standard class in the SFS.
///
/// It represents a collection of [Feature] and is mainly used to
/// deserialized GeoJSON objects with `type: "FeatureCollection"`.
class FeatureCollection extends Object with IterableMixin<Feature> {
  List<Feature> _features;

  /// Creates a new feature collection with [features].
  ///
  /// [features] is option. If present, if may be null or empty, in which
  /// cases an empty feature collection is created.
  FeatureCollection([Iterable<Feature> features]) {
    if (features == null || features.isEmpty) return;
    _features = List.from(features, growable: false);
  }

  Iterator<Feature> get iterator =>
      _features == null ? [].iterator : _features.iterator;
}
