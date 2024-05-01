part of dart_jts;

class GeometryPrecisionReducer {
  ///
  /// Reduces precision of a geometry, ensuring output geometry is valid.
  /// Collapsed linear and polygonal components are removed.
  /// Duplicate vertices are removed.
  /// The geometry precision model is not changed.
  /// <p>
  /// Invalid input geometry may cause an error,
  /// unless the invalidity is below the scale of the precision reduction.
  ///
  ///* @param g the geometry to reduce
  /// @param precModel the precision model to use
  /// @return the reduced geometry
  /// @throws IllegalArgumentException if the reduction fails due to invalid input geometry
  ///
  Geometry? reduce(Geometry g, PrecisionModel precModel) {
    GeometryPrecisionReducer reducer = GeometryPrecisionReducer(precModel);
    return reducer.reduceGeometry(g);
  }

  /**
   * Reduces precision of a geometry, ensuring output polygonal geometry is valid,
   * and preserving collapsed linear elements.
   * Duplicate vertices are removed.
   * The geometry precision model is not changed.
   * <p>
   * Invalid input geometry may cause an error,
   * unless the invalidity is below the scale of the precision reduction.
   *
   * @param g the geometry to reduce
   * @param precModel the precision model to use
   * @return the reduced geometry
   * @throws IllegalArgumentException if the reduction fails due to invalid input geometry
   */
  static Geometry? reduceKeepCollapsed(Geometry geom, PrecisionModel pm) {
    GeometryPrecisionReducer reducer = GeometryPrecisionReducer(pm);
    reducer.setRemoveCollapsedComponents(false);
    return reducer.reduceGeometry(geom);
  }

  /**
   * Reduce precision of a geometry in a pointwise way.
   * All input geometry elements are preserved in the output,
   * including invalid polygons and collapsed polygons and linestrings.
   * The output may not be valid, due to collapse or self-intersection.
   * Duplicate vertices are not removed.
   * The geometry precision model is not changed.
   * <p>
   * Invalid input geometry is allowed.
   *
   * @param g the geometry to reduce
   * @param precModel the precision model to use
   * @return the reduced geometry
   */
  static Geometry? reducePointwise(Geometry g, PrecisionModel precModel) {
    GeometryPrecisionReducer reducer = GeometryPrecisionReducer(precModel);
    reducer.setPointwise(true);
    return reducer.reduceGeometry(g);
  }

  PrecisionModel? targetPM;
  bool removeCollapsed = true;
  bool changePrecisionModel = false;
  bool isPointwise = false;

  GeometryPrecisionReducer(PrecisionModel pm) {
    targetPM = pm;
  }

  /**
   * Sets whether the reduction will result in collapsed components
   * being removed completely, or simply being collapsed to an (invalid)
   * Geometry of the same type.
   * The default is to remove collapsed components.
   *
   * @param removeCollapsed if <code>true</code> collapsed components will be removed
   */
  void setRemoveCollapsedComponents(bool removeCollapsed) {
    this.removeCollapsed = removeCollapsed;
  }

  /**
   * Sets whether the {@link PrecisionModel} of the new reduced Geometry
   * will be changed to be the {@link PrecisionModel} supplied to
   * specify the precision reduction.
   * <p>
   * The default is to <b>not</b> change the precision model
   *
   * @param changePrecisionModel if <code>true</code> the precision model of the created Geometry will be the
   * the precisionModel supplied in the constructor.
   */
  void setChangePrecisionModel(bool changePrecisionModel) {
    this.changePrecisionModel = changePrecisionModel;
  }

  /**
   * Sets whether the precision reduction will be done
   * in pointwise fashion only.
   * Pointwise precision reduction reduces the precision
   * of the individual coordinates only, but does
   * not attempt to recreate valid topology.
   * This is only relevant for geometries containing polygonal components.
   *
   * @param isPointwise if reduction should be done pointwise only
   */
  void setPointwise(bool isPointwise) {
    this.isPointwise = isPointwise;
  }

  /**
   * Reduces the precision of a geometry,
   * according to the specified strategy of this reducer.
   *
   * @param geom the geometry to reduce
   * @return the precision-reduced geometry
   * @throws IllegalArgumentException if the reduction fails due to invalid input geometry is invalid
   */
  Geometry? reduceGeometry(Geometry geom) {
    Geometry reduced;
    if (isPointwise) {
      reduced = PointwisePrecisionReducerTransformer.reduce(geom, targetPM!);
    } else {
      reduced =
          PrecisionReducerTransformer.reduce(geom, targetPM!, removeCollapsed);
    }

    if (changePrecisionModel) {
      return _changePM(reduced, targetPM!);
    }
    return reduced;
  }

  /**
   * Duplicates a geometry to one that uses a different PrecisionModel,
   * without changing any coordinate values.
   *
   * @param geom the geometry to duplicate
   * @param newPM the precision model to use
   * @return the geometry value with a new precision model
   */
  Geometry? _changePM(Geometry geom, PrecisionModel newPM) {
    GeometryEditor geomEditor = _createEditor(geom.getFactory(), newPM);
    // this operation changes the PM for the entire geometry tree
    return geomEditor.edit(geom, NoOpGeometryOperation());
  }

  GeometryEditor _createEditor(
      GeometryFactory geomFactory, PrecisionModel newPM) {
    // no need to change if precision model is the same
    if (geomFactory.getPrecisionModel() == newPM) return GeometryEditor.empty();
    // otherwise create a geometry editor which changes PrecisionModel
    GeometryFactory newFactory = _createFactory(geomFactory, newPM);
    GeometryEditor geomEdit = new GeometryEditor(newFactory);
    return geomEdit;
  }

  GeometryFactory _createFactory(
      GeometryFactory inputFactory, PrecisionModel pm) {
    GeometryFactory newFactory = new GeometryFactory(pm, inputFactory.getSRID(),
        inputFactory.getCoordinateSequenceFactory());
    return newFactory;
  }
}
