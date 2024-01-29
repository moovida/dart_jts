part of dart_jts;

/**
 * A simple elevation model used to populate missing Z values
 * in overlay results.
 * <p>
 * The model divides the extent of the input geometry(s)
 * into an NxM grid.
 * The default grid size is 3x3.
 * If the input has no extent in the X or Y dimension,
 * that dimension is given grid size 1.
 * The elevation of each grid cell is computed as the average of the Z values
 * of the input vertices in that cell (if any).
 * If a cell has no input vertices within it, it is assigned
 * the average elevation over all cells.
 * <p>
 * If no input vertices have Z values, the model does not assign a Z value.
 * <p>
 * The elevation of an arbitrary location is determined as the
 * Z value of the nearest grid cell.
 * <p>
 * An elevation model can be used to populate missing Z values
 * in an overlay result geometry.
 *
 * @author Martin Davis
 *
 */
class ElevationModel {
  static final int DEFAULT_CELL_NUM = 3;

  /**
   * Creates an elevation model from two geometries (which may be null).
   *
   * @param geom1 an input geometry
   * @param geom2 an input geometry, or null
   * @return the elevation model computed from the geometries
   */
  static ElevationModel create(Geometry geom1, Geometry? geom2) {
    Envelope extent = geom1.getEnvelopeInternal().copy();
    if (geom2 != null) {
      extent.expandToIncludeEnvelope(geom2.getEnvelopeInternal());
    }
    ElevationModel model =
        ElevationModel(extent, DEFAULT_CELL_NUM, DEFAULT_CELL_NUM);
    model.addGeometry(geom1);
    if (geom2 != null) model.addGeometry(geom2);
    return model;
  }

  late Envelope extent;
  late int numCellX;
  late int numCellY;
  late double cellSizeX;
  late double cellSizeY;
  late List<List<ElevationCell?>> cells;
  bool isInitialized = false;
  bool hasZValue = false;
  double averageZ = double.nan;

  /**
   * Creates a new elevation model covering an extent by a grid of given dimensions.
   *
   * @param extent the XY extent to cover
   * @param numCellX the number of grid cells in the X dimension
   * @param numCellY the number of grid cells in the Y dimension
   */
  ElevationModel(this.extent, this.numCellX, this.numCellY) {
    cellSizeX = extent.getWidth() / numCellX;
    cellSizeY = extent.getHeight() / numCellY;
    if (cellSizeX <= 0.0) {
      this.numCellX = 1;
    }
    if (cellSizeY <= 0.0) {
      this.numCellY = 1;
    }
    cells = List.empty(growable: true);
    for (var i = 0; i < numCellX; i++) {
      cells.add(List.filled(numCellY, null));
    }
  }

  /**
   * Updates the model using the Z values of a given geometry.
   *
   * @param geom the geometry to scan for Z values
   */
  void addGeometry(Geometry geom) {
    geom.applyCSF(AddGeometryFilter(this));
  }

  void add(double x, double y, double z) {
    if (z.isNaN) return;
    hasZValue = true;
    ElevationCell? cell = getCell(x, y, true);
    cell!.add(z);
  }

  void init() {
    isInitialized = true;
    int numCells = 0;
    double sumZ = 0.0;

    for (int i = 0; i < cells.length; i++) {
      for (int j = 0; j < cells[0].length; j++) {
        ElevationCell? cell = cells[i][j];
        if (cell != null) {
          cell.compute();
          numCells++;
          sumZ += cell.getZ();
        }
      }
    }
    averageZ = double.nan;
    if (numCells > 0) {
      averageZ = sumZ / numCells;
    }
  }

  /**
   * Gets the model Z value at a given location.
   * If the location lies outside the model grid extent,
   * this returns the Z value of the nearest grid cell.
   * If the model has no elevation computed (i.e. due
   * to empty input), the value is returned as {@link Double#NaN}.
   *
   * @param x the x ordinate of the location
   * @param y the y ordinate of the location
   * @return the computed model Z value
   */
  double getZ(double x, double y) {
    if (!isInitialized) init();
    ElevationCell? cell = getCell(x, y, false);
    if (cell == null) return averageZ;
    return cell.getZ();
  }

  /**
   * Computes Z values for any missing Z values in a geometry,
   * using the computed model.
   * If the model has no Z value, or the geometry coordinate dimension
   * does not include Z, the geometry is not updated.
   *
   * @param geom the geometry to populate Z values for
   */
  void populateZ(Geometry geom) {
    // short-circuit if no Zs are present in model
    if (!hasZValue) return;

    if (!isInitialized) init();

    geom.applyCSF(PopulateZFilter(this));
  }

  ElevationCell? getCell(double x, double y, bool isCreateIfMissing) {
    int ix = 0;
    if (numCellX > 1) {
      ix = (x - extent.getMinX()) ~/ cellSizeX;
      ix = ix.clamp(0, numCellX - 1);
    }
    int iy = 0;
    if (numCellY > 1) {
      iy = (y - extent.getMinY()) ~/ cellSizeY;
      iy = iy.clamp(0, numCellY - 1);
    }
    ElevationCell? cell = cells[ix][iy];
    if (isCreateIfMissing && cell == null) {
      cell = ElevationCell();
      cells[ix][iy] = cell;
    }
    return cell;
  }
}

class AddGeometryFilter extends CoordinateSequenceFilter {
  ElevationModel model;
  AddGeometryFilter(this.model);
  bool _hasZ = true;


  @override
  void filter(CoordinateSequence seq, int i) {
    if (!seq.hasZ()) {
      // if no Z then short-circuit evaluation
      _hasZ = false;
      return;
    }
    // if Z not populated then assign using model
    if ((seq.getZ(i)).isNaN) {
      double z = model.getZ(
          seq.getOrdinate(i, Coordinate.X), seq.getOrdinate(i, Coordinate.Y));
      seq.setOrdinate(i, Coordinate.Z, z);
    }
  }

  @override
  bool isDone() {
    return !_hasZ;
  }

  @override
  bool isGeometryChanged() {
    // geometry extent is not changed
    return false;
  }
}

class PopulateZFilter extends CoordinateSequenceFilter {
  bool _isDone = false;
  ElevationModel model;
  PopulateZFilter(this.model);
  @override
  void filter(CoordinateSequence seq, int i) {
    if (!seq.hasZ()) {
      // if no Z then short-circuit evaluation
      _isDone = true;
      return;
    }
    // if Z not populated then assign using model
    if ((seq.getZ(i)).isNaN) {
      double z = model.getZ(seq.getOrdinate(i, Coordinate.X),
          seq.getOrdinate(i, Coordinate.Y));
      seq.setOrdinate(i, Coordinate.Z, z);
    }
  }

  @override
  bool isDone() {
    return _isDone;
  }

  @override
  bool isGeometryChanged() {
    // geometry extent is not changed
    return false;
  }

}
