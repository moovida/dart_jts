part of dart_jts;

/**
 * Builds an edge graph from geometries containing edges.
 *
 * @author mdavis
 *
 */
class EdgeGraphBuilder {
  static EdgeGraph build(List<Geometry> geoms) {
    EdgeGraphBuilder builder = new EdgeGraphBuilder();
    builder.addGeometries(geoms);
    return builder.getGraph();
  }

  EdgeGraph graph = new EdgeGraph();

  EdgeGraphBuilder() {}

  EdgeGraph getGraph() {
    return graph;
  }

  /**
   * Adds the edges of a Geometry to the graph.
   * May be called multiple times.
   * Any dimension of Geometry may be added; the constituent edges are
   * extracted.
   *
   * @param geometry geometry to be added
   */
  void add(Geometry geometry) {
    geometry.applyGCF(EdgeGraphBuilderFilter(this));
  }

  /**
   * Adds the edges in a collection of {@link Geometry}s to the graph.
   * May be called multiple times.
   * Any dimension of Geometry may be added.
   *
   * @param geometries the geometries to be added
   */
  void addGeometries(List<Geometry> geometries) {
    for (var geometry in geometries) {
      add(geometry);
    }
  }

  void addLineString(LineString lineString) {
    CoordinateSequence seq = lineString.getCoordinateSequence();
    for (int i = 1; i < seq.size(); i++) {
      graph.addEdge(seq.getCoordinate(i - 1), seq.getCoordinate(i));
    }
  }
}

class EdgeGraphBuilderFilter extends GeometryComponentFilter {
  EdgeGraphBuilder builder;
  EdgeGraphBuilderFilter(this.builder);
  void filter(Geometry component) {
    if (component is LineString) {
      builder.add(component);
    }
  }
}
