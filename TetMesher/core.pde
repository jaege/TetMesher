//import java.util.Map; //<>// //<>// //<>// //<>// //<>// //<>// //<>// //<>// //<>// //<>// //<>//

// A simplified sphere representation in 3D
class Sphere {
  pt O;
  float r2;

  Sphere() {
    this.O = P();
    this.r2 = 0;
  }

  Sphere(pt O, float r2) {
    this.O = O;
    this.r2 = r2;
  }

  float r() {
    return sqrt(r2);
  }
}

/**
 * use barycentric coordinates to calculate circumcenter
 * assume A, B and C are not on the same line
 */
Sphere circumcircle(pt A, pt B, pt C) {
  Sphere sph = new Sphere();
  // unnormalized version
  //float a2 = d2(B, C), b2 = d2(A, C), c2= d2(A, B);
  //float r = a2 * (b2 + c2 - a2), s = b2 * (c2 + a2 - b2), t = c2 * (a2 + b2 - c2);
  //float n = 1.0 / (r + s + t);
  //sph.O = P(r*n, A, s*n, B, t*n, C);

  // normalized version
  float p = d2(A, B), q = d2(A, C), r = dot(V(A, B), V(A, C));
  float pq = p * q, n = 0.5/(pq - r * r);
  float s = (pq - q * r) * n, t = (pq-p * r) * n;
  sph.O = P(1 - s - t, A, s, B, t, C);

  sph.r2 = d2(A, sph.O);
  return sph;
}

/**
 * use barycentric coordinates to calculate circumcenter
 * assume A, B, C and D are not on the same ???
 */
Sphere circumsphere(pt A, pt B, pt C, pt D) {
  Sphere sph = new Sphere();
  vec AB = V(A, B), AC = V(A, C), AD = V(A, D);
  float a = dot(AB, AB), b = dot(AC, AC), c = dot(AD, AD), 
    d = dot(AB, AC), e = dot(AB, AD), f = dot(AC, AD), 
    n = 0.5/(c*d*d - 2*d*e*f + b*e*e + a*f*f - a*b*c), 
    s = (a*f*f - a*b*c + b*c*d + b*c*e - b*e*f - c*d*f) * n, 
    t = (b*e*e - a*b*c + a*c*d + a*c*f - a*e*f - c*d*e) * n, 
    u = (c*d*d - a*b*c + a*b*e + a*b*f - a*d*f - b*d*e) * n;
  sph.O = P(1 - s - t - u, A, s, B, t, C, u, D);
  sph.r2 = d2(A, sph.O);
  return sph;
}

Sphere circumsphere(Vertex A, Vertex B, Vertex C, float r) {
  Sphere sph = new Sphere();
  Sphere c = circumcircle(A.p, B.p, C.p);
  vec n0 = V(A.n, B.n, C.n);
  vec n = U(N(A.p, B.p, C.p));
  if (dot(n, n0) < 0) n.rev();
  sph.r2 = r * r;
  if (sph.r2 >= c.r2) {
    sph.O = P(c.O, sqrt(sph.r2 - c.r2), n);
  } else {
    sph.r2 = -1;
  }
  return sph;
}

int sign(float f) {
  return f >= 0 ? 1 : -1;
}

float bulge(pt A, pt B, pt C, pt D) {
  Sphere s1 = circumcircle(A, B, C);
  Sphere s = circumsphere(A, B, C, D);
  return s.r() + sign(dot(s1.O, D, s.O)) * d(s.O, s1.O);
}

class Edge {
  int a, b;  // id of point: a to b
  boolean boundary, inner;
  int o;  // opposite
  Edge(int a, int b) {
    this.a = a;
    this.b = b;
    this.boundary = this.inner = false;
  }
  Edge(int a, int b, int o) {
    this.a = a;
    this.b = b;
    this.o = o;
    this.boundary = this.inner = false;
  }
  Boolean isSame(Edge e) {
    return this.a == e.a && this.b == e.b; // || this.a == e.b && this.b == e.a;
  }
}

class Vertex {
  pt p;
  vec n;
  Vertex(pt p, vec n) {
    this.p = p;
    this.n = n;
  }
}

class TriangleMesh2 {
  ArrayList<Vertex> v = new ArrayList<Vertex>();
  ArrayList<Edge> edges = new ArrayList<Edge>();
  ArrayList<int []> triangles = new ArrayList<int []>();  // loop: id1 -> id2 -> id3 -> id1

  TriangleMesh2() {
  }

  /**
   * generate triangle mesh using ball pivoting algorithm
   */
  TriangleMesh2(ArrayList<Vertex> pointCloud, float r) {
    this.v = pointCloud;
    this.edges = new ArrayList<Edge>();
    this.triangles = new ArrayList<int []>();

    ArrayList<Edge> front = new ArrayList<Edge>();
    HashMap<Integer, HashMap<Integer, Edge>> edgesMap = new HashMap<Integer, HashMap<Integer, Edge>>();

    println("TriangleMesh2: Sampled", this.v.size(), "vertices.");
    //for (;; ) {
    int [] seed = findSeed(r);
    println(seed.length);
    if (seed.length != 3) {
      println("TriangleMesh2: Cannot find seed triangle.");
      //break;
    } else {
      println("TriangleMesh2: Found seed triangle", seed[0], seed[1], seed[2]);
      this.triangles.add(seed);
      Edge e1 = addEdge(edgesMap, seed[0], seed[1]);
      e1.o = seed[2];
      front.add(e1);
      Edge e2 = addEdge(edgesMap, seed[1], seed[2]);
      e2.o = seed[0];
      front.add(e2);
      Edge e3 = addEdge(edgesMap, seed[2], seed[0]);
      e3.o = seed[1];
      front.add(e3);   

      while (!front.isEmpty()) {
        Edge e = front.get(front.size()-1);
        front.remove(front.size()-1);
        this.edges.add(e);
        if (e.inner || e.boundary) {
          continue;
        }
        println("TriangleMesh2: Try to find pivot from edge", e.a, e.b);
        int vid = getPivot(e, r);
        if (vid == -1) {
          println("TriangleMesh2: Found boundary edge", e.a, e.b);
          e.boundary = true;
          continue;
        }
        println("TriangleMesh2: Found pivot", vid, "from edge", e.a, e.b);
        this.triangles.add(new int [] {e.a, vid, e.b});

        if (hasEdge(edgesMap, vid, e.a)) {
          println("TriangleMesh2: Found inner edge", vid, e.a);
          Edge es = getEdge(edgesMap, vid, e.a);
          es.inner = true;
        } else if (hasEdge(edgesMap, e.a, vid)) {
          println("TriangleMesh2:Warn: Found duplicated edge", vid, e.a);
          Edge es = getEdge(edgesMap, e.a, vid);
          es.inner = true;
        } else {
          println("TriangleMesh2: Add new edge", e.a, vid);
          Edge es = addEdge(edgesMap, e.a, vid);
          es.o = e.b;
          front.add(es);
        }

        if (hasEdge(edgesMap, e.b, vid)) {
          println("TriangleMesh2: Found inner edge", e.b, vid);
          Edge et = getEdge(edgesMap, e.b, vid);
          et.inner = true;
        } else if (hasEdge(edgesMap, vid, e.b)) {
          println("TriangleMesh2:Warn: Found duplicated edge", vid, e.b);
          Edge es = getEdge(edgesMap, vid, e.b);
          es.inner = true;
        } else {
          println("TriangleMesh2: Add new edge", vid, e.b);
          Edge et = addEdge(edgesMap, vid, e.b);
          et.o = e.a;
          front.add(et);
        }
      }
    }
    //}
    println("TriangleMesh2: Created mesh with ", this.v.size(), "vertices,", 
      this.edges.size(), "edges and", this.triangles.size(), "triangles.");
  }

  int getPivot(Edge e, float r /*, boolean [] visited*/) {
    //println("getPivot: Try to find pivot from edge", e.a, e.b);
    Vertex a = this.v.get(e.a), b = this.v.get(e.b);
    Sphere s = circumsphere(a, b, this.v.get(e.o), r);
    if (s.r2 < 0) {
      println("getPivot:Error: Cannot find circumsphere.");
      return -1;
    }
    pt M = P(a.p, b.p);
    float dd = sq(d(M, s.O) + s.r());
    float thetaMin = TWO_PI;
    int candidate = -1;
    for (int i = 0; i < this.v.size(); ++i) {
      if (i == e.a || i == e.b || i == e.o) continue;
      Vertex vi = this.v.get(i);
      if (d2(vi.p, s.O) > dd) continue;
      //if (v.inner) continue;
      if (!compatible(a, vi, b)) continue;
      Sphere s2 = circumsphere(a, vi, b, r);
      if (s2.r2 < 0) continue;
      if (anyInSphere(s2, this.v, e.a, e.b, i)) continue;
      println("getPivot: Found possible vertex", i);            
      float theta = acos(dot(U(M, s2.O), U(M, s.O)));
      //println(theta);
      if (cw(V(a.p, b.p), V(M, s2.O), V(M, s.O))) theta = TWO_PI - theta;
      //println(theta);
      if (thetaMin > theta) {
        thetaMin = theta;
        candidate = i;
      }
    }
    return candidate;
  }

  int [] findSeed(float r) {
    for (int i = 0; i < this.v.size() - 2; ++i) {
      println("findSeed: Try", i);
      for (int j = i + 1; j < this.v.size() - 1; ++j) {
        //println("findSeed: Try", i, j);
        for (int k = j + 1; k < this.v.size(); ++k) {
          if (!compatible(this.v.get(i), this.v.get(j), this.v.get(k)) &&
            !compatible(this.v.get(i), this.v.get(k), this.v.get(j))) continue;
          //println("findSeed: Try", i, j, k);
          Sphere s = circumsphere(this.v.get(i), this.v.get(j), this.v.get(k), r);
          if (s.r2 < 0) continue;
          if (anyInSphere(s, this.v, i, j, k)) continue;
          if (cw(s.O, this.v.get(i).p, this.v.get(j).p, this.v.get(k).p))
            return new int [] {i, k, j};
          else
            return new int [] {i, j, k};
        }
      }
    }
    return new int [0];
  }
}

boolean hasEdge(HashMap<Integer, HashMap<Integer, Edge>> edges, int a, int b) {
  return edges.containsKey(a) && edges.get(a).containsKey(b);
}

Edge addEdge(HashMap<Integer, HashMap<Integer, Edge>> edges, int a, int b) {
  if (!edges.containsKey(a)) edges.put(a, new HashMap<Integer, Edge>());
  if (edges.get(a).containsKey(b))
    println("addEdge:Warn: Rewite edge", a, b);
  edges.get(a).put(b, new Edge(a, b));
  return edges.get(a).get(b);
}

Edge getEdge(HashMap<Integer, HashMap<Integer, Edge>> edges, int a, int b) {
  //if (!edges.containsKey(a)) edges.put(a, new HashMap<Integer, Edge>());
  //if (!edges.get(a).containsKey(b)) edges.get(a).put(b, new Edge(a, b));
  return edges.get(a).get(b);
}

boolean compatible(Vertex A, Vertex B, Vertex C) {
  vec n = N(A.p, B.p, C.p);
  return dot(n, A.n) > 0 && dot(n, B.n) > 0 && dot(n, C.n) > 0;
  //|| dot(n, A.n) < 0 && dot(n, B.n) < 0 && dot(n, C.n) < 0;
}

class TetMesh {
  ArrayList<Edge> floorEdges = new ArrayList<Edge>();  // floor pt id1 < floor pt id2
  ArrayList<Edge> ceilingEdges = new ArrayList<Edge>();  // ceiling pt id1 < ceiling pt id2
  ArrayList<Edge> mixedEdges = new ArrayList<Edge>();  // floor pt id1, ceiling pt id2
  ArrayList<int []> floorTriangles = new ArrayList<int []>();  // id1-3: floor
  ArrayList<int []> ceilingTriangles = new ArrayList<int []>();  // id1-3: ceiling
  ArrayList<int []> floorTets = new ArrayList<int []>();  // id1-3: floor, id4: ceiling
  ArrayList<int []> ceilingTets = new ArrayList<int []>();  // id1: floor, id2-4: ceiling
  ArrayList<int []> mixedTets = new ArrayList<int []>();  // id1-2: floor, id3-4: ceiling
  pts floor;
  pts ceiling;

  void removeDup() {
    removeDuplicatedEdges(floorEdges);
    removeDuplicatedEdges(ceilingEdges);
    removeDuplicatedEdges(mixedEdges);
  }

  TriangleMesh2 triangulate() {
    int sphereSubdivision = 1;
    int tubeSubdivision = 10;
    float r = rb+5;

    ArrayList<Vertex> v = new ArrayList<Vertex>();
    for (Edge e : this.mixedEdges) {
      v.addAll(v.size(), sampleTube(this.floor.G[e.a], this.ceiling.G[e.b], rt, tubeSubdivision));
    }
    for (Edge e : this.floorEdges) {
      v.addAll(v.size(), sampleTube(this.floor.G[e.a], this.floor.G[e.b], rt, tubeSubdivision));
    }
    for (Edge e : this.ceilingEdges) {
      v.addAll(v.size(), sampleTube(this.ceiling.G[e.a], this.ceiling.G[e.b], rt, tubeSubdivision));
    }
    for (int i = 0; i < this.floor.nv; ++i) {
      v.addAll(v.size(), sampleSphere(this.floor.G[i], rb, sphereSubdivision));
    }
    for (int i = 0; i < this.ceiling.nv; ++i) {
      v.addAll(v.size(), sampleSphere(this.ceiling.G[i], rb, sphereSubdivision));
    }

    return new TriangleMesh2(v, r);
  }
}

ArrayList<Vertex> sampleTube(pt P, pt Q, float r, int s) {
  ArrayList<Vertex> v = new ArrayList<Vertex>();
  vec V = V(P, Q);
  vec I = U(Normal(V));
  vec J = U(N(I, V));
  float n = V.norm();
  float da = TWO_PI/s;
  int kk = int(n / (0.8660254 * da * r));
  float dv = 1.0 / kk;
  for (int k = 0; k <= kk; ++k) {
    float a0 = k % 2 == 0 ? 0 : da / 2;
    pt p0 = P(P, k * dv, V);
    for (float a = a0; a <= TWO_PI+a0; a += da) {
      pt p = P(p0, r*cos(a), I, r*sin(a), J);
      v.add(new Vertex(p, U(p0, p)));
    }
  }
  return v;
}

ArrayList<Vertex> sampleSphere(pt O, float r, int s) {
  if (s < 0) s = 0;
  ArrayList<pt> v = new ArrayList<pt>();
  vec I = V(1, 0, 0), J = V(0, 1, 0), K = V(0, 0, 1);
  pt P1 = P(O, K), P2 = P(O, I), P3 = P(O, J), P4 = P(O, -1, I), P5 = P(O, -1, J), P6 = P(O, -1, K);
  v.add(P1);
  v.add(P2);
  v.add(P3);
  v.add(P4);
  v.add(P5);
  v.add(P6);
  float k = 2 << s;
  for (int i = 1; i < k; ++i) {
    float t = i / k;
    v.add(P(P1, t, V(P1, P2)));
    v.add(P(P1, t, V(P1, P3)));
    v.add(P(P1, t, V(P1, P4)));
    v.add(P(P1, t, V(P1, P5)));
    v.add(P(P6, t, V(P6, P2)));
    v.add(P(P6, t, V(P6, P3)));
    v.add(P(P6, t, V(P6, P4)));
    v.add(P(P6, t, V(P6, P5)));
    v.add(P(P2, t, V(P2, P3)));
    v.add(P(P3, t, V(P3, P4)));
    v.add(P(P4, t, V(P4, P5)));
    v.add(P(P5, t, V(P5, P2)));
  }
  v.addAll(v.size(), sampleTriangle(P1, P2, P3, s));
  v.addAll(v.size(), sampleTriangle(P1, P3, P4, s));
  v.addAll(v.size(), sampleTriangle(P1, P4, P5, s));
  v.addAll(v.size(), sampleTriangle(P1, P5, P2, s));
  v.addAll(v.size(), sampleTriangle(P6, P3, P2, s));
  v.addAll(v.size(), sampleTriangle(P6, P4, P3, s));
  v.addAll(v.size(), sampleTriangle(P6, P5, P4, s));
  v.addAll(v.size(), sampleTriangle(P6, P2, P5, s));

  ArrayList<Vertex> ve = new ArrayList<Vertex>();
  for (int i = 0; i < v.size(); ++i) {
    vec vn = U(O, v.get(i));
    pt p = P(O, r, vn);
    ve.add(new Vertex(p, vn));
  }

  return ve;
}

/**
 * return inner vertices, i.e. not include vertices on the edgek
 */
ArrayList<pt> sampleTriangle(pt A, pt B, pt C, int s) {
  ArrayList<pt> vertices = new ArrayList<pt>();
  s = 2 << s;
  vec di = V(A, B).div(s), dj = V(B, C).div(s);
  for (int i = 1; i < s; ++i) {
    pt P0 = P(A, i, di);
    for (int j = 1; j < i; ++j) {
      vertices.add(P(P0, j, dj));
    }
  }
  return vertices;
}

void removeDuplicatedEdges(ArrayList<Edge> edges) {
  for (int i = edges.size() - 1; i > 0; --i) {
    for (int j = i - 1; j >= 0; --j) {
      if (edges.get(i).isSame(edges.get(j))) {
        edges.remove(i);
        break;
      }
    }
  }
}

class TriangleMesh {
  ArrayList<Edge> edges = new ArrayList<Edge>();
  ArrayList<int []> triangles = new ArrayList<int []>();  // id1 < id2 < id3
  pts vertices;
}

boolean contains(int[] array, int v) {
  for (int i : array) {
    if (i == v) return true;
  }
  return false;
}

/**
 * O(n)
 */
Boolean anyInSphere(Sphere s, pts v, int ... excepts) {
  for (int m = 0; m < v.nv; ++m) {
    if (!contains(excepts, m) && d2(v.G[m], s.O) < s.r2) return true;
  }
  return false;
}

Boolean anyInSphere(Sphere s, ArrayList<Vertex> v, int ... excepts) {
  for (int m = 0; m < v.size(); ++m) {
    if (!contains(excepts, m) && d2(v.get(m).p, s.O) <= s.r2) return true;
  }
  return false;
}

/** O(n^4)
 * assume all vertices are in the same plane
 */
TriangleMesh triangulate(pts vertices) {
  TriangleMesh mesh = new TriangleMesh();
  mesh.vertices = vertices;
  for (int i = 0; i < vertices.nv - 2; ++i) {
    for (int j = i + 1; j < vertices.nv - 1; ++j) {
      for (int k = j + 1; k < vertices.nv; ++k) {
        Sphere s = circumcircle(vertices.G[i], vertices.G[j], vertices.G[k]);
        if (anyInSphere(s, vertices, i, j, k)) continue;
        mesh.edges.add(new Edge(i, j));
        mesh.edges.add(new Edge(i, k));
        mesh.edges.add(new Edge(j, k));
        mesh.triangles.add(new int [] {i, j, k});
      }
    }
  }
  removeDuplicatedEdges(mesh.edges);
  return mesh;
}

/**
 * O(n^4 + nm + nm(n+m))
 */
TetMesh tetrahedralize_v2(pts floor, pts ceiling) {
  TetMesh tetMesh = new TetMesh();
  tetMesh.floor = floor;
  tetMesh.ceiling = ceiling;
  TriangleMesh floorMesh = triangulate(floor);
  TriangleMesh ceilingMesh = triangulate(ceiling);

  if (showFloorTet) {  // case 1: 3 down + 1 up
    tetMesh.floorEdges.addAll(tetMesh.floorEdges.size(), floorMesh.edges);
    tetMesh.floorTriangles.addAll(tetMesh.floorTriangles.size(), floorMesh.triangles);
    for (int [] t : floorMesh.triangles) {
      int i = t[0], j = t[1], k = t[2];
      int m_min = 0;
      float b_min = bulge(floor.G[i], floor.G[j], floor.G[k], ceiling.G[m_min]);
      for (int m = 1; m < ceiling.nv; ++m) {
        float b = bulge(floor.G[i], floor.G[j], floor.G[k], ceiling.G[m]);
        if (b_min > b) {
          b_min = b;
          m_min = m;
        }
      }
      tetMesh.mixedEdges.add(new Edge(i, m_min));
      tetMesh.mixedEdges.add(new Edge(j, m_min));
      tetMesh.mixedEdges.add(new Edge(k, m_min));
      tetMesh.floorTets.add(new int []{i, j, k, m_min});
    }
  }

  if (showCeilingTet) {  // case 2: 1 down + 3 up
    tetMesh.ceilingEdges.addAll(tetMesh.ceilingEdges.size(), ceilingMesh.edges);
    tetMesh.ceilingTriangles.addAll(tetMesh.ceilingTriangles.size(), ceilingMesh.triangles);
    for (int [] t : ceilingMesh.triangles) {
      int i = t[0], j = t[1], k = t[2];
      int m_min = 0;
      float b_min = bulge(ceiling.G[i], ceiling.G[j], ceiling.G[k], floor.G[m_min]);
      for (int m = 1; m < floor.nv; ++m) {
        float b = bulge(ceiling.G[i], ceiling.G[j], ceiling.G[k], floor.G[m]);
        if (b_min > b) {
          b_min = b;
          m_min = m;
        }
      }
      tetMesh.mixedEdges.add(new Edge(m_min, i));
      tetMesh.mixedEdges.add(new Edge(m_min, j));
      tetMesh.mixedEdges.add(new Edge(m_min, k));
      tetMesh.ceilingTets.add(new int []{m_min, i, j, k});
    }
  }

  if (showMixedTet) {  // case 3: 2 down + 2 up
    for (Edge ef : floorMesh.edges) {
      for (Edge ec : ceilingMesh.edges) {
        int i = ef.a, j = ef.b, m = ec.a, n = ec.b;
        Sphere s = circumsphere(floor.G[i], floor.G[j], ceiling.G[m], ceiling.G[n]);
        if (anyInSphere(s, floor, i, j) || anyInSphere(s, ceiling, m, n)) continue;
        tetMesh.floorEdges.add(new Edge(i, j));
        tetMesh.mixedEdges.add(new Edge(i, m));
        tetMesh.mixedEdges.add(new Edge(i, n));
        tetMesh.mixedEdges.add(new Edge(j, m));
        tetMesh.mixedEdges.add(new Edge(j, n));
        tetMesh.ceilingEdges.add(new Edge(m, n));
        tetMesh.mixedTets.add(new int []{i, j, m, n});
        break;
      }
    }
  }

  tetMesh.removeDup();
  return tetMesh;
}


TetMesh tetrahedralize_v1(pts floor, pts ceiling) {
  TetMesh mesh = new TetMesh();
  mesh.floor = floor;
  mesh.ceiling = ceiling;
  // case 1: 3 down + 1 up
  for (int i = 0; i < floor.nv - 2; ++i) {
    for (int j = i + 1; j < floor.nv - 1; ++j) {
      for (int k = j + 1; k < floor.nv; ++k) {
        Boolean skip = false;
        //System.out.printf("Test triangle %d-%d-%d: ", i, j, k);
        Sphere s1 = circumcircle(floor.G[i], floor.G[j], floor.G[k]);
        for (int m = 0; m < floor.nv; ++m) {
          if (m != i && m != j && m != k && d2(floor.G[m], s1.O) < s1.r2) {
            //System.out.printf("floor pt %d inside, try next.\n", m, i, j, k);
            skip = true;
            break;
          }
        }
        if (skip) continue;
        //System.out.printf("Found triangle %d-%d-%d\n", i, j, k);
        mesh.floorEdges.add(new Edge(i, j));
        mesh.floorEdges.add(new Edge(i, k));
        mesh.floorEdges.add(new Edge(j, k));

        // Method 1: naive two loops
        //for (int m = 0; m < ceiling.nv; ++m) {
        //  skip = false;
        //  Sphere s = circumsphere(floor.G[i], floor.G[j], floor.G[k], ceiling.G[m]);
        //  //System.out.printf("  ceiling pt %d: ", m);
        //  for (int n = 0; n < ceiling.nv; ++n) {
        //    if (n != m && d2(ceiling.G[n], s.O) < s.r2) {
        //      //System.out.printf("ceiling pt %d inside, try next.\n", n);
        //      skip = true;
        //      break;
        //    }
        //  }
        //  if (skip) continue;
        //  System.out.printf("found 3-1 tet: %d-%d-%d-%d\n", i, j, k, m);
        //  mesh.mixedEdges.add(new Edge(i, m));
        //  mesh.mixedEdges.add(new Edge(j, m));
        //  mesh.mixedEdges.add(new Edge(k, m));
        //  mesh.floorTet.add(new int []{i, j, k, m});
        //  break;
        //}
        // Method 2: bulge
        int m_min = 0;
        float b_min = bulge(floor.G[i], floor.G[j], floor.G[k], ceiling.G[m_min]);
        //System.out.printf("  ceiling pt %d, bulge=%.2f\n", m_min, b_min);
        for (int m = 1; m < ceiling.nv; ++m) {
          float b = bulge(floor.G[i], floor.G[j], floor.G[k], ceiling.G[m]);
          //System.out.printf("  ceiling pt %d, bulge=%.2f\n", m, b);
          if (b_min > b) {
            b_min = b;
            m_min = m;
          }
        }
        //System.out.printf("found 3-1 tet: %d-%d-%d-%d\n", i, j, k, m_min);
        mesh.mixedEdges.add(new Edge(i, m_min));
        mesh.mixedEdges.add(new Edge(j, m_min));
        mesh.mixedEdges.add(new Edge(k, m_min));
        mesh.floorTets.add(new int []{i, j, k, m_min});
      }
    }
  }


  // case 2: 1 down + 3 up
  for (int i = 0; i < ceiling.nv - 2; ++i) {
    for (int j = i + 1; j < ceiling.nv - 1; ++j) {
      for (int k = j + 1; k < ceiling.nv; ++k) {
        Boolean skip = false;
        Sphere s1 = circumcircle(ceiling.G[i], ceiling.G[j], ceiling.G[k]);
        for (int m = 0; m < ceiling.nv; ++m) {
          if (m != i && m != j && m != k && d2(ceiling.G[m], s1.O) < s1.r2) {
            skip = true;
            break;
          }
        }
        if (skip) continue;
        mesh.ceilingEdges.add(new Edge(i, j));
        mesh.ceilingEdges.add(new Edge(i, k));
        mesh.ceilingEdges.add(new Edge(j, k));

        // Method 2: bulge
        int m_min = 0;
        float b_min = bulge(ceiling.G[i], ceiling.G[j], ceiling.G[k], floor.G[m_min]);
        for (int m = 1; m < floor.nv; ++m) {
          float b = bulge(ceiling.G[i], ceiling.G[j], ceiling.G[k], floor.G[m]);
          if (b_min > b) {
            b_min = b;
            m_min = m;
          }
        }
        //System.out.printf("found 1-3 tet: %d-%d-%d-%d\n", m_min, i, j, k);
        mesh.mixedEdges.add(new Edge(m_min, i));
        mesh.mixedEdges.add(new Edge(m_min, j));
        mesh.mixedEdges.add(new Edge(m_min, k));
        mesh.ceilingTets.add(new int []{m_min, i, j, k});
      }
    }
  }

  // case 3: 2 down + 2 up

  // wrong method
  //for (int i = 0; i < floor.nv - 1; ++i) {
  //  for (int j = i + 1; j < floor.nv; ++j) {
  //    for (int m = 0; m < ceiling.nv - 1; ++m) {
  //      for (int n = m + 1; n < ceiling.nv; ++n) {
  //        Boolean skip = false;
  //        Sphere s = circumsphere(floor.G[i], floor.G[j], ceiling.G[m], ceiling.G[n]);
  //        for (int k = 0; k < floor.nv; ++k) {
  //          if (k != i && k != j && d2(floor.G[k], s.O) < s.r2) {
  //            //System.out.printf("floor pt %d inside, try next.\n", n);
  //            skip = true;
  //            break;
  //          }
  //        }
  //        if (skip) continue;
  //        for (int k = 0; k < ceiling.nv; ++k) {
  //          if (k != m && k != n && d2(ceiling.G[k], s.O) < s.r2) {
  //            //System.out.printf("ceiling pt %d inside, try next.\n", n);
  //            skip = true;
  //            break;
  //          }
  //        }
  //        if (skip) continue;
  //        System.out.printf("found 2-2 tet: %d-%d-%d-%d\n", i, j, m, n);
  //        mesh.floorEdges.add(new Edge(i, j));
  //        mesh.mixedEdges.add(new Edge(i, m));
  //        mesh.mixedEdges.add(new Edge(i, n));
  //        mesh.mixedEdges.add(new Edge(j, m));
  //        mesh.mixedEdges.add(new Edge(j, n));
  //        mesh.ceilingEdges.add(new Edge(m, n));
  //        mesh.mixedTet.add(new int []{i, j, m, n});
  //        break;
  //      }
  //    }
  //  }
  //}

  mesh.removeDup();
  ArrayList<Edge> floorEdges = (ArrayList<Edge>)mesh.floorEdges.clone();
  ArrayList<Edge> ceilingEdges = (ArrayList<Edge>)mesh.ceilingEdges.clone();
  for (Edge ef : floorEdges) {
    for (Edge ec : ceilingEdges) {
      int i = ef.a, j = ef.b, m = ec.a, n = ec.b;
      Boolean skip = false;
      Sphere s = circumsphere(floor.G[i], floor.G[j], ceiling.G[m], ceiling.G[n]);
      for (int k = 0; k < floor.nv; ++k) {
        if (k != i && k != j && d2(floor.G[k], s.O) < s.r2) {
          //System.out.printf("floor pt %d inside, try next.\n", n);
          skip = true;
          break;
        }
      }
      if (skip) continue;
      for (int k = 0; k < ceiling.nv; ++k) {
        if (k != m && k != n && d2(ceiling.G[k], s.O) < s.r2) {
          //System.out.printf("ceiling pt %d inside, try next.\n", n);
          skip = true;
          break;
        }
      }
      if (skip) continue;
      //System.out.printf("found 2-2 tet: %d-%d-%d-%d\n", i, j, m, n);
      mesh.floorEdges.add(new Edge(i, j));
      mesh.mixedEdges.add(new Edge(i, m));
      mesh.mixedEdges.add(new Edge(i, n));
      mesh.mixedEdges.add(new Edge(j, m));
      mesh.mixedEdges.add(new Edge(j, n));
      mesh.ceilingEdges.add(new Edge(m, n));
      mesh.mixedTets.add(new int []{i, j, m, n});
      break;
    }
  }

  mesh.removeDup();
  return mesh;
}