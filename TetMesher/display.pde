void showCircumcircle(pt A, pt B, pt C) {
  Sphere s = circumcircle(A, B, C);
  fill(red, 100);
  show(A, rb+8);
  fill(green, 100);
  show(B, rb+8);
  fill(blue, 100);
  show(C, rb+8);
  fill(black, 100);
  show(s.O, rb+8);
  fill(magenta, 100);
  show(s.O, s.r());
}

void showCircumsphere(pt A, pt B, pt C, pt D) {
  Sphere s = circumsphere(A, B, C, D);
  fill(red, 100);
  show(A, rb+8);
  fill(green, 100);
  show(B, rb+8);
  fill(blue, 100);
  show(C, rb+8);
  fill(black, 100);
  show(D, rb+8);
  show(s.O, rb+8);
  fill(magenta, 100);
  show(s.O, s.r());
}

void showBulge(pt A, pt B, pt C, pt D) {
  Sphere s = circumsphere(A, B, C, D);
  fill(red, 100);
  show(A, rb+8);
  fill(green, 100);
  show(B, rb+8);
  fill(blue, 100);
  show(C, rb+8);
  fill(black, 100);
  show(D, rb+8);

  float b = bulge(A, B, C, D);
  Sphere s2 = circumcircle(A, B, C);
  fill(yellow);
  show(s.O, rb+8);
  show(s2.O, rb);
  //println(-sign(dot(s2.O, s.O, D)), b);
  pt P = P(s2.O, b*sign(dot(s2.O, s.O, D)), U(s2.O, s.O));
  //println(s.O, s2.O, P, b, U(s2.O, s.O));
  show(P, rb);
  beam(s2.O, P, 10);

  fill(magenta, 100);
  show(s.O, s.r());
}

void showTriangle(pt A, pt B, pt C) {
  float r = 0.075;
  pt A2 = P(1-2*r, A, r, B, r, C);
  pt B2 = P(r, A, 1-2*r, B, r, C);
  pt C2 = P(r, A, r, B, 1-2*r, C);
  fill(yellow, 100);
  triangle(A2, B2, C2);
}

void showTetrahedra(pt A, pt B, pt C, pt D) {
  showTriangle(A, B, C);
  showTriangle(A, B, D);
  showTriangle(A, C, D);
  showTriangle(B, C, D);
}

void testTriangleMesh(pts vertices) {
  TriangleMesh mesh = triangulate(vertices);
  println("Triangle mesh edge#=" + mesh.edges.size(), "triangle#=" + mesh.triangles.size());
  for (int [] t : mesh.triangles) {
    showTriangle(mesh.vertices.G[t[0]], mesh.vertices.G[t[1]], mesh.vertices.G[t[2]]);
  }
}

void testSampleSphere(pt P, int s) {
  ArrayList<Vertex> vs = sampleSphere(P, 100, s);
  fill(red);
  for (Vertex v : vs) {
    show(v.p, 3);
  }
  fill(blue, 100);
  show(P, 100);
}

void testSampleTube(pt P, pt Q, float r, int s) {
  ArrayList<Vertex> vs = sampleTube(P, Q, r, s);
  fill(red);
  for (Vertex v : vs) {
    show(v.p, 3);
  }
  fill(blue, 100);
  collar(P, V(P, Q), r, r);
}

void testTriangulateSphere(pt P, float r) {
  ArrayList<Vertex> v = new ArrayList<Vertex>();
  //fill(blue, 100);
  //show(P, r);
  v.addAll(v.size(), sampleSphere(P, r, 1));
  fill(red);
  //int i = 0;
  for (Vertex vv : v) {
    show(vv.p, rm);
    //text(i++, vv.p.x, vv.p.y, vv.p.z);
    //arrow(vv.p, r/2, vv.n, rm);
  }       
  //Vertex A = v.get(0), B = v.get(6), C = v.get(7);
  //arrow(A.p, r/2, A.n, rm);
  //fill(yellow);
  //arrow(B.p, r/2, B.n, rm);
  //fill(blue);
  //arrow(C.p, r/2, C.n, rm);
  
  //Sphere s = circumsphere(A, B, C, r);
  //fill(blue);
  //show(s.O, rm);       
  //fill(red, 100);
  //show(s.O, s.r());   
  //vec n = N(A.p, B.p, C.p);
  //arrow(A.p, B.p, rm);
  //arrow(A.p, C.p, rm);
  //arrow(P(A.p, B.p, C.p), 50, U(n), rm);
  //println(compatible(A, B, C), compatible(A, C, B));
  //println(A.p, B.p, C.p, N(V(A.p, B.p), V(A.p, C.p)));
  
  TriangleMesh2 triMesh = new TriangleMesh2(v, r);
  for (Edge e : triMesh.edges) beam(triMesh.v.get(e.a).p, triMesh.v.get(e.b).p, rm);
}

void testTriangulateTube(pt P, pt Q, float r) {
  ArrayList<Vertex> v = new ArrayList<Vertex>();
  v.addAll(v.size(), sampleTube(P, Q, r, 10));

  fill(red);
  int i = 0;
  for (Vertex vv : v) {
    show(vv.p, rm);
    text(i++, vv.p.x, vv.p.y, vv.p.z);
    //arrow(vv.p, r/2, vv.n, rm);
  }
  Vertex A = v.get(0), B = v.get(1), C = v.get(9);
  arrow(A.p, r/2, A.n, rm);
  fill(yellow);    
  arrow(B.p, r/2, B.n, rm);
  fill(blue);
  arrow(C.p, r/2, C.n, rm);
  
  Sphere s = circumsphere(A, B, C, r);
  fill(blue);
  show(s.O, 40);           
  //fill(red, 100);
  //show(s.O, s.r());   
  vec n = N(A.p, B.p, C.p);
  arrow(A.p, B.p, rm);
  arrow(A.p, C.p, rm);
  arrow(P(A.p, B.p, C.p), 50, U(n), rm);
  println(compatible(A, B, C), compatible(A, C, B));
  println(A.p, B.p, C.p, N(V(A.p, B.p), V(A.p, C.p)));
  println(dot(n, A.n), dot(n, B.n), dot(n, C.n));
  
  fill(blue, 100);
  collar(P, V(P, Q), r, r);
  println("Sampled", v.size(), "vertices.");
  TriangleMesh2 triMesh = new TriangleMesh2(v, r+10);
  for (Edge e : triMesh.edges) beam(triMesh.v.get(e.a).p, triMesh.v.get(e.b).p, rm); //<>//
}

void testTriangulateTubeAndSphere(pt P, pt Q, float r) {
  ArrayList<Vertex> v = new ArrayList<Vertex>();
  v.addAll(v.size(), sampleTube(P, Q, r, 10));
  v.addAll(v.size(), sampleSphere(P, r, 1));
  v.addAll(v.size(), sampleSphere(Q, r, 1));

  fill(red);
  int i = 0;
  for (Vertex vv : v) {
    show(vv.p, rm);
    //text(i++, vv.p.x, vv.p.y, vv.p.z);
    //arrow(vv.p, r/2, vv.n, rm);
  }
  //Vertex A = v.get(0), B = v.get(1), C = v.get(9);
  //arrow(A.p, r/2, A.n, rm);
  //fill(yellow);    
  //arrow(B.p, r/2, B.n, rm);
  //fill(blue);
  //arrow(C.p, r/2, C.n, rm);
  
  //Sphere s = circumsphere(A, B, C, r);
  //fill(blue);
  //show(s.O, 40);           
  ////fill(red, 100);
  ////show(s.O, s.r());   
  //vec n = N(A.p, B.p, C.p);
  //arrow(A.p, B.p, rm);
  //arrow(A.p, C.p, rm);
  //arrow(P(A.p, B.p, C.p), 50, U(n), rm);
  //println(compatible(A, B, C), compatible(A, C, B));
  //println(A.p, B.p, C.p, N(V(A.p, B.p), V(A.p, C.p)));
  //println(dot(n, A.n), dot(n, B.n), dot(n, C.n));
  
  //fill(blue, 100);
  //collar(P, V(P, Q), r, r);
  //show(P, r);
  //show(Q, r);
  println("Sampled", v.size(), "vertices.");
  TriangleMesh2 triMesh = new TriangleMesh2(v, r+10);
  for (Edge e : triMesh.edges) beam(triMesh.v.get(e.a).p, triMesh.v.get(e.b).p, rm);
}

void testTriangulateTubeAndSphere2(pt P, pt Q, pt R, float r) {
  ArrayList<Vertex> v = new ArrayList<Vertex>();
  v.addAll(v.size(), sampleTube(P, Q, r, 10));
  v.addAll(v.size(), sampleTube(P, R, r, 10));
  v.addAll(v.size(), sampleTube(R, Q, r, 10));
  v.addAll(v.size(), sampleSphere(P, r, 1));
  v.addAll(v.size(), sampleSphere(Q, r, 1));
  v.addAll(v.size(), sampleSphere(R, r, 1));

  fill(red);
  int i = 0;
  for (Vertex vv : v) {
    show(vv.p, rm);
    //text(i++, vv.p.x, vv.p.y, vv.p.z);
    //arrow(vv.p, r/2, vv.n, rm);
  }
  println("Sampled", v.size(), "vertices.");
  TriangleMesh2 triMesh = new TriangleMesh2(v, r);
  for (Edge e : triMesh.edges) beam(triMesh.v.get(e.a).p, triMesh.v.get(e.b).p, rm);
}

void testTriangulateTubeAndSphere3(pt P, pt Q, pt R, pt S, float r) {
  ArrayList<Vertex> v = new ArrayList<Vertex>();
  v.addAll(v.size(), sampleTube(P, Q, r, 10));
  v.addAll(v.size(), sampleTube(P, R, r, 10));
  v.addAll(v.size(), sampleTube(R, Q, r, 10));
  v.addAll(v.size(), sampleTube(P, S, r, 10));
  v.addAll(v.size(), sampleTube(Q, S, r, 10));
  v.addAll(v.size(), sampleSphere(P, r, 1));
  v.addAll(v.size(), sampleSphere(Q, r, 1));
  v.addAll(v.size(), sampleSphere(R, r, 1));
  v.addAll(v.size(), sampleSphere(S, r, 1));

  fill(red);
  int i = 0;
  for (Vertex vv : v) {
    show(vv.p, rm);
    //text(i++, vv.p.x, vv.p.y, vv.p.z);
    //arrow(vv.p, r/2, vv.n, rm);
  }
  println("Sampled", v.size(), "vertices.");
  TriangleMesh2 triMesh = new TriangleMesh2(v, r);
  for (Edge e : triMesh.edges) beam(triMesh.v.get(e.a).p, triMesh.v.get(e.b).p, rm);
}

TriangleMesh2 testTriangulateTubeAndSphere4(pt P, pt Q, pt R, pt S, pt T, float r) {
  ArrayList<Vertex> v = new ArrayList<Vertex>();
  v.addAll(v.size(), sampleTube(P, Q, r, 10));
  v.addAll(v.size(), sampleTube(P, R, r, 10));
  v.addAll(v.size(), sampleTube(R, Q, r, 10));
  
  v.addAll(v.size(), sampleTube(P, S, r, 10));
  v.addAll(v.size(), sampleTube(Q, S, r, 10));
  v.addAll(v.size(), sampleTube(P, T, r, 10));
  v.addAll(v.size(), sampleTube(Q, T, r, 10));
  v.addAll(v.size(), sampleTube(R, T, r, 10));
  
  v.addAll(v.size(), sampleTube(S, T, r, 10));
  
  v.addAll(v.size(), sampleSphere(P, r, 1));
  v.addAll(v.size(), sampleSphere(Q, r, 1));
  v.addAll(v.size(), sampleSphere(R, r, 1));
  v.addAll(v.size(), sampleSphere(S, r, 1));
  v.addAll(v.size(), sampleSphere(T, r, 1));

  fill(red);
  int i = 0;
  for (Vertex vv : v) {
    show(vv.p, rm);
    //text(i++, vv.p.x, vv.p.y, vv.p.z);
    //arrow(vv.p, r/2, vv.n, rm);
  }
  println("Sampled", v.size(), "vertices.");
  TriangleMesh2 triMesh = new TriangleMesh2(v, r);
  return triMesh;
}

void showTriangleMesh(TriangleMesh2 mesh) {
  fill(red);
  int i = 0;
  for (Vertex vv : mesh.v) {
    show(vv.p, rm);
    //text(i++, vv.p.x, vv.p.y, vv.p.z);
    //arrow(vv.p, r/2, vv.n, rm);
  }
  for (Edge e : mesh.edges) beam(mesh.v.get(e.a).p, mesh.v.get(e.b).p, rm);
}