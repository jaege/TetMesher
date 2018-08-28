// ******************* LITM: Layer-Interpolating Tet Mesh, 2017 ***********************
Boolean
  animating = true, 
  pickedFocus = false, 
  center = true, 
  //track = false, 
  //showViewer = false, 
  showBalls = true, 
  //showControl = true, 
  //showCurve = true, 
  //showPath = true, 
  //showKeys = true, 
  //showSkater = false, 
  //scene1 = false, 
  //solidBalls = false, 
  //showCorrectedKeys = true, 
  //showQuads = true, 
  //showVecs = true,
  createdMesh = false, 
  showVertices = false, 
  showMesh = true, 
  showTube = true, 
  showFloorTet = true, 
  showCeilingTet = true, 
  showMixedTet = true, 
  showFloorTube = true, 
  showCeilingTube = true, 
  showMixedTube = true, 
  showBulge = false, 
  showByFrame = true, 
  showTriangleFace = false, 
  flipped = false;
float
  h_floor = 0, h_ceiling = 600, h = h_floor, 
  t = 0, dt = 0.1, 
  rb = 30, rt = rb, // radius of the balls and tubes
  rm = 2;

int frm = 0;
//int f = 0, maxf = 2*30, level = 4, method = 5;
//String SDA = "angle";
//float defectAngle = 0;
pts P = new pts();  // polyloop in 3D
pts Q = new pts();  // second polyloop in 3D
pts R, S;
TetMesh tetMesh;
TriangleMesh2 triangleMesh;
TriangleMesh2 testMesh;

void setup() {
  //QuPix = loadImage("data/Qu.jpg");  // load image from file Qu.jpg in folder data
  //YuanPix = loadImage("data/Yuan.jpg");  // load image from file Yuan.jpg in folder data
  textureMode(NORMAL);
  size(1200, 1200, P3D);  // P3D means that we will do 3D graphics
  P.declare();
  Q.declare();  // P is a polyloop in 3D: declared in pts
  //P.resetOnCircle(6,100); Q.copyFrom(P);  // use this to get started if no model exists on file: move points, save to file, comment this line
  P.loadPts("data/pts");
  Q.loadPts("data/pts2");  // loads saved models from file (comment out if they do not exist yet)
  noSmooth();
  frameRate(30);
  R = P;
  S = Q;

  testMesh = new TriangleMesh2();
}

void draw() {
  if (animating) {
    t += dt;
    if (t > 1) {
      t = 0;
    }
  }

  background(255);
  hint(ENABLE_DEPTH_TEST);
  pushMatrix();  // to ensure that we can restore the standard view before writing on the canvas
  setView();  // see pick tab
  showFloor(h);  // draws dance floor as yellow mat
  doPick();  // sets Of and axes for 3D GUI (see pick Tab)
  R.SETppToIDofVertexWithClosestScreenProjectionTo(Mouse());  // for picking (does not set P.pv)


  if (showBalls) {
    fill(orange);
    P.drawBalls(rb);
    fill(green);
    Q.drawBalls(rb);
    fill(red, 100);
    R.showPicked(rb+5);

    //fill(black);
    //textSize(26);
    //P.drawLabels(rb+10);
    //Q.drawLabels(rb+10);
  }

  //mesh = tetrahedralize_v1(P, Q);
  tetMesh = tetrahedralize_v2(P, Q);

  if (showTube) {
    if (showMixedTube) {
      fill(grey);
      for (Edge e : tetMesh.mixedEdges) beam(P.G[e.a], Q.G[e.b], rt);
    }
    if (showFloorTube) {
      fill(orange);
      for (Edge e : tetMesh.floorEdges) beam(P.G[e.a], P.G[e.b], rt);
    }   
    if (showCeilingTube) {
      fill(green);
      for (Edge e : tetMesh.ceilingEdges) beam(Q.G[e.a], Q.G[e.b], rt);
    }

    if (showBulge) {
      if (showByFrame) {
        if (showFloorTet) {
          int[] t = tetMesh.floorTets.get(frm%tetMesh.floorTets.size());
          showBulge(P.G[t[0]], P.G[t[1]], P.G[t[2]], Q.G[t[3]]);
        }
        if (showCeilingTet) {
          int[] t = tetMesh.ceilingTets.get(frm%tetMesh.ceilingTets.size());
          showBulge(Q.G[t[1]], Q.G[t[2]], Q.G[t[3]], P.G[t[0]]);
        }
        if (showMixedTet) {
          int[] t = tetMesh.mixedTets.get(frm%tetMesh.mixedTets.size());
          showCircumsphere(P.G[t[0]], P.G[t[1]], Q.G[t[2]], Q.G[t[3]]);
        }
      } else {
        if (showFloorTet) for (int[] t : tetMesh.floorTets) showBulge(P.G[t[0]], P.G[t[1]], P.G[t[2]], Q.G[t[3]]);
        if (showCeilingTet) for (int[] t : tetMesh.ceilingTets) showBulge(Q.G[t[1]], Q.G[t[2]], Q.G[t[3]], P.G[t[0]]);
        if (showMixedTet) for (int[] t : tetMesh.mixedTets) showCircumsphere(P.G[t[0]], P.G[t[1]], Q.G[t[2]], Q.G[t[3]]);
      }
    }

    if (showTriangleFace) {
      if (showByFrame) {
        if (showFloorTet) {
          int[] t = tetMesh.floorTets.get(frm%tetMesh.floorTets.size());
          showTetrahedra(P.G[t[0]], P.G[t[1]], P.G[t[2]], Q.G[t[3]]);
        }
        if (showCeilingTet) {
          int[] t = tetMesh.ceilingTets.get(frm%tetMesh.ceilingTets.size());
          showTetrahedra(Q.G[t[1]], Q.G[t[2]], Q.G[t[3]], P.G[t[0]]);
        }
        if (showMixedTet) {
          int[] t = tetMesh.mixedTets.get(frm%tetMesh.mixedTets.size());
          showTetrahedra(P.G[t[0]], P.G[t[1]], Q.G[t[2]], Q.G[t[3]]);
        }
      } else {
        if (showFloorTet) for (int[] t : tetMesh.floorTets) showTetrahedra(P.G[t[0]], P.G[t[1]], P.G[t[2]], Q.G[t[3]]);
        if (showCeilingTet) for (int[] t : tetMesh.ceilingTets) showTetrahedra(Q.G[t[1]], Q.G[t[2]], Q.G[t[3]], P.G[t[0]]);
        if (showMixedTet) for (int[] t : tetMesh.mixedTets) showTetrahedra(P.G[t[0]], P.G[t[1]], Q.G[t[2]], Q.G[t[3]]);
      }
    }
  }

  if (createdMesh) {
    if (showVertices) {
      fill(red);
      //int i = 0;
      for (Vertex v : triangleMesh.v) {
        show(v.p, rm);
        //text(i++, v.p.x, v.p.y, v.p.z);
        //arrow(v.p, 50, v.n, rm);
      }
    }
    if (showMesh) {
      fill(red);
      for (Edge e : triangleMesh.edges) beam(triangleMesh.v.get(e.a).p, triangleMesh.v.get(e.b).p, rm);
    }
  }

  //testTriangulateTubeAndSphere(Q.G[1], Q.G[2], rb);
  //testTriangulateTubeAndSphere2(Q.G[1], Q.G[2], Q.G[0], rb);
  //testTriangulateTubeAndSphere3(Q.G[1], Q.G[2], Q.G[0], P.G[0], rb);
  //showTriangleMesh(testMesh);

  popMatrix();  // done with 3D drawing. Restore front view for writing text on canvas
  hint(DISABLE_DEPTH_TEST);  // no z-buffer test to ensure that help text is visible

  scribeHeader("Site count: "+P.nv+" floor + "+Q.nv+" ceiling", 1);
  scribeHeader("Beam count: "+tetMesh.floorEdges.size()+" floor + "+tetMesh.ceilingEdges.size()+" ceiling + "+tetMesh.mixedEdges.size()+" mixed", 2);
  scribeHeader("Triangle count: "+tetMesh.floorTriangles.size()+" floor + "+tetMesh.ceilingTriangles.size()+" ceiling", 3);
  scribeHeader("Tet count: "+tetMesh.floorTets.size()+" floor + "+tetMesh.ceilingTets.size()+" ceiling + "+tetMesh.mixedTets.size()+" mixed", 4);

  // used for demos to show red circle when mouse/key is pressed and
  // what key (disk may be hidden by the 3D model)
  if (mousePressed) {
    pen(cyan, 3);
    noFill();
    ellipse(mouseX, mouseY, 20, 20);
  }
  if (keyPressed) {
    pen(red, 1);
    fill(white);
    ellipse(mouseX+14, mouseY-10, 26, 26);
    fill(red);
    text(key, mouseX+14-4, mouseY-10+4);
  }
  if (scribeText) displayHeader();  // dispalys header on canvas, including my face
  if (scribeText && !filming) displayFooter();  // shows menu at bottom, only if not filming
  if (filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++, 4)+".tif");  // save next frame to make a movie
  change = false;  // to avoid capturing frames when nothing happens (change is set uppn action)
  change = true;
}