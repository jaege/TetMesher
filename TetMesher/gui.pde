void keyPressed() {
  //if(key == '`') picking = true;
  if (key == 'p') R.perturb(1);
  if (key == 'P') R.perturb(200);
  if (key == '?') scribeText = !scribeText;
  if (key == '!') snapPicture();
  if (key == '~') filming = !filming;
  if (key == '.') showBalls = !showBalls;
  if (key == '\\') showTube = !showTube;  // toggle showing the tubes
  if (key == 'f') {
    flipped = !flipped;
    if (flipped) {
      R = Q;
      S = P;
      h = h_ceiling;
    } else {
      R = P;
      S = Q;
      h = h_floor;
    }
  }
  if (key == '1') showFloorTube = !showFloorTube;
  if (key == '2') showCeilingTube = !showCeilingTube;
  if (key == '3') showMixedTube = !showMixedTube;
  if (key == '4') showFloorTet = !showFloorTet;
  if (key == '5') showCeilingTet = !showCeilingTet;
  if (key == '6') showMixedTet = !showMixedTet;
  if (key == '7') showTriangleFace = !showTriangleFace;
  if (key == '9') showBulge = !showBulge;
  if (key == '0') showByFrame = !showByFrame;
  if (key == 'c') {
    println("Calculating triangle mesh:");
    //testTriangulateTube(Q.G[1], Q.G[2], 100);
    //testMesh = testTriangulateTubeAndSphere4(Q.G[1], Q.G[2], Q.G[0], P.G[0], P.G[1], rb);
    triangleMesh = tetMesh.triangulate();
    createdMesh = true;
  }
  if (key == 'j') showVertices = !showVertices;
  if (key == 'k') showMesh = !showMesh;
  if (key == '+') frm += 1;
  if (key == '-') frm = frm > 1 ? frm - 1 : 0;
  if (key == 'h') rt = rb/2;
  if (key == 'H') rt = rb;
  if (key == 'q') Q.copyFrom(P);
  if (key == 'd') {
    R.set_pv_to_pp();
    R.deletePicked();
  }
  if (key == 'W') {
    P.savePts("data/pts");
    Q.savePts("data/pts2");
  }  // save vertices to pts2
  if (key == 'L') {
    P.loadPts("data/pts");
    Q.loadPts("data/pts2");
  }  // loads saved model
  if (key == 'w') P.savePts("data/pts");  // save vertices to pts
  if (key == 'l') P.loadPts("data/pts");
  if (key == 'a') {
    animating = !animating;
  }  // toggle animation
  if (key == '|') {
    P.setZ(h_floor);
    Q.setZ(h_ceiling);
  }  // project all sites on their respective plane
  if (key == '#') exit();
  change = true;  // to save a frame for the movie when user pressed a key
}

void mouseWheel(MouseEvent event) {
  dz -=  event.getAmount();
  change = true;
}

void mousePressed() {
  //if (!keyPressed) picking = true;
  if (!keyPressed) {
    R.set_pv_to_pp();
    println("picked vertex "+R.pp);
  }
  if (keyPressed && key == 'a') {
    R.addPt(Of);
  }
  //if (!keyPressed) P.setPicked();
  change = true;
}

void mouseMoved() {
  if (keyPressed && key == ' ') {
    rx -= PI*(mouseY-pmouseY)/height;
    ry += PI*(mouseX-pmouseX)/width;
  }
  if (keyPressed && key == '`') dz += (float)(mouseY-pmouseY);  // approach view (same as wheel)
  change = true;
}

void mouseDragged() {
  if (!keyPressed) R.setPickedTo(Of);
  //  if (!keyPressed) { Of.add(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); }
  if (keyPressed && key == CODED && keyCode == SHIFT) {
    Of.add(ToK(V((float)(mouseX-pmouseX), (float)(mouseY-pmouseY), 0)));
  }
  if (keyPressed && key == 'x') R.movePicked(ToIJ(V((float)(mouseX-pmouseX), (float)(mouseY-pmouseY), 0)));
  if (keyPressed && key == 'z') R.movePicked(ToK(V((float)(mouseX-pmouseX), (float)(mouseY-pmouseY), 0)));
  if (keyPressed && key == 'X') R.moveAll(ToIJ(V((float)(mouseX-pmouseX), (float)(mouseY-pmouseY), 0)));
  if (keyPressed && key == 'Z') R.moveAll(ToK(V((float)(mouseX-pmouseX), (float)(mouseY-pmouseY), 0)));
  if (keyPressed && key == 't') {  // move focus point on plane
    if (center) F.sub(ToIJ(V((float)(mouseX-pmouseX), (float)(mouseY-pmouseY), 0)));
    else F.add(ToIJ(V((float)(mouseX-pmouseX), (float)(mouseY-pmouseY), 0)));
  }
  if (keyPressed && key == 'T') {  // move focus point vertically
    if (center) F.sub(ToK(V((float)(mouseX-pmouseX), (float)(mouseY-pmouseY), 0)));
    else F.add(ToK(V((float)(mouseX-pmouseX), (float)(mouseY-pmouseY), 0)));
  }
  change = true;
}

// **** Header, footer, help text on canvas
void displayHeader() {  // Displays title and authors face on screen
  pen(black, 1);
  scribeHeader(title, 0);
  scribeHeaderRight(name);
  fill(white);
  //image(YuanPix, width-YuanPix.width/2, 25, YuanPix.width/2, YuanPix.height/2);
  //image(QuPix, width-YuanPix.width*0.6-QuPix.width/2, 25, QuPix.width/2, QuPix.height/2);
}

void displayFooter() {  // Displays help text at the bottom
  pen(black, 1);
  scribeFooter(guide, 1);
  scribeFooter(menu, 0);
}

String
  title = "Delaunay Tetrahedralization Constructions and Smooth Boundary Mesh Rendering", 
  name = "Yu Qu,  Quan Yuan", 
  menu = "?:help, t/T:move view, space:rotate view, `/wheel:zoom, !:picture, ~:(start/stop) filming,  #:quit", 
  guide = "c:calculate mesh, click&drag:pick&slide, f:flip ceiling/floor, x/X:move picked/all, p/P:perturb, X:slide All, |:snap heights, l/L:load, w/W:write";  // user's guide