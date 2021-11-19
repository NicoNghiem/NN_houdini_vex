#ifndef __anisotropy_h__
#define __anisotropy_h__

// Codebase to do helper operations for anisotropy/ellipsoidal computations


//Building anisotropic tensors

matrix3 build_tensor_q(vector N,up,scale){
  //build q, square root of Q = inv(G)
  vector Oz = normalize(N);
  vector Oy = up - dot(Oz,up) * Oz;
  matrix3 O = maketransform(Oz,Oy);
  matrix3 D = ident();
  for(int i = 0 ; i < 3 ; i++){
      setcomp(D,scale[i],i,i);
  }
  return transpose(O) * D * O;
}

matrix3 build_tensor_Q(vector N,up,scale){
  //build Q tensor
  matrix3 q = build_tensor_q(N,up,scale);
  return q * q;
}

matrix3 build_tensor_F(vector N,up,scale){
  //build F one of G square roots
  vector Oz = normalize(N);
  vector Oy = up - dot(Oz,up) * Oz;
  matrix3 O = maketransform(Oz,Oy);
  matrix3 D = ident();
  for(int i = 0 ; i < 3 ; i++){
      setcomp(D,1./scale[i],i,i);
  }
  return transpose(O) * D * O;
}

matrix3 build_tensor_G(vector N,up,scale){
  //build the G tensor such that ellipsoid is xtGx =1
  matrix3 F = build_tensor_F(N,up,scale);
  return F*F;
}

//version using directly geo index and point id
matrix3 build_tensor_q(int geo ; int i){
  vector N = point(geo,'N',i);
  vector up = point(geo,'up',i);
  vector scale = point(geo,'scale',i);

  N = normalize(N);
  up = normalize(up - dot(up,N) * N);
  return build_tensor_q(N,up,scale);
}

matrix3 build_tensor_Q(int geo ; int i){
  vector N = point(geo,'N',i);
  vector up = point(geo,'up',i);
  vector scale = point(geo,'scale',i);

  N = normalize(N);
  up = normalize(up - dot(up,N) * N);
  return build_tensor_Q(N,up,scale);
}

matrix3 build_tensor_F(int geo ; int i){
  vector N = point(geo,'N',i);
  vector up = point(geo,'up',i);
  vector scale = point(geo,'scale',i);

  N = normalize(N);
  up = normalize(up - dot(up,N) * N);
  return build_tensor_F(N,up,scale);
}

matrix3 build_tensor_G(int geo ; int i){
  vector N = point(geo,'N',i);
  vector up = point(geo,'up',i);
  vector scale = point(geo,'scale',i);

  N = normalize(N);
  up = normalize(up - dot(up,N) * N);
  return build_tensor_G(N,up,scale);
}

//distances

float anisotropic_distance(vector P0,P1 ; vector N,up,scale){
  //distance in the space with the tensor from (N,up,scale)
  matrix3 F = build_tensor_F(N,up,scale);
  return length2((P0-P1)*F);
}

//ELLIPSOID / Point interaction

float point_ellipsoid_collision(vector P1,C,N,up,scale ; vector push){
  //compute the anisotropic distance and the vector necessary to push the point to solve collision
  //push it the vector such that P1+push is on the ellipsoid
  matrix3 F = build_tensor_F(N,up,scale);
  matrix3 q = invert(F);
  vector diff = P1-C;
  float dist = anisotropic_distance(C,P1,N,up,scale);
  float invsqrt = 1. / sqrt(dist);
  push = C + diff * invsqrt - P1;
  return dist;
}

float point_ellipsoid_collision(int geoPoints, geoEllipsoid; int id_point, id_ellipsoid ; vector push){
  vector C = point(geoEllipsoid,'P',id_ellipsoid);
  vector N = point(geoEllipsoid,'N',id_ellipsoid);
  vector up = point(geoEllipsoid,'up',id_ellipsoid);
  vector scale = point(geoEllipsoid,'scale',id_ellipsoid);

  vector P1 = point(geoPoints,'P',id_point);
  return point_ellipsoid_collision(P1,C,N,up,scale,push);
}

//ELLIPSOID / ELLIPSOID interaction

matrix3 ellipsoid_q_mink_min_approx(vector P0,P1,N0,N1,up0,up1,scale0,scale1){
  //compute the q tensor of the min Ellipsoidal approximation of the Minkowski sum of both ellipsoids
  vector dir = normalize(P0-P1);

  matrix3 q0 = build_tensor_q(N0,up0,scale0);
  matrix3 q1 = build_tensor_q(N1,up1,scale1);

  matrix3 S0 = dihedral(dir,normalize(q0*dir));
  matrix3 S1 = dihedral(dir,normalize(q1*dir));

  matrix3 qstar = S0 * q0 + S1 * q1;

  return qstar;
}

float ellipsoid_collision_mink(vector P0,P1,N0,N1,up0,up1,scale0,scale1 ; vector push){
  //compute ellipsoid collision via minkowski sum approximators
  // push is the vector such that P0 translated by push solves the collision
  matrix3 qstar = ellipsoid_q_mink_min_approx(P0,P1,N0,N1,up0,up1,scale0,scale1);
  matrix3 Fstar = invert(qstar);

  float dist = length(Fstar*(P0-P1));

  vector dir = qstar * normalize(Fstar*(P0-P1));
  push = (P1 + dir - P0);

  return dist;
}

float ellipsoid_collision_mink(int geo ; int i0,i1 ; vector push){
  vector P0 = point(geo,'P',i0);
  vector up0 = point(geo,'up',i0);
  vector N0 = point(geo,'N',i0);
  vector scale0 = point(geo,'scale',i0);
  vector P1 = point(geo,'P',i1);
  vector up1 = point(geo,'up',i1);
  vector N1 = point(geo,'N',i1);
  vector scale1 = point(geo,'scale',i1);

  return ellipsoid_collision_mink(P0,P1,N0,N1,up0,up1,scale0,scale1,push);
}

float ellipsoid_collision_mink(int geo0,geo1 ; int i0,i1 ; vector push){
  vector P0 = point(geo0,'P',i0);
  vector up0 = point(geo0,'up',i0);
  vector N0 = point(geo0,'N',i0);
  vector scale0 = point(geo0,'scale',i0);
  vector P1 = point(geo1,'P',i1);
  vector up1 = point(geo1,'up',i1);
  vector N1 = point(geo1,'N',i1);
  vector scale1 = point(geo1,'scale',i1);

  return ellipsoid_collision_mink(P0,P1,N0,N1,up0,up1,scale0,scale1,push);
}

//


#endif
