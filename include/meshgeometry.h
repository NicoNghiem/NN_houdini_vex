#ifndef __meshgeometry_h__
#define __meshgeometry_h__

#include <geometryhelpers.h>

//Basic Geometry Operations

//GEOMETRY

//face_area_vector(int geo ; int faceid)
// - INPUT: geo, handle to a geometry ; faceid: id of the face
// - OUTPUT: face normal with length equals to the area of the face

//face_normal(int geo ; int faceid)
// - INPUT: geo, handle to a geometry ; faceid: id of the face
// - OUTPUT: face normal


//CURVATURE:

//pt_gaussian_curvature_tri(int geo ; int ptid) (Gauss Bonnet preserving)
// - INPUT: geo, handle to the geometry ; ptid: id of the point
// - OUTPUT: discrete gaussian curvature at ptid

//pt_mean_curvature_vector_tri(int geo ; int ptid) (Gauss Bonnet preserving)
// - INPUT: geo, handle to the geometry ; ptid: id of the point
// - OUTPUT: discrete mean curvature vector at ptid

//BEGIN CODE

//GEOMETRY

//NORMALS

function vector face_area_vector(int geo ; int id){

    //Give the area vector of prim id

    int point_ids[] = primpoints(geo,id);
    point_ids = reverse(point_ids);
    int N = len(point_ids);

    vector result = 0;
    for(int i = 0 ; i < N ; i++){
        vector vi = point(geo,'P',point_ids[i]);
        vector vinext = point(geo,'P',point_ids[(i+1)%N]);
        result += cross(vi,vinext);
    }

    return 0.5 * result;
}

function vector face_normal(int geo ; int id){
    //return the face normal
    vector f_vec = face_area_vector(geo,id);
    return normalize(f_vec);
}

function float face_area(int geo ; int id){
    //return the area of the face
    vector f_vec = face_area_vector(geo,id);
    return length(f_vec);
}

function vector point_normal(int geo ; int id){
    //return the point normal (weighted by face areas)
    int faces[] = pointprims(geo,id);
    vector normal = 0;
    foreach(int faceid ; faces){
        normal += face_area_vector(geo,faceid);
    }
    return normalize(normal);
}

//FACES AND AREAS

function vector face_center(int geo ; int id){
    // Face center using houdini's native data structure
    int vertIds[] = primpoints(geo,id);
    vector center = 0 ;
    foreach(int ptId ; vertIds){
        center+=point(geo,'P',ptId);
    }
    return center / len(vertIds);
}

function float point_barycentric_area(int geo ; int ptid){
    //barycentric area linked to a vertex (area of the dual polygon associated to the vertex)
    int prims[] = pointprims(geo,ptid);
    float area = 0 ;

    foreach(int primid ; prims){
        area += face_area(geo,primid);
    }
    area /= 3.;
    return area;
}

//CURVATURE COMPUTATION

function float pt_gaussian_curvature_tri(int geo ; int ptid){
    // compute the discrete gaussian curvature at ptid
    // require a triangulation
    // using: http://multires.caltech.edu/pubs/diffGeoOps.pdf  by Meyer et al.

    int prims[] = pointprims(geo,ptid);
    float curv = 2 * 3.1415926535897932384;
    vector P = point(geo, "P" , ptid);

    foreach(int faceid ; prims){
        if(is_triangle(geo,faceid)==0)  error("Need a triangulation. Prim "+itoa(faceid)+" is not a triangle");
        int vtxids[] = primpoints(geo,faceid);
        removevalue(vtxids,ptid); // remove the original pt from the array
        vector P0 = point(geo,"P",vtxids[0]);
        vector P1 = point(geo, "P" , vtxids[1]);
        curv -= compute_angle(P0-P,P1-P);
    }
    return curv;
}

function vector pt_mean_curvature_vector_tri(int geo ; int ptid){
    // compute the mean curvature vector at ptid, using the gradient area vector formula
    // require a triangulation
    // using: http://multires.caltech.edu/pubs/diffGeoOps.pdf  by Meyer et al.

    vector P = point(geo,'P',ptid);
    int neigh[] = neighbours(geo,ptid);

    vector grad = 0;

    foreach(int id ; neigh){
        vector Q = point(geo,'P',id);
        int hedge = pointhedge(geo,ptid,id);
        int hedge_opp = pointhedge(geo,id,ptid);

        int pt0 = hedge_dstpoint(geo,hedge_next(geo,hedge));
        int pt1 = hedge_dstpoint(geo,hedge_next(geo,hedge_opp));

        vector P0 = point(geo,'P',pt0);
        vector P1 = point(geo,'P',pt1);

        grad += (cotan(P-P0,Q-P0) + cotan(P-P1,Q-P1)) * (P - Q);
    }
    grad *= 0.5;
    return grad;
}

#endif
