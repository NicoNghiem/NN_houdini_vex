#ifndef __triangulations_h__
#define __triangulations_h__

#include <geometryhelpers.h>

//Code for working on triangulations
//TODO DOC

// EDGE FLIP OPERATIONS

function vector2 hedge_flip(int geo ; int hedge){
    //given an half edge, perform an edge flip
    //return -1 if failed, 1 if success

    int hedge_opposite = hedge_opposite(geo,hedge);
    int face0 = hedge_prim(geo,hedge);
    int face1 = hedge_prim(geo,hedge_opposite);

    //we need a triangulation
    if(len(primpoints(geo,face0))!=3 || len(primpoints(geo,face1))!=3){
        warning("Not triangles");
        return -1;
    }

    //find pt0 and pt1 which are the unshared pts
    int pt0 = hedge_dstpoint(geo,hedge_next(geo,hedge));
    int pt1 = hedge_dstpoint(geo,hedge_next(geo,hedge_opposite));

    //shared pts
    int spt0 = hedge_srcpoint(geo,hedge);
    int spt1 = hedge_srcpoint(geo,hedge_opposite);

    //create new faces
    int newface0 = addprim(geo,'poly');
    int newface1 = addprim(geo,'poly');

    addvertex(0,newface0,pt1);
    addvertex(0,newface0,pt0);
    addvertex(0,newface0,spt0);

    addvertex(0,newface1,pt0);
    addvertex(0,newface1,pt1);
    addvertex(0,newface1,spt1);

    //remove the old faces
    removeprim(geo,face0,0);
    removeprim(geo,face1,0);

    vector2 result = set(newface0,newface1);
    return result;
}

function vector2 face_edge_flip(int geo ; int face0,face1){
    // given two triangular faces, perform a edge flip if possible
    // return -1,-1 if fails, handle to the faces if success

    //we need a triangulation
    if(len(primpoints(geo,face0))!=3 || len(primpoints(geo,face1))!=3){
        warning("Not triangles");
        return {-1,-1};
    }

    //check if face0 and face1 shares an edge
    int shared_hedge = common_hedge(geo , face0 , face1);
    if(shared_hedge==-1){
        warning("Face "+itoa(face0)+" and face "+itoa(face1)+" do not share an edge");
        return {-1,-1};
    }
    vector2 result = hedge_flip(geo,shared_hedge);
    return result;
}


#endif
