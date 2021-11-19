#ifndef __geometryhelpers_h__
#define __geometryhelpers_h__

// Collection of useful geometry helper functions

//HEDGES

// hedge_opposite(int geo ; int hedge)
// - INPUT: geo, handle to a geometry ; hedge a headge of that geometry
// - OUTPUT: the index of the flipped hedge (add a warning in case the geometry is not manifold)

// common_hedge ( int geo ; int face0 , face1)
// - INPUT: geo, handle to a geometry ; face0 and face1 the index of 2 faces
// - OUTPUT: the hedge of face0 s.t hedge flipped belongs to face1. Return -1 if the faces do not share an edge

// ANGLES

// compute_angle(vector x , y)
// - INPUT: vector x and y
// - OUTPUT: angle (in radians) between x and y

//ALGEBRA

//project(vector x , y)
// - INPUT: vector x and y
// - OUTPUT: x projected on the direction defined by y



//BEGIN CODE

//HEDGES

function int hedge_opposite(int geo ; int hedge){
    //check if the triangulation is manifold (ie hedge has valence exactly 2)
    //and return the opposite hedge
    int hedge_opposite = hedge_nextequiv(geo,hedge);
    if(hedge_opposite==hedge || hedge_nextequiv(geo,hedge_opposite)!=hedge){
        warning("Non manifold mesh");
    }
    return hedge_opposite;
}

function int common_hedge(int geo ; int face0 , face1){
    //Return the hedge of face0 that, once flipped is an hedge of face1
    //Return -1 if the faces do not share an edge
    int initial_hedge = primhedge(geo,face0);
    int c_hedge = initial_hedge;

    do{
        c_hedge = hedge_next(geo,c_hedge); // next hedge
        int eq_hedge = c_hedge; //iterate on the equivalent hedges
        do{
            eq_hedge = hedge_nextequiv(geo,eq_hedge);
            if(hedge_prim(geo,eq_hedge)==face1){
                return c_hedge; //if the opposite hedge belongs to face1 then c_hedge is the hedge we were looking for
            }
        }while(c_hedge!=eq_hedge);

    }while(initial_hedge!=c_hedge); //stop when we circled

    //if we reach this part of the code it means the faces do not share an edge
    warning("Faces do not share an edge");
    return -1;
}

function int is_triangle(int geo ; int faceid){
    int pts[] = primpoints(geo,faceid);
    return (len(pts)==3) ;
}

function int is_quad(int geo ; int faceid){
    int pts[] = primpoints(geo,faceid);
    return (len(pts)==4);
}


//ANGLES

function float compute_angle(vector x,y){
    //Compute the angle in radians between x and y
    return 2. * atan(length(x*length(y) - length(x)*y) / length(x * length(y) + length(x) * y));
}

function float cotan(vector x,y){
    //return the cotan of the angle between x and y
    return dot(x,y) / length(cross(x,y));
}

//ALGEBRA

function vector project(vector x,ydir){
    //compute x projected on direction y
    vector dir = normalize(ydir);
    return dot(x,dir) * dir;
}

function vector reject(vector x , y){
    return x - project(x,y);
}

function vector change_basis(vector x ; matrix3 basis){
    /*
    we assume basis is an orthogonal matrix
    */
    return basis * x * transpose(basis);
}


#endif
