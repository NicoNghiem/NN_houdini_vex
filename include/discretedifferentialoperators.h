#ifndef __discretedifferentialoperators_h__
#define __discretedifferentialoperators_h__

#include <matrixNM.h>
#include <geometryhelpers.h>
#include <meshgeometry.h>

// Implementation of discrete differential operators on polygonal meshes as described in:
// Discrete Differential Operators on Polygonal Meshes, De Goes, Butts, Desbrun, SIGGRAPH 2020
// https://graphics.pixar.com/library/PolyDDG/paper.pdf

//OPERATORS and QUANTITIES
//Some operators have a matrix3 houdini equivalent that allow faster computations for triangle meshes
// mat_cross(p) -> P such that P*q = cross(p,q) for all q

// QUANTITIES:
// Xf(geo,faceid) -> vertex positions (n_vert x 3)_mat
// Ef(geo,faceid) -> edge vectors (n_edge x 3) = (n_vert x 3)_mat
// Bf(geo,faceid) -> edge midpoints (n_vert x 3)_mat
// Cf(geo,faceid) -> face center (barycenter) (3x1)_mat
//TODO: Finish doc

// OPERATORS:
// Af() -> average operator
// Df() -> difference operator
//TODO: Finish doc

/**
-----
-----
Part 3 . Gradient Operator
-----
-----
*/

function matrix3 mat3_cross(vector p){

    //return the matrix P such that Pq = cross(p,q) for all q

    matrix3 result = 0;
    result.xy = -p.z;
    result.xz = p.y;
    result.yx = p.z;
    result.yz = -p.x;
    result.zx = -p.y;
    result.zy = p.x;
    return result;
}

function matrix_NM mat_cross(vector p){
    //matrix_NM version of previous function

    float values[] = {};
    append(values,0);
    append(values,-p.z);
    append(values,p.y);

    append(values,p.z);
    append(values,0);
    append(values,-p.x);

    append(values,-p.y);
    append(values,p.x);
    append(values,0);
    matrix_NM result = matrix_NM(values,3,3);
    return result;
}


function matrix3 Xf3(int geo ; int id){
    /*
    Vertex position operator for a triangle
    */
    int point_ids[] = primpoints(geo,id);
    point_ids = reverse(point_ids);
    if(len(point_ids)!=3){
        error("Expecting a triangulation");
    }
    matrix3 result = 0;
    for(int i = 0 ; i < 3 ; i++){
        vector pt = point(geo,'P',point_ids[i]);
        for(int j = 0 ; j < 3 ; j++){
            setcomp(result,pt[j],j,i);
        }
    }
    return result;
}

function matrix_NM Xf(int geo ; int id){
    //  vertex position operator: each row is the position of a vertex
    int point_ids[] = primpoints(geo,id);
    point_ids = reverse(point_ids);
    int N = len(point_ids);
    int M = 3;
    float values[] = {};
    for(int i = 0 ; i<N ; i++){
        vector pt = point(geo,'P',point_ids[i]);
        for(int j = 0 ; j < M ; j++){
            append(values,pt[j]);
    }}
    return matrix_NM(values,N,M);
}

function matrix_NM Nf(int geo ; int id){
    //  normal operator: each row is the normal of a vertex
    int point_ids[] = primpoints(geo,id);
    point_ids = reverse(point_ids);
    int N = len(point_ids);
    int M = 3;
    float values[] = {};
    for(int i = 0 ; i<N ; i++){
        vector pt = point_normal(geo,point_ids[i]);
        for(int j = 0 ; j < M ; j++){
            append(values,pt[j]);
    }}
    return matrix_NM(values,N,M);
}

function matrix3 Af3(){
    /*
    return the average operator of dimension 3
    */
    matrix3 result = 0;
    result.xx = 0.5;
    result.xy = 0.5;
    result.yy = 0.5;
    result.yz = 0.5;
    result.zz = 0.5;
    result.zx = 0.5;
    return result;
}

function matrix Af4(){

    //return the average operator of dimension 3

    matrix result = 0;
    result.xx = 0.5;
    result.xy = 0.5;
    result.yy = 0.5;
    result.yz = 0.5;
    result.zz = 0.5;
    result.zw = 0.5;
    result.ww = 0.5;
    result.wx = 0.5;
    return result;
}

function matrix_NM Af(int dim){
    //average operator of dimension dim
    float values[] = {};
    int N = dim;
    int M = dim;
    for(int i = 0 ; i<N ; i++){
    for(int j = 0 ; j < M ; j++){
        float val = (i-j==0 || (i-j)%M==M-1) ? 0.5 : 0;
        append(values,val);
    }}
    return matrix_NM(values,N,M);
}

function matrix3 Df3(){
    /*
    return the difference operator of size 3
    */
    matrix3 result = 0;
    result.xx=-1;
    result.yy=-1;
    result.zz=-1;
    result.xy=1;
    result.yz=1;
    result.zx=1;
    return result;
}

function matrix Df4(){
    /*
    return the difference operator of size 4
    */
    matrix result = 0;
    result.xx=-1;
    result.yy=-1;
    result.zz=-1;
    result.ww =-1;
    result.xy=1;
    result.yz=1;
    result.zw=1;
    result.wx=1;
    return result;
}

function matrix_NM Df(int dim){
    // difference operator
    float values[] = {};
    int N = dim;
    int M = dim;
    for(int i = 0 ; i<N ; i++){
    for(int j = 0 ; j < M ; j++){
        float val = i==j ? -1 : (i-j)%M==M-1 ? 1 : 0;
        append(values,val);
    }}
    return matrix_NM(values,N,M);
}

function matrix3 Ef3(int geo ; int id){
    //edge for triangle meshes
    return Xf3(geo,id) * Df3() ;
}

function matrix_NM Ef(int geo ; int id){
    /*
    edge operator: each row is an edge vector
    */
    matrix_NM Xf = Xf(geo,id);
    return mult(Df(Xf.N) , Xf);
}

function matrix3 Bf3(int geo ; int id){
    //edge midpoint
    return Xf3(geo,id) * Af3();
}

function matrix_NM Bf(int geo ; int id){
    /*
    edge midpoint operator: each row is an edge midpoint
    */
    matrix_NM Xf = Xf(geo,id);
    matrix_NM Af = Af(Xf.N);
    return mult(Af,Xf);
}

function matrix_NM Cf(int geo ; int id){
    // Face center - do not use it, it's slower
    // please use face_center()
    matrix_NM Xf = Xf(geo,id);
    int nf = Xf.N;
    Xf = transpose(Xf);
    matrix_NM ones = ones(nf,1);
    matrix_NM cf = mult(Xf,ones);
    cf = mult(1./nf,cf);
    return cf;
}

function vector phif3(int geo ; int id ; string attr){
    int point_ids[] = primpoints(geo,id);
    point_ids = reverse(point_ids);
    if(len(point_ids)!=3){
        error("Expecting a triangulation");
    }

    vector result = 0;

    for(int i = 0 ; i<3 ; i++){
        float val = point(geo,attr,point_ids[i]);
        setcomp(result,val,i);
    }
    return result;
}

function matrix_NM phif(int geo ; int id ; string attr){
    /*
    Return a vector where each row is the value of attribute attr at vertex i (read counter clockwise)
    ONLY for scalar attributes
    */
    int point_ids[] = primpoints(geo,id);
    point_ids = reverse(point_ids);

    float values[] = {};
    int N = len(point_ids);
    int M = 1;
    for(int i=0; i < N ; i++){
        float val = point(geo,attr,point_ids[i]);
        append(values,val);
    }
    return matrix_NM(values,N,M);
}

function matrix_NM phif_vector(int geo ; int faceid ; string attr){
    /*
    Return a matrix where each row is the value of attribute attr at vertex i (read counter clockwise)
    ONLY for vector attributes
    */
    int point_ids[] = primpoints(geo,faceid);
    point_ids = reverse(point_ids);

    float values[] = {};
    int N = len(point_ids);
    int M = 3;
    for(int i=0; i < N ; i++){
        vector val = point(geo,attr,point_ids[i]);
        for(int j = 0 ; j < M ; j++){
            append(values,getcomp(val,j));
        }
    }
    return matrix_NM(values,N,M);
}



//gradient operator

function matrix3 Gf3(int geo ; int id){
    //return (-1./length(face_area_vector(geo,id))) * Af3() * transpose(Ef3(geo,id)) * mat3_cross(face_normal(geo,id));
    return (-1.) * Af3() * transpose(Ef3(geo,id)) * mat3_cross(face_normal(geo,id));
}

function matrix_NM Gf(int geo ;int id){
    int vertex_ids[] = primpoints(geo,id);
    int dim = len(vertex_ids);
    matrix_NM nf = mat_cross(face_normal(geo,id)) ;
    //return mult((-1./face_area(geo,id)) , mult( nf , mult( transpose(Ef(geo,id)) , Af(dim) )));
    return mult(-1. , mult( nf , mult( transpose(Ef(geo,id)) , Af(dim) )));
}

/*
-----
-----
Part 4. Operators on Discrete Forms
-----
-----
*/

//flat operator

function matrix_NM Vf(int geo ; int id){
    //Flat operator: encode vector circulation along the edges
    //size: n_face x 3
    matrix_NM I = ident(3,3);
    vector normal = face_normal(geo,id);
    matrix3 out = - outerproduct(normal,normal);
    matrix_NM outNM = tomatrixNM(out);
    outNM = add(I,outNM);
    matrix_NM Ef = Ef(geo,id);
    return mult(Ef,outNM); // Vf = Ef ( I - ntn)
}

//sharp operator

function matrix_NM Uf(int geo ; int id){
    //Sharp operator: aggregate values of a discrete 1-form into
    //a tangent vector per face
    //size: 3 x n_faces
    float area = face_area(geo,id);

    matrix_NM outer_n = mat_cross(face_normal(geo,id));
    matrix_NM Bf = Bf(geo,id);
    Bf = transpose(Bf);
    matrix_NM cf = Cf(geo,id);
    cf = mult(cf , ones(1,Bf.M));
    cf = mult(-1,cf);
    return mult(1./area,mult(outer_n,add(Bf,cf)));
}

//projection operator

function matrix_NM Pf(int geo ; int id){
    //projection operator:
    //Note: Im Pf = Ker Uf and Ker Pf = Im Vf
    matrix_NM Vf = Vf(geo,id);
    matrix_NM Uf = Uf(geo,id);
    matrix_NM I = ident(Vf.N,Uf.M);
    return add(I,mult(-1.,mult(Vf,Uf)));
}

//inner product operator

function matrix_NM Mf(int geo ; int id ; float lambda){
    //inner product operator
    // the parameter lambda weights the projection term
    float area = face_area(geo,id);
    matrix_NM Uf = Uf(geo,id);
    Uf = mult(area,mult(transpose(Uf),Uf));
    matrix_NM Pf = Pf(geo,id);
    Pf = mult(lambda , mult(transpose(Pf),Pf));
    return add(Uf,Pf);
}

function matrix_NM Mf(int geo ; int id){
    // In 7.1 of Discrete Differential Operators on Polygonal Meshes
    // a value of lambda = 1 is required
    return Mf(geo,id,1.);
}

//divergence operator

function matrix_NM divf(int geo ; int id ; float lambda){
    //divergence operator: computes the divergence of a discrete 1-form
    matrix_NM Mf = Mf(geo,id,lambda);
    matrix_NM Df = Df(Mf.N);

    return mult(transpose(Df), Mf);
}

function matrix_NM divf(int geo ; int id){
    // In 7.1 of Discrete Differential Operators on Polygonal Meshes
    // a value of lambda = 1 is required
    return divf(geo,id,1.);
}

//curl operator

function matrix_NM curlf ( int geo ; int id ){
    /*
    curl operator: computes the curl of a discrete 1-form
    */
    int vertex_ids[] = primpoints(geo,id);
    int dim = len(vertex_ids);
    matrix_NM curl = ones( 1 , dim );
    return curl;
}

//Laplace-Beltrami operator

function matrix_NM Lf(int geo ; int id ; float lambda){
    matrix_NM Mf = Mf(geo,id,lambda);
    matrix_NM Df = Df(Mf.N);

    return mult(transpose(Df), mult(Mf,Df));
}

function matrix_NM Lf(int geo ; int id ){
    // In 7.1 of Discrete Differential Operators on Polygonal Meshes
    // a value of lambda = 1 is required
    return Lf(geo,id,1.);
}

/*
------
------
PART 5. OPERATORS ON DIRECTIONAL FIELDS
------
------


*/

function matrix_NM Tf(int geo ; int faceid){
    //return Tf such that (face_normal,Tf) is orthonormal frame
    //Tf is a 3x2 matrix

    vector normal = face_normal(geo,faceid);

    vector tangent = {0,0,1};
    if(length(cross(tangent,normal))==0)    tangent = {1,0,0};

    tangent = normalize(tangent - dot(tangent,normal)*normal);
    vector frame =  - cross(tangent,normal);
    vector cols[] = {};
    append(cols,tangent);
    append(cols,frame);

    return matNMfromcols(cols);

}

function matrix_NM Tv(int geo ; int ptid){
    //return Tv such that (vertex_normal,Tv) is orthonormal frame
    // Tv is a 3x2 matrix
    vector normal = point_normal(geo,ptid);

    vector tangent = {0,0,1};
    if(length(cross(tangent,normal))==0)    tangent = {1,0,0};

    tangent = normalize(tangent - dot(tangent,normal)*normal);
    vector frame = - cross(tangent,normal);
    vector cols[] = {};
    append(cols,tangent);
    append(cols,frame);
    return matNMfromcols(cols);
}

function float[] project_to_localframe (vector v ; matrix_NM T_frame){
    /*
    Project v to the local frame T_frame (we assume T_frame is indeed 3x2 and orthonormal)
    Should be a vec2 but having it directly as a float[] is easier
    */
    vector N = getcolumnvector(T_frame,0);
    vector up = getcolumnvector(T_frame,1);
    matrix3 basis = maketransform(N,up);
    vector new_vec = change_basis(v,basis);
    //vector2 result = set(new_vec.z , new_vec.y);
    float result[] = {};
    append(result,new_vec.z);
    append(result,new_vec.y);
    return result;
}

function float[] project_to_localframe ( int geo ; string attr ; int ptid){
    /*
    Project attr evaluated at ptid to the local frame T_frame (we assume T_frame is indeed 3x2 and orthonormal)
    Should be a vec2 but having it directly as a float[] is easier
    */
    vector vec = point(geo,attr,ptid);
    // make vec a tangential vector
    vector N = point_normal(geo,ptid);
    vec = reject(vec,N);
    // get local frame
    matrix_NM Tv = Tv(geo,ptid);
    return project_to_localframe(vec,Tv);
}

function matrix_NM project_to_localframe(int geo ; string attr ; int ptid){
    /*
    Project attr evaluated at ptid to the local frame T_frame (we assume T_frame is indeed 3x2 and orthonormal)
    with the matrix_NM format
    */
    int N = 2;
    int M = 1;
    float values[] = project_to_localframe(geo,attr,ptid);
    return matrix_NM(values,N,M);
}

function matrix_NM assemble_localframe_vectors (int geo ; string attr ; int faceid){
    /*
    Assemble a (1 x 2 n_vertex) column where each 2 consecutive values are the coordinates of the tangential component of attr
    projected on the local frame defined as before
    attr MUST be a vector attribute
    */
    int point_ids[] = primpoints(geo,faceid);
    point_ids = reverse(point_ids);

    float colvector[] = {};

    foreach( int ptid ; point_ids){
        float localcoord[] = project_to_localframe(geo,attr,ptid);
        append(colvector,localcoord);
    }
    return matrix_NM(colvector, 1 , 2*len(point_ids));
}

function matrix_NM discrete_levi_citia_connection(int geo ; int faceid ; int ptid){
    /*
    Returns a 2x2 matrix
    */
    matrix_NM Tf = Tf(geo,faceid);
    matrix_NM Tv = Tv(geo,ptid);
    vector nf = face_normal(geo,faceid);
    vector nv = point_normal(geo,ptid);
    matrix3 Q3 = dihedral(nv,nf);
    matrix_NM Q = tomatrixNM(Q3);

    return mult(transpose(Tf) , mult(Q,Tv));
}

function matrix_NM parallel_transported_localframe_vectors(int geo ; string attr ; int faceid){
    /*
    Returns a (n_vertex , 2) matrix containing each local frame vectors (u_v) parallel transported to the face faceid
    u^{\gamma}_f in the paper -> i.e. Rfv * uv
    */
    int point_ids[] = primpoints(geo,faceid);
    point_ids = reverse(point_ids);

    float values[] = {};

    foreach(int ptid ; point_ids){
        matrix_NM Rfv = discrete_levi_citia_connection(geo,faceid,ptid);
        matrix_NM uv = project_to_localframe(geo,attr,faceid);
        matrix_NM transported = mult(Rfv,uv);
        append(values , transported.values);
    }
    return matrix_NM(values,len(point_ids),2);
}



function matrix_NM covariant_derivative(int geo ; string attr ; int faceid){
    // parallel transport the vectors uf (assembled local face vectors) to the local f coordinate
    // return a 2x2 matrix

    matrix_NM Tf = Tf(geo,faceid);
    Tf = transpose(Tf);
    matrix_NM Gf = Gf(geo,faceid);
    matrix_NM uf_gamma = parallel_transported_localframe_vectors(geo,attr,faceid);

    matrix_NM cov_derivative = mult(Tf , mult(Gf , uf_gamma));
    return cov_derivative;
}

function matrix_NM covariant_projection(int geo ; string attr ; int faceid){
    /*
    Linear mapping that parallel-transports the vertex-based vectors to the tangent space of the face
    cf. the projection operator defined in 4.3 of the paper
    */
    matrix_NM Pf = Pf(geo , faceid);
    matrix_NM Df = Df(len(primpoints(geo,faceid)));
    matrix_NM uf_gamma = parallel_transported_localframe_vectors(geo,attr,faceid);

    matrix_NM cov_proj = mult(Pf , mult(Df , uf_gamma));
    return cov_proj;
}

function float vector_laplacian(int geo ; string attr ; int faceid ; float lambda){
    /*

    */
    matrix_NM cov_derivative = covariant_derivative(geo,attr,faceid);
    matrix_NM cov_proj = covariant_projection(geo,attr,faceid);
    float face_area = face_area(geo,faceid);

    return face_area * frobenius_norm2(cov_derivative) + lambda * frobenius_norm2(cov_proj);
}

function float vector_laplacian(int geo ; string attr ; int faceid){
    // Paper suggests to take lambda=1
    return vector_laplacian(geo,attr,faceid,1.);
}

/*
------
------
PART 6. Additional Operators
------
------

*/

//Shape operator

function matrix_NM Sf(int geo ; int faceid){
    // Broken ?? or not ?? TODO

    matrix_NM Tf = Tf(geo,faceid);
    matrix_NM Nf = Nf(geo,faceid);
    matrix_NM Gf = Gf(geo,faceid);

    matrix_NM result = add( mult(Gf,Nf) , mult(transpose(Nf),transpose(Gf)));
    result = mult(result,Tf);
    result = mult(transpose(Tf),result);
    result = mult(0.5,result);
    return result;
}

/*
------
------
CONVENIENCE FUNCTIONS
------
------
Apply operator directly to the attribute matrix


*/

function vector face_gradient3(int geo  ; string attr ; int id){
    return phif3(geo,id,attr) * Gf3(geo,id);
}

function matrix_NM face_gradient(int geo  ; string attr ; int id){
    /*
    Returns the gradient of a scalar field
    */
    matrix_NM Gf = Gf(geo,id);
    matrix_NM phif = phif(geo,id,attr);
    return mult(Gf,phif);
}

function vector face_gradient(int geo ; string attr ; int id ){
    /*
    Returns the gradient of a scalar field
    */
    //if triangle we can optimize computations by using Houdini's matrix3 rather than the not optmized but general matrix_NM
    int N_pts = len(primpoints(geo,id));
    vector grad_vec = 0;
    if(N_pts==3){
        grad_vec = face_gradient3(geo,attr,id);
    }
    else{
        matrix_NM grad = face_gradient(geo, attr , id);
        grad_vec = matNMtohou(grad);
    }
    return grad_vec;
}

function float face_dirichlet_energy(int geo ; string attr ; int id ; float lambda){
    //compute the dirichlet energy: phif_t * Lf * phif
    //i.e the Laplacian

    matrix_NM phif = phif(geo,id,attr);
    matrix_NM Lf = Lf(geo,id,lambda);
    matrix_NM e_dirichlet = mult(transpose(phif),mult(Lf,phif));
    float e_f = matNMtohou(e_dirichlet);
    return e_f;
}

function float face_dirichlet_energy(int geo ; string attr ; int id){
    // In 7.1 of Discrete Differential Operators on Polygonal Meshes
    // a value of lambda = 1 is required
    return face_dirichlet_energy(geo,attr,id,1.);
}


#endif
