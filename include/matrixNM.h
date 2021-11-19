#ifndef __matrixNM_h__
#define __matrixNM_h__

/*
Allow the use of NxM matrices in Houdini. Operations are slower than Houdini native matrices so prefer using those when possible
*/

struct matrix_NM{
    float values[];
    int N;
    int M;

    //helpers
    vector2 shape(){
        return set(N,M);
    }

    string string_dim(){
        // mostly for debug messages
        // OUTPUT-> '(N,M)'
        return '('+itoa(N)+','+itoa(M)+')';
    }

    int get_flattened_index(int i,j){
        /*
        Get the flattened index
        */
        return M*i+j;
    }

    int __set(float val ; int i,j){
        /*
        Set the value i,j. Returns 1 if it worked -1 else
        */
        values[this->get_flattened_index(i,j)] = val;
        return ( i<N && j<M ) ? 1 : -1;
    }

    float __get(int i,j){
        if ( i >= N || j >= M)  error("Index is out of bound. Queried ("+itoa(i)+","+itoa(j)+") ; Dimensions are: "+this->string_dim());
        return values[this->get_flattened_index(i,j)];
    }

    void __const_add(float lambda){
        for(int idx = 0 ; idx < len(values) ; idx++){
            values[idx] += lambda;
        }
    }

    void __add(matrix_NM mat){
        /*
        Element-wise addition
        */
        if ( M != mat.M || N != mat.N)  error("Not same size. self is: "+this->string_dim()+". Other is: "+mat->string_dim());

        for(int i = 0 ; i < N ; i++){
            for(int j = 0 ; i < M ; j++){
                values[this->get_flattened_index(i,j)] += mat.values[this->get_flattened_index(i,j)];
            }
        }
        return;
    }

    void __const_mult(float lambda){
        for(int idx = 0 ; idx < len(values) ; idx++){
            values[idx] *= lambda;
        }
    }

    void __mult(matrix_NM mat){
        /*
        Element-wise multiplication
        */
        if ( M != mat.M || N != mat.N)  error("Not same size. self is: "+this->string_dim()+". Other is: "+mat->string_dim());

        for(int i = 0 ; i < N ; i++){
            for(int j = 0 ; i < M ; j++){
                values[this->get_flattened_index(i,j)] *= mat.values[this->get_flattened_index(i,j)];
            }
        }
        return;
    }


    void __matmul(matrix_NM mat){
        /*
        Matrix multiplication
        */
        // TODO use mult(a,b) in the meantime
    }
}

/*
Constructors
*/

function matrix_NM matNMfromcols(vector cols[]){
    /*
    Build a 3 x len(cols) matrix with the i-th column being the i-th vector of cols
    */
    int N = 3;
    int M = len(cols);
    float values[] = {};
    for(int i = 0 ; i < N ; i++){
    for(int j = 0 ; j < M ; j++){
        append(values, getcomp(cols[j],i));
    }}
    return matrix_NM(values,N,M);
}

function matrix_NM matNMfromrows(vector rows[]){
    /*
    Build a len(rows) x 3 matrix with the i-th row being the i-th vector of rows
    */
    int N = len(rows);
    int M = 3;
    float values[] = {};
    for(int i = 0 ; i < N ; i++){
    for(int j = 0 ; j < M ; j++){
        append(values, getcomp(rows[i],j));
    }}
    return matrix_NM(values,N,M);
}

function matrix_NM ones(int N,M){
    /*
    Return a matrix filled with 1
    */
    float values[] = {};
    for(int r = 0 ; r < N ; r++){
    for(int c = 0 ; c < M ; c++){
        append(values,1);
    }}
    return matrix_NM(values,N,M);
}

function matrix_NM ident(int N,M){
    /*
    Return an identity matrix (i.e. 0 everywhere and 1 where i==j)
    */
    float values[] = {};
    for(int r = 0 ; r < N ; r++){
    for(int c = 0 ; c < M ; c++){
        float val = r==c;
        append(values,val);
    }}
    return matrix_NM(values,N,M);
}


/*
Overloading VEX attribute functions
*/

function void setattrib(int geohandle ; string class ; string name ; int elemnum ; int vtxnum ; matrix_NM mat){
    setattrib(geohandle, class , name+'_values',elemnum,vtxnum,mat.values);
    setattrib(geohandle, class , name+'_N',elemnum,vtxnum,mat.N);
    setattrib(geohandle, class , name+'_M',elemnum,vtxnum,mat.M);
}

function void setpointattrib(int geohandle; string name ; int ptnum ; matrix_NM mat){
    setpointattrib(geohandle,name+'_values',ptnum,mat.values);
    setpointattrib(geohandle,name+'_N',ptnum,mat.N);
    setpointattrib(geohandle,name+'_M',ptnum,mat.M);
}

function void setprimattrib(int geohandle; string name ; int primnum ; matrix_NM mat){
    setprimattrib(geohandle,name+'_values',primnum,mat.values);
    setprimattrib(geohandle,name+'_N',primnum,mat.N);
    setprimattrib(geohandle,name+'_M',primnum,mat.M);
}

function void setdetailattrib(int geohandle; string name ; matrix_NM mat){
    setdetailattrib(geohandle,name+'_values',mat.values);
    setdetailattrib(geohandle,name+'_N',mat.N);
    setdetailattrib(geohandle,name+'_M',mat.M);
}

function matrix_NM point(int geohandle ; string name ; int ptnum){
    float values[] = point(geohandle,name+'_values',ptnum);
    int N = point(geohandle,name+'_N',ptnum);
    int M = point(geohandle,name+'_M',ptnum);

    return matrix_NM(values,N,M);
}

function matrix_NM prim(int geohandle ; string name ; int primnum){
    float values[] = prim(geohandle,name+'_values',primnum);
    int N = prim(geohandle,name+'_N',primnum);
    int M = prim(geohandle,name+'_M',primnum);

    return matrix_NM(values,N,M);
}

function matrix_NM detail(int geohandle ; string name){
    float values[] = detail(geohandle,name+'_values');
    int N = detail(geohandle,name+'_N');
    int M = detail(geohandle,name+'_M');

    return matrix_NM(values,N,M);
}

/*
Useful functions to get values:

Set:
 - setcomp( matrix , val , i , j)
 - TODO setcolumn , setrow

Get:
 - getcomp(matrix , i , j)
 - getcolumn(matrix , c)
 - getrow(matrix , r)

*/

function void setcomp(matrix_NM m ; float val ; int i,j){
    /*
    Set m[i][j]
    */
    m.values[m.M*i+j] = val;
}

function float getcomp(matrix_NM m ; int i,j){
    /*
    Return m[i][j]
    */
    return m.values[m.M * i + j];
}



function float[] getcolumn(matrix_NM m ; int c){
    /*
    Return the c-th column of m (i.e. m[:][c])
    as an Houdini array
    */
    float result[] = {};
    for(int i = 0 ; i<m.N ; i++){
        append(result,getcomp(m,i,c));
    }
    return result;
}

function matrix_NM getcolumn(matrix_NM m ; int c){
    /*
    Return the c-th column of m (i.e. m[:][c])
    as a matrix_NM
    */
    float values[] = getcolumn(m,c);
    return matrix_NM(values,m.N,1);
}

function vector getcolumnvector(matrix_NM m ; int c){
    /*
    Return the c-th column of m (i.e. m[:][c])
    as a vector (N must be equal to 3)
    */
    if(m.N!=3)  error("Column not the right size. Shape is: "+m->string_dim());
    vector result = 0;
    for(int i = 0 ; i<m.N ; i++){
        setcomp(result,getcomp(m,i,c),i);
    }
    return result;
}

function float[] getrow(matrix_NM m ; int r){
    /*
    Return the r-th row of m (i.e. m[r][:])
    as a Houdini array
    */
    float result[] = {};
    for(int i = 0 ; i<m.M ; i++){
        append(result,getcomp(m,r,i));
    }
    return result;
}

function matrix_NM getrow(matrix_NM m ; int r){
    /*
    Return the r-th row of m (i.e. m[r][:])
    as a matrix_NM
    */
    float values[] = getrow(m,r);
    return matrix_NM(values,1,m.M);
}

function vector getrowvector(matrix_NM m ; int r){
    /*
    Return the r-th row of m (i.e. m[r][:])
    as a vector (M must be equal to 3)
    */

    if(m.M!=3)  error("Row not the right size. Shape is: "+m->string_dim());
    vector result = 0;
    for(int i = 0 ; i<m.M ; i++){
        setcomp(result,getcomp(m,r,i),i);
    }
    return result;
}

function vector[] getcolvectors(matrix_NM m){
    /*
    Return an array of vector, each vector being a column of m
    */
    if(m.N != 3)    error("More than 3 rows. Columns are not vectors: dim="+itoa(m.N));
    vector result[] = {};
    for(int col = 0 ; col < m.M ; col++){
        vector val = 0;
        for(int row = 0 ; row < m.N ; row++){
            setcomp(val,getcomp(m,row,col),row);
        }
        append(result,val);
    }
    return result;
}

function vector[] getrowvectors(matrix_NM m){
    /*
    Return an array of vector, each vector being a row of m
    */
    if(m.M != 3)    error("More than 3 rows. Columns are not vectors: dim="+itoa(m.M));
    vector result[] = {};
    for(int row = 0 ; row < m.N ; row++){
        vector val = 0;
        for(int col = 0 ; col < m.M ; col++){
            setcomp(val,getcomp(m,row,col),col);
        }
        append(result,val);
    }
    return result;
}

/*
Conversion Utils
*/

/*
tomatrixNM will convert Houdini native matrices to the corresponding matrix_NM
*/

function matrix_NM tomatrixNM(float m[]){
    return matrix_NM(m,1,len(m));
}

function matrix_NM tomatrixNM(vector v){
    float values[] = {};
    append(values,v.x);
    append(values,v.y);
    append(values,v.z);
    return matrix_NM(values,1,3);
}

function matrix_NM tomatrixNM(vector2 v){
    float values[] = {};
    append(values,v.x);
    append(values,v.y);
    return matrix_NM(values,1,2);
}

function matrix_NM tomatrixNM(matrix3 m){
    float values[] = {};
    matrix3 tm = m;
    for(int i = 0 ; i < 3 ; i++){
    for(int j = 0 ; j < 3 ; j++){
        float val = getcomp(tm,i,j);
        append(values,val);
    }}
    return matrix_NM(values,3,3);
}

function matrix_NM tomatrixNM(matrix m){
    float values[] = {};
    matrix tm = m;
    for(int i = 0 ; i < 4 ; i++){
    for(int j = 0 ; j < 4 ; j++){
        float val = getcomp(tm,i,j);
        append(values,val);
    }}
    return matrix_NM(values,4,4);
}

/*
matNMtohou will convert a native Houdini struct to a matrix_NM
*/

function matrix2 matNMtohou(matrix_NM mat){
    int N = mat.N;
    int M = mat.M;
    matrix2 mat2 = 0;
    if( N!=2 || M!=2)   error("mat is not a matrix2. Dimensions are: "+mat->string_dim());
    for(int i = 0 ; i < N ; i++){
    for(int j = 0 ; j < M ; j++){
        setcomp(mat2,getcomp(mat,i,j),i,j);
    }}
    return mat2;
}

function matrix3 matNMtohou(matrix_NM mat){
    int N = mat.N;
    int M = mat.M;
    matrix3 mat3 = 0;
    if( N!=3 || M!=3)   error("mat is not a matrix3. Dimensions are: "+mat->string_dim());
    for(int i = 0 ; i < N ; i++){
    for(int j = 0 ; j < M ; j++){
        setcomp(mat3,getcomp(mat,i,j),i,j);
    }}
    return mat3;
}

function matrix matNMtohou(matrix_NM mat){
    int N = mat.N;
    int M = mat.M;
    matrix mat4 = 0;
    if( N!=4 || M!=4)   error("mat is not a matrix4. Dimensions are: "+mat->string_dim());
    for(int i = 0 ; i < N ; i++){
    for(int j = 0 ; j < M ; j++){
        setcomp(mat4,getcomp(mat,i,j),i,j);
    }}
    return mat4;
}

function vector matNMtohou(matrix_NM mat){
    int N = mat.N;
    int M = mat.M;
    float values[] = mat.values;
    if( N*M != 3)   error("mat does not have vector dimensions. Dimensions are: "+mat->string_dim());
    return set(values[0],values[1],values[2]);
}

function vector2 matNMtohou(matrix_NM mat){
    float values[] = mat.values;
    if( len(values) != 2)   error("mat does not have vector dimensions. Dimensions are: "+mat->string_dim());
    return set(values[0],values[1]);
}

function float matNMtohou(matrix_NM mat){
    int N = mat.N;
    int M = mat.M;
    if( N != 1 || M != 1)   error("Not a single value. Could not convert to float. Dimensions are: "+mat->string_dim());
    return mat.values[0];
}

// Operation on matrices

function matrix_NM transpose(matrix_NM m){
    /*
    Return the transpose of m
    */
    float values[] = {};
    for(int i = 0 ; i < m.M ; i++){
        append(values,getcolumn(m,i));
    }
    return matrix_NM(values,m.M,m.N);
}

function float trace(matrix_NM m){
    /*
    Return the trace of m
    */
    float trace = 0 ;
    int dim = min(m.M,m.N);
    for(int i = 0 ; i < dim ; i++){
        trace += getcomp(m,i,i);
    }
    return trace;
}

function matrix_NM add(matrix_NM a,b){
    /*
    Return the sum of a and b
    */
    if(a.N!=b.N || a.M!=b.M){
        error('Not same shape: ('+itoa(a.N)+','+itoa(a.M)+') and ('+itoa(b.N)+','+itoa(b.M)+')');
    }
    matrix_NM result = a;
    for(int i = 0 ; i < a.N*a.M ; i++){
        result.values[i] += b.values[i];
    }
    return result;
}

function float dot(matrix_NM a,b){
    /*
    For a and b two vectors, return the euclidean dot product of a and b
    */
    if(( a.N!=1 && a.M!=1) || (b.N!=1 && b.M!=1)){
        error('not two vectors');
    }
    if(max(a.N,a.M)!=max(b.N,b.M)){
        error('not two same size vector');
    }
    float result = 0;
    for(int i=0 ; i < max(a.N,a.M) ; i++){
        result += a.values[i] * b.values[i];
    }
    return result;
}

function matrix_NM mult(matrix_NM a,b){
    /*
    Return the matmul between a and b
    */
    if(a.M!=b.N){
        error('Shape not compatible for mult: '+a->string_dim()+' and '+b->string_dim());
    }
    float values[]={};

    for(int r = 0 ;  r<a.N ; r++){
        //we fill values row by row
        float curr_row[] = {};
        matrix_NM arow = getrow(a,r);
        for(int c = 0; c<b.M ; c++){
            matrix_NM bcol = getcolumn(b,c);
            append(values,dot(arow,bcol));
        }
    }
    return matrix_NM(values,a.N,b.M);
}

function matrix_NM mult(matrix_NM a ; vector b){
    matrix_NM conv = tomatrixNM(b);
    return mult(a,conv);
}

function matrix_NM mult(matrix_NM a ; vector2 b){
    matrix_NM conv = tomatrixNM(b);
    return mult(a,conv);
}

function matrix_NM mult(matrix_NM a ; matrix b){
    matrix_NM conv = tomatrixNM(b);
    return mult(a,conv);
}

function matrix_NM mult(matrix_NM a ; matrix3 b){
    matrix_NM conv = tomatrixNM(b);
    return mult(a,conv);
}

//

function matrix_NM mult(vector b ;matrix_NM a){
    matrix_NM conv = tomatrixNM(b);
    return mult(conv,a);
}

function matrix_NM mult(vector2 b ; matrix_NM a){
    matrix_NM conv = tomatrixNM(b);
    return mult(conv,a);
}

function matrix_NM mult(matrix b ; matrix_NM a){
    matrix_NM conv = tomatrixNM(b);
    return mult(conv,a);
}

function matrix_NM mult(matrix3 b ; matrix_NM a){
    matrix_NM conv = tomatrixNM(b);
    return mult(conv,a);
}


function matrix_NM mult(matrix_NM a ; float f){
    matrix_NM result = a;
    for(int i = 0 ; i < len(a.values) ; i++){
        result.values[i] *=f;
    }
    return result;
}

function matrix_NM mult(float f ;matrix_NM a){
    matrix_NM result = a;
    for(int i = 0 ; i < len(a.values) ; i++){
        result.values[i] *=f;
    }
    return result;
}


/*
----
----
NORMS
----
----
*/

function float L_pq_norm(float p,q ;matrix_NM mat){
    /*
    Returns the L_pq norm of a matrix
    -> for p=q=2 use rather frobenius_norm()
    */
    int N = mat.N;
    int M = mat.M;

    float sum = 0;
    for(int j = 0 ; i < N ; i++){
        float temp_sum = 0;
        for(int i = 0 ; j < M ; j++){
            float value = mat->__get(i,j);
            temp_sum += pow( abs(value) , p );
        }
        sum+= pow(temp_sum , q / p );
    }
    return pow(sum, 1./q);
}

function float frobenius_norm2(matrix_NM mat){
    /*
    Returns the Frobenius norm pow2 of a matrix
    Also known as: trace( transpose(mat) * mat)
    */
    int N = mat.N;
    int M = mat.M;

    float sum = 0;
    for(int i = 0 ; i < N ; i++){
        for(int j = 0 ; j < M ; j++){
            float value = mat->__get(i,j);
            sum += value * value;
        }
    }
    return sum;
}

function float frobenius_norm(matrix_NM mat){
    /*
    Returns the Frobenius norm of a matrix
    Also known as: sqrt(trace( transpose(mat) * mat))
    */
    int N = mat.N;
    int M = mat.M;

    float sum = 0;
    for(int i = 0 ; i < N ; i++){
        for(int j = 0 ; j < M ; j++){
            float value = mat->__get(i,j);
            sum += value * value;
        }
    }
    return sqrt(sum);
}

#endif
