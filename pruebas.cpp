#include "math.cpp"

int main(){

    Matrix4 mat1 = Matrix4(-1,2,3,-4,5,-6,20,8,9.8,10,-11,12,13,14,-150.39,16);

    cout<<" entrada a00 = "<< to_string(mat1.a00)<< " entrada a01 = "<< to_string(mat1.a01)<< " entrada a02 = "<< to_string(mat1.a02)<<" entrada a02 = "<<to_string(mat1.a03)<<"\n"
                 << " entrada a10 = "<< to_string(mat1.a10)<< " entrada a11 = "<< to_string(mat1.a11)<< " entrada a12 = "<< to_string(mat1.a12)<<" entrada a12 = "<<to_string(mat1.a13)<<"\n"
                 << " entrada a20 = "<< to_string(mat1.a20)<< " entrada a21 = "<< to_string(mat1.a21)<< " entrada a22 = "<< to_string(mat1.a22)<<" entrada a22 = "<<to_string(mat1.a23)<<"\n"
                 << " entrada a20 = "<< to_string(mat1.a30)<< " entrada a21 = "<< to_string(mat1.a31)<< " entrada a22 = "<< to_string(mat1.a32)<<" entrada a22 = "<<to_string(mat1.a33)<<"\n";
    //0.0000001,0.0000001,0.0000001,0.0000001,
    Matrix4 mat2 = mat1.multiplyByScalar(-2.53);

    cout<<"\n"<<"\n"<<"\n";

    Matrix4 mat3 = Matrix4::multiply(mat1,mat2);

    mat3 = Matrix4::subtract(mat1,mat2);

    mat3 = mat3.transpose();

    Vector4 vec1 = Vector4( 14,-59,-36,53.47);

    cout<<" entrada a00 = "<< to_string(mat3.a00)<< " entrada a01 = "<< to_string(mat3.a01)<< " entrada a02 = "<< to_string(mat3.a02)<<" entrada a03 = "<<to_string(mat3.a03)<<"\n"
        <<" entrada a10 = "<< to_string(mat3.a10)<< " entrada a11 = "<< to_string(mat3.a11)<< " entrada a12 = "<< to_string(mat3.a12)<<" entrada a13 = "<<to_string(mat3.a13)<<"\n"
        <<" entrada a20 = "<< to_string(mat3.a20)<< " entrada a21 = "<< to_string(mat3.a21)<< " entrada a22 = "<< to_string(mat3.a22)<<" entrada a23 = "<<to_string(mat3.a23)<<"\n"
        <<" entrada a30 = "<< to_string(mat3.a30)<< " entrada a31 = "<< to_string(mat3.a31)<< " entrada a32 = "<< to_string(mat3.a32)<<" entrada a33 = "<<to_string(mat3.a33)<<"\n";
    cout<<"\n"<<"\n"<<"\n";

    Vector4 vec2 = mat1.multiplyVector(vec1);

    cout<<"la coordenada en x es "<< to_string(vec2.x)<<"la coordenada en y es "<< to_string(vec2.y)<<"la coordenada en z es "<< to_string(vec2.z)<<"la coordenada en w es "<< to_string(vec2.w)<< "\n";






}