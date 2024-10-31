#include <SFML/Graphics.hpp>
#include <iostream>
#include "bresenham_and_bezzier.cpp"
#include <fstream>

/**
 * @brief genera un par con dos vectores de coordenas, first:representa a los vertices second:representa a los indices;
 * @param name_document string que contiene el nombre del obj a renderizar 
 * @return par con dos vectores de Vector4
 */
pair<vector<Vector4>,vector<Vector4>> obj_vectores_y_triangulos(string name_document){

    std::ifstream obj(name_document);
    string renglon;

    if(!obj.is_open())
    std::cout<<"error al abrir el archivo"<<"\n";
    //variable para guardar indices
    int aux_index = 0;
    //coordenadas de los vectores
    double x,y,z,w; 
    //vector a agregar
    Vector4 agregado;
    //vector de puntos
    vector<Vector4> points;
    //vector de indices
    vector<Vector4> indexes;
    //string auxiliar
    string txt_storage;

    while(!obj.eof()){
        

        std::getline(obj,renglon);

        if(renglon == "")
        continue;

        if(renglon[0] == 'v'){
            if(renglon[1] == 'n')
            continue;

            renglon.erase(0,2);
            //x
            aux_index = renglon.find(" ");

            x = stod(renglon.substr(0,aux_index));

            renglon.erase(0,aux_index+1);

            //y
            aux_index = renglon.find(" ");

            y = stod(renglon.substr(0,aux_index));

            renglon.erase(0,aux_index+1);

            //z
            aux_index = renglon.find(" ");

            z = stod(renglon.substr(0,aux_index));

            agregado.set(x,y,z,1);

            points.push_back(agregado);

        }

        if(renglon[0] == 'f'){
            renglon.erase(0,2);
            w = -1;

            //x
            aux_index = renglon.find(" ");

            txt_storage = renglon.substr(0,aux_index+1);

            //cout<<txt_storage.substr(0,txt_storage.find("/"))<<"\n";

            x = stod(txt_storage.substr(0,txt_storage.find("/")));

            renglon.erase(0,aux_index+1);

            //y
            aux_index = renglon.find(" ");

            txt_storage = renglon.substr(0,aux_index+1);

            y = stod(txt_storage.substr(0,txt_storage.find("/")));

            renglon.erase(0,aux_index+1);

            //z
            

            if(renglon.find(" ") == -1){
                //cout<<renglon.substr(0,renglon.find("/"))<<"\n";

                z = stod(renglon.substr(0,renglon.find("/")));
                agregado.set(x,y,z,w);

                indexes.push_back(agregado);
                continue;
            }else{

                aux_index = renglon.find(" ");

                txt_storage = renglon.substr(0,aux_index+1);

                z = stod(txt_storage.substr(0,txt_storage.find("/")-1));

                renglon.erase(0,aux_index+1);

            }

            cout<<"mal"<<"\n";

            //w

            if(renglon.find("/") == -1){

                w = stod(renglon);

            }else{

                w = stod(renglon.substr(0,renglon.find("/")-1));

            }           

            agregado.set(x,y,z,w);

            indexes.push_back(agregado);

        }

        //std::cout<<renglon<<"\n";
    }

    pair<vector<Vector4>,vector<Vector4>> result;

    result.first = points;
    result.second = indexes;


    return result;
}

Matrix4 Create_final_matrix(Matrix4 viewport,Matrix4 perspective, Matrix4 lookAt){
    Matrix4 pivot = Matrix4::multiply(perspective,lookAt);
    Matrix4 final_matrix = Matrix4::multiply(viewport,pivot);

    return final_matrix;
}

vector<Vector4> aply_matrix(vector<Vector4> points,Matrix4 matrix)
{   
    Vector4 added = Vector4();
    vector<Vector4> vectors_multipied;

    for(Vector4 point: points){

        added = matrix.multiplyVector(point);
        vectors_multipied.push_back(added);

    }

    return vectors_multipied;
}

Vector3 get_center_of_model(vector<Vector4> points){

    Vector3 center;

    double center_x = 0;

    double center_y = 0;

    double center_z = 0;


    for(Vector4 point:points){

        center_x += point.x;
        center_y += point.y;
        center_z += point.z;

    }

    center_x = center_x/points.size();
    center_y = center_y/points.size();
    center_z = center_z/points.size();

    center.set(center_x,center_y,center_z);

    return center;
}


vector<double> get_model_cage(vector<Vector4> points){
    double min_x;
    double min_y;
    double min_z;
    double max_x;
    double max_y;
    double max_z;

    min_x = points[0].x;
    min_y = points[0].y;
    min_z = points[0].z;
    max_x = points[0].x;
    max_y = points[0].y;
    max_z = points[0].z;

    for(Vector4 point : points){

        if(point.x<min_x)
        min_x = point.x;

        if(point.y<min_y)
        min_y = point.y;

        if(point.z<min_z)
        min_z = point.z;

        if(max_x<point.x)
        max_x = point.x;

        if(max_y<point.y)
        max_y = point.y;

        if(max_z<point.z)
        max_z = point.z;

    }

    vector<double> cage;
    cage.push_back(min_x);
    cage.push_back(max_x);

    cage.push_back(min_y);
    cage.push_back(max_y);

    cage.push_back(min_z);
    cage.push_back(max_z);

    return cage;
}

void paint_vector(Vector4 vec, sf::RenderWindow &window){
    put_pixxel(vec.x/vec.w,vec.y/vec.w,window);
}

void bresenham_vectors(Vector4 vec1,Vector4 vec2,sf::RenderWindow &window){
    draw_line(vec1.x/vec1.w,vec1.y/vec1.w,vec2.x/vec2.w,vec2.y/vec2.w,window);
}

void draw_model_triangles(pair<vector<Vector4>,vector<Vector4>> model,sf::RenderWindow &window){
    int index1,index2,index3,index4;

    for(Vector4 indexes:model.second){
        index1 = (int) indexes.x;
        index2 = (int) indexes.y;
        index3 = (int) indexes.z;
        index4 = (int) indexes.w;

        bresenham_vectors(model.first[index1-1],model.first[index2-1],window);
        bresenham_vectors(model.first[index2-1],model.first[index3-1],window);
        if(indexes.w = -1){
            bresenham_vectors(model.first[index3-1],model.first[index1-1],window);
        }else{
            bresenham_vectors(model.first[index3-1],model.first[index4-1],window);
            bresenham_vectors(model.first[index4-1],model.first[index1-1],window);
        }
    }
}




int main()
{   

    int height = 720;
    int wide = 720;

    string name_document = "Happy_Buddha.obj";

    pair<vector<Vector4>,vector<Vector4>> model = obj_vectores_y_triangulos(name_document);

    vector<double> model_cage = get_model_cage(model.first);

    Vector3 center = Vector3((model_cage[0]+model_cage[1])/2.0,(model_cage[2]+model_cage[3])/2.0,(model_cage[4]+model_cage[5])/2.0);

    Vector3 camera = Vector3(model_cage[1],model_cage[3],model_cage[5]);

    Vector3 min = Vector3(model_cage[0]*1.5,model_cage[2]*1.5,model_cage[4]*1.5);

    double model_distance = Vector3::distance(Vector3::subtract(camera,min),Vector3()) * 1.5;

    sf::RenderWindow window(sf::VideoMode(height,wide),"obj_render",sf::Style::Default);
    window.setFramerateLimit(30);

    double theta = 0;

    while(window.isOpen())
    { 

        sf::Event event;
        while (window.pollEvent(event)){
            if(event.type == sf::Event::EventType::MouseButtonPressed){
                camera = Vector3(camera.x,camera.y+1.0,camera.z);
            }

            if (event.type == sf::Event::Closed)
                window.close();
        }

        Matrix4 look_at = Matrix4::lookAt(camera,center,Vector3(0.0,0.0,1.0));

        Matrix4 perspective = Matrix4::perspective(90,1,0.1,100);

        Matrix4 view_port = Matrix4::viewPort(720,720);

        Matrix4 final_Matrix = Create_final_matrix(view_port,perspective,look_at);

        window.clear();

        Matrix4 RotationMatrix = Matrix4::rotateZ(theta);

        vector<Vector4> transformed_points = aply_matrix(model.first,RotationMatrix);

        vector<Vector4> transformed_points2 = aply_matrix(transformed_points,final_Matrix);

        pair<vector<Vector4>,vector<Vector4>> model_transformed;

        model_transformed.first = transformed_points2;

        model_transformed.second = model.second;

        //draw_model_triangles(model_transformed,window);

        

        for(Vector4 point : transformed_points2){
            paint_vector(point,window);
        }

        theta = theta + 0.1;

        
        window.display();
        
    }
    return 0;
}