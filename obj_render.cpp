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

            x = stod(txt_storage.substr(0,txt_storage.find("/")-1));

            renglon.erase(0,aux_index+1);

            //y
            aux_index = renglon.find(" ");

            txt_storage = renglon.substr(0,aux_index+1);

            y = stod(txt_storage.substr(0,txt_storage.find("/")-1));

            renglon.erase(0,aux_index+1);

            //z
            

            if(renglon.find(" ") == -1){
                z = stod(renglon);
                agregado.set(x,y,z,w);

                indexes.push_back(agregado);
                continue;
            }else{

                aux_index = renglon.find(" ");

                txt_storage = renglon.substr(0,aux_index+1);

                z = stod(txt_storage.substr(0,txt_storage.find("/")-1));

                renglon.erase(0,aux_index+1);

            }

            

            //w

            if(renglon.find("/") == -1){

                w = stod(renglon);

            }else{

                w = stod(renglon.substr(0,renglon.find("/")-1));

            }           

            agregado.set(x,y,z,w);

            indexes.push_back(agregado);

        }

        std::cout<<renglon<<"\n";
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

vector<Vector4> Aplicar_Matriz(vector<Vector4> points,Matrix4 matrix)
{   
    Vector4 added = Vector4();
    vector<Vector4> vectors_multipied;

    for(Vector4 point: points){

        added = matrix.multiplyVector(point);
        vectors_multipied.push_back(added);

    }

    return vectors_multipied;
}





int main()
{   

    int height = 720;
    int wide = 720;

    string name_document = "Cube_Triangles.obj";

    pair<vector<Vector4>,vector<Vector4>> modelo = obj_vectores_y_triangulos(name_document);



    sf::RenderWindow window(sf::VideoMode(height,wide),"obj_render",sf::Style::Default);
    window.setFramerateLimit(30);

    double theta = 0; 
    bool click = false;
    Vector3 camera = Vector3(10,10,5);

    while(window.isOpen())
    { 

        sf::Event event;
        while (window.pollEvent(event)){
            if(event.type == sf::Event::EventType::MouseButtonPressed){
                camera = Vector3(camera.x-1,camera.y-1,camera.y-0.5);
            }

            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();

        std::vector <Vector4> puntos;

        puntos.push_back(Vector4(1,1,1,1));
        puntos.push_back(Vector4(1,1,-1,1));
        puntos.push_back(Vector4(1,-1,-1,1));
        puntos.push_back(Vector4(-1,-1,-1,1));
        puntos.push_back(Vector4(-1,1,1,1));
        puntos.push_back(Vector4(-1,1,-1,1));
        puntos.push_back(Vector4(-1,-1,1,1));
        puntos.push_back(Vector4(1,-1,1,1));

        Matrix4 RotationMatrix = Matrix4::rotateZ(theta);

        for(int i=0; i<8;i++){
            puntos[i] = RotationMatrix.multiplyVector(puntos[i]);
        }
        //Matrix4::orthographic(-2,2,-3,3,2,8)
        //Matrix4::Matrix4::lookAt(camera,Vector3(),Vector3(0,0,1))

        Matrix4 Matriz_perspective_look = Matrix4::multiply(Matrix4::perspective(90,1,2,25),Matrix4::lookAt(camera,Vector3(),Vector3(0,0,1)));

        Matrix4 final_tranformation = Matrix4::multiply(Matrix4::viewPort(720,720),Matriz_perspective_look);

        window.clear();

        Vector4 aux = Vector4();
        std::vector <Vector4> puntos_camara;

        for(Vector4 e : puntos){
            aux = final_tranformation.multiplyVector(e);
            puntos_camara.push_back(aux);
            //put_pixxel(aux.x,aux.y,window);
        }

        draw_line(puntos_camara[0].x/puntos_camara[0].w,puntos_camara[0].y/puntos_camara[0].w,puntos_camara[1].x/puntos_camara[1].w,puntos_camara[1].y/puntos_camara[1].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[1].x/puntos_camara[1].w,puntos_camara[1].y/puntos_camara[1].w,puntos_camara[5].x/puntos_camara[5].w,puntos_camara[5].y/puntos_camara[5].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[5].x/puntos_camara[5].w,puntos_camara[5].y/puntos_camara[5].w,puntos_camara[4].x/puntos_camara[4].w,puntos_camara[4].y/puntos_camara[4].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[4].x/puntos_camara[4].w,puntos_camara[4].y/puntos_camara[4].w,puntos_camara[0].x/puntos_camara[0].w,puntos_camara[0].y/puntos_camara[0].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[7].x/puntos_camara[7].w,puntos_camara[7].y/puntos_camara[7].w,puntos_camara[2].x/puntos_camara[2].w,puntos_camara[2].y/puntos_camara[2].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[2].x/puntos_camara[2].w,puntos_camara[2].y/puntos_camara[2].w,puntos_camara[3].x/puntos_camara[3].w,puntos_camara[3].y/puntos_camara[3].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[3].x/puntos_camara[3].w,puntos_camara[3].y/puntos_camara[3].w,puntos_camara[6].x/puntos_camara[6].w,puntos_camara[6].y/puntos_camara[6].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[6].x/puntos_camara[6].w,puntos_camara[6].y/puntos_camara[6].w,puntos_camara[7].x/puntos_camara[7].w,puntos_camara[7].y/puntos_camara[7].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[7].x/puntos_camara[7].w,puntos_camara[7].y/puntos_camara[7].w,puntos_camara[0].x/puntos_camara[0].w,puntos_camara[0].y/puntos_camara[0].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[2].x/puntos_camara[2].w,puntos_camara[2].y/puntos_camara[2].w,puntos_camara[1].x/puntos_camara[1].w,puntos_camara[1].y/puntos_camara[1].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[3].x/puntos_camara[3].w,puntos_camara[3].y/puntos_camara[3].w,puntos_camara[5].x/puntos_camara[5].w,puntos_camara[5].y/puntos_camara[5].w,window,sf::Color::Yellow);
        draw_line(puntos_camara[6].x/puntos_camara[6].w,puntos_camara[6].y/puntos_camara[6].w,puntos_camara[4].x/puntos_camara[4].w,puntos_camara[4].y/puntos_camara[4].w,window,sf::Color::Yellow);

        theta = theta + 0.1;

        
        window.display();
        
    }
    return 0;
}