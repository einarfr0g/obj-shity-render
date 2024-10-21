#include <SFML/Graphics.hpp>
#include <iostream>
#include "bresenham_and_bezzier.cpp"
#include <fstream>

pair<vector<Vector4>,vector<Vector4>> obj_vector_and_triangles(string name_document){

    std::ifstream obj(name_document);
    string renglon;

    if(!obj.is_open())
    cout<<"error al abrir el archivo"<<"\n";

    while(!obj.eof()){
        std::getline(obj,renglon);
        if(renglon[0] == 'v'){
            
        }

        //cout<<renglon<<"\n";
    }

    pair<vector<int>,vector<int>> result;


    return result;
}



int main()
{   

    int height = 720;
    int wide = 720;

    string name_document = "Cube_Triangles.obj";



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

        draw_line(puntos_camara[0].x/puntos_camara[0].w,puntos_camara[0].y/puntos_camara[0].w,puntos_camara[1].x/puntos_camara[1].w,puntos_camara[1].y/puntos_camara[1].w,window,sf::Color::Blue);
        draw_line(puntos_camara[1].x/puntos_camara[1].w,puntos_camara[1].y/puntos_camara[1].w,puntos_camara[5].x/puntos_camara[5].w,puntos_camara[5].y/puntos_camara[5].w,window,sf::Color::Blue);
        draw_line(puntos_camara[5].x/puntos_camara[5].w,puntos_camara[5].y/puntos_camara[5].w,puntos_camara[4].x/puntos_camara[4].w,puntos_camara[4].y/puntos_camara[4].w,window,sf::Color::Blue);
        draw_line(puntos_camara[4].x/puntos_camara[4].w,puntos_camara[4].y/puntos_camara[4].w,puntos_camara[0].x/puntos_camara[0].w,puntos_camara[0].y/puntos_camara[0].w,window,sf::Color::Blue);
        draw_line(puntos_camara[7].x/puntos_camara[7].w,puntos_camara[7].y/puntos_camara[7].w,puntos_camara[2].x/puntos_camara[2].w,puntos_camara[2].y/puntos_camara[2].w,window,sf::Color::Blue);
        draw_line(puntos_camara[2].x/puntos_camara[2].w,puntos_camara[2].y/puntos_camara[2].w,puntos_camara[3].x/puntos_camara[3].w,puntos_camara[3].y/puntos_camara[3].w,window,sf::Color::Blue);
        draw_line(puntos_camara[3].x/puntos_camara[3].w,puntos_camara[3].y/puntos_camara[3].w,puntos_camara[6].x/puntos_camara[6].w,puntos_camara[6].y/puntos_camara[6].w,window,sf::Color::Blue);
        draw_line(puntos_camara[6].x/puntos_camara[6].w,puntos_camara[6].y/puntos_camara[6].w,puntos_camara[7].x/puntos_camara[7].w,puntos_camara[7].y/puntos_camara[7].w,window,sf::Color::Blue);
        draw_line(puntos_camara[7].x/puntos_camara[7].w,puntos_camara[7].y/puntos_camara[7].w,puntos_camara[0].x/puntos_camara[0].w,puntos_camara[0].y/puntos_camara[0].w,window,sf::Color::Blue);
        draw_line(puntos_camara[2].x/puntos_camara[2].w,puntos_camara[2].y/puntos_camara[2].w,puntos_camara[1].x/puntos_camara[1].w,puntos_camara[1].y/puntos_camara[1].w,window,sf::Color::Blue);
        draw_line(puntos_camara[3].x/puntos_camara[3].w,puntos_camara[3].y/puntos_camara[3].w,puntos_camara[5].x/puntos_camara[5].w,puntos_camara[5].y/puntos_camara[5].w,window,sf::Color::Blue);
        draw_line(puntos_camara[6].x/puntos_camara[6].w,puntos_camara[6].y/puntos_camara[6].w,puntos_camara[4].x/puntos_camara[4].w,puntos_camara[4].y/puntos_camara[4].w,window,sf::Color::Blue);

        theta = theta + 0.1;

        
        window.display();
        
    }
    return 0;
}