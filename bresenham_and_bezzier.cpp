#include <SFML/Graphics.hpp>
#include <utility>
#include <iostream>
#include "math.cpp"

int put_pixxel(int x,int y,sf::RenderWindow& window, sf::Color color = sf::Color::White){

    sf::RectangleShape rectangle(sf::Vector2f(1.f, 1.f));

    if(x<0 || y<0){
        return -1;
    }

    rectangle.move(x,y);

    rectangle.setFillColor(color);
    
    window.draw(rectangle);

    //std::cout<<"Ponemos pixxel en "<<x<<" "<<y<<"\n";

    return 0;
}



void Bresenham_horizontal(int x0,int y0,int x1,int y1,sf::RenderWindow& window , sf::Color color = sf::Color::White){

    if(x0>x1){

        int aux = x0;
        x0 = x1;
        x1 = aux;

        aux = y0;

        y0 = y1;
        y1 = aux;
        
    }

    int current_y=0;
    int y_change_param = 0;

    int change_in_x = x1-x0;
    int change_in_y = y1-y0;

    int direccion_y;

    if(change_in_y<0)
    direccion_y = -1;
    else 
    direccion_y = 1;

    change_in_y *= direccion_y;

    if(change_in_x != 0){

        current_y = y0;
        y_change_param = 2*change_in_y - change_in_x;

        for(int i = 0;i<change_in_x+1;i++){

            put_pixxel(x0 + i, current_y,window, color);

            if(y_change_param >= 0){

                current_y += direccion_y;
                y_change_param = y_change_param -2*change_in_x;

            }

            y_change_param = y_change_param + 2*change_in_y;

        }

    }



}

void Bresenham_vertical(int x0,int y0,int x1,int y1,sf::RenderWindow& window , sf::Color color = sf::Color::White){
    if(y0>y1){
        
        int aux = x0;
        x0 = x1;
        x1 = aux;

        aux = y0;

        y0 = y1;
        y1 = aux;

    }

    int current_x=0;
    int x_change_param = 0;

    int change_in_x = x1-x0;
    int change_in_y = y1-y0;

    int direccion_x;

    if(change_in_x<0)
    direccion_x = -1;
    else 
    direccion_x = 1;

    change_in_x *= direccion_x;

    if(change_in_y != 0){

        current_x = x0;
        x_change_param = 2*change_in_x - change_in_y;

        for(int i = 0;i<change_in_y+1;i++){

            put_pixxel(current_x, y0 + i,window, color);

            if(x_change_param >= 0){

                current_x += direccion_x;
                x_change_param = x_change_param -2*change_in_y;

            }

            x_change_param = x_change_param + 2*change_in_x;

        }

    }
}

void draw_line(int x0,int y0,int x1,int y1,sf::RenderWindow& window , sf::Color color = sf::Color::White){

    if(abs(x1-x0)>=abs(y1-y0))
    Bresenham_horizontal(x0,y0,x1,y1,window,color);
    else
    Bresenham_vertical(x0,y0,x1,y1,window,color);

}

std::pair<int,int> lerp(int x0,int y0,int x1,int y1, double param){

    std::pair<int,int> result;

    double new_x = (x0*((1-param))) + (param*x1);

    double new_y = (y0*(1-param)) + (param*y1);

    std::cout<<new_x<<" "<<new_y<<"\n";

    result.first = new_x;

    result.second = new_y;

    return result;
    
}

void Iterative_lerp(int x0,int y0,int x1,int y1, sf::RenderWindow& window){

    for(int i = 0; i<=1000;i++){

        std::pair<int,int> point_A = lerp(x0,y0,x1,y1,(double)(i)/(double)(1000));
        //std::cout<<"ponermos pixxel en "<<point_A.first<<" "<<point_A.second<<"\n";
        put_pixxel(point_A.first,point_A.second,window);

    }

}

void bezzier_curve(int x0,int y0,int x1,int y1,int x2,int y2,sf::RenderWindow& window , sf::Color color = sf::Color::White){
    
    for(int i = 0; i<=1000;i++){

        std::pair<int,int> point_A = lerp(x0,y0,x1,y1,(double)(i)/(double)(1000));

        std::cout<<"primer punto "<<point_A.first<< " "<<point_A.second<<"\n";

        std::pair<int,int> point_B = lerp(x1,y1,x2,y2,(double)(i)/(double)(1000));

        std::pair<int,int> point_C = lerp(point_A.first,point_A.second,point_B.first,point_B.second,(double)(i)/(double)(1000));

        //std::cout<<"iteracion "<<i<<"\n";
        
        put_pixxel(point_C.first,point_C.second,window,color);
    }

}
