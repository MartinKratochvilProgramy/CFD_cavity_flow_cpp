#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <math.h>
#define PI 3.14159265359

//idea taken from https://nbviewer.org/github/barbagroup/CFDPython/blob/master/lessons/14_Step_11.ipynb
//TODO: create draw function?

const int WIDTH = 820;
const int HEIGHT = 820;
const float dx = 0.02;
const float dy = 0.02;
const int nx = 41;
const int ny = 41;
const float dt = 0.00001;
const float rho = 1.;
const float nu = 0.1;
const int n_it = 50;
const int ARROW_SCALE = 100;


std::vector<std::vector<float>> matrix_zeros(int matrix_width, int matrix_height) {
    //returns matrix of 0s with given width/height

    std::vector<std::vector<float>> zeros;

    for (int y = 0; y < matrix_height; y++)
    {
        std::vector<float> row;
        for (int x = 0; x < matrix_width; x++)
        {    
            row.push_back(0.);
        }
        zeros.push_back(row);
    }
    return zeros;
}

std::vector<std::vector<float>> matrix_ones(int matrix_width, int matrix_height) {
    //returns matrix of 1s with given width/height

    std::vector<std::vector<float>> ones;

    for (int y = 0; y < matrix_height; y++)
    {
        std::vector<float> row;
        for (int x = 0; x < matrix_width; x++)
        {
            row.push_back(1.);
        }
        ones.push_back(row);
    }
    return ones;
}

void build_up_b(std::vector<std::vector<float>>& b, std::vector<std::vector<float>>& u, std::vector<std::vector<float>>& v) {
    for (int y = 1; y < ny-1; y++)
    {
        for (int x = 1; x < nx-1; x++)
        {
            b[y][x] = (rho * (1. / dt *
                ((u[y][x + 1] - u[y][x - 1]) / (2 * dx) 
                    + (v[y + 1][x] - v[y - 1][x]) / (2 * dy)) -
                pow((u[y][x + 1] - u[y][x - 1]) / (2 * dx), 2) -
                2 * ((u[y + 1][x] - u[y - 1][x]) / (2 * dy) *
                    (v[y][x + 1] - v[y][x - 1]) / (2 * dx)) -
                pow((v[y + 1][x] - v[y - 1][x]) / (2 * dy), 2)));
        }
    }
}

void pressure_poisson(std::vector<std::vector<float>>& p, std::vector<std::vector<float>>& pn, std::vector<std::vector<float>>& b) {

    pn = p;
    
    for (int i = 0; i < n_it; i++)
    {
        pn = p;
        for (int y = 1; y < ny - 1; y++)
        {
            for (int x = 1; x < nx - 1; x++)
            {
                p[y][x] = (((pn[y][x + 1] + pn[y][x - 1]) * pow(dy, 2) +
                            (pn[y + 1][x] + pn[y - 1][x]) * pow(dx, 2))/
                            (2 * (pow(dx, 2) + pow(dy, 2))) -
                            rho * pow(dx, 2) * pow(dy, 2) / (2 * (pow(dx, 2) + pow(dy, 2))) *
                            b[y][x]);
                        
            }
        }
        // dp/dx = 0 @right
        for (int y = 0; y < ny; y++)
        {
            p[y][nx - 1] = p[y][nx - 2];
        }
        // dp/dx = 0 @left
        for (int y = 0; y < ny; y++)
        {
            p[y][0] = p[y][1];
        }
        // dp/dy = 0 @bottom
        for (int x = 0; x < nx; x++)
        {
            p[ny - 1][x] = 0; //p[ny - 2][x];
        }
        // dp/dy = 0 @top
        for (int x = 0; x < nx; x++)
        {
            p[0][x] = p[1][x];
        }
    }   
}

void cavity_flow(std::vector<std::vector<float>>& u, std::vector<std::vector<float>>& un, std::vector<std::vector<float>>& v, std::vector<std::vector<float>>& vn, std::vector<std::vector<float>>& p, std::vector<std::vector<float>>& pn, std::vector<std::vector<float>>& b) {
    //this function handles the main calculation

    un = u;
    vn = v;
    b = matrix_zeros(nx, ny);

    build_up_b(b, u, v);
    pressure_poisson(p, pn, b);

    for (int y = 1; y < ny - 1; y++)
    {
        for (int x = 1; x < nx - 1; x++)
        {
            u[y][x] = (un[y][x] -
                un[y][x] * dt / dx *
                (un[y][x] - un[y][x - 1]) -
                vn[y][x] * dt / dy *
                (un[y][x] - un[y - 1][x]) -
                dt / (2 * rho * dx) * (p[y][x + 1] - p[y][x - 1]) +
                nu * (dt / pow(dx, 2) *
                    (un[y][x + 1] - 2 * un[y][x] + un[y][x - 1]) +
                    dt / pow(dy, 2) *
                    (un[y + 1][x] - 2 * un[y][x] + un[y - 1][x])));

            v[y][x] = (vn[y][x] -
                un[y][x] * dt / dx *                //vn ??
                (vn[y][x] - vn[y][x - 1]) -
                vn[y][x] * dt / dy *
                (vn[y][x] - vn[y - 1][x]) -
                dt / (2 * rho * dy) * (p[y + 1][x] - p[y - 1][x]) +
                nu * (dt / pow(dx, 2) *
                    (vn[y][x + 1] - 2 * vn[y][x] + vn[y][x - 1]) +
                    dt / pow(dy, 2) *
                    (vn[y + 1][x] - 2 * vn[y][x] + vn[y - 1][x])));
        }
    }
    //@ right 
    for (int y = 0; y < ny; y++)
    {
        u[y][nx - 1] = 0.;
        v[y][nx - 1] = 0.;
    }
    //@ left 
    for (int y = 0; y < ny; y++)
    {
        u[y][0] = 0.;
        v[y][0] = 0.;
    }
    //@ bottom 
    for (int x = 0; x < nx; x++)
    {
        u[ny - 1][x] = 0.;
        v[ny - 1][x] = 0.;
    }
    //@ top 
    for (int x = 0; x < nx; x++)
    {
        u[0][x] = 1.;
        v[0][x] = 0.;
    }

}

void print_matrix(std::vector<std::vector<float>>& matrix) {
    //prints given matrix to console
    for (int y = 0; y < matrix.size(); y++)
    {
        for (int x = 0; x < matrix[0].size(); x++)
        {
            std::cout << " " << matrix[y][x];
        }
        std::cout << std::endl;
    }
}

void update_min_max_field_values(std::vector<std::vector<float>>& p, float& p_min, float & p_max, std::vector<std::vector<float>>& u, std::vector<std::vector<float>>& v, float& vel_max) {
    //finds and update values pressure min, pressure max and absolute volocity max
    
    p_min = 0;
    p_max = 0;
    vel_max = 0;
    float vel;
    for (int y = 0; y < ny; y++)
    {
        for (int x = 0; x < nx; x++)
        {
            p_min = std::min(p_min, p[y][x]);
            p_max = std::max(p_max, p[y][x]);

            vel = pow(u[x][x] * u[x][x] + v[x][x] * v[x][x], 0.5);
            vel_max = std::max(vel_max, vel);
            
        }
    }
}


void draw_velocty_vector(sf::RenderWindow& window, sf::RectangleShape& velocity_vector, std::vector<std::vector<float>>& u, std::vector<std::vector<float>>& v, float total_speed, int _x, int _y, int scale) {
    //draw line depending on the velocity magnitude and direction @position _x, _y on window class

    velocity_vector.setSize(sf::Vector2f(total_speed * scale, 1));
    velocity_vector.setRotation(atan(-v[_y][_x] / u[_y][_x]) * 180 / PI);
    velocity_vector.setPosition(_x * dx * 1000, _y * dy * 1000);
    window.draw(velocity_vector);
}

void draw_pressure_contour(sf::RenderWindow& window, sf::RectangleShape& pressure_rectangle, std::vector<std::vector<float>>& p, float _p_max, float _p_min, int _x, int _y) {
    //draws contours of pressure on window class

    pressure_rectangle.setSize(sf::Vector2f(dx * 1000, dy * 1000));
    pressure_rectangle.setPosition(_x * dx * 1000, _y * dy * 1000);
    pressure_rectangle.setFillColor(sf::Color(125 * (-_p_min + p[_y][_x]) / _p_max, 125 * (-_p_min + p[_y][_x]) / _p_max, 125 * (-_p_min + p[_y][_x]) / _p_max));
    window.draw(pressure_rectangle);
}


int main()
{
    std::vector<std::vector<float>> u;
    std::vector<std::vector<float>> un;
    std::vector<std::vector<float>> v;
    std::vector<std::vector<float>> vn;
    std::vector<std::vector<float>> p;
    std::vector<std::vector<float>> pn;
    std::vector<std::vector<float>> b;

    float abs_vel;
    float p_max, p_min, vel_max;
    float velocity_arrow_scale_factor;

    u = matrix_zeros(nx, ny);
    un = u;
    v = matrix_zeros(nx, ny);
    vn = v;
    p = matrix_zeros(nx, ny);
    pn = p;
    b = matrix_zeros(nx, ny);

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Cavity");

    sf::RectangleShape pressure_rectangle;
    sf::RectangleShape velocity_vector;

    window.setFramerateLimit(60);

    while (window.isOpen())
    {
        window.clear();
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed) window.close();

            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) window.close();
        }

        cavity_flow(u, un, v, vn, p, pn, b);

        update_min_max_field_values(p, p_min, p_max, u, v, vel_max);

        velocity_arrow_scale_factor = vel_max * 100.;

        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                abs_vel = sqrt(u[y][x] * u[y][x] + v[y][x] * v[y][x]);

                draw_pressure_contour(window, pressure_rectangle, p, p_max, p_min, x, y);

                if (abs(u[y][x]) > 0.01) { 
                    draw_velocty_vector(window, velocity_vector, u, v, abs_vel, x, y, velocity_arrow_scale_factor);
                }

            }
        }

        window.display();
    }
}