#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

struct Vec3f
{
    Vec3f() :x(0.0f), y(0.0f), z(0.0f)
    {}
    Vec3f(float _x, float _y, float _z) : x(_x), y(_y), z(_z)
    {
    }
    float &operator[](int index)
    {
        int a = index % 3;
        switch(a)
        {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                return x;
        }

        return x;
    }
    float x;
    float y;
    float z;
};

struct Sphere
{
    
};

void render()
{
    const int width = 1024;
    const int height = 768;
    std::vector<Vec3f> framebuffer(width * height);

    for (size_t j = 0; j < height; j++)
    {
        for (size_t i = 0; i < width; i++)
        {
            float a = j / (float)height;
            float b = i / (float)width;
            framebuffer[i + j * width] = Vec3f(a, b, 0.0f);
        }
    }

    std::ofstream ofs;
    ofs.open("./out.ppm", std::ofstream::out | std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t j = 0; j < height*width; ++j)
    {
        for (size_t i = 0; i < 3; i++)
        {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[j][i])));
        }
    }
}

int main()
{
    render();
    return 0;
}