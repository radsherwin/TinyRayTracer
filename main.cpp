#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
struct Vec3f
{
    Vec3f() :x(0.0f), y(0.0f), z(0.0f)
    {
    }
    Vec3f(const float &_x, const float &_y, const float &_z) : x(_x), y(_y), z(_z)
    {
    }
    const float &operator[](const int &index) const
    {
        int a = index % 3;
        switch (a)
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

    Vec3f operator-(const Vec3f &v) const
    {
        return Vec3f(x - v.x, y - v.y, z - v.z);
    }
    Vec3f operator+(const Vec3f &v) const
    {
        return Vec3f(x + v.x, y + v.y, z + v.z);
    }

    Vec3f operator*(const Vec3f &v) const
    {
        return Vec3f(x * v.x, y * v.y, z * v.z);
    }

    Vec3f operator*(const float &f) const
    {
        return Vec3f(x * f, y * f, z * f);
    }

    Vec3f operator/=(const Vec3f &v)
    {
        this->x /= v.x;
        this->y /= v.y;
        this->z /= v.z;

        return *this;
    }

    float dot(const Vec3f &v)
    {
        return (this->x * v.x + this->y * v.y + this->z * v.z);
    }

    const float mag()
    {
        return sqrtf(x * x + y * y + z * z);
    }

    Vec3f norm()
    {
        const float mag = this->mag();

        x /= mag;
        y /= mag;
        z /= mag;

        return *this;

    }

    float x;
    float y;
    float z;
};

struct Sphere
{
    Vec3f center;
    float radius;

    Sphere(const Vec3f &c, const float &r) : center(c), radius(r)
    {
    }

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const
    {
        Vec3f L = center - orig;
        float tca = L.dot(dir);
        float d2 = L.dot(L) - tca * tca;
        if (d2 > radius * radius) return false;
        float thc = sqrtf(radius * radius - d2);
        t0 = tca - thc;
        float t1 = tca + thc;
        if(t0 > t1) std::swap(t0, t1);
        if (t0 < 0)
        {
            t0 = t1;
            return false;
        }
        return true;
    }
};

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const Sphere &sphere)
{
    float sphere_dist = std::numeric_limits<float>::max();
    if (!sphere.ray_intersect(orig, dir, sphere_dist))
    {
        return Vec3f(0.2f, 0.7f, 0.3f);
    }

    return Vec3f(0.4f, 0.4f, 0.3f);
}

// Draw a "ray" from center of camera, through exch pixel, and check if that ray intersects with the sphere
void render(const Sphere &sphere)
{
    const int width = 1024;
    const int height = 768;
    const int fov = 3.14159/2;
    std::vector<Vec3f> framebuffer(width * height);

    for (size_t j = 0; j < height; j++)
    {
        for (size_t i = 0; i < width; i++)
        {
            float x = (2 * (i + 0.5f) / (float)width - 1.0f) * tanf(fov / 2.0f) * width / (float)height;
            float y = -(2 * (j + 0.5f) / (float)width - 1.0f) * tanf(fov / 2.0f);
            Vec3f dir = Vec3f(x, y, -1.0f).norm();
            framebuffer[i + j * width] = cast_ray(Vec3f(0,0,0), dir, sphere);
        }
    }

    std::ofstream ofs;
    ofs.open("./out.ppm", std::ofstream::out | std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t j = 0; j < height * width; ++j)
    {
        for (size_t i = 0; i < 3; i++)
        {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[j][i])));
        }
    }
}

int main()
{
    Sphere sphere(Vec3f(-3,0,-16), 2.0f);
    render(sphere);
    return 0;
}