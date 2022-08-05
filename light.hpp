#ifndef LIGHT_H
#define LIGHT_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"
#include "color.hpp"
#include "constants.hpp"
#include <random>
struct Photon
{
    glm::vec3 point;
    glm::vec3 direction;
    Color col;
};

class Light
{
public:
    Color col;
    glm::vec3 position;
    double Intensity; // По сути влияет на количество лучей, выпущенных при делании фотонных карт
    virtual Photon getRandomRay() = 0;
    virtual ~Light(){}
};
class PointLight : public Light
{
    virtual Photon getRandomRay()
    {
        double x;
        double y;
        double z;
        do
        {
            x = 2 * double(std::rand()) / (RAND_MAX - 1) - 1.0;
            y = 2 * double(std::rand()) / (RAND_MAX - 1) - 1.0;
            z = 2 * double(std::rand()) / (RAND_MAX - 1) - 1.0;
            //std::cout << x << std::endl;
        }
        while (x*x + y*y + z*z > 1 || x*x + y*y + z*z < 0.01 );
        Photon Answ = {position, glm::normalize(glm::vec3(x, y, z)), col};
        return Answ;
    }
    virtual ~PointLight(){}
};
class ParallelLight : public Light
{
    glm::vec3 direction;
    virtual ~ParallelLight() {}
};
class BluredLight: public Light
{
public:
    double radius;
    virtual Photon getRandomRay()
    {
        double x;
        double y;
        double z;
        do
        {
            x = 2 * double(std::rand()) / (RAND_MAX - 1.0) - 1.0;
            y = 2 * double(std::rand()) / (RAND_MAX - 1.0) - 1.0;
            z = 2 * double(std::rand()) / (RAND_MAX - 1.0) - 1.0;
        }
        while (x*x + y*y + z*z > 1);
        Photon Answ = {position, glm::normalize(glm::vec3(x, y, z)), col};
        do
        {
            x = 2 * radius * double(std::rand()) / (RAND_MAX - 1.0) - radius;
            y = 2 * radius * double(std::rand()) / (RAND_MAX - 1.0) - radius;
            z = 2 * radius * double(std::rand()) / (RAND_MAX - 1.0) - radius;
        }
        while (x*x + y*y + z*z > radius*radius);
        Answ.point += glm::vec3(x, y, z);
        return Answ;
    }
    virtual ~BluredLight(){}
};

#endif // LIGHT_H
