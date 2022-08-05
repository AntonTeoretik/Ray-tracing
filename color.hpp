#ifndef COLOR_H
#define COLOR_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "glm/glm.hpp"
#include "glm/gtx/projection.hpp"
#include "constants.hpp"
#include <math.h>
//#include <vector>
#include "bitimage/bitmap_image.hpp"
#include <time.h>
#include <chrono>
#include <thread>

double minl = 0.43;
double maxl = 0.73;
double dif = maxl - minl;

std::vector <double> RED =     {1.0, 1.0, 0.5, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1};
//std::vector <double> ORANGE =  {0.1, 0.5, 1.0, 0.5, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//std::vector <double> YELLOW =  {0.0, 0.1, 0.5, 1.0, 0.5, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0};
std::vector <double> GREEN  =  {0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.5, 0.1, 0.0, 0.0, 0.0};
//std::vector <double> BLUE_L =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0, 0.5, 0.1, 0.0};
std::vector <double> BLUE   =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 1.0, 0.3, 0.1};
//std::vector <double> PURPLE =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.5, 1.0};

/*
struct Color
{
    std::vector <double> spectrum;

    Color(glm::vec3 col = glm::vec3(0,0,0), double intensity = 1.0)
    {
        spectrum = std::vector < double >(11);
        //spectrum[1] = col.r; // Магия со спектрами
        //spectrum[2] = col.r; // Магия со спектрами
        //spectrum[4] = col.g; // (т.е. костыли со спектрами)
        //spectrum[5] = col.g; // (т.е. костыли со спектрами)
        //spectrum[8] = col.b; //
        for(int i = 0; i < spectrum.size(); i++)
            spectrum[i] = col.r * RED[i] + col.g * GREEN[i] + col.b * BLUE[i];
        normalize();
        multValue(intensity);
    }
    Color(std::vector <double> _spectrum)
    {
        spectrum = _spectrum;
    }

    double getIntensity()
    {
        double I = 0.0;
        for(double d : spectrum)
            I += d * d / spectrum.size(); //Шаг интегрирования
        return std::sqrt(I);
    }

    Color normalize()
    {
        double I = getIntensity();
        if(I == 0.0) return *this;
        for(double &d : spectrum)
            d /= I;
        return *this;
    }
    Color normalizeMax()
    {
        double max = 0;
        for(double d : spectrum)
            max = std::max(d, max);
        for(double &d : spectrum)
            d /= max;
        return *this;
    }
    rgb_t getColor()
    {
        Color C;
        C.spectrum = spectrum;

        multColors(Color(RED).normalize());
        double R = std::min(getIntensity(), 1.0);
        spectrum = C.spectrum;

        multColors(Color(GREEN).normalize());
        double G = std::min(getIntensity(), 1.0);
        spectrum = C.spectrum;

        multColors(Color(BLUE).normalize());
        double B = std::min(getIntensity(), 1.0);
        spectrum = C.spectrum;

        rgb_t Col = {(unsigned char)(R * 255), (unsigned char)(G * 255), (unsigned char)(B * 255)};
        //rgb_t Col = {0, int(G * 255), 0};
        return Col;
    }

    void multValue(double k)
    {
        for(double &d : spectrum)
            d *= k;
    }
    void sumColors(Color L2, double a = 1.0, double b = 1.0)
    {
        for(int i = 0; i < spectrum.size(); i++)
            spectrum[i] = a * spectrum[i] + b * L2.spectrum[i];
    }
    void multColors(Color L2, double a = 1.0, double b = 1.0)
    {
        for(int i = 0; i < spectrum.size(); i++)
            spectrum[i] = std::sqrt(a * spectrum[i] * b * L2.spectrum[i]);
    }
    void printCol()
    {
        std::cout << std::endl;
        for(double d : spectrum) {
            for(int i = 0; i < d * 10; i++)
                std::cout << "#";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
};
*/
struct Color
{
    std::vector <double> spectrum = std::vector <double> (3);

    Color(glm::vec3 col = glm::vec3(0,0,0), double intensity = 1.0)
    {
        //spectrum = std::vector < double >(3);
        spectrum[0] = col.r; // Магия со спектрами
        spectrum[1] = col.g; // Магия со спектрами
        spectrum[2] = col.b; // (т.е. костыли со спектрами)
        //normalize();
        //multValue(intensity);
    }
    Color(std::vector <double> _spectrum)
    {
        spectrum = _spectrum;
    }
    ~Color(){}
    double getIntensity()
    {
        double I = 0.0;
        for(double d : spectrum)
            I += d * d / spectrum.size(); //Шаг интегрирования
        return std::sqrt(I);
    }

    Color normalize()
    {
        double I = getIntensity();
        if(I == 0.0) return *this;
        for(double &d : spectrum)
            d /= I;
        return *this;
    }
    Color normalizeMax()
    {
        double max = 0;
        for(double d : spectrum)
            max = std::max(d, max);
        for(double &d : spectrum)
            d /= max;
        return *this;
    }
    rgb_t getColor()
    {
        rgb_t Col = {(unsigned char)( std::min(spectrum[0], 1.0) * 255),
                     (unsigned char)( std::min(spectrum[1], 1.0) * 255),
                     (unsigned char)( std::min(spectrum[2], 1.0) * 255)};
        return Col;
    }

    Color multValue(double k)
    {
        for(double &d : spectrum)
            d *= k;
        return *this;
    }
    Color sumColors (Color L2, double a = 1.0, double b = 1.0)
    {
        for(int i = 0; i < spectrum.size(); i++)
            spectrum[i] = a * spectrum[i] + b * L2.spectrum[i];
        return *this;
    }
    Color multColors(Color L2, double a = 1.0, double b = 1.0)
    {
        for(int i = 0; i < spectrum.size(); i++)
            spectrum[i] = std::sqrt(a * spectrum[i] * b * L2.spectrum[i]);
        return *this;
    }
    void printCol()
    {
        std::cout << std::endl;
        for(double d : spectrum) {
            for(int i = 0; i < d * 10; i++);
                std::cout << ".";
                //std::printf(".");

            std::cout << std::endl;
        }
    }
};

Color multValue(Color L1, double k)
{
    Color C;
    for(int i = 0; i < C.spectrum.size(); i++)
        C.spectrum[i] = k * L1.spectrum[i];
    return C;
}
Color sumColors(Color L1, Color L2, double a = 1.0, double b = 1.0)
{
    Color C;
    for(int i = 0; i < C.spectrum.size(); i++)
        C.spectrum[i] = a * L1.spectrum[i] + b * L2.spectrum[i];
    return C;
}
Color multColors(Color L1, Color L2, double a = 1.0, double b = 1.0)
{
    Color C;
    for(int i = 0; i < C.spectrum.size(); i++)
        C.spectrum[i] = std::sqrt(a * L1.spectrum[i] * b * L2.spectrum[i]);
    return C;
}



#endif // COLOR_H
