#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP
#include "glm/glm.hpp"

// Критическое значение вектора - проверка на бесконечность
const glm::vec3 INFVEC = glm::vec3(INFINITY, INFINITY, INFINITY);
bool isInf(glm::vec3 v3) {return v3.x == INFINITY and v3.y == INFINITY and v3.z == INFINITY;}

// Шаг и точность
const double MAX_RADIUS = 200; // Дальше этого радуса сцены нет
const double STEP = 0.01; // Ближе этого шага мы считаем что точки близко
const double PRECISION = 0.001; // Ближе этого шага мы считаем что точки близко

// Радиус Гаусса
const double GAUSS_RAD = 0.1;

// Число вокселей на единицу
int VOX_NUM = 10;

#endif // CONSTANTS_HPP
