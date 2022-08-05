#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#define INFINITY

#include <math.h>

#include <time.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>


#include "bitimage/bitmap_image.hpp" //The library, using for create image, set pixel and save to bmp.

#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"
#include "scene.hpp"
#include "object.hpp"
#include "operation.hpp"
//#include <add_func.hpp>

void print_mat4x4(glm::mat4 M)
{
    std::cout << M[0][0] << " " << M[0][1] << " " << M[0][2] << " " << M[0][3] << std::endl;
    std::cout << M[1][0] << " " << M[1][1] << " " << M[1][2] << " " << M[1][3] << std::endl;
    std::cout << M[2][0] << " " << M[2][1] << " " << M[2][2] << " " << M[2][3] << std::endl;
    std::cout << M[3][0] << " " << M[3][1] << " " << M[3][2] << " " << M[3][3] << std::endl;
}

// Сделана система отрезков. Добить эту чертову хрень!
// Реализовать трассировку унарной и бинарной операции быстрым поиском (если позволяет)
// Подкласс ограниченных объектов (хитбоксы). Аналогично с операциями
// Нормали

/// Починить mode 1
/// Отражение - нечеткое. Зеркальное и диффузное.
// Преломление.
// Освещение - фотонные карты - хэш функции!       //
/// Графика - кусками.
///
/// Постобработка

int main()
{

    BluredLight* L = new BluredLight;
    Ball* B = new Ball(2);

    L->col = Color(glm::vec3(1, 1, 1));
    L->radius = 1;
    L->position = glm::vec3(0.0, -2.0, -2.0);
    L->Intensity = 0.1;


    B->color = Color(glm::vec3(1.0, 1.0, 0.0));
    B->setParams(2.0, 0.0, 0.1, 1, 0.0, 1.0);

    Scene scene;
    scene.Lights.push_back(L);
    scene.Objects.push_back(B);

    scene.cameraPosition = glm::vec3(0, 0, -4);
    scene.cameraDirection = glm::normalize(glm::vec3(0, 0, 1));

    bitmap_image bm = scene.render(0.1, 0.1080, 0.0720, 0);
    bm.save_image("test_img.bmp");

    delete L;
    delete B;

    return 0;
}

