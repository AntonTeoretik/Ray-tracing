#ifndef SCENE_HPP
#define SCENE_HPP

#include <omp.h>

#include "object.hpp"
#include "light.hpp"
#include "constants.hpp"
#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"
#include "glm/gtx/projection.hpp"

#include "bitimage/bitmap_image.hpp"
#include <vector>
#include <random>
#include <time.h>

struct Voxel
{
    int x;
    int y;
    int z;
    Color col;
    int count;
};
struct PhotonMap
{
    int gridCapacity; // Число вокселей на метр

    std::vector <std::vector< Voxel > > Map;

    PhotonMap(int Max = 3951161) : gridCapacity(VOX_NUM)
    {
        Map = std::vector< std::vector <Voxel> >(Max);
    }
    unsigned getHash(glm::vec3 v3)
    {
        int x = v3.x * gridCapacity;
        int y = v3.y * gridCapacity;
        int z = v3.z * gridCapacity;
        unsigned t = unsigned((x * 101) + (y * 36833) + (z*11)) % Map.size();
        return t;
    }
    Color getColor(glm::vec3 v3)
    {
        unsigned i = getHash(v3);
        if(!Map[i].size()) return Color(glm::vec3(0.0, 0.0, 0.0));
        for(Voxel vox : Map[i])
            if(
                vox.x == int(v3.x * gridCapacity) &&
                vox.y == int(v3.y * gridCapacity) &&
                vox.z == int(v3.z * gridCapacity)
              )
            {
                return vox.col;
            }

        return Color(glm::vec3(0.0, 0.0, 0.0));
    }
    int getValue(glm::vec3 v3)
    {
        unsigned i = getHash(v3);
        if(!Map[i].size()) return 0;
        for(Voxel vox : Map[i])
            if(
                vox.x == int(v3.x * gridCapacity) &&
                vox.y == int(v3.y * gridCapacity) &&
                vox.z == int(v3.z * gridCapacity)
              )
            {
                return vox.count;
            }
        return 0;
    }

    void addPhoton(Photon F)
    {
        unsigned i = getHash(F.point);
        for(Voxel vox : Map[i])
            if(
                    vox.x == int(F.point.x * gridCapacity) &&
                    vox.y == int(F.point.y * gridCapacity) &&
                    vox.z == int(F.point.z * gridCapacity)
              ) {
                vox.col.sumColors(F.col);
                vox.count = vox.col.getIntensity();
                return;
            }
        Map[i].push_back({int(F.point.x * gridCapacity), int(F.point.y * gridCapacity), int(F.point.z * gridCapacity),
                         F.col, 1});
    }
    void printPhotonMap()
    {
        for(int i = 0; i < Map.size() - 1; i++)
        {
            if(Map[i].size())
                std::cout << "Hash = " << i << std::endl;
            for(Voxel V : Map[i])
            {
                std::cout << V.col.getColor().red << " "
                          << V.col.getColor().green << " "
                          << V.col.getColor().blue << std::endl;
                std::cout << V.x << " " << V.y << " " << V.z << std::endl;
                std::cout << "   " << V.count << std::endl;
            }
        }
    }
};

// Технические структуры
struct rayInfo
{
    glm::vec3 Pos; //Точка в пространстве (точка падения)
    glm::vec3 Dir; //Направление падения
    Object* Ob; // если null - не привязан ни к чему
};

class Scene
{
public:
    std::vector < Object* > Objects; // Объекты
    std::vector < Light* > Lights; // Источники света
    PhotonMap PhotMap; // Фотонная карта

    //Камера
    glm::vec3 cameraPosition;
    glm::vec3 cameraDirection;
    glm::vec3 cameraHorizont;

    Scene() : cameraPosition(glm::vec3(0, 0, 0)),
              cameraDirection(glm::vec3(0, 0, 1)), //Вдоль по z
              cameraHorizont(glm::vec3(1, 0, 0))  //Параллельно x
    {
        srand(time(NULL));
    }

    Color trace(glm::vec3 startPos, glm::vec3 dir, int recDepth = 4, double intensity = 1.0, int mode = 0)
    {
        if(recDepth == 0 || isInf(dir))
            return Color(glm::vec3(0.0, 0.0, 0.0));
        rayInfo rI = traceRay(startPos, dir);


        if(isInf(rI.Pos) || rI.Ob == nullptr) {
            return Color(glm::vec3(0.0, 0.0, 0.0));
        }

        std::vector <glm::vec3> v = generateRays(rI);
        // расчитываем вероятность (по сути распределение интенсивности по лучам)
        glm::vec3 reflect = v[0];
        glm::vec3 refract = v[1];
        glm::vec3 normal = rI.Ob->getNormal(rI.Pos);

        bool inside = (glm::dot(normal, rI.Dir) > 0);
        // Если мы внутри объекта - мы могли поголотиться
        if(inside)
        {
            double dist = glm::distance(rI.Pos, startPos);
            double prob = 1.0 - std::pow(rI.Ob->absorptionIndex, dist);
            intensity *= prob;
            //if(recDepth < 4)
                //std::cout << recDepth << " " << rI.Ob->absorptionIndex << std::endl;
        }
        //if(intensity < 0.0001) //Энергия меньше определенной
        //    return Color(glm::vec3(0.0, 0.0, 0.0));

        /// ---------------------------------------------------------------------------------------------- */

        //Теперь расчитываем преломление, отражение и рассеяние
        double refractionPart = 0.0;
        if(refract.x != INFINITY) {
            refractionPart = std::fabs (glm::dot (glm::normalize(refract), glm::normalize(normal) ));
            refractionPart = std::pow(refractionPart, rI.Ob->gamma);
            refractionPart *= rI.Ob->transparentPart;
        }

        double reflectionPart = rI.Ob->reflectionPart;
        double diffusePart    = rI.Ob->diffusePart;

        double sum = reflectionPart + refractionPart + diffusePart; // Нормализуем

        reflectionPart /= sum;
        refractionPart /= sum;
        diffusePart    /= sum;

        // Рекурсия!
        Color reflectCol = trace(rI.Pos, reflect, recDepth - 1, reflectionPart, mode);
        Color refractCol = trace(rI.Pos, refract, recDepth - 1, refractionPart, mode);

        Color diffuseCol;
        if(mode == 0 or mode == 2) {
            diffuseCol = getMedColor(rI.Pos, GAUSS_RAD); /// Тут константа! (Радиус Гаусса)
            if(mode == 2)
                diffuseCol.multValue(0.8); /// Константа!
        }

        if(mode == 1 or mode == 2) { // Не фотоны
            diffuseCol.sumColors( multValue(rI.Ob->color, 0.2) ); ///Константа

            for(Light* L : Lights)
            {
                glm::vec3 dirL = rI.Pos - L->position;
                if(glm::dot(-dirL, normal) > 0) { // Свет не внутри объекта
                    rayInfo rIL = traceRay(L->position, dirL);

                    //if(glm::distance(rI.Pos, L->position) > PRECISION)
                    //std::cout << glm::distance(rI.Pos, L->position) << std::endl;

                    if(glm::length(dirL) < glm::distance(rIL.Pos, L->position) + 1 * STEP) // Свет ближе объекта
                    {
                        //std::cout << "Here\n";
                        Color RC = multColors( rI.Ob->color, L->col);
                        double koef = intensity * std::fabs(glm::dot(normal, glm::normalize(rI.Dir)));

                        RC.multValue( koef );
                        diffuseCol.sumColors(RC, 1.0, reflectionPart *
                                             std::fabs(glm::dot(normal, glm::normalize(dirL))));
                    }
                }
            }
        }
        // Операции над цветами

        reflectCol.sumColors(refractCol, reflectionPart, refractionPart); //Складываем результаты отражения и преломления
        reflectCol.multColors(rI.Ob->color); // Проецируем на цвет объекта

        diffuseCol.sumColors(reflectCol, diffusePart, reflectionPart + refractionPart);


        return diffuseCol;
    }

    bitmap_image render(double screenDistance, double xScreen, double yScreen, int mode = 0) //пробная функция (по фонгу)
    {
        if(mode == 0 or mode == 2) {
            std::cout << "Building map... " << std::endl;
            updatePhotonMap(); // Построение карты
            //PhotMap.printPhotonMap();
            std::cout << "Photon map has been built" << std::endl;
        }

        int xResol = 1080; /// Константы
        int yResol = 720;
        int rayCount = 2;

        bitmap_image img(xResol, yResol);
        glm::vec3 vert = glm::cross(cameraHorizont, cameraDirection);
        for (int i = 0; i < xResol; i++)
        {
            std::cout << i << std::endl;

            for (int j = 0; j < yResol; j++)
            {
                Color col;
                for(int ii = 0; ii < rayCount; ii++)
                {
                    for(int jj = 0; jj < rayCount; jj++)
                    {
                        double x = (i*rayCount + ii) * 2 * xScreen / (xResol*rayCount) - xScreen;
                        double y = (j*rayCount + jj) * 2 * yScreen / (yResol*rayCount) - yScreen;

                        glm::vec3 startPos = cameraPosition + float(screenDistance)*cameraDirection +
                                             float(x) * cameraHorizont + float(y) * vert;
                        glm::vec3 dir = startPos - cameraPosition;

                        Color color = trace(startPos, dir, 4, 1.0, mode);
                        col.sumColors(color);
                    }
                }
                //std::cout << "j = " << j << std::endl;
                // Антиалиазинг
                col.multValue(1.0 / (rayCount * rayCount));
                //std::cout << (int)col.getColor().red << " " << (int)col.getColor().green << " " << (int)col.getColor().blue << std::endl;
                img.set_pixel(i, j, col.getColor());
            }
            //img.save_image("TEST.bmp");
        }
        return img;
    }

    void updatePhotonMap () // Перерасчет карты
    {
        srand( time(NULL) );
        for (Light * L: Lights) {//Идем по всем источникам света
            int i = 0;
            int intensity = L->Intensity * 5000000;
            double time = omp_get_wtime();

            std::vector <std::vector <Photon> > Ph(4);
            #pragma omp parallel for num_threads(1)

            for (i = 0; i < intensity; i++) { //В зависимости от яркости источника кидаем фотоны
                if(!(i % 10000 ))
                    std::cout << "   "<< i / 10000 << std::endl;
                Photon photon = L->getRandomRay();
                bool finish = false;
                bool setVal = true;
                rayInfo rI;
                while (!finish) { // Пока фотон не поглотился
                    rI = traceRay(photon.point, photon.direction);
                    if(isInf(rI.Pos)|| rI.Ob == nullptr)
                    {
                        finish = true;
                        setVal = false;
                        break;
                    }

                    std::vector <glm::vec3> v = generateRays(rI);

                    // расчитываем вероятность
                    glm::vec3 reflect = v[0];
                    glm::vec3 refract = v[1];
                    glm::vec3 normal = rI.Ob->getNormal(rI.Pos);
                    bool inside = (glm::dot(normal, rI.Dir) > 0);

                    // Если мы внутри объекта - мы могли поглотиться
                    if(inside)
                    {
                        double dist = glm::distance(rI.Pos, photon.point);
                        double prob = 1.0 - std::pow(rI.Ob->absorptionIndex, dist);
                        double r = std::rand() / (RAND_MAX - 1.0);
                        if(r > prob) // Поглотились!
                        {
                            finish = true;
                            setVal = false;
                            break;
                        }
                    }
                    //Теперь расчитываем преломление, отражение и рассеяние
                    double refractionPart = 0.0;
                    if(refract.x != INFINITY) {
                        refractionPart = std::fabs (glm::dot (glm::normalize(refract), glm::normalize(normal) ));
                        refractionPart = std::pow(refractionPart, rI.Ob->gamma);
                        refractionPart *= rI.Ob->transparentPart;
                    }
                    double reflectionPart = rI.Ob->reflectionPart;
                    double diffusePart = rI.Ob->diffusePart;

                    double sum = reflectionPart + refractionPart + diffusePart;
                    reflectionPart /= sum;
                    refractionPart /= sum;
                    diffusePart    /= sum;

                    double r = std::rand() / (RAND_MAX - 1.0);
                    //std::cout << r << " " << reflectionPart << std::endl;
                    photon.point = rI.Pos;
                    photon.col.multColors(rI.Ob->color);
                    //std::cout << Photon.col.r << " " << Photon.col.g << " " << Photon.col.b << std::endl;
                    //std::cout << " " << rI.Ob->color.r << " " << rI.Ob->color.g << " " << rI.Ob->color.b << std::endl;
                    if(r < refractionPart) {                          // Преломились
                        photon.direction = refract;
                        //std::cout << "Here" << std::endl;
                    } else if (r < reflectionPart + refractionPart) { // Отразились
                        photon.direction = reflect;
                    } else {                                          // Рассеялось
                        //std::cout << Photon.col.r << " " << Photon.col.g << " " << Photon.col.b << std::endl;
                        finish = true;
                        setVal = true;
                        break;
                    }
                }
                if(setVal) {;
                    Ph[omp_get_thread_num()].push_back(photon);
                }
            }
            for(int i = 0; i < 1; i++)
            {
                for(Photon p : Ph[i]) {
                    //std::cout << p.point.x << " " << p.point.y << " " << p.point.z << std::endl;
                    PhotMap.addPhoton(p);
                }
            }

            std::cout << omp_get_wtime() - time << std::endl;
        }

    }
    rayInfo traceRay (glm::vec3 startPos, glm::vec3 dir) // Трассировка одного луча
    {
        rayInfo rI = {startPos, dir, nullptr};
        double dist = INFINITY;

        for(Object * Ob: Objects)
        {
            glm::vec3 hit = Ob->trace(startPos, dir);

            if(glm::dot(hit - startPos, hit - startPos) < dist)
            {
                dist = glm::dot(hit - startPos, hit - startPos);
                rI.Pos = hit;
                rI.Ob = Ob;
            }
        }
        rI.Dir = dir;
        return rI;
    }
    std::vector <glm::vec3> generateRays (rayInfo HI) // Тут все законы преломления и отражения, порождают новые направления
    {
        //std::cout << "      Here1  " << HI.Pos.x << " " << HI.Pos.y << " " << HI.Pos.z << std::endl;
        glm::vec3 normal = HI.Ob->getNormal(HI.Pos);
        //std::cout << "      Here2" << std::endl;
        glm::vec3 reflect = glm::reflect(HI.Dir, normal);

        bool inside = (glm::dot(normal, HI.Dir) > 0);
        float n1 = inside ? HI.Ob->refractionIndex : 1.0; // Показатель преломления (откуда падаем)
        float n2 = inside ? 1.0 : HI.Ob->refractionIndex; // Показатель преломления (куда падаем)

        glm::vec3 refract(0.0, 0.0, 0.0);
        double cosA = std::fabs(glm::dot(glm::normalize(normal), glm::normalize(HI.Dir)) );
        double sinA = std::sqrt (1.0 - cosA * cosA);

        if (sinA * n1 / n2 >= 1.0) {
            //std::cout << "Hi!" << std::endl;
            refract.x = INFINITY;
        }
        else
        {
            glm::vec3 proj = HI.Dir - glm::proj(HI.Dir, normal);
            proj = glm::normalize(proj);

            double sinB = sinA * (n1 / n2);
            float x = sinB / std::sqrt(1.0 - sinB * sinB);
            //if (sinB != sinA)
                //std::cout << x << std::endl;
            refract = (inside ? normal : -normal ) + x * proj;
        }

        std::vector <glm::vec3> V;
        V.push_back(reflect);
        V.push_back(refract);
        return V;

    }
    Color getSelfColor (rayInfo HI, double rad = 0.1) //Расчет собственного цвета (без учета остальных!)
    {
        Color col = HI.Ob->color;
        col.multValue(getEnergy(HI.Pos, rad));
        return col;
    }
    double getEnergy(glm::vec3 pos, double rad)
    {
        double d = 0.0;
        for(double x = -rad; x <= rad; x += 1.0 / PhotMap.gridCapacity)
            for(double y =  std::sqrt(rad*rad - x*x);
                       y <= std::sqrt(rad*rad - x*x);
                       y += 1.0 / PhotMap.gridCapacity)
                for(double z =  rad*rad - x*x - y*y;
                           z <= rad*rad - x*x - y*y;
                           z += 1.0 / PhotMap.gridCapacity)
                    d += PhotMap.getValue(pos + glm::vec3(x, y, z));

        d /= M_PI * rad * rad;
        return d;
    }
    Color getMedColor(glm::vec3 pos, double rad, glm::vec3 n = INFVEC) //Сложение цветов!
    {
        Color col;
        double step = 1.0 / PhotMap.gridCapacity;

        for(double x = -rad; x <= rad; x += step)
            for(double y =  -rad; y <= rad; y += step)
                for(double z = -rad; z <= rad; z += step)
                {
                    if(x*x + y*y + z*z < rad*rad) {
                        Color c1 = PhotMap.getColor(pos + glm::vec3(x, y, z));
                        double koef = std::exp(-glm::dot(glm::vec3(x, y, z), glm::vec3(x, y, z)) );
                    //koef = 1.0 + koef;
                        col.sumColors(c1, 1.0, koef);
                    }
                }


        //col.multValue(1 / (M_PI * rad * rad * PhotMap.gridCapacity * PhotMap.gridCapacity));
        //col.multValue(1 / (4 * rad * rad * PhotMap.gridCapacity * PhotMap.gridCapacity));
        //col.multValue(1.0);

        return col;
        //return Color(0, 0, 0);
    }
};




#endif // SCENE_HPP
