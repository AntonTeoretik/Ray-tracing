#ifndef OBJECT_H
#define OBJECT_H

#include <functional>
#include <algorithm>
#include <vector>
#include <iostream>
#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"
#include "segment.hpp"
#include "color.hpp"
#include "constants.hpp"
#include "bitimage/bitmap_image.hpp"

/*!
    XXX Hitbox XXX
    Hitsphere
    Object - множество с кучей параметров.
                        Центр
                        Функция множества
                        Хитбокс (базово равен бесконечности)
        Operation - абстрактная операциями над множествами или другими операциями
            UnaryOperation
                Minus
            BinaryOperation
                Union
                Intersection
                Difference
                SymmDifferense
        SimpleObject - множество с картой нормалей,
                            параметризацией и простым уравнением пересечения
                            и конечным хитбоксом
            Plane
            Ball


    Scene
*/

void print_vec3(glm::vec3 v, int precision = 3)
{
    std::cout << std::endl
              << std::round(v.x * std::pow(10, precision))/std::pow(10, precision) << " "
              << std::round(v.y * std::pow(10, precision))/std::pow(10, precision) << " "
              << std::round(v.z * std::pow(10, precision))/std::pow(10, precision) << std::endl;
}


struct Hitspere
{

    double rad;
    glm::vec3 center;

    Hitspere(double rad = INFINITY, glm::vec3 center = glm::vec3()) : rad(rad), center(center) {}

    glm::vec3 hit(glm::vec3 start, glm::vec3 dir, bool inside = false)
    {
        //Проверки
        if(inside and glm::distance(start, center) <= rad)
            return start;

        if(!inside and rad == INFINITY)
            return INFVEC;
        //решаем квадратное уравнение
        double a = glm::dot(dir, dir);
        double b = -2 * glm::dot(dir, center - start);
        double c = glm::dot(start - center, start - center) - rad*rad;
        double D = b*b - 4*a*c;
        if(D < 0) {
            return INFVEC; //Нет пересечения
        }
        D = std::sqrt(D);
        double t1 = (-b - D) / (2 * a);
        double t2 = (-b + D) / (2 * a);
        if(t2 <= 0)
            return INFVEC; //Нет пересечения
        if(t1 <= 0 && t2 / glm::length(dir) < PRECISION)
            return INFVEC; //На самом деле мы на границе, ну чуть-чуть не дошли
        if(t1 <= 0) {
            return start + float(t2) * dir; //
        }
        if(t1 > 0 && t1 / glm::length(dir) < PRECISION)
            return start + float(t2) * dir; // У границы
        return start + float(t1) * dir;
    }

    static Hitspere UnionSphere(Hitspere A, Hitspere B)
    {
        if(A.rad < B.rad) // Теперь у А радиус больше
        {
            Hitspere C = A;
            A = B;
            B = C;
        }
        // Проверка критических значений
        if(A.rad == INFINITY && B.rad == INFINITY)
            return Hitspere(INFINITY, (A.center + B.center) * 0.5f );
        if(A.rad == INFINITY)
            return Hitspere(INFINITY, A.center);
        // Теперь радиусы только конечные
        double dist = glm::length(A.center - B.center);

        if(dist + B.rad <= A.rad) // B внутри А
            return A;

        double R = (A.rad + dist + B.rad) * 0.5;
        float t = (dist + B.rad - A.rad) * 0.5;
        t /= dist;

        glm::vec3 cent = t * A.center + (1 - t) * B.center;

        return Hitspere(R, cent);

    }
    static Hitspere IntersectSphere(Hitspere A, Hitspere B)
    {
        if(A.rad < B.rad) // Теперь у А радиус больше
        {
            Hitspere C = A;
            A = B;
            B = C;
        }
        // Проверка критических значений
        if(A.rad == INFINITY && B.rad == INFINITY)
            return Hitspere(INFINITY, (A.center + B.center) * 0.5f );
        if(A.rad == INFINITY)
            return Hitspere(B.rad, B.center);

        // Теперь радиусы только конечные
        double dist = glm::length(A.center - B.center);

        if(dist + B.rad <= A.rad) // B внутри А
            return B;
        // B не внутри А
        if(A.rad + B.rad < dist) // Пустое множество
            return Hitspere(-1, INFVEC); // Критическое значение -
        // удаленный на бесконечность шар с отрицательным радиусом

        // Шары пересекаются
        // Вычисляем радиус
        double R = (dist + A.rad + B.rad) / 2; // Полупериметр

        R = std::sqrt( R * (R - A.rad) * (R - B.rad) * (R - dist) ); // Площадь треугольника
        R = (2 / dist) * R; // Высота - она же радиус

        float t = std::sqrt(A.rad*A.rad - R*R);
        t /= dist;
        glm::vec3 cent = t * A.center + (1 - t) * B.center;

        return Hitspere(R, cent);
    }
};

///Максимально абстрактный класс объекта
class Object
{
protected:
public:
    // Физические параметры
    Color color;
    double refractionIndex;    // Коэффициэнт преломления (будет еще зависеть от угла)
    double absorptionIndex;    // Коэффициэнт поглощения (часть, которая поглотится при преломлении)
    // Коэффициэнт отражения - (остальное (reflection + refraction = 1.0)

    double reflectionPart;    // Зеркальное отражение.
    double diffusePart;       // Диффузное отражение - это учитывается в фотонных картах. (Просто рассеяние)

    double transparentPart;   // Часть, которая преломится
    double gamma;             // Часть степенного закона

    //
    bool _isSimple; //редкой разновидности костыль!

    void setParams(double _refractionIndex = 2.0,
                   double _absorptionIndex = INFINITY,
                   double _reflectionPart = 0.0,
                   double _diffusePart = 1.0,
                   double _transparentPart = 0.0,
                   double _gamma = 1.0)
    {
        refractionIndex = _refractionIndex;
        absorptionIndex = _absorptionIndex;
        reflectionPart = _reflectionPart;
        diffusePart = _diffusePart;
        transparentPart = _transparentPart;
        gamma = _gamma;
    }

    glm::vec3 center = glm::vec3(0, 0, 0);
    std::function < double(const glm::vec3&) > setFunc; //>0 если принадлежит множеству, <0 если нет
    glm::mat4 transformeMatrix = glm::mat4(); //Матрица положения в пространстве.
    //проверка на простоту объекта.
    bool isSimple(){ return _isSimple; }

    Object(std::function<double(const glm::vec3&)> Pos = [](const glm::vec3&){return -1;}) :
        setFunc(Pos) { _isSimple = false; }
    virtual ~Object(){}
    ///Матрицы преобразования
    virtual Object* trasform(const glm::mat4& M) //Преобразовать наш объект
    {
        transformeMatrix = M * transformeMatrix;
        inverseTransformeMatrix *= glm::inverse(M);
        glm::vec4 v4 = M * glm::vec4(center, 1.0);
        v4 /= v4.w;
        center = glm::vec3(v4);
        return this;
    }
    virtual Object* rotate_origin(double fi, const glm::vec3& v3) // Отдельная реализация, чтобы не терять точность и время на обратных матрицах
    {
        glm::mat4 M = glm::rotate<double, glm::highp>(fi, v3);
        transformeMatrix = M * transformeMatrix;
        inverseTransformeMatrix *= glm::transpose(M);
        glm::vec4 v4 = M * glm::vec4(center, 1.0);
        v4 /= v4.w;
        center = glm::vec3(v4);
        return this;
    }
    virtual Object* rotate_origin(double fi, double x, double y, double z)
    {
        return rotate_origin(fi, glm::vec3(x, y, z));
    }
    virtual Object* rotate_center(double fi, const glm::vec3& v3)
    {
        glm::vec3 cur_center = center;
        translate(-cur_center);
        rotate_origin(fi, v3);
        translate(cur_center);
        return this;
    }
    virtual Object* rotate_center(double fi, double x, double y, double z)
    {
        return rotate_center(fi, glm::vec3(x, y, z));
    }
    virtual Object* rotate_point(double fi, const glm::vec3& v3, const glm::vec3& point)
    {
        translate(-point);
        rotate_origin(fi, v3);
        translate(point);
        return this;
    }
    virtual Object* translate(const glm::vec3& v3)
    {
        glm::mat4 M = glm::translate(v3);
        transformeMatrix = M * transformeMatrix;
        inverseTransformeMatrix *= glm::translate(-v3);
        glm::vec4 v4 = M * glm::vec4(center, 1.0);
        v4 /= v4.w;
        center = glm::vec3(v4);
        return this;
    }
    virtual Object* translate(double x, double y, double z)
    {
        return translate(glm::vec3(x, y, z));
    }

    ///Значение функции
    virtual double operator ()(const glm::vec3& v3)
    {
        //print_vec3(v3);
        glm::vec4 v4 = glm::vec4(v3, 1.0);
        v4 = inverseTransformeMatrix*v4;
        //std::cout << v4.x << " " << v4.y << " " << v4.z << std::endl;
        v4 /= v4.w;
        return setFunc(glm::vec3(v4));
    }
    ///Трассировка луча (ищем первое пересечение, перешагиваем через границу)
    virtual glm::vec3 trace(const glm::vec3& startPosition, const glm::vec3& direction)
    {
        /*! Алгоритм:
         * Проверяем на столкновение с хитбоксом, если не попали - (inf, 0, 0)
         * Идем шажками STEP, пока не поменяется знак.
         * Используем бинпоиск, чтобы увеличить точность
        */
        glm::vec3 dir = glm::normalize(direction);

        glm::vec3 currentPos = startPosition; // начало в первом пересечении



        //Если нет пересечения - вернем вектор (INFINITY, 0, 0)

        /*if(glm::distance(currentPos, startPosition) < PRECISION)
        {
            currentPos += float(STEP) * dir;
        }*/
        /*if(std::isnan(currentPos.x) || currentPos.x == -INFINITY)
        {
            print_vec3(currentPos);
            //print_vec3(previosPos);
            print_vec3(dir);
            print_vec3(direction);
            system("Pause");
        }*/
        bool sign = (operator ()(currentPos) >= 0);
        bool newsign = (operator ()(currentPos) >= 0);
        //Первый проход - идем шагами step, пока не меняется знак множества
        //Если мы изначально не попали в границу - выходим, когда улетаем за пределы сцены

        //Проверяем пробную точку!

        glm::vec3 previosPos = currentPos;

        int t = 0; // !
        //std::cout << "A" << std::endl;
        while(sign == newsign)
        {
            if(std::isnan(currentPos.x) || currentPos.x == -INFINITY)
            {
                std::cout << t << " " << operator()(previosPos) << std::endl;
                print_vec3(currentPos);
                print_vec3(previosPos);
                print_vec3(dir);
                print_vec3(direction);
                system("Pause");
            }

            if(glm::dot(currentPos, currentPos) >= MAX_RADIUS*MAX_RADIUS) {
                return INFVEC; //Костыль на проверку за выход за пределы карты. Мб по умному сделать? Хотя зачем...
            }

            sign = newsign;

            previosPos = currentPos;
            double koef = std::fabs( 5 * operator()(currentPos) ); // Можем делать шаги не больше n * STEP
            koef *= std::sqrt(std::fabs(glm::dot(getNormal(currentPos), dir)));
            koef = std::max(1.0, koef);
            //double koef = 1;

            currentPos += float(STEP * koef) * dir;
            newsign = (operator ()(currentPos) >= 0);


            if(t > 1000) {
                std::cout << operator ()(currentPos) << " " << std::endl;
                print_vec3(currentPos);
            }
            t += 1; // !
        }
        //std::cout << "B" << std::endl;
        //std::cout << "  " << t << std::endl; // !
        //Итак, теперь у нас есть два разных по знаку значения. Приближаемся бинпоиском.
        double currentValue = operator ()(currentPos);
        //glm::vec3 previosPos = currentPos - stepVec;
        glm::vec3 middle = currentPos;
        double precQ = PRECISION*PRECISION;
        t = 0;

        while( glm::dot(previosPos - currentPos, previosPos - currentPos) >= precQ) //Пока расстояние между точками справа и слева больше точности
        {

            middle = float(0.5) * (previosPos + currentPos);
            if(operator ()(middle) * operator ()(previosPos) > 0) {// знак центра совпадает со знаком предыдущего
                previosPos = middle;
            } else {
                currentPos = middle;
            }
            currentValue = operator ()(currentPos); //Эта строчка означает, что мы сто процентов пересечем границу
            //std::cout << " "<< operator ()(currentPos)<< " ";
            t += 1;
        }

        //std::cout << "  " << t << std::endl;
        //std::cout << middle.x << std::endl;
        return currentPos;
    }
    virtual glm::vec3 getNormal(const glm::vec3& Position)
    {
        //std::cout << "          Here1 " << std::endl;
        //std::cout << "          Here2 " << currentValue << std::endl;
        // Ищем градиент функции - он же нормаль?
        double dFx  = operator ()(Position + glm::vec3(STEP, 0, 0));
        double dFxM = operator ()(Position - glm::vec3(STEP, 0, 0));

        double dFy  = operator ()(Position + glm::vec3(0, STEP, 0));
        double dFyM = operator ()(Position - glm::vec3(0, STEP, 0));

        double dFz  = operator ()(Position + glm::vec3(0, 0, STEP));
        double dFzM = operator ()(Position - glm::vec3(0, 0, STEP));

        return -glm::normalize(glm::vec3(dFx - dFxM, dFy - dFyM, dFz - dFzM));
    }

private:
    glm::mat4 inverseTransformeMatrix = glm::mat4(); //Обратная матрица (нужна для проверки значения функции)
};
class SimpleObject : public Object
{
public:
    // Новые переменные - карта нормалей!
    std::function <glm::vec3(glm::vec3&)> normal;

    //Конструктор!
    SimpleObject(std::function<double(const glm::vec3&)> Pos = [](const glm::vec3&){return -1;},
                 std::function<glm::vec3(glm::vec3&)> n = [](const glm::vec3&){return glm::vec3(0,0,0);}) :
        normal(n)
    {
        setFunc = Pos; _isSimple = true;
    }

    ~SimpleObject(){}

    //
};

//Шары
class Ball : public SimpleObject
{
    double r;
public:
    Ball(double r = 1, const glm::vec3& pos = glm::vec3(0, 0, 0)): r(r)
    {
        setFunc = [r, pos](const glm::vec3& v3)
        {

            return r*r - ((v3.x - pos.x)*(v3.x - pos.x)
                       + (v3.y - pos.y)*(v3.y - pos.y)
                       + (v3.z - pos.z)*(v3.z - pos.z) );

        };
        //std::cout << setFunc(glm::vec3(0, 0, 0)) << std::endl;
        center = pos;
    }
    virtual glm::vec3 trace(const glm::vec3& startPosition, const glm::vec3& direction)
    {
        //решаем квадратное уравнение
        double a = glm::dot(direction, direction);
        double b = -2 * glm::dot(direction, center - startPosition);
        double c = glm::dot(startPosition - center, startPosition - center) - r*r;
        double D = b*b - 4*a*c;
        if(D < 0) {
            return INFVEC; //Нет пересечения
        }
        D = std::sqrt(D);
        double t1 = (-b - D) / (2 * a);
        double t2 = (-b + D) / (2 * a);
        if(t2 <= 0)
            return INFVEC; //Нет пересечения
        if(t1 <= 0 && (t2*t2)/glm::dot(direction, direction) < STEP * STEP)
            return INFVEC; //На самом деле мы на границе, ну чуть-чуть не дошли
        if(t1 <= 0) {
            return startPosition + float(t2) * direction; //
        }
        if(t1 > 0 && (t1*t1)/glm::dot(direction, direction) < STEP * STEP)
            return startPosition + float(t2) * direction; // У границы
        return startPosition + float(t1) * direction;
    }

    ~Ball(){}
};
//Плоскость
class Plane : public SimpleObject
{
    glm::vec3 n;
public:
    Plane(const glm::vec3& normal = glm::vec3(0,0,1), const glm::vec3& pos = glm::vec3(0,0,0)): n(normal)
    {
        n = glm::normalize(normal);
        center = pos;
        setFunc = [this](const glm::vec3 & v3)
        {
            return glm::dot(this->n, this->center - v3);
        };
    }
    virtual double operator ()(const glm::vec3& v3)
    {
        return setFunc(v3);
    }

    virtual Object* trasform(const glm::mat4& M) //Преобразовать наш объект
    {
        glm::vec4 n4 = M * glm::vec4(n, 1.0);
        n4 /= n4.w;
        n = glm::normalize(glm::vec3(n4));
        return this->Object::trasform(M);
    }
    virtual Object* rotate_origin(double fi, const glm::vec3& v3) // Отдельная реализация, чтобы не терять точность и время на обратных матрицах
    {
        glm::mat4 M = glm::rotate<double, glm::highp>(fi, v3);
        glm::vec4 n4 = M * glm::vec4(n, 1.0);
        n4 /= n4.w;
        n = glm::normalize(glm::vec3(n4));
        return this->Object::rotate_origin(fi, v3);
    }
    virtual Object* rotate_center(double fi, const glm::vec3 &v3)
    {
        glm::mat4 M = glm::rotate<double, glm::highp>(fi, v3);
        glm::vec4 n4 = M * glm::vec4(n, 1.0);
        n4 /= n4.w;
        n = glm::normalize(glm::vec3(n4));
        return this->Object::rotate_center(fi, v3);
    }
    virtual Object* rotate_point(double fi, const glm::vec3 &v3, const glm::vec3 &point)
    {
        glm::mat4 M = glm::rotate<double, glm::highp>(fi, v3);
        glm::vec4 n4 = M * glm::vec4(n, 1.0);
        n4 /= n4.w;
        n = glm::normalize(glm::vec3(n4));
        return this->Object::rotate_point(fi, v3, point);
    }

    virtual ~Plane(){}

    virtual glm::vec3 trace(const glm::vec3& startPosition, const glm::vec3& direction)
    {
        glm::vec3 dir = glm::normalize(direction);
        if(std::fabs(glm::dot(n, center - startPosition)) <= STEP)
            return INFVEC; //Слишком близко
        //решаем линейное уравнение
        double a = glm::dot(n, dir);
        double b = glm::dot(n, center - startPosition);

        if(a == 0) {
            //std::cout << "Here" << std::endl;
            return INFVEC; //Не попали
        }
        double t = b / a;
        if(t <= STEP)
            return INFVEC; //Не попали
        return startPosition + float(t) * dir;

    }
};


#endif // OBJECT_H
