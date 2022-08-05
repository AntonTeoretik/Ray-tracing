#ifndef OPERATION_H
#define OPERATION_H

#include "object.hpp"
#include "constants.hpp"

/// Класс для операции с множествами (AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA)
class Operation : public Object
{
public:
    std::vector< Object* > Args;
    virtual double operator () (const glm::vec3&){return -1;}
    ///Матрицы преобразования
    virtual Object* trasform(const glm::mat4& M) //Преобразовать наш объект
    {
        for(Object* ob : Args)
            ob->trasform(M);

        return this;
    }
    virtual Object* rotate_origin(double fi, const glm::vec3& v3) // Отдельная реализация, чтобы не терять точность и время на обратных матрицах
    {
        for(Object* ob : Args)
            ob->rotate_origin(fi, v3);

        return this;
    }
    virtual Object* rotate_origin(double fi, double x, double y, double z)
    {
        return rotate_origin(fi, glm::vec3(x, y, z));
    }
    virtual Object* rotate_center(double fi, const glm::vec3& v3)
    {
        for(Object* ob : Args)
            ob->rotate_point(fi, v3, center);

        return this;
    }
    virtual Object* rotate_center(double fi, double x, double y, double z)
    {
        return rotate_center(fi, glm::vec3(x, y, z));
    }
    virtual Object* rotate_point(double fi, const glm::vec3& v3, const glm::vec3& point)
    {
        for(Object* ob : Args)
            ob->rotate_point(fi, v3, point + center);

        return this;
    }
    virtual Object* translate(const glm::vec3& v3)
    {
        for(Object* ob : Args)
            ob->translate(v3);
        return this;
    }
    virtual Object* translate(double x, double y, double z)
    {
        return translate(glm::vec3(x, y, z));
    }

    virtual ~Operation(){}
};

// Унарные операции
class UnaryOperation : public Operation
{
public:
    UnaryOperation()
    {
        Args = std::vector <Object*>(1);
    }

    virtual double operator () (const glm::vec3&){return -1;}
    virtual ~UnaryOperation();
};
class UMinus : public UnaryOperation
{
public:
    UMinus(Object* A)
    {
        Args[0] = A;
        _isSimple = A->isSimple();
    }
    virtual double operator ()(const glm::vec3& v3){return -1 * (*Args[0])(v3);}
    virtual ~UMinus();

    virtual glm::vec3 trace(const glm::vec3 &startPosition, const glm::vec3 &direction)
    {
        if (Args[0]->isSimple()) {
            return Args[0]->trace(startPosition, direction);
        } else {
            return Object::trace(startPosition, direction);
        }
    }
};

// Бинарные операции
class BinaryOperation : public Operation
{
public:
    BinaryOperation()
    {
        Args = std::vector <Object*>(2);
    }
    virtual double operator () (const glm::vec3&){return -1;}
    virtual ~BinaryOperation(){}
};
class Union: public BinaryOperation
{
public:
    Union(Object* A, Object* B)
    {
       Args[0] = A;
       Args[1] = B;
       _isSimple = A->isSimple() && B->isSimple();
       center = (A->center + B->center) / 2.0f;
    }
    virtual double operator () (const glm::vec3& v3)
    {
        return std::max((*Args[0])(v3), (*Args[1])(v3));
    }
    virtual ~Union(){}

    virtual glm::vec3 trace(const glm::vec3 &startPosition, const glm::vec3 &direction)
    {
        if (Args[0]->isSimple() && Args[1]->isSimple()) {
            //std::cout << "Simple U" << std::endl;
            Object* A = Args[0];
            Object* B = Args[1];

            std::vector < double > vertA;
            std::vector < double > vertB;
            glm::vec3 interA = startPosition;
            glm::vec3 interB = startPosition;
            //Берем все точки пересечения, которые есть
            while(!isInf(interA)) { // Ищем все точки пересечения для A
                interA = A->trace(interA, direction);
                double t = std::sqrt( glm::dot(interA - startPosition, interA - startPosition) /
                                      glm::dot(direction, direction) );
                vertA.push_back(t);
            }
            vertA.pop_back(); //Удаляем INFINITY
            //std::cout << vertA.size() << std::endl;
            while(!isInf(interB)) { // Ищем все точки пересечения для В
                interB = B->trace(interB, direction);
                //std::cout << "HERE  "<< interB.x << " " << interB.y << " " << interB.z << std::endl;

                double t = std::sqrt( glm::dot(interB - startPosition, interB - startPosition) /
                                      glm::dot(direction, direction) );
                vertB.push_back(t);
            }
            vertB.pop_back(); //Удаляем INFINITY
            //if(vertB.size() != 0) std::cout << vertA.size() << " " << vertB.size() << std::endl;
            // Немного костылей -
            // делаем систему отрезков. если у нас было нечетное число точке пересечения, то добавим INFINITY или -1

            if(vertA.size() % 2 && (*Args[0])(startPosition) > 0) vertA.insert(vertA.begin(), -1); // Внутри A
            if(vertA.size() % 2 && (*Args[0])(startPosition) < 0) vertA.insert(vertA.begin(), INFINITY); // Вне A
            if(vertB.size() % 2 && (*Args[1])(startPosition) > 0) vertB.insert(vertB.begin(), -1); // Внутри B
            if(vertB.size() % 2 && (*Args[1])(startPosition) < 0) vertB.insert(vertB.begin(), INFINITY); // Вне B

            // Костыыыыль!!!
            // Надо найти объединение всех отрезков - это тоже система отрезков.
            SegmentSystem SA, SB;
            for (int i = 0; i < vertA.size(); i += 2)
                SA.Seg.push_back(Segment(vertA[i], vertA[i+1]));
            for (int i = 0; i < vertB.size(); i += 2)
                SB.Seg.push_back(Segment(vertB[i], vertB[i+1]));
            SA = SA.UnionS(SB);///Отличие только тут
            if(SA.Seg.size() == 0)
                return INFVEC;
            if(SA.Seg[0].a < 0) { // == -1
                //std::cout << SA.Seg[0].a << " " <<  SA.Seg[0].b << std::endl;
                return startPosition + float(SA.Seg[0].b)*direction;
            }
            return startPosition + float(SA.Seg[0].a)*direction;

        } else {
            return this->Object::trace(startPosition, direction);
        }
    }

};
class Intersection: public BinaryOperation
{
public:
    Intersection(Object* A, Object* B)
    {
       Args[0] = A;
       Args[1] = B;
       _isSimple = A->isSimple() && B->isSimple();

       center = (A->center + B->center) / 2.0f;
    }
    virtual double operator () (const glm::vec3& v3)
    {
        return std::min((*Args[0])(v3), (*Args[1])(v3));
    }
    virtual ~Intersection(){}

    virtual glm::vec3 trace(const glm::vec3 &startPosition, const glm::vec3 &direction)
    {
        if (Args[0]->isSimple() && Args[1]->isSimple()) {
            //std::cout << "Simple I" << std::endl;
            Object* A = Args[0];
            Object* B = Args[1];

            std::vector < double > vertA;
            std::vector < double > vertB;
            glm::vec3 interA = startPosition;
            glm::vec3 interB = startPosition;
            //Берем все точки пересечения, которые есть
            while(!isInf(interA)) { // Ищем все точки пересечения для A
                interA = A->trace(interA, direction);
                double t = std::sqrt( glm::dot(interA - startPosition, interA - startPosition) /
                                      glm::dot(direction, direction) );
                //std::cout << "  Simple A   "
                //          << startPosition.x << " " << startPosition.y << " " << startPosition.z << std::endl;
                vertA.push_back(t);
            }
            vertA.pop_back(); //Удаляем INFINITY
            while(!isInf(interB)) { // Ищем все точки пересечения для В
                interB = B->trace(interB, direction);
                //std::cout << "  Simple B" << std::endl;

                double t = std::sqrt( glm::dot(interB - startPosition, interB - startPosition) /
                                      glm::dot(direction, direction) );
                vertB.push_back(t);
            }
            vertB.pop_back(); //Удаляем INFINITY
            // Немного костылей -
            // делаем систему отрезков. если у нас было нечетное число точке пересечения, то добавим -1


            if(vertA.size() % 2 && (*Args[0])(startPosition) > 0) vertA.insert(vertA.begin(), -1); // Внутри A
            if(vertA.size() % 2 && (*Args[0])(startPosition) <= 0) vertA.insert(vertA.begin(), INFINITY); // Вне A
            if(vertA.size() == 0 && (*Args[0])(startPosition) > 0)
            {
                vertA.insert(vertA.begin(), INFINITY);
                vertA.insert(vertA.begin(), -1);
            }

            if(vertB.size() % 2 && (*Args[1])(startPosition) > 0) vertB.insert(vertB.begin(), -1); // Внутри B
            if(vertB.size() % 2 && (*Args[1])(startPosition) <= 0) vertB.insert(vertB.begin(), INFINITY); // Вне B
            if(vertB.size() == 0 && (*Args[1])(startPosition) > 0)
            {
                vertB.insert(vertB.begin(), INFINITY);
                vertB.insert(vertB.begin(), -1);
            }

            // Костыыыыль!!!
            // Надо найти пересечение всех отрезков - это тоже система отрезков.
            SegmentSystem SA, SB;
            for (int i = 0; i < vertA.size(); i += 2)
                SA.Seg.push_back(Segment(vertA[i], vertA[i+1]));
            for (int i = 0; i < vertB.size(); i += 2)
                SB.Seg.push_back(Segment(vertB[i], vertB[i+1]));

//            if(SA.Seg.size()) {
//                std::cout << "Plane 1" << std::endl;
//                SA.printSegments();
//                std::cout << "Plane 2" << std::endl;
//                SB.printSegments();
//            }

            //std::cout << "HI" << std::endl;
            int size = SA.Seg.size();
            SA = SA.IntersectS(SB);///Отличие только тут

            //std::cout << "INT" << std::endl;
            //SA.printSegments();

            if(SA.Seg.size() == 0) {
                //if(size != 0 && SB.Seg.size() != 0)
                    //std::cout << "    Hi" << std::endl;
                return INFVEC;
            }
            if(SA.Seg[0].a < 0) // == -1
            {
                if(SA.Seg[0].b != INFINITY)
                    return startPosition + float(SA.Seg[0].b)*direction;
                return INFVEC;
            }

            glm::vec3 Answ = startPosition + float(SA.Seg[0].a)*direction;
            //std::cout << Answ.x << " " << Answ.y << " " << Answ.z << "  ---  " << SA.Seg[0].a << std::endl;
            return Answ;

        } else {
            return this->Object::trace(startPosition, direction);
        }
    }
};
class Difference: public BinaryOperation // A \ B
{
public:
    Difference(Object* A, Object* B)
    {
       Args[0] = A;
       Args[1] = B;
       center = A->center;
       _isSimple = A->isSimple() && B->isSimple();
    }
    virtual double operator () (const glm::vec3& v3)
    {
        return std::min((*Args[0])(v3), -(*Args[1])(v3));
    }
    virtual ~Difference(){}

    virtual glm::vec3 trace(const glm::vec3 &startPosition, const glm::vec3 &direction)
    {
        if (Args[0]->isSimple() && Args[1]->isSimple()) {
            //std::cout << "Simple I" << std::endl;
            Object* A = Args[0];
            Object* B = Args[1];

            std::vector < double > vertA;
            std::vector < double > vertB;
            glm::vec3 interA = startPosition;
            glm::vec3 interB = startPosition;
            //Берем все точки пересечения, которые есть
            while(!isInf(interA)) { // Ищем все точки пересечения для A
                interA = A->trace(interA, direction);
                double t = std::sqrt( glm::dot(interA - startPosition, interA - startPosition) /
                                      glm::dot(direction, direction) );
                vertA.push_back(t);
            }
            vertA.pop_back(); //Удаляем INFINITY
            while(!isInf(interB)) { // Ищем все точки пересечения для В
                interB = B->trace(interB, direction);
                //std::cout << "HERE  "<< interB.x << " " << interB.y << " " << interB.z << std::endl;

                double t = std::sqrt( glm::dot(interB - startPosition, interB - startPosition) /
                                      glm::dot(direction, direction) );
                vertB.push_back(t);
            }
            vertB.pop_back(); //Удаляем INFINITY
            // Немного костылей -
            // делаем систему отрезков. если у нас было нечетное число точке пересечения, то добавим -1
            if(vertA.size() % 2 && (*Args[0])(startPosition) > 0) vertA.insert(vertA.begin(), -1); // Внутри A
            if(vertA.size() % 2 && (*Args[0])(startPosition) < 0) vertA.insert(vertA.begin(), INFINITY); // Вне A
            if(vertB.size() % 2 && (*Args[1])(startPosition) > 0) vertB.insert(vertB.begin(), -1); // Внутри B
            if(vertB.size() % 2 && (*Args[1])(startPosition) < 0) vertB.insert(vertB.begin(), INFINITY); // Вне B
            // Костыыыыль!!!
            // Надо найти разность всех отрезков - это тоже система отрезков.
            SegmentSystem SA, SB;
            for (int i = 0; i < vertA.size(); i += 2)
                SA.Seg.push_back(Segment(vertA[i], vertA[i+1]));
            for (int i = 0; i < vertB.size(); i += 2)
                SB.Seg.push_back(Segment(vertB[i], vertB[i+1]));
            SA = SA.DifferenceS(SB);  ///Отличие только тут
            if(SA.Seg.size() == 0)
                return INFVEC;
            if(SA.Seg[0].a < 0)
                return startPosition + float(SA.Seg[0].b)*direction;
            return startPosition + float(SA.Seg[0].a)*direction;

        } else {
            return this->Object::trace(startPosition, direction);
        }
    }
};
class SymDifference: public BinaryOperation
{
public:
    SymDifference(Object* A, Object* B)
    {
       Args[0] = A;
       Args[1] = B;
       _isSimple = A->isSimple() && B->isSimple();

       center = (A->center + B->center) / 2.0f;

    }
    virtual double operator () (const glm::vec3& v3)
    {
        double a1 = (*Args[0])(v3);
        double a2 = (*Args[1])(v3);
        return std::min(  std::max( a1,  a2 ),
                         -std::min( a1,  a2 ));
    }
    virtual ~SymDifference(){}

    virtual glm::vec3 trace(const glm::vec3 &startPosition, const glm::vec3 &direction)
    {
        if (Args[0]->isSimple() && Args[1]->isSimple()) {
            //std::cout << "Simple I" << std::endl;
            Object* A = Args[0];
            Object* B = Args[1];

            std::vector < double > vertA;
            std::vector < double > vertB;
            glm::vec3 interA = startPosition;
            glm::vec3 interB = startPosition;
            //Берем все точки пересечения, которые есть
            while(isInf(interA)) { // Ищем все точки пересечения для A
                interA = A->trace(interA, direction);
                double t = std::sqrt( glm::dot(interA - startPosition, interA - startPosition) /
                                      glm::dot(direction, direction) );
                vertA.push_back(t);
            }
            vertA.pop_back(); //Удаляем INFINITY
            while(isInf(interB)) { // Ищем все точки пересечения для В
                interB = B->trace(interB, direction);
                //std::cout << "HERE  "<< interB.x << " " << interB.y << " " << interB.z << std::endl;

                double t = std::sqrt( glm::dot(interB - startPosition, interB - startPosition) /
                                      glm::dot(direction, direction) );
                vertB.push_back(t);
            }
            vertB.pop_back(); //Удаляем INFINITY
            // Немного костылей -
            // делаем систему отрезков. если у нас было нечетное число точке пересечения, то добавим -1
            if(vertA.size() % 2 && (*Args[0])(startPosition) > 0) vertA.insert(vertA.begin(), -1); // Внутри A
            if(vertA.size() % 2 && (*Args[0])(startPosition) < 0) vertA.insert(vertA.begin(), INFINITY); // Вне A
            if(vertB.size() % 2 && (*Args[1])(startPosition) > 0) vertB.insert(vertB.begin(), -1); // Внутри B
            if(vertB.size() % 2 && (*Args[1])(startPosition) < 0) vertB.insert(vertB.begin(), INFINITY); // Вне B
            // Костыыыыль!!!
            // Надо найти разность всех отрезков - это тоже система отрезков.
            SegmentSystem SA, SB;
            for (int i = 0; i < vertA.size(); i += 2)
                SA.Seg.push_back(Segment(vertA[i], vertA[i+1]));
            for (int i = 0; i < vertB.size(); i += 2)
                SB.Seg.push_back(Segment(vertB[i], vertB[i+1]));
            SA = SA.SymDifferenceS(SB);  ///Отличие только тут
            if(SA.Seg.size() == 0)
                return INFVEC;
            if(SA.Seg[0].a < 0)
                return startPosition + float(SA.Seg[0].b)*direction;
            return startPosition + float(SA.Seg[0].a)*direction;

        } else {
            return this->Object::trace(startPosition, direction);
        }
    }

};



#endif // OPERATION_H
