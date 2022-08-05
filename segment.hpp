#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>

struct Segment
{
    double a;
    double b;
    Segment(double a = -INFINITY, double b = INFINITY, bool sort = true) :
        a(sort ? std::min(a, b): a),
        b(sort ? std::max(a, b): b)
    {} // сортируем если true
    void printSeg()
    {
        std::cout << a << " " << b << std::endl;
    }
};

struct SegmentSystem
{
    std::vector <Segment> Seg;
    SegmentSystem(std::vector <Segment> segments = std::vector <Segment>()) : Seg(segments) {}
    void sortSegs()
    {
        std::sort(Seg.begin(), Seg.end(), [](Segment A, Segment B){
            return A.a < B.a;});
    }
    SegmentSystem Normalyze() // Делаем адекватный вид
    {
        for(int i = 0; i < Seg.size(); i++) {
            if(Seg[i].b - Seg[i].a < 0.0001)
            {
                Seg.erase(Seg.begin() + i);
            }
        }
        if (Seg.size() == 0) return SegmentSystem();
        // Сначала удаляем отрезки-точки


        //Идея со скобочками
        std::vector < Segment > brackets; // Segment.a = координата. Segment.b = скобочка (число 1 или -1)
        for(Segment i : Seg) {
            brackets.push_back( Segment (i.a,  1, false) );
            brackets.push_back( Segment (i.b, -1, false) );
        }
        // Кладем вершины в массив с соответсвующими скобочками

        std::sort(brackets.begin(), brackets.end(), [](Segment A, Segment B){return A.a <= B.a;} );
        // Сортируем - получаем скобочную последовательность (дожны получить)

        double counter = 0.0;
        SegmentSystem S;
        Segment A;

        for(Segment i : brackets)
        {
            counter += i.b;
            if(i.b > 0 && std::fabs(counter - 1.0) < 0.1 )
                A.a = i.a; // Кладем левый конец
            if(i.b < 0 && std::fabs(counter) < 0.1 ) {
                A.b = i.a; // Кладем правый конец и сохраняем
                S.Seg.push_back(A);
            }
        }
        // проверяем на совпадающие точки
        for (int i = 0; i < S.Seg.size()-1; i++)
        {
            if( S.Seg[i].b >= S.Seg[i+1].a) {
                S.Seg[i].b = S.Seg[i+1].b;
                S.Seg.erase(S.Seg.begin() + i + 1); // Удаляем i+1-й элемент
            }
        }
        return S;
    }
    SegmentSystem Minus()
    {
        Seg = this->Normalyze().Seg; // уже отсортирован
        Segment A;
        SegmentSystem S;
        if(Seg.size() == 0) {
            A.a = -INFINITY;
            A.b = INFINITY;
            S.Seg.push_back(A);
            return S;
        }
        // есть хотя бы один отрезок
        if(Seg[0].a != -INFINITY) {
            A.b = Seg[0].a;
            S.Seg.push_back(A);
        }
        for(int i = 0; i < Seg.size() - 1; i++) {
            A.a = Seg[i].b;
            A.b = Seg[i + 1].a;
            S.Seg.push_back(A);
        }
        if(Seg[Seg.size() - 1].b != INFINITY) {
            A.a = Seg[Seg.size() - 1].b;
            A.b = INFINITY;
            S.Seg.push_back(A);
        }
        return S;
    }
    SegmentSystem UnionS(SegmentSystem Sys)
    {
        SegmentSystem Res;
        Res.Seg = Seg;
        for(Segment i : Sys.Seg)
            Res.Seg.push_back(i);
        return Res.Normalyze();
    }
    SegmentSystem IntersectS(SegmentSystem Sys)
    {
        Sys = Sys.Normalyze();

        for(Segment i : (this->Normalyze()).Seg)
            Sys.Seg.push_back(i);
        //Идея со скобочками
        std::vector < Segment > brackets; // Segment.a = координата. Segment.b = скобочка (число 1 или -1)
        for(Segment i : Sys.Seg) {
            brackets.push_back( Segment (i.a,  1, false) );
            brackets.push_back( Segment (i.b, -1, false) );
        }
        // Кладем вершины в массив с соответсвующими скобочками

        std::sort(brackets.begin(), brackets.end(), [](Segment A, Segment B){return A.a < B.a;} );
        // Сортируем - получаем скобочную последовательность (дожны получить)
        SegmentSystem S;
        Segment A;
        double count = 0.0;
        if(brackets.size() == 0) {
            return S;
        }
        for (int i = 0; i < brackets.size()-1; i++)
        {
            count += brackets[i].b;
            //std::cout << "HI3" << std::endl;
            if(brackets[i].b > 0 && brackets[i + 1].b < 0 && count > 1) {

                A.a = brackets[i].a; // Кладем левый конец
                A.b = brackets[i+1].a; // Кладем правый конец
                S.Seg.push_back(A); // и сохраняем
            }
        }
        //std::cout << "before" << std::endl;
        //S.printSegments();
        S = S.Normalyze();
        //std::cout << "after "<< std::endl;
        return S;
    }
    SegmentSystem DifferenceS(SegmentSystem Sys) {
        return this->IntersectS(Sys.Minus());
    }
    SegmentSystem SymDifferenceS(SegmentSystem Sys) {
        return this->UnionS(Sys).DifferenceS(this->IntersectS(Sys));
    }

    void printSegments()
    {
        if(Seg.size() != 0)
            for (Segment i : Seg)
                i.printSeg();
        std::cout << std::endl;
    }
};

#endif // SEGMENT_HPP
