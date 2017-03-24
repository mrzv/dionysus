#ifndef DIONYSUS_DIAGRAM_H
#define DIONYSUS_DIAGRAM_H

#include <vector>
#include <tuple>

namespace dionysus
{

template<class Value_, class Data_>
class Diagram
{
    public:
        using Value = Value_;
        using Data  = Data_;
        struct Point
        {
                  Point(Value b, Value d, Data dd):
                      birth(b), death(d), data(dd)      {}

            Value birth, death;
            Data  data;
        };

        using Points         = std::vector<Point>;
        using iterator       = typename Points::iterator;
        using const_iterator = typename Points::const_iterator;

    public:
        const_iterator  begin() const           { return points.begin(); }
        const_iterator  end() const             { return points.end(); }
        iterator        begin()                 { return points.begin(); }
        iterator        end()                   { return points.end(); }

        size_t  size() const                    { return points.size(); }
        void    push_back(const Point& p)       { points.push_back(p); }
        template<class... Args>
        void    emplace_back(Args&&... args)    { points.emplace_back(std::forward<Args>(args)...); }

    private:
        std::vector<Point>      points;
};

}

#endif
