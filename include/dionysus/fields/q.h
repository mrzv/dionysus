#ifndef DIONYSUS_Q_H
#define DIONYSUS_Q_H

#include <iostream>

// TODO: eventually need to be able to adaptively switch to arbitrary precision arithmetic

namespace dionysus
{

template<typename Element_ = long>
class Q
{
    public:
        using           BaseElement = Element_;
        struct Element
        {
            BaseElement numerator, denominator;

            bool        operator==(Element o) const                     { return numerator == o.numerator && denominator == o.denominator; }
            bool        operator!=(Element o) const                     { return !((*this) == o); }

            friend
            std::ostream&   operator<<(std::ostream& out, Element e)    { out << e.numerator << '/' << e.denominator; return out; }
        };

        Element         id()  const                         { return { 1,1 }; }
        Element         zero()  const                       { return { 0,1 }; }
        Element         init(BaseElement a) const           { return { a,1 }; }

        Element         neg(Element a) const                { return { -a.numerator, a.denominator }; }
        Element         add(Element a, Element b) const     { return { a.numerator*b.denominator + b.numerator*a.denominator, a.denominator*b.denominator }; }

        Element         inv(Element a) const                { return { a.denominator, a.numerator }; }
        Element         mul(Element a, Element b) const     { return { a.numerator*b.numerator, a.denominator*b.denominator }; }
        Element         div(Element a, Element b) const     { return mul(a, inv(b)); }

        bool            is_zero(Element a) const            { return a.numerator == 0; }
};

}

#endif
