#ifndef HERA_COMMON_DIAGRAM_TRAITS_H
#define HERA_COMMON_DIAGRAM_TRAITS_H

namespace hera {

template<class PairContainer_, class PointType_ = typename std::remove_const_t<std::remove_reference_t< decltype(*std::declval<PairContainer_>().begin())>> >
struct DiagramTraits
{
    using Container = PairContainer_;
    using PointType = PointType_;
    using RealType  = typename std::remove_reference< decltype(std::declval<PointType>()[0]) >::type;

    static RealType get_x(const PointType& p)       { return p[0]; }
    static RealType get_y(const PointType& p)       { return p[1]; }
    static int     get_id(const PointType& p)       { return p.id; }
};

template<class PairContainer_, class RealType_>
struct DiagramTraits<PairContainer_, std::pair<RealType_, RealType_>>
{
    using RealType  = RealType_;
    using PointType = std::pair<RealType, RealType>;

    static RealType get_x(const PointType& p)       { return p.first; }
    static RealType get_y(const PointType& p)       { return p.second; }
    static int     get_id(const PointType&)         { return 0; }
};

} // end namespace hera


#endif // HERA_COMMON_DIAGRAM_TRAITS_H
