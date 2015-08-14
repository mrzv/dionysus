#include <cassert>

template<class K, class T, class F>
std::tuple<bool,bool,bool>
RelativeLZZ<K,T,F>::
lower_upper(const Simplex& s, const ValueVertex& vval) const
{
    Vertex v = std::get<1>(vval);
    bool lower = false, upper = false, has_v = false;
    for (auto& u : s)
    {
        if (u == v)
        {
            has_v = true;
            continue;
        }

        ValueVertex uval(function(u), u);
        if (uval < vval)
            lower = true;
        else if (uval > vval)
            upper = true;
    }

    return std::make_tuple(lower,upper,has_v);
}

template<class K, class T, class F>
void
RelativeLZZ<K,T,F>::
add_both(const Simplex& s)
{
    if (cmplx.find(s) != cmplx.end())       // already in the complex, skip
        return;

    //fmt::print("Adding to both: {}\n", s);

    zz.add_both(s.boundary(zz.field()) |
                           ba::transformed([this](const CellChainEntry& e)
                           {
                               auto it = cmplx.find(e.index());
                               std::get<0>(it->second)++;      // increment cofaces
                               Index idx = std::get<1>(it->second);
                               return ChainEntry(e.element(), idx);
                           }));
    op++;

    cmplx.emplace(s, CofacesIndex {0,idx++});
}

template<class K, class T, class F>
void
RelativeLZZ<K,T,F>::
remove_both(const Simplex& s)
{
    auto it = cmplx.find(s);

    if (std::get<0>(it->second) > 0)        // cofaces left, keep the simplex
        return;

    //fmt::print("Removing from both: {}\n", s);

    Index idx = std::get<1>(it->second);
    zz.remove_both(idx);
    op++;

    cmplx.erase(s);

    for (const auto& sb : s.boundary())
    {
        auto bit = cmplx.find(sb);
        std::get<0>(bit->second)--;
    }
}

template<class K, class T, class F>
template<class ReportPair>
void
RelativeLZZ<K,T,F>::
operator()(const ReportPair& report_pair)
{
    Vertices vertices;
    for(const Vertex& v : topology.vertices())
        vertices.push_back(ValueVertex(function(v), v));
    std::sort(vertices.begin(), vertices.end());

    for(auto& vval : vertices)
    {
        Value   val;
        Vertex  v;
        std::tie(val,v) = vval;
        //fmt::print("Vertex: {} {}\n", v, val);

        // generate the closures of upper and lower stars
        std::vector<Simplex>    lower_star, upper_star, lower_link, upper_link;
        for (auto& s : topology.closed_star(v))
        {
            bool lower, upper, has_v;
            std::tie(lower,upper,has_v) = lower_upper(s, vval);
            if (lower && !upper)
            {
                if (has_v)
                    lower_star.emplace_back(std::move(s));
                else
                    lower_link.emplace_back(std::move(s));
            }
            else if (!lower && upper)
            {
                if (has_v)
                    upper_star.emplace_back(std::move(s));
                else
                    upper_link.emplace_back(std::move(s));
            } else if (!lower && !upper)
            {
                // s = { v }
                lower_star.push_back(s);
                upper_star.push_back(s);
            }
        }

        // add upper link and star to zz
        std::sort(upper_star.begin(), upper_star.end());        // order by dimension
        std::sort(upper_link.begin(), upper_link.end());        // order by dimension
        for (auto& s : upper_link)
            add_both(s);
        for (auto& s : upper_star)
            add_both(s);

        for (auto it = upper_star.rbegin(); it != upper_star.rend(); ++it)
        {
            Index pair = remove_relative(*it);
            if (pair != zz.unpaired)
            {
                Index pair_op = std::upper_bound(ops.begin(), ops.end(), pair) - ops.begin();
                if (pair_op != ops.size())
                    report_pair(it->dimension(), vertices[pair_op], vval, birth_type_[pair], false);

                birth_type_.erase(pair);
            }
        }

        std::sort(lower_star.begin(), lower_star.end());        // order by dimension
        std::sort(lower_link.begin(), lower_link.end());        // order by dimension
        for (auto& s : lower_star)
        {
            Index pair = add_relative(s);
            if (pair != zz.unpaired)
            {
                Index pair_op = std::upper_bound(ops.begin(), ops.end(), pair) - ops.begin();
                if (pair_op != ops.size())
                    report_pair(s.dimension(), vertices[pair_op], vval, birth_type_[pair], true);

                birth_type_.erase(pair);
            }
        }

        // go through both lower_star and lower_link
        auto sit = lower_star.rbegin();
        auto lit = lower_link.rbegin();
        while(sit != lower_star.rend() || lit != lower_link.rend())
        {
            decltype(sit) it;
            if (sit == lower_star.rend())
                it = lit++;
            else if (lit == lower_link.rend())
                it = sit++;
            else if (sit->dimension() > lit->dimension())
                it = sit++;
            else
                it = lit++;

            remove_both(*it);
        }

        ops.push_back(op);
    }

    //fmt::print("Left alive: {}\n", zz.size_alive());
    assert(zz.size_alive() == 0);
}

template<class K, class T, class F>
typename RelativeLZZ<K,T,F>::Index
RelativeLZZ<K,T,F>::
remove_relative(const Simplex& s)
{
    //fmt::print("Removing from relative: {}\n", s);

    Index idx = std::get<1>(cmplx.find(s)->second);
    op++;
    Index pair = zz.remove(idx);
    if (pair == zz.unpaired)
        birth_type_[idx] = false;

    return pair;
}

template<class K, class T, class F>
typename RelativeLZZ<K,T,F>::Index
RelativeLZZ<K,T,F>::
add_relative(const Simplex& s)
{
    //fmt::print("Adding to relative: {}\n", s);

    Index idx = std::get<1>(cmplx.find(s)->second);
    op++;
    Index pair = zz.add(idx, s.boundary(zz.field()) |
                               ba::transformed([this](const CellChainEntry& e)
                               {
                                   auto it = cmplx.find(e.index());
                                   Index idx = std::get<1>(it->second);
                                   return ChainEntry(e.element(), idx);
                               }));
    if (pair == zz.unpaired)
        birth_type_[idx] = true;

    return pair;
}
