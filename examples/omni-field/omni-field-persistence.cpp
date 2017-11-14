#include <iostream>
#include <vector>

#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
#include <dionysus/omni-field-persistence.h>
#include <dionysus/diagram.h>

namespace d = dionysus;

#include <format.h>

using Simplex       = d::Simplex<>;
using Filtration    = d::Filtration<Simplex>;
using Persistence   = d::OmniFieldPersistence<>;

int main()
{
    // Klein bottle
    Filtration filtration
    {
      {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8},
      {0,1}, {1,2}, {2,0}, {0,3}, {3,4}, {4,0},
      {1,5}, {5,6}, {6,2}, {2,7}, {7,8}, {8,1},
      {3,5}, {5,7}, {7,3}, {4,6}, {6,8}, {8,4},
      {0,5}, {1,7}, {2,3}, {3,6}, {5,8}, {7,4},
      {4,2}, {6,1}, {8,0},
      {0,3,5}, {0,1,5}, {1,5,7}, {1,2,7}, {2,7,3}, {2,3,0},
      {3,4,6}, {3,5,6}, {5,6,8}, {5,7,8}, {7,8,4}, {7,3,4},
      {4,0,2}, {4,6,2}, {6,2,1}, {6,8,1}, {8,1,0}, {8,4,0}
    };

    fmt::print("Boundary matrix over Q\n");
    d::Q<> q;
    for (auto& s : filtration)
    {
        fmt::print("{} at {}\n", s, filtration.index(s));
        for (auto sb : s.boundary(q))
            fmt::print("   {} * {} at {}\n", sb.element(), sb.index(), filtration.index(sb.index()));
    }

    Persistence     persistence;
    for(auto& s : filtration)
    {
        using SimplexChainEntry = d::ChainEntry<Persistence::Field, Simplex>;
        using ChainEntry        = d::ChainEntry<Persistence::Field, Persistence::Index>;
        persistence.add(s.boundary(persistence.field()) |
                                                 ba::transformed([&filtration](const SimplexChainEntry& e)
                                                 { return ChainEntry(e.element(), filtration.index(e.index())); }));
    }
    fmt::print("Reduction finished\n");

    fmt::print("Special primes:");
    for (auto x : persistence.primes())
        fmt::print(" {}", x);
    fmt::print("\n");

    unsigned i = 0;
    fmt::print("Q chains finished\n");
    for (auto& c : persistence.q_chains())
    {
        fmt::print("{}: ", i);
        for (auto& ce : c)
            fmt::print(" + {} * {}", ce.element(), ce.index());
        fmt::print("\n");
        ++i;
    }

    fmt::print("Zp chains finished\n");
    for (auto& x : persistence.zp_chains())
    {
        unsigned i = x.first;
        fmt::print("{}:\n", i);
        for (auto& ec : x.second)
        {
            auto& e = ec.first;
            auto& c = ec.second;

            fmt::print("  mod {}:", e);

            for (auto& ce : c)
                fmt::print(" + {} * {}", ce.element(), ce.index());
        }
        fmt::print("\n");
    }

    auto primes = persistence.primes();
    primes.emplace(primes.begin(), 1);
    for (auto& p : primes)
    {
        if (p == 1)
            fmt::print("Over Z_p (for all p, except those specified explicitly)\n");
        else
            fmt::print("Over Z_{}:\n", p);
        auto diagrams = init_diagrams(prime_adapter(persistence, p), filtration,
                                      [&](const Simplex& s) -> float  { return filtration.index(s); },        // inefficient, but works
                                      [](Persistence::Index i)        { return i; });
        i = 0;
        for (auto& dgm : diagrams)
        {
            fmt::print("  Dimension {}:\n", i++);
            for (auto& pt : dgm)
                fmt::print("    {} {}\n", pt.birth(), pt.death());
        }
    }
}
