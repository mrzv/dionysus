#ifndef DIONYSUS_OMNI_FIELD_REDUCTION_H
#define DIONYSUS_OMNI_FIELD_REDUCTION_H

#include <vector>
#include <unordered_map>

#include "fields/q.h"
#include "fields/zp.h"
#include "chain.h"

namespace dionysus
{

template<typename Index_ = unsigned, class Comparison_ = std::less<Index_>>
class OmniFieldPersistence
{
    public:
        using   Index       = Index_;
        using   Q           = ::dionysus::Q<>;
        using   Field       = Q;
        using   Comparison  = Comparison_;

        using   BaseElement = Q::BaseElement;
        using   Zp          = ::dionysus::ZpField<BaseElement>;
        using   Zps         = std::unordered_map<BaseElement, Zp>;

        using   QElement    = Q::Element;
        using   QEntry      = ChainEntry<Q,Index>;
        using   QChain      = std::vector<QEntry>;

        using   ZpElement   = Zp::Element;
        using   ZpEntry     = ChainEntry<Zp, Index>;
        using   ZpChain     = std::vector<ZpEntry>;

        using   QChains     = std::vector<QChain>;
        using   ZpChains    = std::unordered_map<Index, std::unordered_map<BaseElement, ZpChain>>;

        using   QLows       = std::unordered_map<Index, Index>;
        using   ZpLows      = std::unordered_map<Index, std::unordered_map<BaseElement, Index>>;

        using   Factors     = std::vector<BaseElement>;

        const Field&        field() const                       { return q_; }

        void                sort(QChain& c)                     { std::sort(c.begin(), c.end(),
                                                                  [this](const QEntry& e1, const QEntry& e2)
                                                                  { return this->cmp_(e1.index(), e2.index()); }); }

        template<class ChainRange>
        void                add(const ChainRange& chain)        { return add(QChain(std::begin(chain), std::end(chain))); }
        void                add(QChain&& chain);

        void                reserve(size_t s)                   { q_chains_.reserve(s); }

        void                reduce(ZpChain& zp_chain, BaseElement p);
        ZpChain             convert(const QChain& c, const Zp& field);
        bool                special(Index i, BaseElement p) const   { auto it = zp_chains_.find(i); if (it == zp_chains_.end()) return false; if (it->second.find(p) == it->second.end()) return false; return true; }

        const Zp&           zp(BaseElement p)                   { auto it = zps_.find(p); if (it != zps_.end()) return it->second; return zps_.emplace(p, Zp(p)).first->second; }

        static Factors      factor(BaseElement x);

        const QChains&      q_chains() const                    { return q_chains_; }
        const ZpChains&     zp_chains() const                   { return zp_chains_; }

        // This is a bit of a hack; it takes advantage of the fact that zp(p)
        // generates field on-demand and memoizes them. So there is an entry in
        // zps_ only if something special happened over the prime.
        Factors             primes() const                      { Factors result; result.reserve(zps_.size()); for (auto& x : zps_) result.push_back(x.first); return result; }

    private:
        QChains     q_chains_;
        ZpChains    zp_chains_;

        QLows       q_lows_;
        ZpLows      zp_lows_;

        Q           q_;
        Zps         zps_;

        Comparison  cmp_;
};

} // dionysus

#include "omni-field-persistence.hpp"

#endif
