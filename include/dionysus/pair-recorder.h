#ifndef DIONYSUS_PAIR_RECORDER_H
#define DIONYSUS_PAIR_RECORDER_H

namespace dionysus
{

template<class Persistence_>
struct PairRecorder: public Persistence_
{
    typedef             Persistence_                    Persistence;
    typedef             typename Persistence::Index     Index;


    using Persistence::Persistence;

    template<class ChainRange>
    Index               add(const ChainRange& chain)
    {
        Index p = Persistence::add(chain);
        pairs_.push_back(p);
        if (p != unpaired())
            pairs_[p] = pairs_.size() - 1;

        return p;
    }

    Index               pair(Index i) const             { return pairs_[i]; }

    void                resize(size_t s)                { Persistence::resize(s); pairs_.resize(s, unpaired()); }
    size_t              size() const                    { return pairs_.size(); }
    static const Index  unpaired()                      { return Reduction<Index>::unpaired; }

    std::vector<Index>  pairs_;
};

}

#endif
