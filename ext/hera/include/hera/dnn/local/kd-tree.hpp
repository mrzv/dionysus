#include <boost/range/counting_range.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/range.hpp>

#include <queue>
#include <stack>

#include "../parallel/tbb.h" // for task_group

template<class T>
hera::ws::dnn::KDTree<T>::
KDTree(const Traits& traits, HandleContainer&& handles, double _wassersteinPower):
    traits_(traits), tree_(std::move(handles)), wassersteinPower(_wassersteinPower)
{ assert(wassersteinPower >= 1.0); init(); }

template<class T>
template<class Range>
hera::ws::dnn::KDTree<T>::
KDTree(const Traits& traits, const Range& range, double _wassersteinPower):
    traits_(traits), wassersteinPower(_wassersteinPower)
{
    assert( wassersteinPower >= 1.0);
    init(range);
}

template<class T>
template<class Range>
void
hera::ws::dnn::KDTree<T>::
init(const Range& range)
{
    size_t sz = std::distance(std::begin(range), std::end(range));
    tree_.reserve(sz);
    weights_.resize(sz, 0);
    subtree_weights_.resize(sz, 0);
    for (PointHandle h : range)
        tree_.push_back(h);
    init();
}

template<class T>
void
hera::ws::dnn::KDTree<T>::
init()
{
    if (tree_.empty())
        return;

#if defined(TBB)
    hera::dnn::task_group g;
    g.run(OrderTree(tree_.begin(), tree_.end(), 0, traits()));
    g.wait();
#else
    OrderTree(tree_.begin(), tree_.end(), 0, traits()).serial();
#endif

    for (size_t i = 0; i < tree_.size(); ++i)
        indices_[tree_[i]] = i;
}

template<class T>
struct
hera::ws::dnn::KDTree<T>::OrderTree
{
                OrderTree(HCIterator b_, HCIterator e_, size_t i_, const Traits& traits_):
                    b(b_), e(e_), i(i_), traits(traits_)            {}

    void        operator()() const
    {
        if (e - b < 1000)
        {
            serial();
            return;
        }

        HCIterator m = b + (e - b)/2;
        CoordinateComparison cmp(i, traits);
        std::nth_element(b,m,e, cmp);
        size_t next_i = (i + 1) % traits.dimension();

        hera::dnn::task_group g;
        if (b < m - 1)  g.run(OrderTree(b,   m, next_i, traits));
        if (e > m + 2)  g.run(OrderTree(m+1, e, next_i, traits));
        g.wait();
    }

    void        serial() const
    {
        std::queue<KDTreeNode> q;
        q.push(KDTreeNode(b,e,i));
        while (!q.empty())
        {
            HCIterator b, e; size_t i;
            std::tie(b,e,i) = q.front();
            q.pop();
            HCIterator m = b + (e - b)/2;

            CoordinateComparison cmp(i, traits);
            std::nth_element(b,m,e, cmp);
            size_t next_i = (i + 1) % traits.dimension();

            // Replace with a size condition instead?
            if (m - b > 1)  q.push(KDTreeNode(b,   m, next_i));
            if (e - m > 2)  q.push(KDTreeNode(m+1, e, next_i));
        }
    }

    HCIterator      b, e;
    size_t          i;
    const Traits&   traits;
};

template<class T>
template<class ResultsFunctor>
void
hera::ws::dnn::KDTree<T>::
search(PointHandle q, ResultsFunctor& rf) const
{
    typedef         typename HandleContainer::const_iterator        HCIterator;
    typedef         std::tuple<HCIterator, HCIterator, size_t>      KDTreeNode;

    if (tree_.empty())
        return;

    DistanceType    D  = std::numeric_limits<DistanceType>::max();

    // TODO: use tbb::scalable_allocator for the queue
    std::queue<KDTreeNode>  nodes;

    nodes.push(KDTreeNode(tree_.begin(), tree_.end(), 0));

    while (!nodes.empty())
    {
        HCIterator b, e; size_t i;
        std::tie(b,e,i) = nodes.front();
        nodes.pop();

        CoordinateComparison cmp(i, traits());
        i = (i + 1) % traits().dimension();

        HCIterator   m    = b + (e - b)/2;

        DistanceType dist = (wassersteinPower == 1.0) ?  traits().distance(q, *m) + weights_[m - tree_.begin()] : std::pow(traits().distance(q, *m), wassersteinPower) + weights_[m - tree_.begin()];


        D = rf(*m, dist);

        // we are really searching w.r.t L_\infty ball; could prune better with an L_2 ball
        Coordinate diff = cmp.diff(q, *m);     // diff returns signed distance

        DistanceType diffToWasserPower = (wassersteinPower == 1.0) ? diff : ((diff > 0 ? 1.0 : -1.0) * std::pow(fabs(diff), wassersteinPower));

        size_t       lm   = m + 1 + (e - (m+1))/2 - tree_.begin();
        if (e > m + 1 && diffToWasserPower - subtree_weights_[lm] >= -D) {
            nodes.push(KDTreeNode(m+1, e, i));
        }

        size_t       rm   = b + (m - b) / 2 - tree_.begin();
        if (b < m && diffToWasserPower + subtree_weights_[rm] <= D) {
            nodes.push(KDTreeNode(b,   m, i));
        }
    }
}

template<class T>
void
hera::ws::dnn::KDTree<T>::
adjust_weights(DistanceType delta)
{
    for(auto& w : weights_)
        w -= delta;

    for(auto& sw : subtree_weights_)
        sw -= delta;
}


template<class T>
void
hera::ws::dnn::KDTree<T>::
change_weight(PointHandle p, DistanceType w)
{
    size_t idx = indices_[p];

    if ( weights_[idx] == w ) {
        return;
    }

    bool weight_increases = ( weights_[idx] < w );
    weights_[idx] = w;

    typedef     std::tuple<HCIterator, HCIterator>      KDTreeNode;

    // find the path down the tree to this node
    // not an ideal strategy, but // it's not clear how to move up from the node in general
    std::stack<KDTreeNode>  s;
    s.push(KDTreeNode(tree_.begin(),tree_.end()));

    do
    {
        HCIterator b,e;
        std::tie(b,e) = s.top();

        size_t im = b + (e - b)/2 - tree_.begin();

        if (idx == im)
            break;
        else if (idx < im)
            s.push(KDTreeNode(b, tree_.begin() + im));
        else    // idx > im
            s.push(KDTreeNode(tree_.begin() + im + 1, e));
    } while(1);

    // update subtree_weights_ on the path to the root
    DistanceType min_w = w;
    while (!s.empty())
    {
        HCIterator b,e;
        std::tie(b,e) = s.top();
        HCIterator m  = b + (e - b)/2;
        size_t     im = m - tree_.begin();
        s.pop();


        // left and right children
        if (b < m)
        {
            size_t lm = b + (m - b)/2 - tree_.begin();
            if (subtree_weights_[lm] < min_w)
                min_w = subtree_weights_[lm];
        }

        if (e > m + 1)
        {
            size_t rm = m + 1 + (e - (m+1))/2 - tree_.begin();
            if (subtree_weights_[rm] < min_w)
                min_w = subtree_weights_[rm];
        }

        if (weights_[im] < min_w) {
            min_w = weights_[im];
        }

        if (weight_increases) {

            if (subtree_weights_[im] < min_w )   // increase weight
                subtree_weights_[im] = min_w;
            else
                break;

        } else {

            if (subtree_weights_[im] > min_w )   // decrease weight
                subtree_weights_[im] = min_w;
            else
                break;

        }
    }
}

template<class T>
typename hera::ws::dnn::KDTree<T>::HandleDistance
hera::ws::dnn::KDTree<T>::
find(PointHandle q) const
{
    hera::dnn::NNRecord<HandleDistance> nn;
    search(q, nn);
    return nn.result;
}

template<class T>
typename hera::ws::dnn::KDTree<T>::Result
hera::ws::dnn::KDTree<T>::
findR(PointHandle q, DistanceType r) const
{
    hera::dnn::rNNRecord<HandleDistance> rnn(r);
    search(q, rnn);
    std::sort(rnn.result.begin(), rnn.result.end());
    return rnn.result;
}

template<class T>
typename hera::ws::dnn::KDTree<T>::Result
hera::ws::dnn::KDTree<T>::
findK(PointHandle q, size_t k) const
{
    hera::dnn::kNNRecord<HandleDistance> knn(k);
    search(q, knn);
    std::sort(knn.result.begin(), knn.result.end());
    return knn.result;
}


template<class T>
struct hera::ws::dnn::KDTree<T>::CoordinateComparison
{
                CoordinateComparison(size_t i, const Traits& traits):
                    i_(i), traits_(traits)                              {}

    bool        operator()(PointHandle p1, PointHandle p2) const        { return coordinate(p1) < coordinate(p2); }
    Coordinate  diff(PointHandle p1, PointHandle p2) const              { return coordinate(p1) - coordinate(p2); }

    Coordinate  coordinate(PointHandle p) const                         { return traits_.coordinate(p, i_); }
    size_t      axis() const                                            { return i_; }

    private:
        size_t          i_;
        const Traits&   traits_;
};

template<class T>
void
hera::ws::dnn::KDTree<T>::
printWeights(void)
{
#ifndef FOR_R_TDA
    std::cout << "weights_:" << std::endl;
    for(const auto ph : indices_) {
        std::cout << "idx = " << ph.second << ": (" << (ph.first)->at(0) << ", " << (ph.first)->at(1) << ") weight = " << weights_[ph.second] << std::endl;
    }
    std::cout << "subtree_weights_:" << std::endl;
    for(size_t idx = 0; idx < subtree_weights_.size(); ++idx) {
        std::cout << idx << " : " << subtree_weights_[idx] << std::endl;
    }
#endif
}

template<class T>
hera::bt::dnn::KDTree<T>::KDTree(const Traits& traits, HandleContainer&& handles):
    traits_(traits),
    tree_(std::move(handles)),
    delete_flags_(handles.size(), static_cast<char>(0) ),
    subtree_n_elems(handles.size(), static_cast<size_t>(0)),
    num_points_(handles.size())
{
    init();
}

template<class T>
template<class Range>
hera::bt::dnn::KDTree<T>::KDTree(const Traits& traits, const Range& range):
    traits_(traits)
{
    init(range);
}

template<class T>
template<class Range>
void hera::bt::dnn::KDTree<T>::init(const Range& range)
{
    size_t sz = std::distance(std::begin(range), std::end(range));
    subtree_n_elems = std::vector<int>(sz, 0);
    delete_flags_ = std::vector<char>(sz, 0);
    num_points_ = sz;
    tree_.reserve(sz);
    for (PointHandle h : range)
        tree_.push_back(h);
    parents_.resize(sz, -1);
    init();
}

template<class T>
void hera::bt::dnn::KDTree<T>::init()
{
    if (tree_.empty())
        return;

#if defined(TBB)
    hera::dnn::task_group g;
    g.run(OrderTree(this, tree_.begin(), tree_.end(), -1, 0, traits()));
    g.wait();
#else
    OrderTree(this, tree_.begin(), tree_.end(), -1, 0, traits()).serial();
#endif

    for (size_t i = 0; i < tree_.size(); ++i)
        indices_[tree_[i]] = i;
    init_n_elems();
}

template<class T>
struct
hera::bt::dnn::KDTree<T>::OrderTree
{
                OrderTree(KDTree* tree_, HCIterator b_, HCIterator e_, ssize_t p_, size_t i_, const Traits& traits_):
                    tree(tree_), b(b_), e(e_), p(p_), i(i_), traits(traits_)            {}

    void        operator()() const
    {
        if (e - b < 1000)
        {
            serial();
            return;
        }

        HCIterator m  = b + (e - b)/2;
        ssize_t    im = m - tree->tree_.begin();
        tree->parents_[im]  = p;

        CoordinateComparison cmp(i, traits);
        std::nth_element(b,m,e, cmp);
        size_t next_i = (i + 1) % traits.dimension();

        hera::dnn::task_group g;
        if (b < m - 1)  g.run(OrderTree(tree, b,   m, im, next_i, traits));
        if (e > m + 2)  g.run(OrderTree(tree, m+1, e, im, next_i, traits));
        g.wait();
    }

    void        serial() const
    {
        std::queue<KDTreeNode> q;
        q.push(KDTreeNode(b,e,p,i));
        while (!q.empty())
        {
            HCIterator b, e; ssize_t p; size_t i;
            std::tie(b,e,p,i) = q.front();
            q.pop();
            HCIterator m  = b + (e - b)/2;
            ssize_t    im = m - tree->tree_.begin();
            tree->parents_[im]  = p;

            CoordinateComparison cmp(i, traits);
            std::nth_element(b,m,e, cmp);
            size_t next_i = (i + 1) % traits.dimension();

            // Replace with a size condition instead?
            if (m - b > 1)
                q.push(KDTreeNode(b,   m, im, next_i));
            else if (b < m)
                tree->parents_[im - 1] = im;
            if (e - m > 2)
                q.push(KDTreeNode(m+1, e, im, next_i));
            else if (e > m + 1)
                tree->parents_[im + 1] = im;
        }
    }

    KDTree*         tree;
    HCIterator      b, e;
    ssize_t         p;
    size_t          i;
    const Traits&   traits;
};

template<class T>
void hera::bt::dnn::KDTree<T>::update_n_elems(ssize_t idx, const int delta)
// add delta to the number of points in node idx and update subtree_n_elems
// for all parents of the node idx
{
    //std::cout << "subtree_n_elems.size = " << subtree_n_elems.size() << std::endl;
    // update the node itself
    while (idx != -1)
    {
        //std::cout << idx << std::endl;
        subtree_n_elems[idx] += delta;
        idx = parents_[idx];
    }
}

template<class T>
void hera::bt::dnn::KDTree<T>::increase_n_elems(const ssize_t idx)
{
    update_n_elems(idx, static_cast<ssize_t>(1));
}

template<class T>
void hera::bt::dnn::KDTree<T>::decrease_n_elems(const ssize_t idx)
{
    update_n_elems(idx, static_cast<ssize_t>(-1));
}

template<class T>
void hera::bt::dnn::KDTree<T>::init_n_elems()
{
    for(size_t idx = 0; idx < tree_.size(); ++idx) {
        increase_n_elems(idx);
    }
}


template<class T>
template<class ResultsFunctor>
void hera::bt::dnn::KDTree<T>::search(PointHandle q, ResultsFunctor& rf) const
{
    typedef         typename HandleContainer::const_iterator        HCIterator;
    typedef         std::tuple<HCIterator, HCIterator, size_t>      KDTreeNode;

    if (tree_.empty())
        return;

    DistanceType    D  = std::numeric_limits<DistanceType>::infinity();

    // TODO: use tbb::scalable_allocator for the queue
    std::queue<KDTreeNode>  nodes;

    nodes.push(KDTreeNode(tree_.begin(), tree_.end(), 0));

    //std::cout << "started kdtree::search" << std::endl;

    while (!nodes.empty())
    {
        HCIterator b, e; size_t i;
        std::tie(b,e,i) = nodes.front();
        nodes.pop();

        CoordinateComparison cmp(i, traits());
        i = (i + 1) % traits().dimension();

        HCIterator   m    = b + (e - b)/2;
        size_t m_idx = m - tree_.begin();
        // ignore deleted points
        if ( delete_flags_[m_idx] == 0 ) {
            DistanceType dist = traits().distance(q, *m);
            // + weights_[m - tree_.begin()];
            //std::cout << "Supplied to functor: m : ";
            //std::cout << "(" << (*(*m))[0] << ", " << (*(*m))[1] << ")";
            //std::cout << " and q : ";
            //std::cout << "(" << (*q)[0] << ", " << (*q)[1] << ")" << std::endl;
            //std::cout << "dist^q + weight = " << dist << std::endl;
            //std::cout << "weight = " << weights_[m - tree_.begin()] << std::endl;
            //std::cout << "dist = " << traits().distance(q, *m) << std::endl;
            //std::cout << "dist^q = " << pow(traits().distance(q, *m), wassersteinPower) << std::endl;

            D = rf(*m, dist);
        }
        // we are really searching w.r.t L_\infty ball; could prune better with an L_2 ball
        Coordinate diff = cmp.diff(q, *m);     // diff returns signed distance
        DistanceType diffToWasserPower = (diff > 0 ? 1.0 : -1.0) * fabs(diff);

        size_t       lm   = m + 1 + (e - (m+1))/2 - tree_.begin();
        if ( e > m + 1 and subtree_n_elems[lm] > 0 ) {
            if (e > m + 1 && diffToWasserPower  >= -D) {
                nodes.push(KDTreeNode(m+1, e, i));
            }
        }

        size_t       rm   = b + (m - b) / 2 - tree_.begin();
        if ( subtree_n_elems[rm] > 0 ) {
            if (b < m && diffToWasserPower  <= D) {
                nodes.push(KDTreeNode(b,   m, i));
            }
        }
    }
    //std::cout << "exited kdtree::search" << std::endl;
}

template<class T>
typename hera::bt::dnn::KDTree<T>::HandleDistance hera::bt::dnn::KDTree<T>::find(PointHandle q) const
{
    hera::dnn::NNRecord<HandleDistance> nn;
    search(q, nn);
    return nn.result;
}

template<class T>
typename hera::bt::dnn::KDTree<T>::Result hera::bt::dnn::KDTree<T>::findR(PointHandle q, DistanceType r) const
{
    hera::dnn::rNNRecord<HandleDistance> rnn(r);
    search(q, rnn);
    //std::sort(rnn.result.begin(), rnn.result.end());
    return rnn.result;
}

template<class T>
typename hera::bt::dnn::KDTree<T>::Result hera::bt::dnn::KDTree<T>::findFirstR(PointHandle q, DistanceType r) const
{
    hera::dnn::firstrNNRecord<HandleDistance> rnn(r);
    search(q, rnn);
    return rnn.result;
}

template<class T>
typename hera::bt::dnn::KDTree<T>::Result hera::bt::dnn::KDTree<T>::findK(PointHandle q, size_t k) const
{
    hera::dnn::kNNRecord<HandleDistance> knn(k);
    search(q, knn);
    // do we need this???
    std::sort(knn.result.begin(), knn.result.end());
    return knn.result;
}

template<class T>
struct hera::bt::dnn::KDTree<T>::CoordinateComparison
{
                CoordinateComparison(size_t i, const Traits& traits):
                    i_(i), traits_(traits)                              {}

    bool        operator()(PointHandle p1, PointHandle p2) const        { return coordinate(p1) < coordinate(p2); }
    Coordinate  diff(PointHandle p1, PointHandle p2) const              { return coordinate(p1) - coordinate(p2); }

    Coordinate  coordinate(PointHandle p) const                         { return traits_.coordinate(p, i_); }
    size_t      axis() const                                            { return i_; }

    private:
        size_t          i_;
        const Traits&   traits_;
};

template<class T>
void hera::bt::dnn::KDTree<T>::delete_point(const size_t idx)
{
    // prevent double deletion
    assert(delete_flags_[idx] == 0);
    delete_flags_[idx] = 1;
    decrease_n_elems(idx);
    --num_points_;
}

template<class T>
void hera::bt::dnn::KDTree<T>::delete_point(PointHandle p)
{
    delete_point(indices_[p]);
}

