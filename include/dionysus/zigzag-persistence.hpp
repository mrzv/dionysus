template<class F, class I, class C>
template<class ChainRange>
typename dionysus::ZigzagPersistence<F,I,C>::Index
dionysus::ZigzagPersistence<F,I,C>::
add(const ChainRange& chain_)
{
    Index op = operations++;

    IndexChain cycles;      // chain_ -> Z*cycles
    Column     z_remainder = Z.reduce(chain_, cycles);
    assert(z_remainder.empty());

    IndexChain  boundaries;
    DequeColumn b_remainder = B.reduce(cycles, boundaries);

    // add up columns of D indexed by boundaries
    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(row(e1), row(e2)); };
    Column      chain;
    for (auto& x : boundaries)
        Chain<Column>::addto(chain, x.element(), C.col(x.index()), field(), row_cmp);
    chain.push_back(Entry(field().id(), IndexPair(cell_indices++,0)));

    if (b_remainder.empty())        // birth
    {
        Index z_col = z_indicies_last++;
        Z.set(z_col, std::move(chain));
        birth_index[z_col] = op;
        return unpaired;
    }
    else                            // death
    {
        Index b_col = b_indices++;
        Index pair  = row(b_remainder.back());
        B.set(b_col, std::move(b_remainder));
        C.set(b_col, std::move(chain));
        return birth_index[pair];
    }
}

template<class F, class I, class C>
typename dionysus::ZigzagPersistence<F,I,C>::Index
dionysus::ZigzagPersistence<F,I,C>::
remove(Index cell)
{
    Index   op    = operations++;

    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(row(e1), row(e2)); };
    typedef     typename DequeColumn::value_type        DequeEntry;
    auto        b_row_cmp = [this](const DequeEntry& e1, const DequeEntry& e2)
                            { return this->cmp()(row(e1), row(e2)); };

    Row&    z_row = Z.row(cell);
    if (z_row.empty())              // birth
    {
        Row&    c_row   = C.row(cell);
        // c_row.front() may not be the first column in order, but that doesn't really matter, does it? (TODO)
        auto&   c_front = c_row.front();

        Index   j     = col(c_front);
        Index   l     = row(B.col(j).back());

        //std::cout << j << ' ' << l << std::endl;

        // cycle = ZB[j] = DC[j]
        Column  cycle;
        for (auto& x : B.col(j))
            Chain<Column>::addto(cycle, x.element(), Z.col(row(x)), field(), row_cmp);

        //std::cout << "Cycle:" << std::endl;
        //for (auto& x : cycle)
        //    std::cout << x.element() << ' ' << row(x) << std::endl;

        // 1: prepend the cycle
        Index   znew        = z_indicies_first--;
        Index   oth         = Z.set(znew, std::move(cycle));        // oth records our collision (used in step 6)
        birth_index[znew]   = op;

        //std::cout << "znew oth: " << znew << ' ' << oth << std::endl;
        //std::cout << "oth column:" << std::endl;
        //for (auto& x : Z.col(oth))
        //    std::cout << x.element() << ' ' << row(x) << std::endl;

        // 2: prepend the row to B
        FieldElement    m = field().neg(field().inv(c_front.element()));        // m = -1/c
        B.prepend_row(znew, m, c_row);
        //std::cout << "Prepended row with multiplier: " << m << std::endl;

        // 3: subtract C[j] from every C[k]
        Column&         Cj    = C.col(j);
        const DequeRow& b_row = B.row(znew);        // use the copy of c_row in B, since c_row will be modified in the following loop

        for (auto it = ++b_row.begin(); it != b_row.end(); ++it)
        {
            Index c = col(*it);
            Chain<Column>::addto(C.col(c), it->element(), Cj, field(), row_cmp);    // using it->element() since b_row = m*c_row
            C.fix(c);                                                               // old elements got removed via auto_unlink_hook
            // we don't need lows in C, so not updating them
        }
        //std::cout << "Done with step 3" << std::endl;

        // 4: subtract B[j] from every B[k] that has l
        //    (we don't need to update C because ZB[j] = 0 after step 2)
        DequeColumn&    Bj      = B.col(j);
        FieldElement    bm      = field().neg(field().inv(Bj.back().element()));    // bm = -1/B[l,j]
        IndexChain      Bl_row; // make a copy of Bl_row, since it will be changing
        for (auto& x : B.row(l))
        {
            if (col(x) == j)
                continue;
            Bl_row.emplace_back(x.element(), col(x));
        }
        for (auto& x : Bl_row)
        {
            Index c = x.index();
            Chain<DequeColumn>::addto(B.col(c), field().mul(bm, x.element()), Bj, field(), b_row_cmp);
            B.fix(c);                                                               // old elements got removed via auto_unlink_hook
            // l cannot be the low in c, so no need to update lows
        }
        //std::cout << "Done with step 4" << std::endl;

        // 5: drop row l and column j from B; drop column l from Z; drop column j from C
        B.drop_col(j);
        B.drop_row(l);
        Z.drop_col(l);
        C.drop_col(j);
        C.drop_row(cell);
        Z.drop_row(cell);
        //std::cout << "Done with step 5" << std::endl;
        if (oth == l)       // we just dropped our collision in Z
            oth = znew;

        // 6: reduce Z
        std::unordered_map<Index, DequeColumn>  b_changes;  // the columns to add in B to apply row changes
        Index cur = znew;
        while (oth != cur)
        {
            Column& cur_col = Z.col(cur);
            Column& oth_col = Z.col(oth);
            std::cout << "--- " << cur_col.size() << ' ' << oth_col.size() << std::endl;
            FieldElement m1 = cur_col.back().element();
            FieldElement m2 = oth_col.back().element();
            FieldElement m2_div_m1 = field().div(m2, m1);
            Chain<Column>::addto(oth_col, field().neg(m2_div_m1), cur_col, field(), row_cmp);
            Z.fix(oth, oth_col);

            // record the changes we need to make in B
            for (auto& x : this->B.row(oth))
                b_changes[col(x)].emplace_front(field().neg(m2_div_m1), cur, col(x));

            cur = oth;
            Index low = row(oth_col.back());
            if (Z.is_low(low))
                oth = Z.low(low);
            Z.update_low(cur);
        }
        // apply changes in B (the complexity here could get ugly)
        for (auto& bx : b_changes)
        {
            Chain<DequeColumn>::addto(B.col(bx.first), field().id(), bx.second, field(), b_row_cmp);
            B.fix(bx.first);
            // no need to update low (additions from bottom up)
        }
        //std::cout << "Done with step 6" << std::endl;

        return unpaired;
    }
    else                            // death
    {
        // 1: change basis to clear z_row
        z_row.sort([this](const Entry& e1, const Entry& e2)
                   { return this->cmp()(col(e1), col(e2)); });      // this adds a log factor, but it makes life easier
        Index        j = col(z_row.front());
        FieldElement e = z_row.front().element();

        // figure out the columns we use for reduction
        typedef     typename Row::reverse_iterator      RowRIterator;
        std::vector<RowRIterator>                       reducers;
        reducers.push_back(z_row.rbegin());
        for (RowRIterator it = ++z_row.rbegin(); it != z_row.rend(); ++it)
        {
            Index c = col(*it);
            if (cmp()(Z.low(col(*reducers.back())), Z.low(c)))
                reducers.push_back(it);
        }
        if (reducers.back() != --z_row.rend())
            reducers.push_back(--z_row.rend());


        std::unordered_map<Index, DequeColumn>  b_changes;  // the columns to add to columns of B with non-zero entries in row l
        auto add_in_z = [this,&b_changes,&row_cmp](Index to, Index from, FieldElement m, FieldElement e)
                        {
                            FieldElement    mult = this->field().mul(m, e);
                            Chain<Column>::addto(Z.col(to), mult, Z.col(from), this->field(), row_cmp);
                            this->Z.fix(to);       // NB: rows will be linked in the back, so the iterators are Ok
                            this->Z.update_low(to);

                            for (auto& x : this->B.row(to))
                                b_changes[col(x)].emplace_front(this->field().neg(mult), from, col(x));
                        };
        for (size_t i = 1; i < reducers.size(); ++i)
        {
            auto rit = reducers[i];
            FieldElement m = field().neg(field().inv(rit->element()));

            // NB: since we are moving it from right to left, z_row will be losing elements at the end,
            //     so the iterators should remain intact
            for (auto it  = std::prev(rit); it != reducers[i-1]; --it)
                add_in_z(col(*it), col(*rit), m, it->element());

            add_in_z(col(*reducers[i-1]), col(*rit), m, reducers[i-1]->element());
        }

        // apply changes in b (the complexity here could get ugly)
        for (auto& bx : b_changes)
        {
            Chain<DequeColumn>::addto(B.col(bx.first), field().id(), bx.second, field(), b_row_cmp);
            B.fix(bx.first);
            // no need to update low (additions from bottom up)
        }

        // 2: subtract cycle from every chain in C
        const Column& Zj = Z.col(j);

        IndexChain Ccols;       // save the columns in C, we'll be modifying C.row(cell)
        for (auto& x : C.row(cell))
            Ccols.emplace_back(x.element(), col(x));

        for (auto& x : Ccols)
        {
            Index           c = x.index();
            FieldElement    m = field().neg(field().mul(x.element(), field().inv(e)));      // m = -C[k][cell]/Z[j][cell]
            Chain<Column>::addto(C.col(c), m, Zj, field(), row_cmp);
            C.fix(c);
            // we don't care about lows in C, so don't update them
        }

        // 3: drop
        Z.drop_col(col(z_row.front()));
        Z.drop_row(cell);
        C.drop_row(cell);
        B.drop_row(j);

        Index birth = birth_index[j];
        birth_index.erase(j);

        return birth;
    }
}
