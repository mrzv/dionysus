template<class Chain_>
bool chain_eq(const Chain_& c1, const Chain_& c2)
{
    if (c1.size() != c2.size())
        return false;

    for (size_t i = 0; i < c1.size(); ++i)
    {
        if (c1[i].index() != c2[i].index() || c1[i].element() != c2[i].element())
            return false;
    }
    return true;
}

template<class Chain_>
bool chain_ne(const Chain_& c1, const Chain_& c2)
{
    return !chain_eq(c1, c2);
}

template<class Chain_>
void init_chain(py::module& m, std::string prefix = "")
{
    py::class_<Chain_>(m, (prefix + "Chain").c_str(), "chain of indices (formal sum with coefficients in Zp)")
        .def("__len__",     &PyReducedMatrix::Chain::size,                                  "size of the chain")
        .def("__getitem__", [](const PyReducedMatrix::Chain& c, size_t i) { return c[i]; }, "access the entry at a given index")
        .def("__iter__",    [](const PyReducedMatrix::Chain& c) { return py::make_iterator(c.begin(), c.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the entries of the chain")
        .def("__eq__",      &chain_eq<Chain_>, "equality comparison")
        .def("__ne__",      &chain_ne<Chain_>, "nonequal comparison")
        .def("__repr__",    [](const PyReducedMatrix::Chain& c)
                            {
                                std::ostringstream oss;
                                auto it = c.begin();
                                if (it == c.end())
                                    return oss.str();
                                oss << it->e << '*' << it->i;
                                while (++it != c.end())
                                    oss << " + " << it->e << '*' << it->i;
                                return oss.str();
                            })
    ;

    using Entry = typename Chain_::value_type;
    py::class_<Entry>(m, (prefix + "ChainEntry").c_str(), "(coefficient, index) entry in a chain)")
        .def_property_readonly("element",   [](const Entry& x) { return x.element(); },
                                            "coefficient of the chain element")
        .def_readonly("index",              &Entry::i,   "index of the chain element")
        .def("__repr__",                    [](const Entry& e)
                                            { std::ostringstream oss; oss << e.e << '*' << e.i; return oss.str(); })
    ;
}
