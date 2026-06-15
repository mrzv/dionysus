#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/linked-multi-filtration.h>
#include <dionysus/multi-filtration.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/row-reduction.h>
#include <dionysus/simplex.h>
#include <dionysus/standard-reduction.h>

namespace
{

using Simplex = dionysus::Simplex<int, float>;
using Filtration = dionysus::Filtration<Simplex>;
using CheckedFiltration = dionysus::Filtration<Simplex,
                                               boost::multi_index::hashed_unique<boost::multi_index::identity<Simplex>>,
                                               true>;
using MultiFiltration = dionysus::MultiFiltration<Simplex, true>;
using LinkedMultiFiltration = dionysus::LinkedMultiFiltration<Simplex, true>;
using Field = dionysus::ZpField<>;
using Persistence = dionysus::OrdinaryPersistence<Field>;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

std::string str(const Simplex& simplex)
{
    std::ostringstream out;
    out << simplex;
    return out.str();
}

std::vector<std::string> boundary_strings(const Simplex& simplex)
{
    std::vector<std::string> faces;
    for (const auto& face : simplex.boundary())
        faces.push_back(str(face));
    return faces;
}

void test_simplex_empty_and_sequence_behavior()
{
    const Simplex empty;
    require(empty.dimension() == -1, "empty simplex dimension should be -1");
    require(empty.size() == 0, "empty simplex size should be zero");
    require(str(empty) == "<>", "empty simplex should stream as <>");
    require(std::vector<int>(empty.begin(), empty.end()).empty(), "empty simplex range should be empty");

    Simplex triangle({2, 0, 1}, 7.0f);
    require(triangle.dimension() == 2, "triangle dimension should be 2");
    require(triangle.size() == 3, "triangle size should be 3");
    require(triangle[0] == 0 && triangle[1] == 1 && triangle[2] == 2,
            "simplex vertices should be sorted");
    require(triangle.contains(1), "simplex should contain existing vertex");
    require(!triangle.contains(3), "simplex should not contain missing vertex");
    require(str(triangle) == "<0,1,2>", "triangle should stream sorted vertices");
}

void test_simplex_boundary_join_and_identity()
{
    Simplex triangle({0, 1, 2}, 7.0f);
    require(boundary_strings(triangle) == std::vector<std::string>({"<1,2>", "<0,2>", "<0,1>"}),
            "triangle boundary should preserve alternating face order");

    Simplex vertex({1}, 2.0f);
    Simplex edge = vertex.join(0);
    require(str(vertex) == "<1>", "join should not mutate the original simplex");
    require(str(edge) == "<0,1>", "join should return a sorted simplex");
    require(edge.data() == 2.0f, "join should preserve simplex data");

    Simplex first({0}, 0.0f);
    Simplex second({0}, 1.0f);
    require(first == second, "simplex identity should ignore data");
    first.data() = 3.0f;
    require(first == second, "mutating data should not change simplex identity");
}

void test_filtration_unique_lookup_and_rearrange()
{
    Filtration filtration({Simplex({0, 1}, 2.0f), Simplex({0}, 0.0f), Simplex({1}, 1.0f)});
    require(filtration.size() == 3, "filtration should contain three unique cells");
    require(filtration.contains(Simplex({0})), "filtration should contain vertex 0");
    require(!filtration.contains(Simplex({2})), "filtration should not contain vertex 2");

    filtration.sort();
    require(str(filtration[0]) == "<0>", "sorted filtration should place vertex 0 first");
    require(str(filtration[1]) == "<1>", "sorted filtration should place vertex 1 second");
    require(str(filtration[2]) == "<0,1>", "sorted filtration should place edge last");
    require(filtration.index(Simplex({0}), 0) == 0, "filtration index should find vertex 0");
    require(filtration.index(Simplex({1}), 0) == 1, "filtration index should find vertex 1");

    filtration.rearrange({2, 0, 1});
    require(str(filtration[0]) == "<0,1>", "rearranged filtration should move edge first");
    require(filtration.index(Simplex({0, 1}), 0) == 0, "rearranged edge index should update");
}

void test_checked_filtration_rejects_missing_cell()
{
    CheckedFiltration filtration({Simplex({0})});
    bool threw = false;
    try
    {
        (void) filtration.index(Simplex({1}), 0);
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }
    require(threw, "checked filtration should reject missing cells");
}

void test_multi_filtration_duplicate_lookup()
{
    MultiFiltration filtration({Simplex({0}, 0.0f), Simplex({0}, 1.0f), Simplex({1}, 0.0f), Simplex({0, 1}, 2.0f)});
    require(filtration.size() == 4, "multi-filtration should preserve duplicate simplices");
    require(filtration.index(Simplex({0}, 99.0f), 1) == 1,
            "duplicate lookup should choose the latest matching simplex before the bound");

    std::vector<size_t> boundary_indices;
    size_t edge_index = 3;
    for (const auto& face : filtration[edge_index].boundary())
        boundary_indices.push_back(filtration.index(face, edge_index));

    require(boundary_indices == std::vector<size_t>({2, 1}),
            "multi-filtration boundary lookup should use latest duplicate before the coface");
}

void test_linked_multi_filtration_linked_lookup_and_rearrange()
{
    LinkedMultiFiltration filtration({Simplex({0}, 0.0f), Simplex({0}, 1.0f), Simplex({1}, 0.0f), Simplex({0, 1}, 2.0f)});

    require(filtration.size() == 4, "linked multi-filtration should preserve duplicate simplices");
    require(filtration.index(Simplex({0}, 99.0f), 1) == 1,
            "linked lookup should use the stored linked index when it matches");

    std::vector<size_t> boundary_indices;
    size_t edge_index = 3;
    for (const auto& face : filtration[edge_index].boundary())
        boundary_indices.push_back(filtration.index(face, edge_index));

    require(boundary_indices == std::vector<size_t>({2, 1}),
            "linked multi-filtration boundary lookup should use latest duplicate before the coface");

    filtration.rearrange({3, 0, 1, 2});
    require(str(filtration[0]) == "<0,1>", "linked rearrange should move the edge first");
    require(filtration.index(Simplex({0, 1}), 0) == 0, "linked rearrange should update cell indices");
}

void test_linked_multi_filtration_sort_updates_linked_indices()
{
    LinkedMultiFiltration filtration;
    filtration.push_back(Simplex({0, 1}, 2.0f), 0);
    filtration.push_back(Simplex({0}, 0.0f), 1);
    filtration.push_back(Simplex({1}, 1.0f), 2);

    filtration.sort();

    require(str(filtration[0]) == "<0>", "linked sort should place vertex 0 first");
    require(str(filtration[1]) == "<1>", "linked sort should place vertex 1 second");
    require(str(filtration[2]) == "<0,1>", "linked sort should place edge last");
    require(filtration.index(Simplex({0, 1}), 2) == 2, "linked sort should update the edge linked index");
}

Filtration edge_filtration()
{
    return Filtration({Simplex({0}, 0.0f), Simplex({1}, 1.0f), Simplex({0, 1}, 2.0f)});
}

void require_edge_pairs(const Persistence& persistence, const std::string& label)
{
    require(persistence.size() == 3, label + " should produce three reduced columns");
    require(persistence.pair(0) == Persistence::unpaired(), label + " should leave vertex 0 unpaired");
    require(persistence.pair(1) == 2, label + " should pair vertex 1 with the edge");
    require(persistence.pair(2) == 1, label + " should pair the edge with vertex 1");
}

void test_standard_reduction_edge_filtration()
{
    Field field(5);
    Persistence persistence(field);
    dionysus::StandardReduction<Persistence> reduce(persistence);

    reduce(edge_filtration());

    require_edge_pairs(persistence, "standard reduction");
}

void test_row_reduction_edge_filtration()
{
    Field field(5);
    dionysus::RowReduction<Field> reduce(field);

    reduce(edge_filtration());

    require_edge_pairs(reduce.persistence(), "row reduction");
}

} // namespace

int main()
{
    test_simplex_empty_and_sequence_behavior();
    test_simplex_boundary_join_and_identity();
    test_filtration_unique_lookup_and_rearrange();
    test_checked_filtration_rejects_missing_cell();
    test_multi_filtration_duplicate_lookup();
    test_linked_multi_filtration_linked_lookup_and_rearrange();
    test_linked_multi_filtration_sort_updates_linked_indices();
    test_standard_reduction_edge_filtration();
    test_row_reduction_edge_filtration();
    return 0;
}
