#include "classes/iface_state_bits.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"

#include <cstdio>
#include <cstdlib>

namespace iface_state_bits {

unsigned layoutGeneration { 0 };

namespace {

std::vector<Layout> layouts {};

int bits_for(std::size_t maxCode)
{
    int bits { 1 };
    while ((std::size_t { 1 } << bits) <= maxCode)
        ++bits;
    return bits;
}

} // namespace

void build_layouts(const std::vector<MolTemplate>& molTemplateList)
{
    layouts.assign(molTemplateList.size(), Layout {});
    for (std::size_t t { 0 }; t < molTemplateList.size(); ++t) {
        Layout& layout = layouts[t];
        int shift { 0 };
        for (const Interface& iface : molTemplateList[t].interfaceList) {
            Field field;
            for (const Interface::State& state : iface.stateList) {
                bool seen { false };
                for (char iden : field.idens)
                    seen = seen || iden == state.iden;
                if (!seen)
                    field.idens.push_back(state.iden);
            }
            field.shift = shift;
            field.width = 1 + bits_for(field.idens.size());
            shift += field.width;
            layout.fields.push_back(field);
        }
        layout.usable = shift <= 63;
    }
    ++layoutGeneration;
}

const Layout* layout_for(int molTypeIndex)
{
    if (molTypeIndex < 0 || static_cast<std::size_t>(molTypeIndex) >= layouts.size()
        || !layouts[molTypeIndex].usable)
        return nullptr;
    return &layouts[molTypeIndex];
}

std::uint64_t encode(const Molecule& mol)
{
    const Layout* layout { layout_for(mol.molTypeIndex) };
    if (layout == nullptr || mol.interfaceList.size() != layout->fields.size())
        return invalid;
    std::uint64_t word { 0 };
    for (std::size_t i { 0 }; i < mol.interfaceList.size(); ++i) {
        const Molecule::Iface& iface = mol.interfaceList[i];
        // hasIntangibles() reads interface relIfaceIndex at that position only
        // when relIndex agrees, and compares the interface's own molTypeIndex.
        if (iface.relIndex != static_cast<int>(i) || iface.molTypeIndex != mol.molTypeIndex)
            return invalid;
        const Field& field = layout->fields[i];
        std::uint64_t code { 0 };
        for (std::size_t k { 0 }; k < field.idens.size(); ++k) {
            if (field.idens[k] == iface.stateIden) {
                code = k + 1;
                break;
            }
        }
        word |= ((code << 1) | static_cast<std::uint64_t>(iface.isBound)) << field.shift;
    }
    return word;
}

Requirement compile_requirement(
    int molTypeIndex, int relIfaceIndex, bool isBound, char state, std::uint64_t& mask, std::uint64_t& value)
{
    const Layout* layout { layout_for(molTypeIndex) };
    if (layout == nullptr)
        return Requirement::unlisted;
    if (relIfaceIndex < 0 || static_cast<std::size_t>(relIfaceIndex) >= layout->fields.size())
        return Requirement::unsatisfiable;
    const Field& field = layout->fields[relIfaceIndex];
    std::uint64_t code { 0 };
    for (std::size_t k { 0 }; k < field.idens.size(); ++k) {
        if (field.idens[k] == state) {
            code = k + 1;
            break;
        }
    }
    if (code == 0)
        return Requirement::unlisted;
    mask = ((std::uint64_t { 1 } << field.width) - 1) << field.shift;
    value = ((code << 1) | static_cast<std::uint64_t>(isBound)) << field.shift;
    return Requirement::ok;
}

void refresh_all(std::vector<Molecule>& moleculeList)
{
    for (Molecule& mol : moleculeList)
        mol.refresh_iface_state_bits();
}

void verify(const Molecule& mol, const char* where)
{
    const std::uint64_t expected { encode(mol) };
    if (mol.ifaceStateBits == expected)
        return;
    std::fprintf(stderr, "STATE BITS MISMATCH at %s: molecule %d (type %d): word %016llx, interfaces give %016llx\n",
        where, mol.index, mol.molTypeIndex, static_cast<unsigned long long>(mol.ifaceStateBits),
        static_cast<unsigned long long>(expected));
    for (const Molecule::Iface& iface : mol.interfaceList)
        std::fprintf(stderr, "  iface %d: isBound %d state '%c' (%d) molType %d\n", iface.relIndex, iface.isBound,
            iface.stateIden ? iface.stateIden : '0', iface.stateIden, iface.molTypeIndex);
    std::abort();
}

void verify_all(const std::vector<Molecule>& moleculeList, const char* where)
{
    for (const Molecule& mol : moleculeList)
        verify(mol, where);
}

} // namespace iface_state_bits

void Molecule::refresh_iface_state_bits() { ifaceStateBits = iface_state_bits::encode(*this); }
