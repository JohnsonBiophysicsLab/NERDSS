#pragma once

#include <cstdint>
#include <vector>

struct MolTemplate;
struct Molecule;

/*!
 * \brief Each Molecule's interface states -- isBound and stateIden of every
 * interface -- packed into one word, Molecule::ifaceStateBits, so that rate
 * selection can test a rate state's ancillary requirements with a few word
 * operations instead of reading the molecule's interfaces.
 *
 * The layout belongs to the MolTemplate: interface i owns a field at
 * Layout::fields[i].shift holding (code << 1) | isBound, where code is 1 + the
 * position of stateIden among the interface's distinct state names, and 0 for
 * a name the template does not list.  The fields are exact, not hashed, so a
 * requirement (interface, isBound, state) holds exactly when its field equals
 * the requirement's value.  A template whose fields need more than 63 bits has
 * no usable layout.
 *
 * The word is derived data.  It is never written to restart files or sent
 * between MPI ranks: every site that changes an interface's isBound or
 * stateIden calls Molecule::refresh_iface_state_bits() afterwards, and the
 * start-up paths and Molecule::deserialize() recompute it.  A word that cannot
 * describe its molecule -- no layout yet, a molecule with the wrong number of
 * interfaces, an interface whose relIndex or molTypeIndex disagrees with its
 * position or molecule -- is `invalid`, and rate selection then reads the
 * interfaces as before.  Building with -DNERDSS_CHECK_STATE_BITS checks every
 * word against its interfaces at every selection and at the top of every step.
 */
namespace iface_state_bits {

//! Set in a word that does not describe its molecule; never set by a field.
constexpr std::uint64_t invalid { std::uint64_t { 1 } << 63 };

struct Field {
    int shift { 0 };
    int width { 0 };
    std::vector<char> idens {}; //!< distinct state names in stateList order; code = position + 1
};

struct Layout {
    bool usable { false };
    std::vector<Field> fields {}; //!< per interface, in interfaceList order
};

/*!
 * \brief Lays out every template's fields and starts a new generation.
 *
 * Call once the templates are final, then refresh_all(); compiled rate lists
 * built under an earlier generation rebuild themselves on their next use.
 */
void build_layouts(const std::vector<MolTemplate>& molTemplateList);

extern unsigned layoutGeneration;
//! 0 until build_layouts() has run.  Read on every compiled rate selection, so inline.
inline unsigned generation() { return layoutGeneration; }

//! The layout of this template, or nullptr when there is none.
const Layout* layout_for(int molTypeIndex);

//! The word for this molecule as its interfaces stand now.
std::uint64_t encode(const Molecule& mol);

//! What a requirement compiles to: its field's mask and the value the field must hold.
enum class Requirement { ok, unsatisfiable, unlisted };

/*!
 * \brief Compiles one requirement on an interface of a molecule of this template.
 *
 * `unsatisfiable`: no molecule of the template can meet it (the interface does
 * not exist).  `unlisted`: it names a state the template does not list, which
 * the word cannot tell apart from other unlisted names, so the caller must not
 * use the word for it.
 */
Requirement compile_requirement(
    int molTypeIndex, int relIfaceIndex, bool isBound, char state, std::uint64_t& mask, std::uint64_t& value);

void refresh_all(std::vector<Molecule>& moleculeList);

/*!
 * \brief Aborts unless mol.ifaceStateBits == encode(mol).  Used by the
 * NERDSS_CHECK_STATE_BITS build; `where` names the caller in the report.
 */
void verify(const Molecule& mol, const char* where);
void verify_all(const std::vector<Molecule>& moleculeList, const char* where);

} // namespace iface_state_bits
