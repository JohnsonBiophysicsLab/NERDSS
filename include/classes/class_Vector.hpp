/*! \file class_vector.hpp
 * \ingroup Associate
 * Created on 5/20/18 by Matthew Varga
 * Purpose: Vector class for association
 * Notes:
 */
#pragma once

#include "class_Coord.hpp"

#include <cmath>

/*! \ingroup Associate
 * \brief Holds a vector with the origin as the start point
 *
 * A 2D vector will always be an x and y coordinate
 *
 * The trivial constructors and the coordinate arithmetic live in this header on
 * purpose.  They are called from the innermost association and propagation
 * loops, so leaving them in class_Vector.cpp turned each one into an
 * out-of-line call: profiling a clathrin run showed Vector::Vector(Coord) alone
 * holding 3.3% of the total samples for what is a three-double copy.  This
 * mirrors what was already done for Coord.
 */
struct Vector : public Coord {
    double magnitude { 0.0 };

    // operators
    Vector operator*(const double val) const { return { x * val, y * val, z * val }; }
    Vector operator-(const double val) const { return { x - val, y - val, z - val }; };
    Vector operator-(const Vector& vec) const { return { this->x - vec.x, this->y - vec.y, this->z - vec.z }; }
    Vector operator-() const { return { -x, -y, -z }; }
    /*! \brief Divides the Vector in place.
     *
     * Returns a reference, like every other compound assignment.  This used to
     * return a fresh Vector and leave *this untouched, so `vec /= val` silently
     * did nothing; the only caller, normalize(), compensated by writing
     * `*this = *this /= magnitude`.  magnitude itself is deliberately not
     * rescaled here: callers that need it call calc_magnitude() afterwards.
     */
    Vector& operator/=(double val)
    {
        x /= val;
        y /= val;
        z /= val;
        return *this;
    }
    Vector operator/(double val) const { return { x / val, y / val, z / val }; }
    /*! \warning Deliberately left non-const, unlike the operators above.
     *
     * The free operator+(const Vector&, const Coord&) below is an equally exact
     * match for a const Vector but returns Coord rather than Vector, so which
     * overload a call site gets — and therefore the type of the result — is
     * currently decided by whether the left operand is const.  Adding const
     * here makes the two ambiguous.  Untangling that pair changes result types
     * at call sites, so it does not belong in a result-preserving change.
     */
    Vector operator+(const Coord& crd) { return { this->x + crd.x, this->y + crd.y, this->z + crd.z }; }
    friend std::ostream& operator<<(std::ostream& os, Vector& vec);
    friend Coord operator+(const Vector& vec, const Coord& crd)
    {
        return { vec.x + crd.x, vec.y + crd.y, vec.z + crd.z };
    }

    /*!
     * \brief Calculates the magnitude of the Vector. Note that this is no longer done when the Vector is constructed.
     */
    void calc_magnitude() { this->magnitude = sqrt(x * x + y * y + z * z); }

    /*!
     * \brief Normalizes the vector to unity.
     */
    void normalize()
    {
        calc_magnitude();
        if (magnitude == 0) {
            // dividing by zero results in NANs
            magnitude = 1;
        }
        *this /= magnitude;
        calc_magnitude();
    }

    /*!
     * \brief Returns a Vectorwhich is the crossproduct of the Vector with another Vector.
     * \params[in] this Vector 1.
     * \params[in] vec Vector 2.
     */
    Vector cross(const Vector& vec) const
    {
        double u0 = this->y * vec.z - this->z * vec.y;
        double u1 = this->z * vec.x - this->x * vec.z;
        double u2 = this->x * vec.y - this->y * vec.x;

        Vector test { u0, u1, u2 };
        test.normalize();
        return test;
    }

    /*!
     * \brief Returns the dot product of the Vector with another Vector.
     * \params[in] this Vector 1.
     * \params[in] vec Vector 2.
     */
    double dot(const Vector& vec) const
    {
        return ((this->x * vec.x) + (this->y * vec.y) + (this->z * vec.z));
    }

    /*!
     * \brief Returns the angle between the Vector and another Vector.
     * \params[in] this Vector 1.
     * \params[in] vec Vector 2.
     */
    double dot_theta(const Vector& vec) const; //!< returns the angle between two vectors

    /*!
     * \brief returns the projection of the vector onto the normal.
     *
     * \$f v2 = v1 - \frac{v1 \dot v2}{v2 \dot v2} \times v2 \$f
     */
    Vector vector_projection(Vector normal);

    Vector() = default;
    Vector(double x, double y)
        : Coord(x, y, 0)
    {
    }
    Vector(const double& x, const double& y, const double& z)
        : Coord(x, y, z)
    {
    }
    explicit Vector(std::array<double, 3>& arr);
    explicit Vector(const std::vector<double>& arr);
    /// constructor where the start point is implied to be the origin
    explicit Vector(const Coord& coord)
        : Coord(coord.x, coord.y, coord.z)
    {
    }
    /// constructor which takes in two points and outputs a vector with the start point at the origin
    Vector(const Coord& coordEnd, const Coord& coordStart)
        : Coord(coordEnd.x - coordStart.x, coordEnd.y - coordStart.y, coordEnd.z - coordStart.z)
    {
    }

    /*
    Function serialize serializes the Vector into array of bytes.
    */
    void serialize(unsigned char* arrayRank, int& nArrayRank) {
        PUSH(x);
        PUSH(y);
        PUSH(z);
        PUSH(magnitude);
    }
    /*
    Function deserialize deserializes the Vector from array of bytes.
    */
    void deserialize(unsigned char* arrayRank, int& nArrayRank) {
        POP(x);
        POP(y);
        POP(z);
        POP(magnitude);
    }
};
