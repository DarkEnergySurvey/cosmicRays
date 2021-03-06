/*
 * LSST Data Management System
 * Copyright 2014-2015 AURA/LSST.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

/// \file
/// \brief This file contains the Ellipse class implementation.

#include "lsst/sphgeom/Ellipse.h"

#include <cmath>
#include <ostream>
#include <stdexcept>

#include "lsst/sphgeom/Box.h"
#include "lsst/sphgeom/Box3d.h"
#include "lsst/sphgeom/Circle.h"
#include "lsst/sphgeom/ConvexPolygon.h"
#include "lsst/sphgeom/codec.h"


namespace lsst {
namespace sphgeom {

Ellipse::Ellipse(UnitVector3d const & f1, UnitVector3d const & f2, Angle alpha) :
    _a(alpha.asRadians() - 0.5 * PI)
{
    if (alpha.isNan()) {
        throw std::invalid_argument("Invalid ellipse semi-axis angle");
    }
    if (f1 == f2) {
        _gamma = Angle(0.0);
    } else if (f1 == -f2) {
        _gamma = Angle(0.5 * PI);
    } else {
        _gamma = 0.5 * NormalizedAngle(f1, f2);
    }
    if (isEmpty()) {
        *this = empty();
        return;
    } else if (isFull()) {
        *this = full();
        return;
    }
    if (_gamma.asRadians() == 0.0) {
        // The foci are identical, so this ellipse is a circle centered at
        // the common focal point. Pick an orthonormal basis containing f1
        // and use it to construct an orthogonal matrix that maps f1 to
        // (0, 0, 1).
        UnitVector3d b0 = UnitVector3d::orthogonalTo(f1);
        UnitVector3d b1 = UnitVector3d(f1.cross(b0));
        _S = Matrix3d(b0.x(), b0.y(), b0.z(),
                      b1.x(), b1.y(), b1.z(),
                      f1.x(), f1.y(), f1.z());
        _b = _a;
        _tana = std::fabs(tan(_a));
        _tanb = _tana;
        return;
    }
    // _gamma != 0 implies that f1 - f2 != 0. Also, if f1 = -f2 then
    // _gamma = PI/2, and the ellipse must either empty or full. So
    // at this stage f1 + f2 != 0.
    Vector3d b0 = f1 - f2;
    Vector3d b2 = f1 + f2;
    Vector3d b1 = b0.cross(b2);
    b0.normalize();
    b1.normalize();
    b2.normalize();
    _S = Matrix3d(b0.x(), b0.y(), b0.z(),
                  b1.x(), b1.y(), b1.z(),
                  b2.x(), b2.y(), b2.z());
    // Compute _b.
    double r = std::min(1.0, std::max(-1.0, cos(alpha) / cos(_gamma)));
    _b = Angle(std::acos(r) - 0.5 * PI);
    if (_a.asRadians() <= 0.0 && _b > _a) {
        _b = _a;
    } else if (_a.asRadians() > 0.0 && _b < _a) {
        _b = _a;
    }
    _tana = std::fabs(tan(_a));
    _tanb = std::fabs(tan(_b));
    return;
}

Ellipse::Ellipse(UnitVector3d const & center,
                 Angle alpha,
                 Angle beta,
                 Angle orientation)
{
    if (!std::isfinite(orientation.asRadians())) {
        throw std::invalid_argument("Invalid ellipse orientation");
    }
    if (alpha.isNan() ||
        beta.isNan() ||
        (alpha.asRadians() <  0.5 * PI && beta.asRadians() >= 0.5 * PI) ||
        (alpha.asRadians() >  0.5 * PI && beta.asRadians() <= 0.5 * PI) ||
        (alpha.asRadians() == 0.5 * PI && beta.asRadians() != 0.5 * PI)) {
        throw std::invalid_argument("Invalid ellipse semi-axis angle(s)");
    }
    if (alpha.asRadians() < 0.0 || beta.asRadians() < 0.0) {
        *this = empty();
        return;
    } else if (alpha.asRadians() > PI || beta.asRadians() > PI ||
               (alpha.asRadians() == PI && beta.asRadians() == PI)) {
        *this = full();
        return;
    }
    if (alpha == beta) {
        // The ellipse is a circle. Pick an orthonormal basis containing the
        // center and use it to construct some orthogonal matrix that maps the
        // center to (0, 0, 1).
        UnitVector3d b0 = UnitVector3d::orthogonalTo(center);
        UnitVector3d b1 = UnitVector3d(center.cross(b0));
        _S = Matrix3d(b0.x(), b0.y(), b0.z(),
                      b1.x(), b1.y(), b1.z(),
                      center.x(), center.y(), center.z());
        _a = alpha - Angle(0.5 * PI);
        _b = _a;
        _gamma = Angle(0.0);
        _tana = std::fabs(tan(_a));
        _tanb = _tana;
        return;
    }
    if ((alpha.asRadians() < 0.5 * PI && alpha < beta) ||
        (alpha.asRadians() > 0.5 * PI && alpha > beta)) {
        std::swap(alpha, beta);
        orientation = orientation + Angle(0.5 * PI);
    }
    UnitVector3d b0 =
        UnitVector3d::northFrom(center).rotatedAround(center, -orientation);
    UnitVector3d b1 = UnitVector3d(b0.cross(center));
    _S = Matrix3d(b0.x(), b0.y(), b0.z(),
                  b1.x(), b1.y(), b1.z(),
                  center.x(), center.y(), center.z());
    _a = alpha - Angle(0.5 * PI);
    _b = beta - Angle(0.5 * PI);
    double d = std::min(1.0, std::max(-1.0, cos(alpha) / cos(beta)));
    _gamma = Angle(std::acos(d));
    _tana = std::fabs(tan(_a));
    _tanb = std::fabs(tan(_b));
}

bool Ellipse::contains(UnitVector3d const & v) const {
    UnitVector3d const c = getCenter();
    double vdotc = v.dot(c);
    Vector3d u;
    double scz;
    // To maintain high accuracy for very small and very large ellipses,
    // decompose v as v = u ?? c near ??c. Then S v = S u ?? S c, and
    // S c = (0, 0, 1).
    if (vdotc > 0.5) {
        u = v - c;
        scz = 1.0;
    } else if (vdotc < -0.5) {
        u = v + c;
        scz = -1.0;
    } else {
        u = v;
        scz = 0.0;
    }
    u = _S * u;
    double x = u.x() * _tana;
    double y = u.y() * _tanb;
    double z = u.z() + scz;
    double d = (x * x + y * y) - z * z;
    if (_a.asRadians() > 0.0) {
        return z >= 0.0 || d >= 0.0;
    } else {
        return z >= 0.0 && d <= 0.0;
    }
}

Box Ellipse::getBoundingBox() const {
    // For now, simply return the bounding box of the ellipse bounding circle.
    //
    // Improving on this seems difficult, mainly because error bounds must be
    // computed to guarantee that the resulting box tightly bounds the ellipse.
    // In case this ends up being important and someone wants to go down
    // this route in the future, here are some thoughts on how to proceed.
    //
    // First of all, if the ellipse contains a pole, its bounding box must
    // contain all longitudes. Otherwise, consider the plane spanned by the
    // basis vectors u = (0, 0, 1) and v = (cos ??, sin ??, 0). To obtain
    // longitude bounds for the ellipse, we must find values of ?? for which
    // this plane is tangent to the elliptical cone that defines the spherical
    // ellipse boundary. This is the case when the plane intersects the cone
    // in a single line.
    //
    // Let ?? v + ?? u be the direction of the line of intersection, where
    // ??, ?? ??? ???. We know that ?? ??? 0 because u is not on the ellipse boundary,
    // so fix ?? = 1. Let Q be the symmetric matrix representation of the
    // ellipse. Then:
    //
    //    (v + ?? u)??? Q (v + ?? u) = 0
    //
    // Expanding gives:
    //
    //    (v??? Q v) + ?? (u??? Q v) + ?? (v??? Q u) + ???? (u??? Q u) = 0
    //
    // By the symmetry of Q:
    //
    //    ???? (u??? Q u) + 2?? (u??? Q v) + (v??? Q v) = 0
    //
    // This is a quadratic equation which has one solution exactly when its
    // discriminant is 0, i.e. when:
    //
    //    (u??? Q v)?? - (u??? Q u) (v??? Q v) = 0
    //
    // Substituting for u and v and simplifying yields:
    //
    //    (Q???????? - Q??????Q??????)tan???? + 2(Q??????Q?????? - Q??????Q??????)tan ?? + (Q???????? - Q??????Q??????) = 0
    //
    // Finding the latitude bounds can be done by parameterizing the ellipse
    // boundary and then solving for the zeros of the derivative of z with
    // respect to the parameter. This looks to be more involved than the
    // longitude bound calculation, and I haven't worked through the details.
    return getBoundingCircle().getBoundingBox();
}

Box3d Ellipse::getBoundingBox3d() const {
    return getBoundingCircle().getBoundingBox3d();
}

Circle Ellipse::getBoundingCircle() const {
    Angle r = std::max(getAlpha(), getBeta()) + 2.0 * Angle(MAX_ASIN_ERROR);
    return Circle(getCenter(), r);
}

Relationship Ellipse::relate(Box const & b) const {
    return getBoundingCircle().relate(b) & (DISJOINT | WITHIN);
}

// For now, implement ellipse-circle and ellipse-ellipse relation
// computation by approximating ellipses via their bounding circles.
//
// It should be possible to improve on this using the following algorithm to
// compute ellipse-ellipse intersection points.
//
// Ellipses that are neither empty nor full have quadratic forms with
// symmetric matrix representations P, Q of full rank. Consider the matrix
// pencil ?? P + ?? Q, where neither ??, ?? ??? ??? is zero. This is a family of
// quadratic forms, where every member passes through the intersection points
// of P and Q. The scalars ??, ?? are homogeneous, meaning that the quadratic
// forms for (??, ??) and (k ??, k ??) are identical for k ??? 0, so we can fix ?? = 1.
//
// If we can find ?? such that R = P - ?? Q is rank deficient, then the resulting
// quadratic form corresponds to a line or to a pair of (possibly identical)
// planes. The intersection of a plane or a line with either P or Q is
// then easily computed (see below). Finding ?? is a generalized eigenvalue
// problem, and one way to solve it is to note that:
//
//    det(P - ?? Q) = 0   ???   det(PQ????? - ?? I) = 0
//
// so that the values of ?? we are interested in are the eigenvalues of PQ?????.
// Use one of these eigenvalues to construct R. Then:
//
// 1) If R has rank 0, P and Q are equivalent up to scale and the corresponding
//    quadratic forms are identical.
//
// 2) If R has rank 1, then the quadratic form R can be factorized as
//    (ax + by + cz)?? = 0. Why? There is an eigen-decomposition of R,
//    R = M D M???, where D is diagonal, M is orthogonal, D?????? is the only
//    non-zero eigenvalue and the first column of M is the unit eigenvector
//    n for D??????. Now
//
//    v??? R v = (M??? v)??? D (M??? v) = (n??v)?? = 0
//
//    This just says that there is some basis in which the quadratic form
//    is x?? = 0, which is the plane defined by x = 0.
//
// 3) Otherwise, R has rank 2, and corresponds to a line or to a pair of
//    planes. To see why, again consider its eigen-decomposition R = M D M???,
//    where M is orthogonal, D is diagonal with non-zero eigenvalues D?????? and
//    D??????, and the first two columns of M are the corresponding eigenvectors.
//
//    If D?????? and D?????? have the same sign, then there is a basis in which
//    the quadratic form is a?? x?? + b?? y?? = 0, i.e. the z-axis. In the
//    original basis this is a line perpendicular to the two eigenvectors of
//    R where P and Q are tangent.
//
//    If D?????? and D?????? have opposite signs, then there is some basis in which
//    the quadratic form reduces to a?? x?? - b?? y?? = (a x - b y)(a x + b y) = 0;
//    that is, solutions lie on two planes.
//
// To compute the intersection of an ellipse with quadratic form Q and a
// plane with unit normal n, we must solve v??? Q v = 0 subject to v??n = 0.
// Pick two unit vectors b???, b??? orthogonal to n, write v = ?? b??? + ?? b???, and
// substitute into the equation of the quadratic form to obtain:
//
//    ???? (b?????? Q b???) + 2???? (b?????? Q b???) + ???? (b?????? Q b???) = 0
//
// This only has solutions when:
//
//    (b?????? Q b???)?? ??? (b?????? Q b???) (b?????? Q b???)
//
// As in the case of the ellipse bounding box computation, actually using the
// above in an implementation of relate() involves a non-trivial error analysis.
//
// TODO(smm): investigate the theory of matrix pencils and generalized
// eigenvalue problems. A couple unanswered questions this could shed some
// light on are:
//
// - What algorithm should be used to solve for the generalized eigenvalues?
//   Note that very small and very large ellipses will have matrix
//   representations with very large condition numbers.
// - Which of the generalized eigenvalues should be chosen / leads to the most
//   accurate computation? Is there some usefully exploitable relationship
//   between them and the degenerate quadratic forms they engender?

Relationship Ellipse::relate(Circle const & c) const {
    return getBoundingCircle().relate(c) & (DISJOINT | WITHIN);
}

Relationship Ellipse::relate(ConvexPolygon const & p) const {
    return getBoundingCircle().relate(p) & (DISJOINT | WITHIN);
}

Relationship Ellipse::relate(Ellipse const & e) const {
    return getBoundingCircle().relate(e.getBoundingCircle()) & DISJOINT;
}

std::vector<uint8_t> Ellipse::encode() const {
    std::vector<uint8_t> buffer;
    uint8_t tc = TYPE_CODE;
    buffer.reserve(ENCODED_SIZE);
    buffer.push_back(tc);
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            encodeDouble(_S(r, c), buffer);
        }
    }
    encodeDouble(_a.asRadians(), buffer);
    encodeDouble(_b.asRadians(), buffer);
    encodeDouble(_gamma.asRadians(), buffer);
    encodeDouble(_tana, buffer);
    encodeDouble(_tanb, buffer);
    return buffer;
}

std::unique_ptr<Ellipse> Ellipse::decode(uint8_t const * buffer, size_t n) {
    if (buffer == nullptr || n != ENCODED_SIZE || buffer[0] != TYPE_CODE) {
        throw std::runtime_error("Byte-string is not an encoded Ellipse");
    }
    std::unique_ptr<Ellipse> ellipse(new Ellipse);
    ++buffer;
    double m00 = decodeDouble(buffer); buffer += 8;
    double m01 = decodeDouble(buffer); buffer += 8;
    double m02 = decodeDouble(buffer); buffer += 8;
    double m10 = decodeDouble(buffer); buffer += 8;
    double m11 = decodeDouble(buffer); buffer += 8;
    double m12 = decodeDouble(buffer); buffer += 8;
    double m20 = decodeDouble(buffer); buffer += 8;
    double m21 = decodeDouble(buffer); buffer += 8;
    double m22 = decodeDouble(buffer); buffer += 8;
    ellipse->_S = Matrix3d(m00, m01, m02,
                           m10, m11, m12,
                           m20, m21, m22);
    double a = decodeDouble(buffer); buffer += 8;
    double b = decodeDouble(buffer); buffer += 8;
    double gamma = decodeDouble(buffer); buffer += 8;
    ellipse->_a = Angle(a);
    ellipse->_b = Angle(b);
    ellipse->_gamma = Angle(gamma);
    double tana = decodeDouble(buffer); buffer += 8;
    double tanb = decodeDouble(buffer); buffer += 8;
    ellipse->_tana = tana;
    ellipse->_tanb = tanb;
    return ellipse;
}

std::ostream & operator<<(std::ostream & os, Ellipse const & e) {
    os << "{\"Ellipse\": [" << e.getTransformMatrix() << ", "
       << e.getAlpha() << ", " << e.getBeta() << "]}";
    return os;
}

}} // namespace lsst::sphgeom
