/*
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef LSST_GEOM_LINEAR_TRANSFORM_H
#define LSST_GEOM_LINEAR_TRANSFORM_H

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/geom/Point.h"
#include "lsst/geom/Angle.h"

namespace lsst {
namespace geom {

LSST_EXCEPTION_TYPE(SingularTransformException, lsst::pex::exceptions::RuntimeError,
                    lsst::geom::SingularTransformException)

/**
 *  A 2D linear coordinate transformation.
 *
 *  The transform is represented by a matrix @f$ \mathbf{M} @f$ such that
 *  @f[
 *     \left[\begin{array}{ c }
 *     x_f \\
 *     y_f
 *     \end{array}\right]
 *     =
 *     \mathbf{M}
 *     \left[\begin{array}{ c }
 *     x_i \\
 *     y_i
 *     \end{array}\right]
 *  @f]
 *  where @f$(x_i,y_i)@f$ are the input coordinates and @f$(x_f,y_f)@f$ are
 *  the output coordinates.
 *
 *  If @f$ x_f(x_i,y_i) @f$ and @f$ y_f(x_i,y_i) @f$ are continuous
 *  differentiable functions, then
 *  @f[
 *     \mathbf{M} = \left[\begin{array}{ c c }
 *     \displaystyle\frac{\partial x_f}{\partial x_i} &
 *     \displaystyle\frac{\partial x_f}{\partial y_i} \\
 *     \displaystyle\frac{\partial y_f}{\partial x_i} &
 *     \displaystyle\frac{\partial y_f}{\partial y_i}
 *     \end{array}\right]
 *  @f]
 *  evaluated at @f$(x_i,y_i)@f$.
 */
class LinearTransform final {
public:
    enum Parameters { XX = 0, YX = 1, XY = 2, YY = 3 };

    typedef Eigen::Matrix<double, 4, 1> ParameterVector;
    typedef Eigen::Matrix<double, 2, 4> TransformDerivativeMatrix;
    typedef Eigen::Matrix<double, 4, 4> ProductDerivativeMatrix;

    typedef Eigen::Matrix<double, 2, 2, Eigen::DontAlign> Matrix;

    /** Construct an empty (identity) LinearTransform. */
    LinearTransform() noexcept : _matrix(Matrix::Identity()) {}

    /** Construct an LinearTransform from an Eigen::Matrix. */
    explicit LinearTransform(Matrix const& matrix) noexcept : _matrix(matrix) {}

    // Can't use both =default and noexcept until Eigen supports noexcept
    LinearTransform(LinearTransform const& other) noexcept : _matrix(other._matrix){};
    LinearTransform(LinearTransform&& other) noexcept : _matrix(std::move(other._matrix)){};
    ~LinearTransform() noexcept = default;

    LinearTransform operator*(LinearTransform const& other) const noexcept {
        return LinearTransform(getMatrix() * other.getMatrix());
    }

    static LinearTransform makeScaling(double s) noexcept {
        return LinearTransform((Matrix() << s, 0.0, 0.0, s).finished());
    }

    static LinearTransform makeScaling(double s, double t) noexcept {
        return LinearTransform((Matrix() << s, 0.0, 0.0, t).finished());
    }

    static LinearTransform makeRotation(Angle t) noexcept {
        return LinearTransform(Matrix(Eigen::Rotation2D<double>(t.asRadians())));
    }

    // Can't use both =default and noexcept until Eigen supports noexcept
    LinearTransform& operator=(LinearTransform const& other) noexcept {
        _matrix = other._matrix;
        return *this;
    }
    LinearTransform& operator=(LinearTransform&& other) noexcept {
        _matrix = std::move(other._matrix);
        return *this;
    }

    LinearTransform& operator+=(LinearTransform const& other) noexcept {
        _matrix += other._matrix;
        return *this;
    }

    LinearTransform operator+(LinearTransform const& other) noexcept {
        LinearTransform tmp(*this);
        tmp += other;
        return tmp;
    }

    LinearTransform& operator-=(LinearTransform const& other) noexcept {
        _matrix -= other._matrix;
        return *this;
    }

    LinearTransform operator-(LinearTransform const& other) noexcept {
        LinearTransform tmp(*this);
        tmp -= other;
        return tmp;
    }

    /**
     * Return the transform matrix elements as a parameter vector
     *
     * The elements will be ordered XX, YX, XY, YY
     */
    ParameterVector const getParameterVector() const noexcept;
    /**
     * Set the transform matrix elements from a parameter vector
     *
     * The parameter vector is ordered XX, YX, XY, YY
     */
    void setParameterVector(ParameterVector const& vector) noexcept;

    Matrix const& getMatrix() const noexcept { return _matrix; }
    Matrix& getMatrix() noexcept { return _matrix; }

    double& operator[](int i) { return _matrix(i % 2, i / 2); }
    double const& operator[](int i) const { return const_cast<Matrix&>(_matrix)(i % 2, i / 2); }

    /**
     * Return the inverse transform.
     *
     * @throws lsst::geom::SingularTransformException if not invertible
     */
    LinearTransform const inverted() const;

    /**
     * Return the determinant of the 2x2 matrix
     */
    double computeDeterminant() const noexcept;

    /** Whether the transform is a no-op. */
    bool isIdentity() const noexcept { return getMatrix().isIdentity(); }

    /**
     *  Transform a Point2D object.
     *
     *  This operation is equivalent to applying the LinearTransform to an
     *  lsst::geom::Extent
     */
    Point2D operator()(Point2D const& p) const noexcept { return Point2D(getMatrix() * p.asEigen()); }

    /**
     *  Transform a Extent2D object.
     *
     *  This operation is equivalent to applying the LinearTransform to an
     *  lsst::geom::Point
     */
    Extent2D operator()(Extent2D const& p) const noexcept { return Extent2D(getMatrix() * p.asEigen()); }

    //@{
    /**
     *  Transform a point given and returned as separate double values.
     *
     *  This interface is intended primarily for use in Python (where it is
     *  vectorized to support NumPy array arguments).
     */
    double applyX(double x, double y) const noexcept { return _matrix(0, 0)*x + _matrix(0, 1)*y; }
    double applyY(double x, double y) const noexcept { return _matrix(1, 0)*x + _matrix(1, 1)*y; }
    //@}

    /**
     * Derivative of (*this)(input) with respect to the transform elements (for Point).
     */
    TransformDerivativeMatrix dTransform(Point2D const& input) const noexcept;

    /// Derivative of (*this)(input) with respect to the transform elements (for Extent);
    TransformDerivativeMatrix dTransform(Extent2D const& input) const noexcept {
        return dTransform(Point2D(input));
    }

private:
    Matrix _matrix;
};

std::ostream& operator<<(std::ostream& os, lsst::geom::LinearTransform const& t);

}  // namespace geom
}  // namespace lsst

#endif  // !LSST_GEOM_LINEAR_TRANSFORM_H
