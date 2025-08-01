#include "catch.hpp"
#include "../cppambergromacs/CppLibAmberGromacs.hpp"

#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <thread>

TEST_CASE("Vector - Constructor and Basic Operations", "[Vector]") {
    SECTION("Default constructor initializes to zero") {
        Vector v;
        REQUIRE(v.x == 0.0);
        REQUIRE(v.y == 0.0);
        REQUIRE(v.z == 0.0);
    }

    SECTION("Parameterized constructor sets components") {
        Vector v(1.0, 2.0, 3.0);
        REQUIRE(v.x == 1.0);
        REQUIRE(v.y == 2.0);
        REQUIRE(v.z == 3.0);
    }

    SECTION("Copy constructor copies components") {
        Vector v1(1.0, 2.0, 3.0);
        Vector v2(v1);
        REQUIRE(v2.x == v1.x);
        REQUIRE(v2.y == v1.y);
        REQUIRE(v2.z == v1.z);
    }

    SECTION("Parameterized constructor with extreme values") {
        Vector v(std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), std::numeric_limits<float>::min());
        REQUIRE(v.x == std::numeric_limits<float>::max());
        REQUIRE(v.y == -std::numeric_limits<float>::max());
        REQUIRE(v.z == std::numeric_limits<float>::min());
    }

    SECTION("Parameterized constructor with small values") {
        Vector v(Vector::EPSILON, -Vector::EPSILON, 0.0);
        REQUIRE(v.x == Approx(Vector::EPSILON));
        REQUIRE(v.y == Approx(-Vector::EPSILON));
        REQUIRE(v.z == 0.0);
    }

    SECTION("Constructor with denormalized numbers") {
        Vector v(std::numeric_limits<float>::denorm_min(), std::numeric_limits<float>::denorm_min(), 0.0);
        REQUIRE(v.x == std::numeric_limits<float>::denorm_min());
        REQUIRE(v.y == std::numeric_limits<float>::denorm_min());
        REQUIRE(v.z == 0.0);
    }
}

TEST_CASE("Vector - Magnitude and Normalization", "[Vector]") {
    Vector v(3.0, 4.0, 0.0);

    SECTION("Magnitude calculation") {
        REQUIRE(v.magnitude() == Approx(5.0));
    }

    SECTION("Has equal magnitude") {
        Vector v2(0.0, 5.0, 0.0);
        REQUIRE(v.hasEqualMagnitude(v2));
        Vector v3(6.0, 8.0, 0.0);
        REQUIRE_FALSE(v.hasEqualMagnitude(v3));
    }

    SECTION("Normalization modifies vector in place") {
        v.normalize();
        REQUIRE(v.magnitude() == Approx(1.0));
        REQUIRE(v.x == Approx(3.0 / 5.0));
        REQUIRE(v.y == Approx(4.0 / 5.0));
        REQUIRE(v.z == 0.0);
    }

    SECTION("Get normalized returns new normalized vector") {
        Vector norm = v.getNormalized();
        REQUIRE(norm.magnitude() == Approx(1.0));
        REQUIRE(norm.x == Approx(3.0 / 5.0));
        REQUIRE(norm.y == Approx(4.0 / 5.0));
        REQUIRE(norm.z == 0.0);
        // Original vector unchanged
        REQUIRE(v.x == 3.0);
        REQUIRE(v.y == 4.0);
        REQUIRE(v.z == 0.0);
    }

    SECTION("Normalization of zero vector") {
        Vector zero(0.0, 0.0, 0.0);
        zero.normalize();
        REQUIRE(zero.x == 0.0);
        REQUIRE(zero.y == 0.0);
        REQUIRE(zero.z == 0.0);
        Vector norm_zero = zero.getNormalized();
        REQUIRE(norm_zero.x == 0.0);
        REQUIRE(norm_zero.y == 0.0);
        REQUIRE(norm_zero.z == 0.0);
    }

    SECTION("Magnitude of very small vector") {
        Vector small(Vector::EPSILON, Vector::EPSILON, Vector::EPSILON);
        REQUIRE(small.magnitude() == Approx(std::sqrt(3 * Vector::EPSILON * Vector::EPSILON)));
    }

    SECTION("Normalization of large vector") {
        Vector large(1e10, 1e10, 1e10);
        Vector norm = large.getNormalized();
        REQUIRE(norm.magnitude() == Approx(1.0));
        REQUIRE(norm.x == Approx(1.0 / std::sqrt(3.0)));
        REQUIRE(norm.y == Approx(1.0 / std::sqrt(3.0)));
        REQUIRE(norm.z == Approx(1.0 / std::sqrt(3.0)));
    }

    SECTION("Has equal magnitude with near-equal vectors") {
        Vector v1(1.0, 0.0, 0.0);
        Vector v2(1.0 + Vector::EPSILON / 2.0, 0.0, 0.0);
        REQUIRE(v1.hasEqualMagnitude(v2));
        Vector v3(1.0 + Vector::EPSILON * 2.0, 0.0, 0.0);
        REQUIRE_FALSE(v1.hasEqualMagnitude(v3));
    }

    SECTION("Magnitude with negative components") {
        Vector v_neg(-3.0, -4.0, 0.0);
        REQUIRE(v_neg.magnitude() == Approx(5.0));
    }
}

TEST_CASE("Vector - Operator Overloads", "[Vector]") {
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(4.0, 5.0, 6.0);

    SECTION("Addition operator") {
        Vector sum = v1 + v2;
        REQUIRE(sum.x == 5.0);
        REQUIRE(sum.y == 7.0);
        REQUIRE(sum.z == 9.0);
    }

    SECTION("Subtraction operator") {
        Vector diff = v1 - v2;
        REQUIRE(diff.x == -3.0);
        REQUIRE(diff.y == -3.0);
        REQUIRE(diff.z == -3.0);
    }

    SECTION("Dot product operator") {
        Real dot = v1 * v2;
        REQUIRE(dot == Approx(1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0));
    }

    SECTION("Scalar multiplication (vector * scalar)") {
        Vector scaled = v1 * 2.0;
        REQUIRE(scaled.x == 2.0);
        REQUIRE(scaled.y == 4.0);
        REQUIRE(scaled.z == 6.0);
    }

    SECTION("Scalar multiplication (scalar * vector)") {
        Vector scaled = 2.0 * v1;
        REQUIRE(scaled.x == 2.0);
        REQUIRE(scaled.y == 4.0);
        REQUIRE(scaled.z == 6.0);
    }

    SECTION("Scalar division") {
        Vector divided = v1 / 2.0;
        REQUIRE(divided.x == 0.5);
        REQUIRE(divided.y == 1.0);
        REQUIRE(divided.z == 1.5);
    }

    SECTION("Cross product operator") {
        Vector cross = v1 % v2;
        REQUIRE(cross.x == Approx(2.0 * 6.0 - 3.0 * 5.0));
        REQUIRE(cross.y == Approx(3.0 * 4.0 - 1.0 * 6.0));
        REQUIRE(cross.z == Approx(1.0 * 5.0 - 2.0 * 4.0));
    }

    SECTION("Compound assignment operators") {
        Vector v = v1;
        v += v2;
        REQUIRE(v.x == 5.0);
        REQUIRE(v.y == 7.0);
        REQUIRE(v.z == 9.0);

        v = v1;
        v -= v2;
        REQUIRE(v.x == -3.0);
        REQUIRE(v.y == -3.0);
        REQUIRE(v.z == -3.0);

        v = v1;
        v *= 2.0;
        REQUIRE(v.x == 2.0);
        REQUIRE(v.y == 4.0);
        REQUIRE(v.z == 6.0);

        v = v1;
        v /= 2.0;
        REQUIRE(v.x == 0.5);
        REQUIRE(v.y == 1.0);
        REQUIRE(v.z == 1.5);

        v = v1;
        v %= v2;
        REQUIRE(v.x == Approx(2.0 * 6.0 - 3.0 * 5.0));
        REQUIRE(v.y == Approx(3.0 * 4.0 - 1.0 * 6.0));
        REQUIRE(v.z == Approx(1.0 * 5.0 - 2.0 * 4.0));
    }

    SECTION("Comparison operators") {
        Vector v3(1.0, 2.0, 3.0);
        REQUIRE(v1 == v3);
        REQUIRE_FALSE(v1 == v2);
        REQUIRE(v1 != v2);
        REQUIRE_FALSE(v1 != v3);

        Vector v4(2.0, 0.0, 0.0); // Magnitude sqrt(4) = 2
        Vector v5(0.0, 2.0, 0.0); // Magnitude sqrt(4) = 2
        Vector v6(3.0, 4.0, 0.0); // Magnitude sqrt(25) = 5
        REQUIRE(v4 < v6);
        REQUIRE(v6 > v4);
        REQUIRE(v4 <= v5);
        REQUIRE(v4 >= v5);
    }

    SECTION("Operators with zero vector") {
        Vector zero(0.0, 0.0, 0.0);
        Vector sum = v1 + zero;
        REQUIRE(sum == v1);
        Vector diff = v1 - zero;
        REQUIRE(diff == v1);
        Real dot = v1 * zero;
        REQUIRE(dot == Approx(0.0));
        Vector cross = v1 % zero;
        REQUIRE(cross == Vector(0.0, 0.0, 0.0));
    }

    SECTION("Scalar division by near-zero") {
        Vector divided = v1 / Vector::EPSILON;
        REQUIRE(divided.x == Approx(1.0 / Vector::EPSILON));
        REQUIRE(divided.y == Approx(2.0 / Vector::EPSILON));
        REQUIRE(divided.z == Approx(3.0 / Vector::EPSILON));
    }

    SECTION("Comparison with near-equal vectors") {
        Vector v_near(1.0 + Vector::EPSILON / 2.0, 2.0, 3.0);
        REQUIRE(v1 == v_near); // Within EPSILON tolerance
        Vector v_diff(1.0 + Vector::EPSILON * 2.0, 2.0, 3.0);
        REQUIRE(v1 != v_diff); // Outside EPSILON tolerance
    }

    SECTION("Compound assignment with self") {
        Vector v = v1;
        v += v;
        REQUIRE(v.x == 2.0);
        REQUIRE(v.y == 4.0);
        REQUIRE(v.z == 6.0);

        v = v1;
        v -= v;
        REQUIRE(v == Vector(0.0, 0.0, 0.0));

        v = v1;
        v %= v;
        REQUIRE(v == Vector(0.0, 0.0, 0.0)); // Cross product with self is zero
    }

    SECTION("Chained operations") {
        Vector result = (v1 + v2) * 2.0 - v1;
        REQUIRE(result.x == Approx(9.0));
        REQUIRE(result.y == Approx(12.0));
        REQUIRE(result.z == Approx(15.0));
    }

    SECTION("Scalar multiplication with negative scalar") {
        Vector scaled = v1 * -2.0;
        REQUIRE(scaled.x == -2.0);
        REQUIRE(scaled.y == -4.0);
        REQUIRE(scaled.z == -6.0);
    }

    SECTION("Comparison with zero magnitude") {
        Vector zero(0.0, 0.0, 0.0);
        REQUIRE(zero <= v1);
        REQUIRE(v1 > zero);
        REQUIRE(zero == zero);
    }
}

TEST_CASE("Vector - Additional Methods", "[Vector]") {
    Vector v1(1.0, 0.0, 0.0);
    Vector v2(0.0, 1.0, 0.0);

    SECTION("Dot product method") {
        REQUIRE(v1.dot(v2) == Approx(0.0));
    }

    SECTION("Cross product method") {
        Vector cross = v1.cross(v2);
        REQUIRE(cross.x == 0.0);
        REQUIRE(cross.y == 0.0);
        REQUIRE(cross.z == Approx(1.0));
    }

    SECTION("Volume box") {
        Vector v(2.0, 3.0, 4.0);
        REQUIRE(v.volumeBox() == Approx(2.0 * 3.0 * 4.0));
    }

    SECTION("Rotation around axis") {
        Vector axis(0.0, 0.0, 1.0); // Rotate around z-axis
        Vector v(1.0, 0.0, 0.0);
        Vector rotated = v.rotatedAround(axis, Vector::PI / 2.0); // 90 degrees
        REQUIRE(rotated.x == Approx(0.0).margin(Vector::EPSILON));
        REQUIRE(rotated.y == Approx(1.0).margin(Vector::EPSILON));
        REQUIRE(rotated.z == Approx(0.0).margin(Vector::EPSILON));
    }

    SECTION("VMD string representation") {
        Vector v(1.0, 2.0, 3.0);
        REQUIRE(v.toVMD() == "{1.000000 2.000000 3.000000}");
    }

    SECTION("Mixed product") {
        Vector a(1.0, 0.0, 0.0);
        Vector b(0.0, 1.0, 0.0);
        Vector c(0.0, 0.0, 1.0);
        REQUIRE(Vector::mixedProduct(a, b, c) == Approx(1.0));
    }

    SECTION("Dot product with parallel vectors") {
        Vector v3(2.0, 0.0, 0.0);
        REQUIRE(v1.dot(v3) == Approx(2.0));
    }

    SECTION("Cross product with parallel vectors") {
        Vector v3(2.0, 0.0, 0.0);
        Vector cross = v1.cross(v3);
        REQUIRE(cross == Vector(0.0, 0.0, 0.0));
    }

    SECTION("Volume box with negative components") {
        Vector v(-2.0, 3.0, -4.0);
        REQUIRE(v.volumeBox() == Approx(24.0));
    }

    SECTION("Rotation with zero angle") {
        Vector axis(0.0, 0.0, 1.0);
        Vector v(1.0, 0.0, 0.0);
        Vector rotated = v.rotatedAround(axis, 0.0);
        REQUIRE(rotated == v);
    }

    SECTION("Rotation with 180 degrees") {
        Vector axis(0.0, 0.0, 1.0);
        Vector v(1.0, 0.0, 0.0);
        Vector rotated = v.rotatedAround(axis, Vector::PI);
        REQUIRE(rotated.x == Approx(-1.0).margin(Vector::EPSILON));
        REQUIRE(rotated.y == Approx(0.0).margin(Vector::EPSILON));
        REQUIRE(rotated.z == Approx(0.0).margin(Vector::EPSILON));
    }

    SECTION("VMD string with negative and small values") {
        Vector v(-1.0, Vector::EPSILON, -Vector::EPSILON);
        std::string expected = "{-1.000000 " + std::to_string(Vector::EPSILON) + " " + std::to_string(-Vector::EPSILON) + "}";
        REQUIRE(v.toVMD() == expected);
    }

    SECTION("Mixed product with coplanar vectors") {
        Vector a(1.0, 0.0, 0.0);
        Vector b(0.0, 1.0, 0.0);
        Vector c(2.0, 2.0, 0.0);
        REQUIRE(Vector::mixedProduct(a, b, c) == Approx(0.0));
    }

    SECTION("Rotation with 360 degrees") {
        Vector axis(0.0, 0.0, 1.0);
        Vector v(1.0, 0.0, 0.0);
        Vector rotated = v.rotatedAround(axis, 2.0 * Vector::PI);
        REQUIRE(rotated.x == Approx(1.0).margin(Vector::EPSILON));
        REQUIRE(rotated.y == Approx(0.0).margin(Vector::EPSILON));
        REQUIRE(rotated.z == Approx(0.0).margin(Vector::EPSILON));
    }

    SECTION("Dot product with negative vectors") {
        Vector v_neg(-1.0, 0.0, 0.0);
        REQUIRE(v1.dot(v_neg) == Approx(-1.0));
    }

    SECTION("Mixed product with zero vector") {
        Vector zero(0.0, 0.0, 0.0);
        REQUIRE(Vector::mixedProduct(v1, v2, zero) == Approx(0.0));
    }
}

TEST_CASE("Vector - Distance and Angle Functions", "[Vector]") {
    Vector a(0.0, 0.0, 0.0);
    Vector b(1.0, 0.0, 0.0);
    Vector box(10.0, 10.0, 10.0);

    SECTION("Distance without PBC") {
        REQUIRE(distanceWithoutPBC(a, b) == Approx(1.0));
    }

    SECTION("Distance with PBC") {
        Vector b_pbc(11.0, 0.0, 0.0);
        REQUIRE(distancePBC(a, b_pbc, box) == Approx(1.0));
    }

    SECTION("2D distance") {
        REQUIRE(distance2D(a, b) == Approx(1.0));
    }

    SECTION("Displacement with PBC") {
        Vector b_pbc(11.0, 0.0, 0.0);
        Vector disp = displacementPBC(a, b_pbc, box);
        REQUIRE(disp.x == Approx(-1.0));
        REQUIRE(disp.y == 0.0);
        REQUIRE(disp.z == 0.0);
    }

    SECTION("Squared distance with PBC") {
        Vector b_pbc(11.0, 0.0, 0.0);
        REQUIRE(squaredDistancePBC(a, b_pbc, box) == Approx(1.0));
    }

    SECTION("Angle between vectors in radians") {
        Vector v1(1.0, 0.0, 0.0);
        Vector v2(0.0, 1.0, 0.0);
        REQUIRE(angleBetweenRadians(v1, v2) == Approx(Vector::PI / 2.0));
    }

    SECTION("Angle between vectors in degrees") {
        Vector v1(1.0, 0.0, 0.0);
        Vector v2(0.0, 1.0, 0.0);
        REQUIRE(angleBetweenDegrees(v1, v2) == Approx(90.0));
    }

    SECTION("Angle between three points with PBC") {
        Vector c1(0.0, 0.0, 0.0);
        Vector c2(1.0, 0.0, 0.0);
        Vector c3(1.0, 1.0, 0.0);
        REQUIRE(getAngle(c1, c2, c3, box) == Approx(Vector::PI / 2.0));
        REQUIRE(getAngle(c2, c1, c3, box) == Approx(Vector::PI / 4.0));
    }

    SECTION("VMD line representation") {
        REQUIRE(lineVMD(a, b) == "draw line {0.000000 0.000000 0.000000} {1.000000 0.000000 0.000000}");
    }

    SECTION("Distance with PBC across multiple box lengths") {
        Vector b_pbc(31.0, 0.0, 0.0); // 3 box lengths + 1
        REQUIRE(distancePBC(a, b_pbc, box) == Approx(1.0));
    }

    SECTION("2D distance ignoring z-component") {
        Vector b_z(1.0, 0.0, 100.0);
        REQUIRE(distance2D(a, b_z) == Approx(1.0));
    }

    SECTION("Angle between parallel vectors") {
        Vector v1(1.0, 0.0, 0.0);
        Vector v2(2.0, 0.0, 0.0);
        REQUIRE(angleBetweenRadians(v1, v2) == Approx(0.0));
        REQUIRE(angleBetweenDegrees(v1, v2) == Approx(0.0));
    }

    SECTION("Angle between opposite vectors") {
        Vector v1(1.0, 0.0, 0.0);
        Vector v2(-1.0, 0.0, 0.0);
        REQUIRE(angleBetweenRadians(v1, v2) == Approx(Vector::PI));
        REQUIRE(angleBetweenDegrees(v1, v2) == Approx(180.0));
    }

    SECTION("Angle with zero vector") {
        Vector v1(1.0, 0.0, 0.0);
        Vector v2(-1.0, 0.0, 0.0);
        Vector zero(0.0, 0.0, 0.0);
        REQUIRE(angleBetweenRadians(v1, zero) == Approx(0.0));
        REQUIRE(angleBetweenDegrees(v1, zero) == Approx(0.0));
    }

    SECTION("Angle between three collinear points with PBC") {
        Vector c1(0.0, 0.0, 0.0);
        Vector c2(1.0, 0.0, 0.0);
        Vector c3(2.0, 0.0, 0.0);
        REQUIRE(getAngle(c1, c2, c3, box) == Approx(Vector::PI));
    }

    SECTION("VMD line with negative coordinates") {
        Vector a_neg(-1.0, -2.0, -3.0);
        Vector b_neg(-4.0, -5.0, -6.0);
        REQUIRE(lineVMD(a_neg, b_neg) == "draw line {-1.000000 -2.000000 -3.000000} {-4.000000 -5.000000 -6.000000}");
    }

    SECTION("Distance with PBC in all dimensions") {
        Vector b_pbc(11.0, 12.0, 13.0);
        REQUIRE(distancePBC(a, b_pbc, box) == Approx(std::sqrt(1+4+9)));
    }

    SECTION("Displacement with PBC in all dimensions") {
        Vector b_pbc(11.0, 12.0, 13.0);
        Vector disp = displacementPBC(a, b_pbc, box);
        REQUIRE(disp.x == Approx(-1.0));
        REQUIRE(disp.y == Approx(-2.0));
        REQUIRE(disp.z == Approx(-3.0));
    }

    SECTION("Angle between three points with same position") {
        Vector c(0.0, 0.0, 0.0);
        REQUIRE(std::isnan(getAngle(c, c, c, box))); // Degenerate case
    }
}

TEST_CASE("Toolkit - String and Array Utilities", "[Toolkit]") {
    SECTION("Strip removes spaces") {
        std::string s = "  hello world  ";
        REQUIRE(ToolKit::strip(s) == "helloworld");
    }

    SECTION("Get bin position") {
        REQUIRE(ToolKit::getBinPosition(5.0, 0.0, 10.0, 10) == 5);
        REQUIRE(ToolKit::getBinPosition(0.0, 0.0, 10.0, 10, true) == 0);
        REQUIRE(ToolKit::getBinPosition(11.0, 0.0, 10.0, 10, true) == 9);
        REQUIRE_THROWS_AS(ToolKit::getBinPosition(11.0, 0.0, 10.0, 10, false), std::out_of_range);
    }

    SECTION("Take time measures execution time") {
        auto func = []() { std::this_thread::sleep_for(std::chrono::milliseconds(1000)); };
        long int elapsed = ToolKit::takeTime(func);
        REQUIRE(elapsed >= 1);
    }

    SECTION("Parallel execution") {
        std::vector<int> list = {1, 2, 3};
        std::vector<int> res(3, 0);
        auto func = [](int x, int& r) { r = x * 2; };
        ToolKit::parallel(func, list, res);
        REQUIRE(res[0] == 2);
        REQUIRE(res[1] == 4);
        REQUIRE(res[2] == 6);
    }

    SECTION("Serial execution") {
        std::vector<int> list = {1, 2, 3};
        std::vector<int> res(3, 0);
        auto func = [](int x, int& r) { r = x * 2; };
        ToolKit::serialExecution(func, list, res);
        REQUIRE(res[0] == 2);
        REQUIRE(res[1] == 4);
        REQUIRE(res[2] == 6);
    }

    SECTION("Strip empty string") {
        std::string s = "";
        REQUIRE(ToolKit::strip(s) == "");
    }

    SECTION("Strip string with only spaces") {
        std::string s = "   ";
        REQUIRE(ToolKit::strip(s) == "");
    }

    SECTION("Get bin position with negative values") {
        REQUIRE(ToolKit::getBinPosition(-5.0, -10.0, 0.0, 10) == 5);
        REQUIRE(ToolKit::getBinPosition(-11.0, -10.0, 0.0, 10, true) == 0);
        REQUIRE_THROWS_AS(ToolKit::getBinPosition(-11.0, -10.0, 0.0, 10, false), std::out_of_range);
    }

    SECTION("Get bin position with zero bins") {
        REQUIRE_THROWS_AS(ToolKit::getBinPosition(5.0, 0.0, 10.0, 0), std::out_of_range);
    }

    SECTION("Parallel execution with empty list") {
        std::vector<int> list;
        std::vector<int> res;
        auto func = [](int x, int& r) { r = x * 2; };
        ToolKit::parallel(func, list, res);
        REQUIRE(res.empty());
    }

    SECTION("Serial execution with single element") {
        std::vector<int> list = {1};
        std::vector<int> res(1, 0);
        auto func = [](int x, int& r) { r = x * 2; };
        ToolKit::serialExecution(func, list, res);
        REQUIRE(res[0] == 2);
    }

    SECTION("Strip string with special characters") {
        std::string s = "  hello\n\t world  ";
        REQUIRE(ToolKit::strip(s) == "hello\n\tworld");
    }

    SECTION("Get bin position with equal min and max") {
        REQUIRE_THROWS_AS(ToolKit::getBinPosition(5.0, 10.0, 10.0, 10), std::out_of_range);
    }

    SECTION("Parallel execution with large list") {
        std::vector<int> list(1000);
        std::vector<int> res(1000, 0);
        for (size_t i = 0; i < list.size(); ++i) list[i] = i;
        auto func = [](int x, int& r) { r = x * 2; };
        ToolKit::parallel(func, list, res);
        for (size_t i = 0; i < res.size(); ++i) {
            REQUIRE(res[i] == 2 * i);
        }
    }

    SECTION("Serial execution with large list") {
        std::vector<int> list(1000);
        std::vector<int> res(1000, 0);
        for (size_t i = 0; i < list.size(); ++i) list[i] = i;
        auto func = [](int x, int& r) { r = x * 2; };
        ToolKit::serialExecution(func, list, res);
        for (size_t i = 0; i < res.size(); ++i) {
            REQUIRE(res[i] == 2 * i);
        }
    }
}

TEST_CASE("Sorter - Sorting Functions", "[Sorter]") {
    SECTION("Sort vector ascending") {
        std::vector<int> vec = {3, 1, 4, 1, 5};
        Sorter::sort(vec);
        REQUIRE(vec == std::vector<int>{1, 1, 3, 4, 5});
    }

    SECTION("Sort vector descending") {
        std::vector<int> vec = {3, 1, 4, 1, 5};
        Sorter::sort(vec, Sorter::Order::Descending);
        REQUIRE(vec == std::vector<int>{5, 4, 3, 1, 1});
    }

    SECTION("Sort array ascending") {
        int arr[] = {3, 1, 4, 1, 5};
        Sorter::sort(arr, 5);
        REQUIRE(arr[0] == 1);
        REQUIRE(arr[1] == 1);
        REQUIRE(arr[2] == 3);
        REQUIRE(arr[3] == 4);
        REQUIRE(arr[4] == 5);
    }

    SECTION("Sort array descending") {
        int arr[] = {3, 1, 4, 1, 5};
        Sorter::sort(arr, 5, Sorter::Order::Descending);
        REQUIRE(arr[0] == 5);
        REQUIRE(arr[1] == 4);
        REQUIRE(arr[2] == 3);
        REQUIRE(arr[3] == 1);
        REQUIRE(arr[4] == 1);
    }

    SECTION("Cosort two vectors") {
        std::vector<double> values = {3.0, 1.0, 2.0};
        std::vector<int> indexes = {0, 1, 2};
        Sorter::cosort(values, indexes);
        REQUIRE(values == std::vector<double>{1.0, 2.0, 3.0});
        REQUIRE(indexes == std::vector<int>{1, 2, 0});
    }

    SECTION("Cosort two vectors descending") {
        std::vector<double> values = {3.0, 1.0, 2.0};
        std::vector<int> indexes = {0, 1, 2};
        Sorter::cosort(values, indexes, Sorter::Order::Descending);
        REQUIRE(values == std::vector<double>{3.0, 2.0, 1.0});
        REQUIRE(indexes == std::vector<int>{0, 2, 1});
    }

    SECTION("Cosort throws on unequal sizes") {
        std::vector<double> values = {1.0, 2.0};
        std::vector<int> indexes = {0, 1, 2};
        REQUIRE_THROWS_AS(Sorter::cosort(values, indexes), std::invalid_argument);
    }

    SECTION("Sort empty vector") {
        std::vector<int> vec;
        Sorter::sort(vec);
        REQUIRE(vec.empty());
    }

    SECTION("Sort array with zero size") {
        int arr[0];
        Sorter::sort(arr, 0);
        // No crash, no change needed
    }

    SECTION("Sort vector with single element") {
        std::vector<int> vec = {42};
        Sorter::sort(vec);
        REQUIRE(vec == std::vector<int>{42});
    }

    SECTION("Cosort with equal values") {
        std::vector<double> values = {2.0, 2.0, 2.0};
        std::vector<int> indexes = {0, 1, 2};
        Sorter::cosort(values, indexes);
        REQUIRE(values == std::vector<double>{2.0, 2.0, 2.0});
        REQUIRE(indexes == std::vector<int>{0, 1, 2}); // Stable sort
    }

    SECTION("Cosort empty vectors") {
        std::vector<double> values;
        std::vector<int> indexes;
        Sorter::cosort(values, indexes);
        REQUIRE(values.empty());
        REQUIRE(indexes.empty());
    }

    SECTION("Sort vector with negative values") {
        std::vector<int> vec = {-3, -1, -4, -1, -5};
        Sorter::sort(vec);
        REQUIRE(vec == std::vector<int>{-5, -4, -3, -1, -1});
    }

    SECTION("Sort large array") {
        std::vector<int> arr(1000);
        for (size_t i = 0; i < arr.size(); ++i) arr[i] = arr.size() - i - 1;
        Sorter::sort(arr.data(), arr.size());
        for (size_t i = 0; i < arr.size(); ++i) {
            REQUIRE(arr[i] == static_cast<int>(i));
        }
    }

    SECTION("Cosort with floating-point edge cases") {
        std::vector<double> values = {std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), 0.0};
        std::vector<int> indexes = {0, 1, 2};
        Sorter::cosort(values, indexes);
        REQUIRE((std::isinf(values[0]) && values[0] < 0));
        REQUIRE(values[1] == Approx(0.0));
        REQUIRE((std::isinf(values[2]) && values[2] > 0));
    }
}

TEST_CASE("Geometrics - Tetrahedron and Plane Calculations", "[Geometrics]") {
    Vector O(0.0, 0.0, 0.0);
    Vector H1(1.0, 0.0, 0.0);
    Vector H2(0.0, 1.0, 0.0);
    Vector bounds(10.0, 10.0, 10.0);

    SECTION("Check PBC adjusts coordinates") {
        Vector H1_pbc(11.0, 0.0, 0.0);
        Vector* adjusted = Geometrics::checkPBC(O, H1_pbc, H2, bounds);
        REQUIRE(adjusted[0] == O);
        REQUIRE(adjusted[1].x == Approx(1.0));
        REQUIRE(adjusted[2] == H2);
        delete[] adjusted;
    }

    SECTION("Determinant 3x3") {
        Real matrix[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        REQUIRE(Geometrics::determinant_3x3(matrix) == Approx(0.0));
    }

    SECTION("Cramer's Rule") {
        Real matrix[3][3] = {{2, 0, 0}, {0, 3, 0}, {0, 0, 4}};
        Real independents[3] = {2, 3, 4};
        Vector result = Geometrics::CramersRule(matrix, independents);
        REQUIRE(result.x == Approx(1.0));
        REQUIRE(result.y == Approx(1.0));
        REQUIRE(result.z == Approx(1.0));
    }

    SECTION("Perfect tetrahedron calculation") {
        Geometrics::TetrahedronVertices tet = Geometrics::getPerfectTetrahedron(O, H1, H2, bounds);
        Vector* vertices = tet.toArray();
        REQUIRE(vertices[0].magnitude() > 0.0);
        REQUIRE(vertices[1].magnitude() > 0.0);
        REQUIRE(vertices[2].magnitude() > 0.0);
        REQUIRE(vertices[3].magnitude() > 0.0);
        delete[] vertices;
    }

    SECTION("Plane from points") {
        Vector a(0, 0, 0), b(1, 0, 0), c(0, 1, 0);
        Geometrics::Plane plane = Geometrics::getPlaneFromPoints(a, b, c, bounds);
        REQUIRE(plane.normal.z == Approx(1.0));
        REQUIRE(plane.independent_term == Approx(0.0));
    }

    SECTION("Distance to plane") {
        Geometrics::Plane plane(Vector(0, 0, 1), 0);
        Vector p(0, 0, 5);
        REQUIRE(Geometrics::distanceToPlane(plane, p) == Approx(5.0));
    }

    SECTION("Distance to line") {
        Geometrics::Line line(Vector(1, 0, 0), Vector(0, 0, 0));
        Vector p(0, 1, 0);
        REQUIRE(Geometrics::distanceToLine(line, p, bounds) == Approx(1.0));
    }

    SECTION("Is in box") {
        Geometrics::Box box(Vector(0, 0, 0), Vector(1, 1, 1));
        Vector p(0.5, 0.5, 0.5);
        REQUIRE(Geometrics::isInBox(p, box));
        Vector p_out(2, 2, 2);
        REQUIRE_FALSE(Geometrics::isInBox(p_out, box));
    }

    SECTION("SphereList management") {
        Geometrics::SphereList spheres(3);
        spheres.centers[0] = Vector(1, 0, 0);
        spheres.centers[1] = Vector(0, 1, 0);
        spheres.centers[2] = Vector(0, 0, 1);
        spheres.size = 3;

        Geometrics::SphereList copy(spheres);
        REQUIRE(copy.size == 3);
        REQUIRE(copy.centers[0] == Vector(1, 0, 0));
        REQUIRE(copy.centers[1] == Vector(0, 1, 0));
        REQUIRE(copy.centers[2] == Vector(0, 0, 1));

        Geometrics::SphereList assigned;
        assigned = spheres;
        REQUIRE(assigned.size == 3);
        REQUIRE(assigned.centers[0] == Vector(1, 0, 0));
    }

    SECTION("Check PBC with all vectors outside box") {
        Vector O_pbc(5.0, 5.0, 5.0);
        Vector H1_pbc(16.0, 15.0, 15.0);
        Vector H2_pbc(15.0, 16.0, 15.0);
        Vector* adjusted = Geometrics::checkPBC(O_pbc, H1_pbc, H2_pbc, bounds);
        REQUIRE(adjusted[0] == O_pbc);
        REQUIRE(adjusted[1].x == Approx(6.0));
        REQUIRE(adjusted[1].y == Approx(5.0));
        REQUIRE(adjusted[1].z == Approx(5.0));
        REQUIRE(adjusted[2].x == Approx(5.0));
        REQUIRE(adjusted[2].y == Approx(6.0));
        REQUIRE(adjusted[2].z == Approx(5.0));
        delete[] adjusted;
    }

    SECTION("Determinant of identity matrix") {
        Real matrix[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        REQUIRE(Geometrics::determinant_3x3(matrix) == Approx(1.0));
    }

    SECTION("Cramer's Rule with zero determinant") {
        Real matrix[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        Real independents[3] = {1, 2, 3};
        Vector result = Geometrics::CramersRule(matrix, independents);
        REQUIRE((std::isnan(result.x) || std::isinf(result.x))); // Expect NaN or inf due to division by zero
    }

    SECTION("Perfect tetrahedron with near-collinear points") {
        Vector H2_near(1.0, Vector::EPSILON, 0.0);
        Geometrics::TetrahedronVertices tet = Geometrics::getPerfectTetrahedron(O, H1, H2_near, bounds);
        Vector* vertices = tet.toArray();
        REQUIRE(vertices[0].magnitude() > 0.0);
        REQUIRE(vertices[1].magnitude() > 0.0);
        REQUIRE(vertices[2].magnitude() > 0.0);
        REQUIRE(vertices[3].magnitude() > 0.0);
        delete[] vertices;
    }

    SECTION("Distance to plane with point on plane") {
        Geometrics::Plane plane(Vector(0, 0, 1), 0);
        Vector p(0, 0, 0);
        REQUIRE(Geometrics::distanceToPlane(plane, p) == Approx(0.0));
    }

    SECTION("Distance to line with point on line") {
        Geometrics::Line line(Vector(1, 0, 0), Vector(0, 0, 0));
        Vector p(2, 0, 0);
        REQUIRE(Geometrics::distanceToLine(line, p, bounds) == Approx(0.0));
    }

    SECTION("Is in box with boundary points") {
        Geometrics::Box box(Vector(0, 0, 0), Vector(1, 1, 1));
        Vector p_min(0, 0, 0);
        Vector p_max(1, 1, 1);
        REQUIRE(Geometrics::isInBox(p_min, box));
        REQUIRE(Geometrics::isInBox(p_max, box));
    }

    SECTION("SphereList with zero size") {
        Geometrics::SphereList spheres(0);
        REQUIRE(spheres.size == 0);
        REQUIRE(spheres.centers == nullptr);
        Geometrics::SphereList copy(spheres);
        REQUIRE(copy.size == 0);
        REQUIRE(copy.centers == nullptr);
    }

    SECTION("Check PBC with same point") {
        Vector* adjusted = Geometrics::checkPBC(O, O, O, bounds);
        REQUIRE(adjusted[0] == O);
        REQUIRE(adjusted[1] == O);
        REQUIRE(adjusted[2] == O);
        delete[] adjusted;
    }

    SECTION("Plane with zero normal") {
        Vector a(0, 0, 0), b(0, 0, 0), c(0, 0, 0);
        Geometrics::Plane plane = Geometrics::getPlaneFromPoints(a, b, c, bounds);
        REQUIRE(plane.normal == Vector(0, 0, 0));
        REQUIRE(plane.independent_term == Approx(0.0));
    }

    SECTION("Distance to plane with zero normal") {
        Geometrics::Plane plane(Vector(0, 0, 0), 0);
        Vector p(1, 1, 1);
        REQUIRE(Geometrics::distanceToPlane(plane, p) == Approx(0.0));
    }

    SECTION("Is in box with negative bounds") {
        Geometrics::Box box(Vector(-1, -1, -1), Vector(0, 0, 0));
        Vector p(-0.5, -0.5, -0.5);
        REQUIRE(Geometrics::isInBox(p, box));
        Vector p_out(-2, -2, -2);
        REQUIRE_FALSE(Geometrics::isInBox(p_out, box));
    }
}

TEST_CASE("Vector - Constants and Edge Cases", "[Vector]") {
    SECTION("Verify PI and RAD2DEG constants") {
        REQUIRE(Vector::PI == Approx(3.14159265358979323846));
        REQUIRE(Vector::RAD2DEG == Approx(180.0 / Vector::PI));
    }

    SECTION("Verify EPSILON constant") {
        Vector v1(1.0, 0.0, 0.0);
        Vector v2(1.0 + Vector::EPSILON / 2.0, 0.0, 0.0);
        REQUIRE(v1 == v2); // Within EPSILON tolerance
    }

    SECTION("Division by zero scalar") {
        Vector v(1.0, 2.0, 3.0);
        Vector divided = v / 0.0;
        REQUIRE((std::isinf(divided.x) || std::isnan(divided.x)));
        REQUIRE((std::isinf(divided.y) || std::isnan(divided.y)));
        REQUIRE((std::isinf(divided.z) || std::isnan(divided.z)));
    }

    SECTION("Rotation with zero axis") {
        Vector zero_axis(0.0, 0.0, 0.0);
        Vector v(1.0, 0.0, 0.0);
        Vector rotated = v.rotatedAround(zero_axis, Vector::PI / 2.0);
        REQUIRE(rotated == v); // No rotation due to zero axis
    }

    SECTION("Compound assignment with zero scalar") {
        Vector v(1.0, 2.0, 3.0);
        v *= 0.0;
        REQUIRE(v == Vector(0.0, 0.0, 0.0));
        v = Vector(1.0, 2.0, 3.0);
        v /= std::numeric_limits<float>::infinity();
        REQUIRE(v.x == Approx(0.0).margin(Vector::EPSILON));
        REQUIRE(v.y == Approx(0.0).margin(Vector::EPSILON));
        REQUIRE(v.z == Approx(0.0).margin(Vector::EPSILON));
    }
}

TEST_CASE("Toolkit - Structs", "[Toolkit]") {
    SECTION("ArrInt initialization and access") {
        ToolKit::ArrInt arr;
        arr.size = 3;
        arr.arr = new int[3]{1, 2, 3};
        REQUIRE(arr.arr[0] == 1);
        REQUIRE(arr.arr[1] == 2);
        REQUIRE(arr.arr[2] == 3);
        delete[] arr.arr;
    }

    SECTION("ArrFloat initialization and access") {
        ToolKit::ArrFloat arr;
        arr.size = 2;
        arr.arr = new Real[2]{1.5, 2.5};
        REQUIRE(arr.arr[0] == Approx(1.5));
        REQUIRE(arr.arr[1] == Approx(2.5));
        delete[] arr.arr;
    }

    SECTION("FlaggedArrFloat construction") {
        Real* data = new Real[2]{1.0, 2.0};
        ToolKit::FlaggedArrFloat arr(2, data, true);
        REQUIRE(arr.size == 2);
        REQUIRE(arr.arr[0] == Approx(1.0));
        REQUIRE(arr.arr[1] == Approx(2.0));
        REQUIRE(arr.flag == true);
        delete[] data;
    }

    SECTION("ArrInt with zero size") {
        ToolKit::ArrInt arr;
        arr.size = 0;
        arr.arr = nullptr;
        REQUIRE(arr.size == 0);
        REQUIRE(arr.arr == nullptr);
    }

    SECTION("FlaggedArrFloat with negative values") {
        Real* data = new Real[2]{-1.0, -2.0};
        ToolKit::FlaggedArrFloat arr(2, data, false);
        REQUIRE(arr.size == 2);
        REQUIRE(arr.arr[0] == Approx(-1.0));
        REQUIRE(arr.arr[1] == Approx(-2.0));
        REQUIRE(arr.flag == false);
        delete[] data;
    }
}
