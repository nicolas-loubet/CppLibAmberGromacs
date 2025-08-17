#include "catch.hpp"
#include "../cppambergromacs/CppLibAmberGromacs.hpp"

#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>

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
    }

    SECTION("VMD line representation") {
        REQUIRE(lineVMD(a, b) == "draw line {0.000000 0.000000 0.000000} {1.000000 0.000000 0.000000}");
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
}

/// Helper: reads a file and returns its contents
string readFile(const string& fname) {
    ifstream in(fname);
    REQUIRE(in.is_open());
    stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

TEST_CASE("CSVWriter", "[CSVWriter]") {
    SECTION("writeHeader and writeRow") {
        string fname = "test_output_basic.csv";
        {
            CSVWriter csv(fname);
            csv.writeHeader({"A","B","C"});
            csv.writeRow(1, 2.5, "hola");
            csv.writeRow(2, 3.1416, "chau");
        }
        string expected =
            "A,B,C\n"
            "1,2.5,hola\n"
            "2,3.1416,chau\n";
        REQUIRE(readFile(fname) == expected);
        REQUIRE(std::remove(fname.c_str()) == 0);
    }

    SECTION("writeDistribution without totals") {
        string fname = "test_output_dist.csv";
        int N_BINS = 3;
        Real MIN_BINS = 0.0, MAX_BINS = 3.0;

        Real* col1 = new Real[N_BINS]{1.0, 1.0, 2.0}; // sum = 4
        Real* col2 = new Real[N_BINS]{0.0, 2.0, 2.0}; // sum = 4

        {
            CSVWriter csv(fname);
            csv.writeDistribution({col1, col2}, N_BINS, MIN_BINS, MAX_BINS, {"x","dist1","dist2"});
        }

        string expected =
            "x,dist1,dist2\n"
            "0.5,0.25,0\n"
            "1.5,0.25,0.5\n"
            "2.5,0.5,0.5\n";
        REQUIRE(readFile(fname) == expected);

        REQUIRE(std::remove(fname.c_str()) == 0);
        delete[] col1;
        delete[] col2;
    }

    SECTION("writeDistribution with explicit totals") {
        string fname = "test_output_dist_totals.csv";
        int N_BINS = 2;
        Real MIN_BINS = 0.0, MAX_BINS = 2.0;

        Real* col1 = new Real[N_BINS]{1.0, 3.0};
        Real* col2 = new Real[N_BINS]{2.0, 2.0};

        vector<Real> totals = {10.0, 5.0}; // normalizadores externos

        {
            CSVWriter csv(fname);
            csv.writeDistribution({col1,col2}, N_BINS, MIN_BINS, MAX_BINS, {"x","dist1","dist2"}, totals);
        }

        string expected =
            "x,dist1,dist2\n"
            "0.5,0.1,0.4\n"
            "1.5,0.3,0.4\n";
        REQUIRE(readFile(fname) == expected);

        REQUIRE(std::remove(fname.c_str()) == 0);
        delete[] col1;
        delete[] col2;
    }

    SECTION("writeTable generic") {
        string fname = "test_output_table.csv";

        int nRows = 2, nCols = 2;
        double** table = new double*[nRows];
        for(int i=0; i<nRows; i++) {
            table[i] = new double[nCols];
            for(int j=0; j<nCols; j++) table[i][j] = i*10+j;
        }
        vector<string> rows = {"row1","row2"};
        vector<string> cols = {"Name","C1","C2"};

        {
            CSVWriter csv(fname);
            csv.writeTable(rows, table, nRows, nCols, cols);
        }

        string expected =
            "Name,C1,C2\n"
            "row1,0,1\n"
            "row2,10,11\n";
        REQUIRE(readFile(fname) == expected);

        for(int i=0; i<nRows; i++) delete[] table[i];
        delete[] table;
        REQUIRE(std::remove(fname.c_str()) == 0);
    }
}
