#include <iostream>
#include <vector>
#include <array>
#include <unordered_map>
#include <random>
#include <cmath>
#include <chrono>

using namespace std;

template <size_t Dim>
inline double dot_product(const array<double, Dim>& a, const array<double, Dim>& b) {
    double result = 0.0;
    for (size_t i = 0; i < Dim; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

template <size_t Dim>
inline int hash_function(const array<double, Dim>& point, const array<double, Dim>& a, double b, double w) {
    double dot_product_result = dot_product(point, a);
    return static_cast<int>((dot_product_result + b) / w);
}

template <size_t Dim>
double compute_avg_distance(const vector<array<double, Dim>>& points) {
    double sum = 0.0;
    int count = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            double dist = 0.0;
            for (size_t k = 0; k < Dim; ++k) {
                dist += (points[i][k] - points[j][k]) * (points[i][k] - points[j][k]);
            }
            sum += sqrt(dist);
            count++;
        }
    }
    return count > 0 ? sum / count : 1.0;
}

template <size_t Dim>
class LSH {
public:
    LSH(int num_tables, int hash_size, double w, mt19937& rng)
            : num_tables(num_tables), hash_size(hash_size), w(w), rng(rng) {
        for (int i = 0; i < num_tables; ++i) {
            vector<array<double, Dim>> table_hashes;
            vector<double> table_offsets;
            for (int j = 0; j < hash_size; ++j) {
                array<double, Dim> hash;
                for (size_t k = 0; k < Dim; ++k) {
                    uniform_real_distribution<double> dist(-1.0, 1.0);
                    hash[k] = dist(rng);
                }
                table_hashes.push_back(hash);
                uniform_real_distribution<double> dist(0.0, w);
                table_offsets.push_back(dist(rng));
            }
            hash_tables.push_back(table_hashes);
            hash_offsets.push_back(table_offsets);
            tables.push_back(unordered_map<int, vector<array<double, Dim>>>());
        }
        cout << "[ОТЛАДКА] LSH инициализирован с " << num_tables << " таблицами и " << hash_size << " хеш-функциями на таблицу, w=" << w << endl;
    }

    void add_point(const array<double, Dim>& point) {
        for (int i = 0; i < num_tables; ++i) {
            int hash_code = 0;
            for (int j = 0; j < hash_size; ++j) {
                hash_code = hash_code * 31 + hash_function(point, hash_tables[i][j], hash_offsets[i][j], w);
            }
            tables[i][hash_code].push_back(point);
        }
        cout << "[ОТЛАДКА] Точка добавлена в LSH." << endl;
    }

    array<double, Dim> find_nearest(const array<double, Dim>& query, const vector<array<double, Dim>>& points) {
        array<double, Dim> nearest = points[0];
        double min_dist = numeric_limits<double>::max();

        cout << "[ОТЛАДКА] Начало поиска ближайшего соседа..." << endl;

        for (int i = 0; i < num_tables; ++i) {
            int hash_code = 0;
            for (int j = 0; j < hash_size; ++j) {
                hash_code = hash_code * 31 + hash_function(query, hash_tables[i][j], hash_offsets[i][j], w);
            }
            cout << "[ОТЛАДКА] Таблица " << i << ", хеш-код: " << hash_code << endl;

            for (int offset = -1; offset <= 1; ++offset) {
                int check_hash = hash_code + offset;
                if (tables[i].find(check_hash) != tables[i].end()) {
                    cout << "[ОТЛАДКА] Найдено " << tables[i][check_hash].size() << " точек в ведре с хеш-кодом " << check_hash << endl;
                    for (const auto& point : tables[i][check_hash]) {
                        double dist = 0.0;
                        for (size_t k = 0; k < Dim; ++k) {
                            dist += (point[k] - query[k]) * (point[k] - query[k]);
                        }
                        cout << "[ОТЛАДКА] Расстояние до точки (";
                        for (size_t k = 0; k < Dim; ++k) {
                            cout << point[k];
                            if (k < Dim - 1) cout << ", ";
                        }
                        cout << "): " << dist << endl;

                        if (dist < min_dist) {
                            min_dist = dist;
                            nearest = point;
                        }
                    }
                } else {
                    cout << "[ОТЛАДКА] В ведре с хеш-кодом " << check_hash << " не найдено точек." << endl;
                }
            }
        }

        cout << "[ОТЛАДКА] Поиск ближайшего соседа завершен." << endl;
        return nearest;
    }

private:
    int num_tables;
    int hash_size;
    double w;
    mt19937& rng;
    vector<vector<array<double, Dim>>> hash_tables;
    vector<vector<double>> hash_offsets;
    vector<unordered_map<int, vector<array<double, Dim>>>> tables;
};

// Функция для 3D тестов
void run_test_3d(const string& test_name, const vector<array<double, 3>>& points, const array<double, 3>& query, const array<double, 3>& expected) {
    mt19937 rng(42);
    double w = compute_avg_distance(points);

    cout << "\n=== Тест: " << test_name << " ===\n";
    LSH<3> lsh(5, 4, w, rng);

    for (const auto& point : points) {
        lsh.add_point(point);
    }

    auto start = chrono::high_resolution_clock::now();
    array<double, 3> nearest = lsh.find_nearest(query, points);
    auto end = chrono::high_resolution_clock::now();

    cout << "Точка запроса: (";
    for (size_t i = 0; i < 3; ++i) {
        cout << query[i];
        if (i < 2) cout << ", ";
    }
    cout << ")" << endl;

    cout << "Ожидаемый результат: (";
    for (size_t i = 0; i < 3; ++i) {
        cout << expected[i];
        if (i < 2) cout << ", ";
    }
    cout << ")" << endl;

    cout << "Найденная точка: (";
    for (size_t i = 0; i < 3; ++i) {
        cout << nearest[i];
        if (i < 2) cout << ", ";
    }
    cout << ")" << endl;

    cout << "Время поиска: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << " микросекунд." << endl;

    if (test_name == "Точки в разных плоскостях (3D)") {
        vector<array<double, 3>> valid_results = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
        bool passed = false;
        for (const auto& valid : valid_results) {
            if (nearest == valid) {
                passed = true;
                break;
            }
        }
        if (passed) {
            cout << "Тест пройден! (Найдена одна из равноудалённых точек)" << endl;
        } else {
            cout << "Тест не пройден!" << endl;
        }
    } else {
        if (nearest == expected) {
            cout << "Тест пройден!" << endl;
        } else {
            cout << "Тест не пройден!" << endl;
        }
    }
}

// Функция для 2D тестов
void run_test_2d(const string& test_name, const vector<array<double, 2>>& points, const array<double, 2>& query, const array<double, 2>& expected) {
    mt19937 rng(42);
    double w = compute_avg_distance(points);

    cout << "\n=== Тест: " << test_name << " ===\n";
    LSH<2> lsh(5, 4, w, rng);

    for (const auto& point : points) {
        lsh.add_point(point);
    }

    auto start = chrono::high_resolution_clock::now();
    array<double, 2> nearest = lsh.find_nearest(query, points);
    auto end = chrono::high_resolution_clock::now();

    cout << "Точка запроса: (";
    for (size_t i = 0; i < 2; ++i) {
        cout << query[i];
        if (i < 1) cout << ", ";
    }
    cout << ")" << endl;

    cout << "Ожидаемый результат: (";
    for (size_t i = 0; i < 2; ++i) {
        cout << expected[i];
        if (i < 1) cout << ", ";
    }
    cout << ")" << endl;

    cout << "Найденная точка: (";
    for (size_t i = 0; i < 2; ++i) {
        cout << nearest[i];
        if (i < 1) cout << ", ";
    }
    cout << ")" << endl;

    cout << "Время поиска: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << " микросекунд." << endl;

    if (nearest == expected) {
        cout << "Тест пройден!" << endl;
    } else {
        cout << "Тест не пройден!" << endl;
    }
}

int main() {
    // 3D тесты
    vector<array<double, 3>> points1 = {{{1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}, {3.0, 3.0, 3.0}}};
    array<double, 3> query1 = {1.5, 1.5, 1.5};
    array<double, 3> expected1 = {1.0, 1.0, 1.0};
    run_test_3d("Точки на одной прямой (3D)", points1, query1, expected1);

    vector<array<double, 3>> points2 = {{{-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}}};
    array<double, 3> query2 = {0.0, 0.0, 0.0};
    array<double, 3> expected2 = {1.0, 1.0, 1.0};
    run_test_3d("Точки в разных квадрантах (3D)", points2, query2, expected2);

    vector<array<double, 3>> points3 = {{{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}}};
    array<double, 3> query3 = {1.5, 1.5, 1.5};
    array<double, 3> expected3 = {1.0, 1.0, 1.0};
    run_test_3d("Точки с одинаковыми координатами (3D)", points3, query3, expected3);

    vector<array<double, 3>> points4 = {{{0.5, 0.5, 0.5}, {1.5, 1.5, 1.5}, {2.5, 2.5, 2.5}}};
    array<double, 3> query4 = {1.0, 1.0, 1.0};
    array<double, 3> expected4 = {0.5, 0.5, 0.5};
    run_test_3d("Точки в случайных позициях (3D)", points4, query4, expected4);

    vector<array<double, 3>> points5 = {{{-1.0, -1.0, -1.0}, {-2.0, -2.0, -2.0}, {-3.0, -3.0, -3.0}}};
    array<double, 3> query5 = {-1.5, -1.5, -1.5};
    array<double, 3> expected5 = {-1.0, -1.0, -1.0};
    run_test_3d("Точки с отрицательными координатами (3D)", points5, query5, expected5);

    vector<array<double, 3>> points6 = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    array<double, 3> query6 = {0.5, 0.5, 0.5};
    array<double, 3> expected6 = {1.0, 0.0, 0.0};
    run_test_3d("Точки в разных плоскостях (3D)", points6, query6, expected6);

    // 2D тесты
    vector<array<double, 2>> points2d1 = {{{1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}};
    array<double, 2> query2d1 = {1.5, 1.5};
    array<double, 2> expected2d1 = {1.0, 1.0};
    run_test_2d("Точки на одной прямой (2D)", points2d1, query2d1, expected2d1);

    vector<array<double, 2>> points2d2 = {{{-1.0, -1.0}, {1.0, 1.0}, {2.0, 2.0}}};
    array<double, 2> query2d2 = {0.0, 0.0};
    array<double, 2> expected2d2 = {-1.0, -1.0};
    run_test_2d("Точки в разных квадрантах (2D)", points2d2, query2d2, expected2d2);

    vector<array<double, 2>> points2d3 = {{{1.0, 1.0}, {1.0, 1.0}, {2.0, 2.0}}};
    array<double, 2> query2d3 = {1.5, 1.5};
    array<double, 2> expected2d3 = {1.0, 1.0};
    run_test_2d("Точки с одинаковыми координатами (2D)", points2d3, query2d3, expected2d3);

    vector<array<double, 2>> points2d4 = {{{0.5, 0.5}, {1.5, 1.5}, {2.5, 2.5}}};
    array<double, 2> query2d4 = {1.0, 1.0};
    array<double, 2> expected2d4 = {0.5, 0.5};
    run_test_2d("Точки в случайных позициях (2D)", points2d4, query2d4, expected2d4);

    vector<array<double, 2>> points2d5 = {{{-1.0, -1.0}, {-2.0, -2.0}, {-3.0, -3.0}}};
    array<double, 2> query2d5 = {-1.5, -1.5};
    array<double, 2> expected2d5 = {-1.0, -1.0};
    run_test_2d("Точки с отрицательными координатами (2D)", points2d5, query2d5, expected2d5);

    return 0;
}