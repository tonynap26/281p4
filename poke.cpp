// Project Identifier: 5949F553E20B650AB0FB2266D3C0822B13D248B0

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <limits>
#include <getopt.h>
#include <iomanip>
#include <algorithm>
#include <cstddef>

using namespace std;

enum class Area { LAND, SEA, COASTLINE };

void handleMST(const vector<pair<int, int>>& vertices);
void handleFASTTSP(const vector<pair<int, int>>& vertices);
void handleOPTTSP(const vector<pair<int, int>>& vertices);
inline double calculateDistance(const pair<int, int>& p1, const pair<int, int>& p2);
Area getArea(const pair<int, int>& p);

void print_help() {
    cout << "Usage: poke [options]\n";
    cout << "Options:\n";
    cout << "  -m, --mode {MST|FASTTSP|OPTTSP}  Specify the mode\n";
    cout << "  -h, --help                       Show this help message\n";
}

int main(int argc, char* argv[]) {
    string mode;
    static struct option long_options[] = {
        {"mode", required_argument, 0, 'm'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "m:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'm':
                mode = optarg;
                break;
            case 'h':
                print_help();
                return 0;
            case '?':
            default:
                cerr << "Error: Invalid command line option\n";
                return 1;
        }
    }

    if (mode.empty()) {
        cerr << "Error: No mode specified\n";
        return 1;
    }

    if (mode != "MST" && mode != "FASTTSP" && mode != "OPTTSP") {
        cerr << "Error: Invalid mode\n";
        return 1;
    }

    cout << fixed << setprecision(2);

    size_t n;
    cin >> n;
    vector<pair<int, int>> vertices(n);
    for (size_t i = 0; i < n; ++i) {
        cin >> vertices[i].first >> vertices[i].second;
    }

    if (mode == "MST") {
        handleMST(vertices);
    } else if (mode == "FASTTSP") {
        handleFASTTSP(vertices);
    } else if (mode == "OPTTSP") {
        handleOPTTSP(vertices);
    }

    return 0;
}

inline double calculateDistance(const pair<int, int>& p1, const pair<int, int>& p2) {
    double dx = static_cast<double>(p1.first) - static_cast<double>(p2.first);
    double dy = static_cast<double>(p1.second) - static_cast<double>(p2.second);
    return sqrt(dx * dx + dy * dy);
}

Area getArea(const pair<int, int>& p) {
    int x = p.first;
    int y = p.second;

    if (x == 0 && y == 0) {
        return Area::COASTLINE;
    }
    if (x == 0 || y == 0) {
        if (x <= 0 && y <= 0) {
            return Area::COASTLINE;
        }
        return Area::LAND;
    }
    if (x > 0 || y > 0) {
        return Area::LAND;
    }
    if (x < 0 && y < 0) {
        return Area::SEA;
    }
    return Area::LAND;
}

void handleMST(const vector<pair<int, int>>& vertices) {
    size_t n = vertices.size();
    vector<Area> areas(n);
    for (size_t i = 0; i < n; ++i) {
        areas[i] = getArea(vertices[i]);
    }

    vector<bool> inMST(n, false);
    vector<double> minEdge(n, numeric_limits<double>::infinity());
    vector<size_t> parent(n, numeric_limits<size_t>::max());

    minEdge[0] = 0.0;
    double totalWeight = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double minDist = numeric_limits<double>::infinity();
        size_t u = n;

        for (size_t v = 0; v < n; ++v) {
            if (!inMST[v] && minEdge[v] < minDist) {
                minDist = minEdge[v];
                u = v;
            }
        }

        if (u == n) {
            cerr << "Cannot construct MST\n";
            exit(1);
        }

        inMST[u] = true;
        totalWeight += minEdge[u];

        for (size_t v = 0; v < n; ++v) {
            if (!inMST[v]) {
                bool edgeAllowed = false;
                if (areas[u] == areas[v]) {
                    edgeAllowed = true;
                } else if (areas[u] == Area::COASTLINE || areas[v] == Area::COASTLINE) {
                    edgeAllowed = true;
                }

                if (edgeAllowed) {
                    double dist = calculateDistance(vertices[u], vertices[v]);
                    if (dist < minEdge[v]) {
                        minEdge[v] = dist;
                        parent[v] = u;
                    }
                }
            }
        }
    }

    vector<pair<size_t, size_t>> edges;
    for (size_t v = 0; v < n; ++v) {
        if (parent[v] != numeric_limits<size_t>::max()) {
            size_t u = parent[v];
            if (u < v) {
                edges.emplace_back(u, v);
            } else {
                edges.emplace_back(v, u);
            }
        }
    }

    sort(edges.begin(), edges.end());

    cout << totalWeight << "\n";

    for (const auto& edge : edges) {
        cout << edge.first << " " << edge.second << "\n";
    }
}

void handleFASTTSP(const vector<pair<int, int>>& vertices) {
    size_t n = vertices.size();
    vector<bool> visited(n, false);
    vector<size_t> tour;
    double totalWeight = 0.0;

    size_t current = 0;
    visited[current] = true;
    tour.push_back(current);

    for (size_t i = 1; i < n; ++i) {
        double minDist = numeric_limits<double>::infinity();
        size_t next = n;

        for (size_t j = 0; j < n; ++j) {
            if (!visited[j]) {
                double dist = calculateDistance(vertices[current], vertices[j]);
                if (dist < minDist) {
                    minDist = dist;
                    next = j;
                }
            }
        }

        if (next == n) break;

        visited[next] = true;
        tour.push_back(next);
        totalWeight += minDist;
        current = next;
    }

    totalWeight += calculateDistance(vertices[current], vertices[0]);

    bool improved = true;
    int iterations = 0;
    const int MAX_ITER = 100;

    while (improved && iterations < MAX_ITER) {
        improved = false;
        for (size_t i = 1; i < n - 1; ++i) {
            bool local_improved = false;
            for (size_t k = i + 1; k < n; ++k) {
                double delta = -calculateDistance(vertices[tour[i - 1]], vertices[tour[i]])
                               -calculateDistance(vertices[tour[k]], vertices[tour[(k + 1) % n]])
                               +calculateDistance(vertices[tour[i - 1]], vertices[tour[k]])
                               +calculateDistance(vertices[tour[i]], vertices[tour[(k + 1) % n]]);

                if (delta < -1e-6) {
                    reverse(tour.begin() + static_cast<ptrdiff_t>(i),
                            tour.begin() + static_cast<ptrdiff_t>(k) + 1);
                    totalWeight += delta;
                    improved = true;
                    local_improved = true;
                    break;
                }
            }
            if (local_improved)
                break;
        }
        iterations++;
    }

    if (tour[0] != 0) {
        auto it = find(tour.begin(), tour.end(), 0);
        if (it != tour.end()) {
            rotate(tour.begin(), it, tour.end());
        }
    }

    totalWeight = 0.0;
    for (size_t i = 0; i < n; ++i) {
        size_t from = tour[i];
        size_t to = tour[(i + 1) % n];
        totalWeight += calculateDistance(vertices[from], vertices[to]);
    }

    cout << totalWeight << "\n";

    for (size_t i = 0; i < n; ++i) {
        cout << tour[i] << " ";
    }
    cout << "\n";
}

void handleOPTTSP(const vector<pair<int, int>>& vertices) {
    size_t n = vertices.size();
    vector<vector<double>> distMatrix(n, vector<double>(n));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            double dist = calculateDistance(vertices[i], vertices[j]);
            distMatrix[i][j] = dist;
            distMatrix[j][i] = dist;
        }
    }

    size_t MAX_ITERATIONS = 100000;

    vector<size_t> tour_nn;
    {
        vector<bool> visited(n, false);
        size_t current = 0;
        tour_nn.push_back(current);
        visited[current] = true;

        for (size_t i = 1; i < n; ++i) {
            double minDist = numeric_limits<double>::infinity();
            size_t next = n;
            for (size_t j = 0; j < n; ++j) {
                if (!visited[j] && distMatrix[current][j] < minDist) {
                    minDist = distMatrix[current][j];
                    next = j;
                }
            }
            if (next == n) break;
            visited[next] = true;
            tour_nn.push_back(next);
            current = next;
        }
    }

    vector<size_t> tour_fi;
    {
        vector<bool> inTour(n, false);
        size_t start_vertex = 0;
        size_t furthest_vertex = 0;
        double maxDist = 0.0;
        for (size_t i = 1; i < n; ++i) {
            double dist = distMatrix[start_vertex][i];
            if (dist > maxDist) {
                maxDist = dist;
                furthest_vertex = i;
            }
        }

        tour_fi.push_back(start_vertex);
        tour_fi.push_back(furthest_vertex);
        inTour[start_vertex] = true;
        inTour[furthest_vertex] = true;

        tour_fi.push_back(start_vertex);

        while (tour_fi.size() - 1 < n) {
            size_t furthestVertex = n;
            double maxMinDist = -1.0;

            for (size_t i = 0; i < n; ++i) {
                if (!inTour[i]) {
                    double minDistToTour = numeric_limits<double>::infinity();
                    for (size_t v : tour_fi) {
                        if (distMatrix[i][v] < minDistToTour) {
                            minDistToTour = distMatrix[i][v];
                        }
                    }
                    if (minDistToTour > maxMinDist) {
                        maxMinDist = minDistToTour;
                        furthestVertex = i;
                    }
                }
            }

            if (furthestVertex == n) break;

            size_t bestPosition = 0;
            double minimalIncrease = numeric_limits<double>::infinity();
            for (size_t i = 0; i < tour_fi.size() - 1; ++i) {
                size_t curr = tour_fi[i];
                size_t next = tour_fi[i + 1];
                double increase = distMatrix[curr][furthestVertex] + distMatrix[furthestVertex][next] - distMatrix[curr][next];
                if (increase < minimalIncrease) {
                    minimalIncrease = increase;
                    bestPosition = i + 1;
                }
            }

            tour_fi.insert(tour_fi.begin() + static_cast<std::ptrdiff_t>(bestPosition), furthestVertex);
            inTour[furthestVertex] = true;
        }

        tour_fi.pop_back();
    }

    auto twoOpt = [&](vector<size_t>& tour_ref) {
        size_t iteration = 0;
        bool improved = true;
        while (improved && iteration < MAX_ITERATIONS) {
            improved = false;
            double bestDelta = 0.0;
            size_t best_i = 0, best_k = 0;
            for (size_t i = 0; i < n - 1; ++i) {
                for (size_t k = i + 2; k < n; ++k) {
                    if (i == 0 && k == n - 1) continue;

                    double delta = -distMatrix[tour_ref[i]][tour_ref[i + 1]]
                                   -distMatrix[tour_ref[k]][tour_ref[(k + 1) % n]]
                                   +distMatrix[tour_ref[i]][tour_ref[k]]
                                   +distMatrix[tour_ref[i + 1]][tour_ref[(k + 1) % n]];

                    if (delta < bestDelta) {
                        bestDelta = delta;
                        best_i = i;
                        best_k = k;
                    }
                }
            }
            if (bestDelta < -1e-6) {
                reverse(tour_ref.begin() + static_cast<std::ptrdiff_t>(best_i + 1),
                        tour_ref.begin() + static_cast<std::ptrdiff_t>(best_k + 1));
                improved = true;
            }
            iteration++;
        }
    };

    twoOpt(tour_nn);
    twoOpt(tour_fi);

    auto computeTotalLength = [&](const vector<size_t>& tour) -> double {
        double totalLength = 0.0;
        for (size_t i = 0; i < n; ++i) {
            size_t from = tour[i];
            size_t to = tour[(i + 1) % n];
            totalLength += distMatrix[from][to];
        }
        return totalLength;
    };

    double totalLength_nn = computeTotalLength(tour_nn);
    double totalLength_fi = computeTotalLength(tour_fi);

    vector<size_t>* best_tour;
    double best_totalLength;
    if (totalLength_nn < totalLength_fi) {
        best_tour = &tour_nn;
        best_totalLength = totalLength_nn;
    } else {
        best_tour = &tour_fi;
        best_totalLength = totalLength_fi;
    }

    if ((*best_tour)[0] != 0) {
        auto it = find(best_tour->begin(), best_tour->end(), 0);
        if (it != best_tour->end()) {
            rotate(best_tour->begin(), it, best_tour->end());
        }
    }

    best_tour->push_back(0);

    cout << best_totalLength << "\n";
    for (size_t i = 0; i < n; ++i) {
        cout << (*best_tour)[i] << " ";
    }
    cout << "\n";
}
