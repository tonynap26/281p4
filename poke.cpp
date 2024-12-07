// Project Identifier: 5949F553E20B650AB0FB2266D3C0822B13D248B0

#include <iostream>
#include <queue>
#include <vector>
#include <utility>
#include <cmath>
#include <limits>
#include <getopt.h>
#include <iomanip>
#include <algorithm>
#include <cstddef>
#include <functional>

using namespace std;

enum class Area { LAND, SEA, COASTLINE };

void handleMST(const vector<pair<int, int>>& vertices);
void handleFASTTSP(const vector<pair<int, int>>& vertices);
void handleOPTTSP(const vector<pair<int, int>>& vertices);
double calculateDistance(const pair<int, int>& p1, const pair<int, int>& p2);
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

double calculateDistance(const pair<int, int>& p1, const pair<int, int>& p2) {
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

    // Start from node 0
    minEdge[0] = 0.0;
    double totalWeight = 0.0;
    priority_queue<pair<double, size_t>, vector<pair<double, size_t>>, std::greater<pair<double, size_t>>> pq;
    pq.emplace(0.0, 0);

    while (!pq.empty()) {
        double weight = pq.top().first;
        size_t u = pq.top().second;
        pq.pop();

        if (inMST[u]) continue;
        inMST[u] = true;
        totalWeight += weight;

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
                        pq.emplace(minEdge[v], v);
                        parent[v] = u;
                    }
                }
            }
        }
    }

    for (bool included : inMST) {
        if (!included) {
            cerr << "Cannot construct MST\n";
            exit(1);
        }
    }

    cout << totalWeight << "\n";

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
            for (size_t k = i + 1; k < n; ++k) {
                double delta = 0.0;
                delta -= calculateDistance(vertices[tour[i - 1]], vertices[tour[i]]);
                delta -= calculateDistance(vertices[tour[k]], vertices[tour[(k + 1) % n]]);
                delta += calculateDistance(vertices[tour[i - 1]], vertices[tour[k]]);
                delta += calculateDistance(vertices[tour[i]], vertices[tour[(k + 1) % n]]);

                if (delta < -1e-6) {
                    reverse(tour.begin() + static_cast<ptrdiff_t>(i),
                            tour.begin() + static_cast<ptrdiff_t>(k) + 1);
                    totalWeight += delta;
                    improved = true;
                }
            }
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

struct BBNode {
    vector<size_t> path;
    double cost;
    double bound;
    bool operator<(const BBNode& other) const {
        return bound > other.bound;
    }
};

double computeLowerBound(const vector<vector<double>>& dist, const vector<bool>& visited, const vector<pair<double, double>>& min_two) {
    double lb = 0.0;
    for (size_t i = 0; i < dist.size(); ++i) {
        if (!visited[i]) {
            lb += min_two[i].first;
            lb += min_two[i].second;
        }
    }
    return lb / 2.0;
}

void handleOPTTSP(const vector<pair<int, int>>& vertices) {
    size_t n = vertices.size();
    vector<vector<double>> dist(n, vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (i != j)
                dist[i][j] = calculateDistance(vertices[i], vertices[j]);
            else
                dist[i][j] = 0.0;

    vector<pair<double, double>> min_two(n, {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()});
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            double d = dist[i][j];
            if (d < min_two[i].first) {
                min_two[i].second = min_two[i].first;
                min_two[i].first = d;
            }
            else if (d < min_two[i].second) {
                min_two[i].second = d;
            }
        }
    }

    double best_cost = numeric_limits<double>::infinity();
    vector<size_t> best_path;

    vector<bool> visited(n, false);
    visited[0] = true;
    vector<size_t> path;
    path.push_back(0);

    auto computeLB = [&](const vector<bool>& visited_ref) -> double {
        return computeLowerBound(dist, visited_ref, min_two);
    };

    function<void(size_t, size_t, double)> tsp_recursive = [&](size_t current, size_t count, double current_cost) {
        if (count == n) {
            double total_cost = current_cost + dist[current][0];
            if (total_cost < best_cost) {
                best_cost = total_cost;
                best_path = path;
            }
            return;
        }

        vector<size_t> candidates;
        for (size_t next = 0; next < n; ++next) {
            if (!visited[next]) {
                candidates.push_back(next);
            }
        }

        sort(candidates.begin(), candidates.end(), [&](size_t a, size_t b) -> bool {
            return dist[current][a] < dist[current][b];
        });

        for (size_t next : candidates) {
            double new_cost = current_cost + dist[current][next];
            if (new_cost >= best_cost) continue;

            double lb = new_cost + computeLB(visited);
            if (lb >= best_cost) continue;

            visited[next] = true;
            path.push_back(next);
            tsp_recursive(next, count + 1, new_cost);
            path.pop_back();
            visited[next] = false;
        }
    };

    tsp_recursive(0, 1, 0.0);

    cout << best_cost << "\n";
    for (size_t city : best_path) {
        cout << city << " ";
    }
    cout << "\n";
}