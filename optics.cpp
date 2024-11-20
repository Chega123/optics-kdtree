#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ctime>
#include <map>
#include <iomanip> 
#include<algorithm>
#include <set>
using namespace std;

vector<double*> load_points_from_csv(string& filename, int dimension) {
    ifstream file(filename);
    vector<double*> points;
    
    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir el archivo " << filename << endl;
        return points;
    }
    
    string line;
    getline(file, line); // Ignorar encabezado
    
    while (getline(file, line)) {
        istringstream ss(line);
        double* point = new double[dimension];
        
        for (int i = 0; i < dimension; ++i) {
            string value;
            if (!getline(ss, value, ',')) break;
            point[i] = stod(value);
        }
        points.push_back(point);
    }
    
    file.close();
    return points;
}

class KDtree{
   

    public:
    int dimension;
    struct Node{
        Node* left;
        Node* right;
        double* point;
        Node(int dimension): left(nullptr), right(nullptr){point = new double[dimension];}
    };
    Node* root;

    KDtree(int dimension) : dimension(dimension), root(nullptr) {}

    void insert(double* point){
        root=insert_node(root,point,0);
    }

    double* nearest(double* point) {
        double best_dist = 10000000;
        return search_near(root, point, 0, best_dist);
    }
        
    double euclidian_dist(double* p1, double* p2){
        double sum=0;
        for(int i=0;i<dimension;i++){
            sum+=(p1[i]-p2[i])*(p1[i]-p2[i]);
        }
        return sqrt(sum);
    }

    void print_node(Node* node, int depth) {
        if (node == nullptr) return;

        std::cout << std::setw(depth * 2) << ""; 
        std::cout << "Nivel " << depth << " | Punto: (";
        for (int i = 0; i < dimension; i++) {
            std::cout << node->point[i];
            if (i < dimension - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;

        // Recursivamente imprimir los hijos
        if (node->left != nullptr) {
            std::cout << std::setw(depth * 2) << "" << "Izquierda -> ";
            print_node(node->left, depth + 1);
        }
        if (node->right != nullptr) {
            std::cout << std::setw(depth * 2) << "" << "Derecha -> ";
            print_node(node->right, depth + 1);
        }
    }

    vector<double*> range_query(double* point, double rad){
        vector<double*>neighbors;
        range_query_function(root,point,rad,0,neighbors);
        return neighbors;
    }

    private:

    Node* insert_node(Node* node,double* point,int dim){
        if(node==nullptr){
            Node* temp= new Node (dimension);
            for(int i=0;i<dimension;i++){
                temp->point[i]=point[i];
            }
            return temp;
        }

        int div=dim%dimension; // se refiere a en que dimension se hara la division

        if(point[div]<node->point[div]){
            node->left=insert_node(node->left,point,dim+1);
        }
        else{
            node->right=insert_node(node->right,point,dim+1);
        }
        return node;
    }

    double* search_near(Node* node,double* point, int dim,double& best_dist){
        if(node==nullptr){return nullptr;}
        double dist=euclidian_dist(node->point,point);
        double* best_point=nullptr;

        if(dist<best_dist){best_dist=dist; best_point=node->point;}

        int div=dim%dimension;
        Node* prolly_branch=(point[div]<node->point[div]) ? node->left:node->right;
        Node* other_branch=(point[div]<node->point[div]) ? node->right:node->left;

        double* next_best_point=search_near(prolly_branch,point,dim+1,best_dist);

        if(next_best_point!=nullptr&&euclidian_dist(next_best_point, point) < best_dist){best_point=next_best_point;}
        
        if(abs(point[div]-node->point[div])<best_dist){
            double* other_best_point=search_near(other_branch,point,dim+1,best_dist);
            if (other_best_point != nullptr&&euclidian_dist(other_best_point, point) < best_dist) {
                best_point = other_best_point;
            }
        }
        return best_point;
    }

    void range_query_function(Node* node, double*point, double rad, int dim, vector<double*>&results){
        if(node == nullptr)return;
        double dist=euclidian_dist(node->point,point);
        if(dist<=rad){
            results.push_back(node->point);
        }
        int div=dim%dimension;
        if(point[div]-rad<=node->point[div]){
            range_query_function(node->left,point,rad,dim+1,results);
        }
        if(point[div]+rad>=node->point[div]){
            range_query_function(node->right,point,rad,dim+1,results);
        }

    }


};


class OPTICS {
public:
    struct PointInfo {
        bool procesado = false;
        double core_distance = -1;
        double reachability_distance = -1;
    };

    map<double*, PointInfo> pointInfo;
    OPTICS(vector<double*>& dataset, KDtree& kdtree, int dimension, double epsilon, int minPts)
        : dataset(dataset), kdtree(kdtree), dimension(dimension), epsilon(epsilon), minPts(minPts) {
        initializePointInfo();
    }

    vector<double*> run() {
        vector<double*> ordered;
        for (double* point : dataset) {
            if (!pointInfo[point].procesado) {
                vector<double*> neighbors = kdtree.range_query(point, epsilon);
                markAsProcessed(point);
                ordered.push_back(point);

                if (coreDistance(point, neighbors) != -1) {
                    set<pair<double, double*>> seeds;
                    update(neighbors, point, seeds);

                    while (!seeds.empty()) {
                        double* next = seeds.begin()->second;
                        seeds.erase(seeds.begin());

                        vector<double*> neighborsNext = kdtree.range_query(next, epsilon);
                        markAsProcessed(next);
                        ordered.push_back(next);

                        if (coreDistance(next, neighborsNext) != -1) {
                            update(neighborsNext, next, seeds);
                        }
                    }
                }
            }
        }
        return ordered;
    }

    vector<vector<double*>> extractClusters(double reachabilityThreshold, vector<double*> ordered) {
        set<double*> uniquePoints(ordered.begin(), ordered.end());
        vector<double*> filteredOrdered(uniquePoints.begin(), uniquePoints.end());

        vector<vector<double*>> clusters;
        vector<double*> currentCluster;

        for (double* point : filteredOrdered) {
            double reachabilityDistance = pointInfo[point].reachability_distance;
            cout << "Point reachability: " << reachabilityDistance << " - Threshold: " << reachabilityThreshold << endl;
            if (reachabilityDistance > reachabilityThreshold || reachabilityDistance == -1) {
                if (!currentCluster.empty()) {
                    cout<<"cluster nuevo: "<<endl;
                    clusters.push_back(currentCluster);
                    currentCluster.clear();
                }
            } else {
                currentCluster.push_back(point);
            }
        }
        if (!currentCluster.empty()) {
            clusters.push_back(currentCluster);
        }

        return clusters;
    }


private:
    
    vector<double*>& dataset;
    KDtree& kdtree;
    int dimension;
    double epsilon;
    int minPts;
    
    void initializePointInfo() {
        for (double* point : dataset) {
            pointInfo[point] = PointInfo();
        }
    }

    void markAsProcessed(double* point) {
        pointInfo[point].procesado = true;
    }

    double coreDistance(double* point, const vector<double*>& neighbors) {
        if (neighbors.size() < minPts) {
            return -1; // No cumple con los puntos mínimos
        }

        vector<double*> sortedNeighbors = neighbors;
        sort(sortedNeighbors.begin(), sortedNeighbors.end(),
             [this, point](double* a, double* b) {
                 return kdtree.euclidian_dist(point, a) < kdtree.euclidian_dist(point, b);
             });

        double coreDist = kdtree.euclidian_dist(point, sortedNeighbors[minPts - 1]);
        pointInfo[point].core_distance = coreDist;
        return coreDist;
    }

    void update(const vector<double*>& neighbors, double* point, set<pair<double, double*>>& seeds) {
        double coreDist = pointInfo[point].core_distance;

        for (double* neighbor : neighbors) {
            if (!pointInfo[neighbor].procesado) {
                double newReachDist = max(coreDist, kdtree.euclidian_dist(point, neighbor));

                if (pointInfo[neighbor].reachability_distance == -1) {
                    pointInfo[neighbor].reachability_distance = newReachDist;
                    seeds.insert({newReachDist, neighbor});
                } else if (newReachDist < pointInfo[neighbor].reachability_distance) {
                    seeds.erase({pointInfo[neighbor].reachability_distance, neighbor});
                    pointInfo[neighbor].reachability_distance = newReachDist;
                    seeds.insert({newReachDist, neighbor});
                }
            }
        }
    }
};


void save_clusters_to_csv(const vector<double*>& dataset, const map<vector<double>, int>& clusterMap, int dimension, const string& outputFilename) {
    ofstream file(outputFilename);
    for (double* point : dataset) {
        vector<double> pointVec(point, point + dimension);
        for (int i = 0; i < dimension; ++i) {
            file << point[i];
            if (i < dimension - 1) {
                file << ", ";
            }
        }
        file << ", " << clusterMap.at(pointVec) << "\n";
    }

    file.close();
    cout << "Clusters guardados en " << outputFilename << endl;
}




int main() {
    string filename = "abalone_real.csv";       
    string outputFilename = "clusters.csv"; 
    int dimension = 6;
    double epsilon = 2;
    int minPts = 8;
    double reachabilityThreshold = 0.18;    

    vector<double*> points = load_points_from_csv(filename, dimension);

    KDtree kdtree(dimension);
    for (double* point : points) {
        kdtree.insert(point);
    }


    OPTICS optics(points, kdtree, dimension, epsilon, minPts);
    vector<double*> ordered = optics.run();

    vector<vector<double*>> clusters = optics.extractClusters(reachabilityThreshold, ordered);
    cout<<clusters.size()<<endl;

    map<vector<double>, int> clusterMap; 
    int clusterId = 1;

    for (const auto& cluster : clusters) {
        for (double* point : cluster) {
            vector<double> pointVec(point, point + dimension); 
            clusterMap[pointVec] = clusterId;
        }
        ++clusterId;
    }
    for (double* point : points) {
        vector<double> pointVec(point, point + dimension); 
        if (clusterMap.count(pointVec) == 0) { 
            clusterMap[pointVec] = 0; 
        }
    }

    // Imprimir resultados
    for (double* point : points) {
        vector<double> pointVec(point, point + dimension);
        cout << "Punto " << pointVec[0] << " -> ClusterID: " << clusterMap[pointVec] << endl;
    }

    save_clusters_to_csv(points, clusterMap, dimension, outputFilename);

    for (double* point : points) {
        delete[] point;
    }

    return 0;
}
