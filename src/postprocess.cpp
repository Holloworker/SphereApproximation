#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <sstream>

struct Sphere{
    float x,y,z;
    float radius;
};

struct SpherePair{
    float index1, index2;
    float distance;
    
    bool operator<(const SpherePair& other) const{
        return distance > other.distance;
    };
};

float ComputeDistance(const Sphere& a, const Sphere& b){
    float centerdistance = std::sqrt(
        std::pow(a.x - b.x, 2) +
        std::pow(a.y - b.y, 2) +
        std::pow(a.z - b.z, 2)
    );
    return centerdistance + a.radius + b.radius;
}

Sphere mergeSpheres(const Sphere& a, const Sphere& b) {
    Sphere newSphere;

    float centerDistance = std::sqrt(
        std::pow(a.x - b.x, 2) +
        std::pow(a.y - b.y, 2) +
        std::pow(a.z - b.z, 2)
    );

    float totalRadius = a.radius + b.radius;
    newSphere.x = (a.x * b.radius + b.x * a.radius) / totalRadius;
    newSphere.y = (a.y * b.radius + b.y * a.radius) / totalRadius;
    newSphere.z = (a.z * b.radius + b.z * a.radius) / totalRadius;

    newSphere.radius = (centerDistance + a.radius + b.radius) / 2;

    return newSphere;
}

void writeActiveSpheresToFile(const std::vector<Sphere>& spheres, const std::unordered_set<int>& activeIndices, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    for (int index : activeIndices) {
        const Sphere& sphere = spheres[index];
        outfile << sphere.x << " " 
                << sphere.y << " " 
                << sphere.z << " " 
                << sphere.radius << " " 
                << 1.000000 << std::endl;
    }

    outfile.close();
    std::cout << "Active spheres written to file: " << filename << std::endl;
}

int main(int argc, const char *argv[]){
    if (argc != 4){
        std::cerr << "Usage: " << argv[0] << "<input_filename> <output_filename> <target_count>" << std::endl;
        return 1;
    }
    std::string inputfile = argv[1];
    std::string outputfile = argv[2];
    int targetcount = std::stoi(argv[3]);
    if (targetcount <= 0) {
        std::cerr << "Target count must be greater than 0." << std::endl;
        return 1;
    }
    std::ifstream file(inputfile);
    if (!file.is_open()){
        std::cerr << "fail to open the file" << std::endl;
    }
    std::vector<Sphere> spheres;
    std::string line;
    while (std::getline(file, line)){
        std::istringstream iss(line);
        Sphere sphere;
        float extra;
        if(iss >> sphere.x >> sphere.y >> sphere.z >> sphere.radius >> extra){
            spheres.push_back(sphere);
        }
    }
    file.close();
    std::unordered_set<int> activeIndices;   
    std::priority_queue<SpherePair> pq;
    for (size_t i = 0; i < spheres.size(); ++i) {
        activeIndices.insert(i);
    }
    for (size_t i = 0; i < spheres.size(); ++i) {
        for (size_t j = i + 1; j < spheres.size(); ++j) {
            float distance = ComputeDistance(spheres[i], spheres[j]);
            pq.push({static_cast<int>(i), static_cast<int>(j), distance});
        }
    }

    while (activeIndices.size() > targetcount) {
        SpherePair closestPair = pq.top();
        pq.pop();

        if (activeIndices.find(closestPair.index1) == activeIndices.end() ||
            activeIndices.find(closestPair.index2) == activeIndices.end()) {
            continue;
        }

        std::cout << "Merging spheres " << closestPair.index1 << " and " << closestPair.index2 << std::endl;

        Sphere newSphere = mergeSpheres(spheres[closestPair.index1], spheres[closestPair.index2]);
        int newIndex = spheres.size();
        spheres.push_back(newSphere);

        activeIndices.erase(closestPair.index1);
        activeIndices.erase(closestPair.index2);
        activeIndices.insert(newIndex);

        for (int index : activeIndices) {
            if (index == newIndex) continue;
            float distance = ComputeDistance(spheres[index], newSphere);
            pq.push({index, newIndex, distance});
        }
    }
    writeActiveSpheresToFile(spheres, activeIndices, outputfile);
}