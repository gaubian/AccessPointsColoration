#include <map>
#include <set>
#include <iostream>
#include <cmath>
#include "attrib_color.hpp"

const double earth_radius = 6379.009 * 1000;
const double pi = 4*atan(1);
const double floor_size = 3;

// Converts degrees to radians.
double to_radians(double x) {
    return pi * x / 180;
}

// Distance in meters between two positions.
double dist(coord x, coord y) {
    double latx = std::get<0>(x),
	   lonx = std::get<1>(x),
	   altx = std::get<2>(x),
	   laty = std::get<0>(y),
	   lony = std::get<1>(y),
	   alty = std::get<2>(y);
    double latdist = to_radians(laty - latx);
    double londist = to_radians(lony - lonx);
    double a = sin(latdist / 2) * sin(latdist / 2)
	    + cos(to_radians(latx)) * cos(to_radians(laty))
	    * sin(londist / 2) * sin(londist / 2);
    double c = 2 * atan2(sqrt(a),sqrt(1-a));
    double d = earth_radius * c;
    double height = altx - alty;
    d = d * d + height * height;
    return sqrt(d);
}

// Sum the weight of an access points' neighbours which have same color.
double score(id name, graph &grph, map_id_int &assigner) {
    double ans = 0;
    neighbours nei = grph[name];
    for(std::pair<double,id> p : nei)
	if(assigner[name] == assigner[p.second]) ans += p.first;
    return ans;
}

// Returns whether switching to a given color an access point improves
// our solution, and do it if it does.
bool improve(id name, color c, graph &grph, map_id_int &assigner) {
    int init_c = assigner[name];
    double init_score = score(name,grph,assigner);
    assigner[name] = c;
    double new_score = score(name,grph,assigner);
    if(new_score >= init_score) assigner[name] = init_c;
    return new_score < init_score;
}

// Select a random element from a set.
id random_elt(std::set<id> se) {
    int n = se.size();
    std::set<id>::iterator it = se.begin();
    advance(it,rand() % n);
    return *it;
}

// Colors trying to minimize sum of weights of same-color neighboors.
map_id_int solver(graph &grph, map_id_int &auth_colors) {
    map_id_int assigner;
    std::set<id> to_see;
    for(std::pair<id,int> p : auth_colors) {
	 assigner[p.first] = 0;
	 to_see.insert(p.first);
    }
    while(!to_see.empty()) {
        id access_point = random_elt(to_see);
        int k = auth_colors[access_point];
        bool flag = false;
	for(color i = 0; i < k && !flag; ++i)
            flag = flag || improve(access_point,i,grph,assigner);
	if(flag)
            for(std::pair<double,id> e : grph[access_point])
		to_see.insert(e.second);
	to_see.erase(access_point);
    }
    return assigner;
}

// Computes volume of two spheres' intersection.
double volume(coord x, double rx, coord y, double ry) {
    double d = dist(x,y);
    if(d >= rx + ry) return 0;
    double rxplusry = rx + ry;
    double rxlessry = rx - ry;
    double rxryd = rxplusry - d;
    double rhs = d * d + 2 * d * rxplusry - 3 * rxlessry * rxlessry;
    return pi * rxryd * rxryd * rhs / (12 * d);
}

// Converts coordinates and radii to a weighted graph of collision.
graph grph_of_coord(std::map<id,coord> &coord_map,
		    std::map<id,double> &radius_map) {
    graph ans;
    for(std::pair<id,coord> p : coord_map) {
	neighbours name_neighbours;
	coord pos = p.second;
        id name = p.first;
	for(std::pair<id,coord> nei : coord_map) {
	    id nei_name = nei.first;
	    if(nei_name == name) continue;
	    coord nei_pos = nei.second;
	    double vol = volume(pos,radius_map[name],
	        nei_pos,radius_map[nei_name]);
	    if(vol > 0) {
	       name_neighbours.insert({vol,nei_name});
	       name_neighbours.insert({vol,name});
	    }
	}
	ans[name] = name_neighbours;
    }
    return ans;
}

int main() {
    id name;
    std::string fl;
    double lon, lat, alt, rad;
    int lim_col = 4;
    map_id_int auth_color;
    std::map<id,coord> coord_map;
    std::map<id,double> radius_map;
    while(std::cin >> name >> lon >> lat>> alt >> fl >> rad) {
	coord c = std::make_tuple(lon,lat,floor_size * alt);
	auth_color[name] = lim_col;
        coord_map[name] = c;
	radius_map[name] = rad;
    }
    graph grph = grph_of_coord(coord_map, radius_map);
    map_id_int assigner = solver(grph,auth_color);
    for(std::pair<id,color> p : assigner)
	std::cout << p.first << " " << p.second << std::endl;
}

