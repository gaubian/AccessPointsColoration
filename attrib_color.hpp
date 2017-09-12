#ifndef ATTRIB_COLOR_HPP
#define ATTRIB_COLOR_HPP

using id = std::string;
using color = int;
using map_id_int = std::map<id,int>;
using neighbours = std::set<std::pair<double,id>>;
using graph = std::map<id,neighbours>;
using coord = std::tuple<double,double,double>;

double to_radians(double x);
double dist(coord x, coord y);
double score(id name, graph &grph, map_id_int &assigner);
bool improve(id name, color c, graph &grph, map_id_int &assigner);
id random_elt(std::set<id> se);
map_id_int solver(graph &grph, map_id_int &auth_colors);
double volume(coord x, double rx, coord y, double ry);
graph grph_of_coord(std::map<id,coord> &coord_map,
		    std::map<id,double> &radius_map);
#endif
