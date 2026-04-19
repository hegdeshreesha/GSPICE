#ifndef GSPICE_NODE_HPP
#define GSPICE_NODE_HPP

#include <string>

namespace gspice {

class Node {
public:
    Node(const std::string& name, int id) : name_(name), id_(id) {}
    
    std::string getName() const { return name_; }
    int getId() const { return id_; }

private:
    std::string name_;
    int id_; // Index in the MNA matrix
};

} // namespace gspice

#endif // GSPICE_NODE_HPP
