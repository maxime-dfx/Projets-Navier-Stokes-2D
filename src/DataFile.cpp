#include "DataFile.h"
#include <fstream>
#include <iostream>
#include <cmath>
// Assurez-vous que le chemin vers toml est correct selon votre installation
#include <toml/toml.hpp> 

using namespace std;

DataFile::DataFile(std::string file_name) : _file_name(file_name)
{
    cout << "Lecture de " << file_name << "..." << endl;
    auto config = toml::parse(file_name);

    const auto& space = toml::find(config, "space");
    _xmin = toml::find<double>(space, "xmin");
    _xmax = toml::find<double>(space, "xmax");
    _ymin = toml::find<double>(space, "ymin");
    _ymax = toml::find<double>(space, "ymax");
    _hx   = toml::find<double>(space, "hx");
    _hy   = toml::find<double>(space, "hy");
    
    // Calcul robuste des dimensions
    _Nx = static_cast<int>(round((_xmax - _xmin) / _hx));
    _Ny = static_cast<int>(round((_ymax - _ymin) / _hy));

    const auto& time = toml::find(config, "time");
    _t0     = toml::find<double>(time, "t0");
    _tfinal = toml::find<double>(time, "tfinal");
    _dt     = toml::find<double>(time, "dt");
    _scheme = toml::find<std::string>(time, "scheme");

    const auto& physics = toml::find(config, "physics");
    _nu  = toml::find<double>(physics, "nu");
    _rho = toml::find<double>(physics, "rho");
    
    const auto& output = toml::find(config, "output");
    _results = toml::find<std::string>(output, "results");

    const auto& boundary = toml::find(config, "boundary");
    _bc_left   = toml::find<std::string>(boundary, "left");
    _bc_right  = toml::find<std::string>(boundary, "right");
    _bc_bottom = toml::find<std::string>(boundary, "bottom");
    _bc_top    = toml::find<std::string>(boundary, "top");
    
    cout << "Config: Grid " << _Nx << "x" << _Ny << " | BC: " 
         << _bc_left << "/" << _bc_right << "/" << _bc_bottom << "/" << _bc_top << endl;
}