#include "DataFile.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <toml/toml.hpp> // Assure-toi que c'est bien toml11

using namespace std;

DataFile::DataFile(std::string file_name) : _file_name(file_name)
{
    cout << "Lecture de " << file_name << "..." << endl;
    
    // 1. Parsing du fichier
    auto config = toml::parse(file_name);

    // --- SPACE ---
    const auto& space = toml::find(config, "space");
    _xmin = toml::find<double>(space, "xmin");
    _xmax = toml::find<double>(space, "xmax");
    _ymin = toml::find<double>(space, "ymin");
    _ymax = toml::find<double>(space, "ymax");
    _hx   = toml::find<double>(space, "hx");
    _hy   = toml::find<double>(space, "hy");
    
    _Nx = static_cast<int>(round((_xmax - _xmin) / _hx));
    _Ny = static_cast<int>(round((_ymax - _ymin) / _hy));

    // --- TIME ---
    const auto& time = toml::find(config, "time");
    _t0     = toml::find<double>(time, "t0");
    _tfinal = toml::find<double>(time, "tfinal");
    _dt     = toml::find<double>(time, "dt");
    _scheme = toml::find<std::string>(time, "scheme");

    // --- PHYSICS ---
    const auto& physics = toml::find(config, "physics");
    _nu  = toml::find<double>(physics, "nu");
    _rho = toml::find<double>(physics, "rho");
    
    // --- OUTPUT ---
    const auto& output = toml::find(config, "output");
    _results = toml::find<std::string>(output, "results");

    // --- BOUNDARY (Corrigé) ---
    // On utilise toml::find_or pour gérer les valeurs par défaut
    // Syntaxe: toml::find_or<Type>(Table, "Clé", Valeur_Par_Defaut)
    
    const auto& boundary = toml::find(config, "boundary");

    _bc_left   = toml::find_or<std::string>(boundary, "left",   "Dirichlet");
    _bc_right  = toml::find_or<std::string>(boundary, "right",  "Dirichlet");
    _bc_bottom = toml::find_or<std::string>(boundary, "bottom", "Dirichlet");
    _bc_top    = toml::find_or<std::string>(boundary, "top",    "Dirichlet");

    _bc_left_dir   = toml::find_or<double>(boundary, "left_val",   0.0);
    _bc_right_dir  = toml::find_or<double>(boundary, "right_val",  0.0);
    _bc_bottom_dir = toml::find_or<double>(boundary, "bottom_val", 0.0);
    _bc_top_dir    = toml::find_or<double>(boundary, "top_val",    0.0);
    
    cout << "Config: Grid " << _Nx << "x" << _Ny << " | BC Types: " 
         << _bc_left << "/" << _bc_right << "/" << _bc_bottom << "/" << _bc_top << endl;
}