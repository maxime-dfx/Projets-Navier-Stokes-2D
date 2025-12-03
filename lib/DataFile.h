#ifndef _DATA_FILE_H
#define _DATA_FILE_H

#include <string>
#include <iostream>
#include <vector>

class DataFile {
private:
   double _xmin, _xmax, _ymin, _ymax, _hx, _hy;
   int _Nx, _Ny;
   double _nu, _rho;
   double _t0, _tfinal, _dt;
   std::string _scheme;
   std::string _results;
   const std::string _file_name;

   std::string _bc_left;
   std::string _bc_right;
   std::string _bc_bottom;
   std::string _bc_top;

   double _bc_left_dir;
   double _bc_right_dir;
   double _bc_bottom_dir;
   double _bc_top_dir;

public:
   DataFile(std::string file_name);

   double Get_xmin() const { return _xmin; }
   double Get_xmax() const { return _xmax; }
   double Get_ymin() const { return _ymin; }
   double Get_ymax() const { return _ymax; }
   double Get_hx() const { return _hx; }
   double Get_hy() const { return _hy; }
   int Get_Nx() const { return _Nx; }
   int Get_Ny() const { return _Ny; }
   double Get_nu() const { return _nu; }
   double Get_rho() const { return _rho; }
   double Get_t0() const { return _t0; }
   double Get_tfinal() const { return _tfinal; }
   double Get_dt() const { return _dt; }
   std::string Get_scheme() const { return _scheme; }
   std::string Get_results() const { return _results; }

   std::string Get_BC_Left()   const { return _bc_left; }
   std::string Get_BC_Right()  const { return _bc_right; }
   std::string Get_BC_Bottom() const { return _bc_bottom; }
   std::string Get_BC_Top()    const { return _bc_top; }
   double Get_BC_Left_dir() const { return _bc_left_dir; }
   double Get_BC_Right_dir() const { return _bc_right_dir; }
   double Get_BC_Bottom_dir() const { return _bc_bottom_dir; }
   double Get_BC_Top_dir() const { return _bc_top_dir; }
};

#endif