// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#define HERMES_REPORT_WARN

#include "../hermes2d.h"

// Essential BC.
EssentialBC::EssentialBC(Hermes::vector<std::string> markers) : markers(markers) {
  current_time = 0.0;
  value_const = 0.0;
};

EssentialBC::EssentialBC(std::string marker) {
  markers.push_back(marker);
  current_time = 0.0;
  value_const = 0.0;
};

EssentialBC::~EssentialBC() {};

scalar EssentialBC::value(double x, double y) const {
  warn("EssentialBC::Function used either for a constant condition, or not redefined for nonconstant condition.");
  return 0.0;
};

void EssentialBC::set_current_time(double time) {
  this->current_time = time;
}

double EssentialBC::get_current_time() const {
  return current_time;
}

// Essential BC Constant
DefaultEssentialBCConst::DefaultEssentialBCConst(Hermes::vector<std::string> markers, scalar value_const) : EssentialBC(markers) {
  this->value_const = value_const;
}

DefaultEssentialBCConst::DefaultEssentialBCConst(std::string marker, scalar value_const) : EssentialBC(Hermes::vector<std::string>()) {
  this->value_const = value_const;
  markers.push_back(marker);
}

// EssentialBCs.
EssentialBCs::EssentialBCs() {
};

EssentialBCs::EssentialBCs(Hermes::vector<EssentialBC *> essential_bcs) {
  add_boundary_conditions(essential_bcs);
};

EssentialBCs::EssentialBCs(EssentialBC * boundary_condition) {
  Hermes::vector<EssentialBC *> boundary_conditions;
  boundary_conditions.push_back(boundary_condition);
  add_boundary_conditions(boundary_conditions);
};

void EssentialBCs::add_boundary_conditions(Hermes::vector<EssentialBC *> boundary_conditions) {
  for(Hermes::vector<EssentialBC *>::iterator it = boundary_conditions.begin(); it != boundary_conditions.end(); it++)
        all.push_back(*it);

  markers.clear();
  create_marker_cache();
};

void EssentialBCs::add_boundary_condition(EssentialBC * boundary_condition) {
  Hermes::vector<EssentialBC *> boundary_conditions;
  boundary_conditions.push_back(boundary_condition);
  add_boundary_conditions(boundary_conditions);
};

Hermes::vector<EssentialBC *>::const_iterator EssentialBCs::begin() const {
  return all.begin();
}

Hermes::vector<EssentialBC *>::const_iterator EssentialBCs::end() const {
  return all.end();
}

EssentialBCs::~EssentialBCs() {
};

void EssentialBCs::create_marker_cache() {
  for(this->iterator = begin(); iterator != end(); iterator++)
    for(Hermes::vector<std::string>::const_iterator it = (*iterator)->markers.begin(); it != (*iterator)->markers.end(); it++) {
      if (markers[*it] != NULL)
        error("Attempt to define more than one description of the BC on the same part of the boundary with marker '%s'.", it->c_str());
      markers[*it] = *iterator;
    }
}


EssentialBC* EssentialBCs::get_boundary_condition(std::string marker) {
  if(markers.find(marker) == markers.end())
    return NULL;
  else
    return markers[marker];
}

void EssentialBCs::set_current_time(double time) {
  for(iterator = begin(); iterator != end(); iterator++)
    (*iterator)->set_current_time(time);
};
