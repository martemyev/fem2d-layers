#include "block_of_layers.h"
#include "fem/rectangle.h"
#include "fem/auxiliary_functions.h"



void BlockOfLayers::init(const Point &min_point, const Point &max_point,
                         unsigned int n_layers, double angle,
                         const std::vector<double> &layers_thickness,
                         const std::vector<double> &layers_coef_alpha,
                         const std::vector<double> &layers_coef_beta)
{
  require(n_layers == layers_thickness.size() &&
          n_layers == layers_coef_alpha.size() &&
          n_layers == layers_coef_beta.size(),
          "The number of layers doesn't coinside with the size of suppoting vectors (thickness, coef_alpha, coef_beta)");

  _min_point = min_point;
  _max_point = max_point;
  _n_layers = n_layers;
  _angle = angle;
  _layers_thickness = layers_thickness;
  _layers_coef_alpha = layers_coef_alpha;
  _layers_coef_beta = layers_coef_beta;

  _layers.resize(_n_layers);
  for (unsigned int i = 0; i < _n_layers; ++i)
    _layers[i].init(i, _layers_thickness, _min_point, _max_point, _angle);
}



bool BlockOfLayers::contains_element(const Rectangle &cell,
                                     const std::vector<Point> &points) const
{
  // we think that a cell belongs to a block if a center of the cell belongs to the block
  double xc = 0., yc = 0.; // center of the cell

  for (unsigned int i = 0; i < cell.n_vertices; ++i)
  {
    xc += points[cell.vertex(i)].coord(0);
    yc += points[cell.vertex(i)].coord(1);
  }
  xc /= cell.n_vertices;
  yc /= cell.n_vertices;

  if (xc <= _max_point.coord(0) &&
      xc >= _min_point.coord(0) &&
      yc <= _max_point.coord(1) &&
      yc >= _min_point.coord(1))
    return true;

  return false;
}



void BlockOfLayers::get_coefs(const Rectangle &cell,
                              const std::vector<Point> &points,
                              double &coef_alpha,
                              double &coef_beta) const
{
  expect(contains_element(cell, points), "This block doesn't have this cell");

  for (unsigned int i = 0; i < _n_layers; ++i)
  {
    if (_layers[i].contains_element(cell, points))
    {
      coef_alpha = _layers_coef_alpha[i];
      coef_beta  = _layers_coef_beta[i];
      return;
    }
  }

  require(false, "The element cannot be found in layers");
}
