#ifndef BLOCK_OF_LAYERS_H
#define BLOCK_OF_LAYERS_H

#include "fem/point.h"
#include "layer.h"
#include <vector>


class Rectangle;



class BlockOfLayers
{
public:
  void init(const Point &min_point, const Point &max_point,
            unsigned int n_layers, double angle,
            const std::vector<double> &layers_thickness,
            const std::vector<double> &layers_coef_alpha,
            const std::vector<double> &layers_coef_beta);

  bool contains_element(const Rectangle &rect,
                        const std::vector<Point> &points) const;

  void get_coefs(const Rectangle &cell,
                 const std::vector<Point> &points,
                 double &coef_alpha,
                 double &coef_beta) const;


private:
  double _angle;
  unsigned int _n_layers;
  Point _min_point; // limits of the block
  Point _max_point;

  std::vector<double> _layers_thickness;
  std::vector<double> _layers_coef_alpha;
  std::vector<double> _layers_coef_beta;

  std::vector<Layer> _layers;
};




#endif // BLOCK_OF_LAYERS_H
