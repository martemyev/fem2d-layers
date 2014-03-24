#ifndef BLOCK_OF_LAYERS_H
#define BLOCK_OF_LAYERS_H

#include "fem/point.h"
#include "fem/rectangle.h"
#include "layer.h"
#include <vector>



class BlockOfLayers
{
public:
  void init(const fem::Point &min_point, const fem::Point &max_point,
            unsigned int n_layers, double angle,
            const std::vector<double> &layers_thickness,
            const std::vector<double> &layers_coef_alpha,
            const std::vector<double> &layers_coef_beta);

  bool contains_element(const fem::Rectangle &rect,
                        const std::vector<fem::Point> &points) const;

  void get_coefs(const fem::Rectangle &cell,
                 const std::vector<fem::Point> &points,
                 double &coef_alpha,
                 double &coef_beta) const;

  const Layer& layer_which_contains(const fem::Rectangle &cell,
                                    const std::vector<fem::Point> &points) const;

  unsigned int n_layers() const;


private:
  double _angle;
  unsigned int _n_layers;
  fem::Point _min_point; // limits of the block
  fem::Point _max_point;

  std::vector<double> _layers_thickness;
  std::vector<double> _layers_coef_alpha;
  std::vector<double> _layers_coef_beta;

  std::vector<Layer> _layers;
};




#endif // BLOCK_OF_LAYERS_H
