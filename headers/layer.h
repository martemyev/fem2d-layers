#ifndef LAYER_H
#define LAYER_H

#include "fem/quadrangle.h"


class Point;
class Rectangle;


/**
 * A layer is a parallelogram, which is in turn a quadrangle.
 * In case of horizontal distribution a layer is a rectangle.
 */
class Layer : public Quadrangle
{
public:
  Layer();
  Layer(const Layer &layer);
  Layer& operator =(const Layer &layer);
  virtual ~Layer();

  void init(unsigned int number, const std::vector<double> &thickness_percent,
            const Point &min_point, const Point &max_point,
            double angle);

            /**
             * Check if this layer contains a rectangle
             * @param rect - a rectangle that we want to check
             */
  bool contains_element(const Rectangle &cell, const std::vector<Point> &points) const;

private:
  unsigned int _number;
  //double _thickness;
  double _thickness_percent;
  double _angle;
  double _angle_rad_abs;
  Point _min_point;
  Point _max_point;

  double _ab, _at;
  double _bb, _bt;

  //double _x_pos_right, _x_pos_left; // some valuable xs

  //double _y_prev, _y_cur;

            /**
             * Transformation matrix for y-coordinates of the points.
             * We don't consider transformation of x-coordinates,
             * because the layers are distributed (nearly) horizontally,
             * and therefore the x-coordinates for the stretched rectangle (parallelogram = layer)
             * stay the same as for the straight rectangle
             */
  //double _transform_y[n_vertices];


            /**
             * The beginning point of the layer.
             * The choice of this point is based not only on the number of the layer,
             * but also on the angle, because in case of positive or negative angle
             * the x-coordinate of the beginning point will be different
             */
  //Point _beg_point;

  //void calc_transformation();
};

#endif // LAYER_H
