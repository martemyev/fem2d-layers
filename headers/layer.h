#ifndef LAYER_H
#define LAYER_H

#include "fem/quadrangle.h"
#include "fem/point.h"
#include "fem/rectangle.h"

//class fem::Point;
//class fem::Rectangle;


/**
 * A layer is a parallelogram, which is in turn a quadrangle.
 * In case of horizontal distribution a layer is a rectangle.
 * These layers are supposed to be nearly horizontally,
 * which means that the angle Pi/2 is prohibited -
 * so we can not get vertical aligned layers.
 */
class Layer : public fem::Quadrangle
{
public:
//  Layer();
//  Layer(const Layer &layer);
//  Layer& operator =(const Layer &layer);
//  virtual ~Layer();

            /**
             * Initialization of the (sloped) layer
             * @param number - the number of the layer
             * @param thickness_percent - the vector of all thicknesses of all layers in the domain bounded by 2 points
             * @param min_point - the min point of the domain where the layers are supposed to be
             * @param max_point - the max point of the domain where the layers are supposed to be
             * @param angle - the slope angle in degrees (in respect of OX axis - so the angle 0 means horizontal layer),
             *                the possible range of the angle values is (-90, 90), not including right angles
             */
  void init(unsigned int number, const std::vector<double> &thickness_percent,
            const fem::Point &min_point, const fem::Point &max_point,
            double angle);

            /**
             * Check if this layer contains a rectangle
             * @param rect - a rectangle that we want to check
             */
  bool contains_element(const fem::Rectangle &cell, const std::vector<fem::Point> &points) const;

  double thickness() const;

private:
  unsigned int _number;
  double _thickness_percent;
  double _thickness;
  double _angle;
  double _angle_rad_abs;
  fem::Point _min_point;
  fem::Point _max_point;

            /**
             * The layer has 2 sides parallel to OY axis, and 2 sides with possible slope.
             * These sloped sides are described by 2 equations:
             * y = _a_bottom * x + _b_bottom,
             * y = _a_top * x + _b_top,
             * Both equations are used to determine whether some point belongs to the layer or not.
             */
  double _a_bottom, _a_top;
  double _b_bottom, _b_top;

//  //double _x_pos_right, _x_pos_left; // some valuable xs

//  //double _y_prev, _y_cur;

//            /**
//             * Transformation matrix for y-coordinates of the points.
//             * We don't consider transformation of x-coordinates,
//             * because the layers are distributed (nearly) horizontally,
//             * and therefore the x-coordinates for the stretched rectangle (parallelogram = layer)
//             * stay the same as for the straight rectangle
//             */
//  //double _transform_y[n_vertices];


//            /**
//             * The beginning point of the layer.
//             * The choice of this point is based not only on the number of the layer,
//             * but also on the angle, because in case of positive or negative angle
//             * the x-coordinate of the beginning point will be different
//             */
//  //Point _beg_point;

//  //void calc_transformation();
};

#endif // LAYER_H
