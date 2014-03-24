#include "layer.h"
#include "fem/rectangle.h"
#include "fem/math_functions.h"
#include "fem/auxiliary_functions.h"


//Layer::Layer() { }

//Layer::Layer(const Layer &layer)
//  : Quadrangle(layer),
//    _number(layer._number)
//{ }

//Layer& Layer::operator =(const Layer &layer)
//{
//  Quadrangle::operator =(layer);
//  _number = layer._number;
//}

//Layer::~Layer() { }



void Layer::init(unsigned int number, const std::vector<double> &thickness_percent,
                 const fem::Point &min_point, const fem::Point &max_point,
                 double angle)
{
  const double right_angle = 90.; // right angle in degrees
  expect(fabs(angle) < right_angle, "Angle doesn't belong to correct range: (-90, 90), its value : " + d2s(angle));
  expect(fabs(fabs(angle) - right_angle) > fem::math::FLOAT_NUMBERS_EQUALITY_TOLERANCE,
         "Angle is equal to right angle (90), what is prohibited");

  const unsigned int n_layers = thickness_percent.size(); // the total amount of layers
  expect(number < n_layers, "The number of current layer (" + d2s(number + 1) +
         ") is bigger than the total number of layers (" + d2s(n_layers) + ")");

  _number = number;
  _thickness_percent = thickness_percent[number];
  _min_point = min_point;
  _max_point = max_point;

  _angle = angle; // the angle between horizontal (OX) axis
  _angle_rad_abs = fabs(_angle * fem::math::PI / 180.); // abs angle in radians

  // the numeration of vertices is the same as in case of rectangle
  // 2 --- 3   Y
  // |     |   |
  // 0 --- 1   ----> X
  // and since the layers are distributed (nearly) horizontally
  // the x-components of the (0-th and 2-nd) and (1-st and 3-rd) points coincide
  // and they coincide to x-limits of the domain
  _X[0] = _X[2] = _min_point.coord(0);
  _X[1] = _X[3] = _max_point.coord(0);

  const double Hy = _max_point.coord(1) - _min_point.coord(1); // the length (y-direction) of the domain
  const double Hx = _max_point.coord(0) - _min_point.coord(0); // the length (x-direction) of the domain

  // define the limits that the layer occupies when it is distributed horizontally
  double y_bottom = _min_point.coord(1) - 1; // it's wrong by default
  double y_top = _min_point.coord(1);
  for (unsigned int i = 0; i < _number + 1; ++i)
  {
    y_bottom = y_top;
    if (i == n_layers - 1) // if it is the last layer, we don't calculate the limit to eliminate the error
      y_top = _max_point.coord(1);
    else
      y_top = y_bottom + 0.01*thickness_percent[i]*Hy;
  }
  expect(y_bottom >= _min_point.coord(1), "Please check y_bottom. It's wrong.");
  expect(y_top > y_bottom, "y_top and/or y_bottom are wrong");

  _thickness = y_top - y_bottom; // real thickness of an unslopped layer

  // now we define y-coordinates of the real layer vertices (taking the angle into account)
  if (fabs(_angle) < fem::math::FLOAT_NUMBERS_EQUALITY_TOLERANCE) // if angle is zero, there is no stretching, and transformation is not required
  {
    _Y[0] = _Y[1] = y_bottom;
    _Y[2] = _Y[3] = y_top;
  }
  else // there is a slope
  {
    // just for brevity
    const double x0 = _min_point.coord(0);
    const double x1 = _max_point.coord(0);
    const double y0 = _min_point.coord(1);

    if (_angle > 0) // if the angle is positive
    {
      // we find the x-coordinates corresponding to y-coordinates of the layer
      // using a diagonal line going through the left top corner of the domain with layers
      // to the bottom right corner
      //         ************
      //         | \        |
      //         |  \       |
      // y_top   -----      |
      //         |    \     |
      //         |     \    |
      // y_bottom--------   |
      //         |       \  |
      //         |        \ |
      //         ************
      //             x_left
      //                x_right
      //
      double x_right = (Hx/Hy) * (-y_bottom + y0) + x1;
      double x_left  = (Hx/Hy) * (-y_top    + y0) + x1;

      const double tg = tan(_angle_rad_abs); // for brevity
      _Y[0] = y_bottom - (x_right - x0) * tg;
      _Y[1] = y_bottom + (x1 - x_right) * tg;
      _Y[2] = y_top - (x_left - x0) * tg;
      _Y[3] = y_top + (x1 - x_left) * tg;
    }
    else // if the angle is negative
    {
      // we find the x-coordinates corresponding to y-coordinates of the layer
      // using a diagonal line going through the bottom left corner of the domain with layers
      // to the top right corner
      //         ***********
      //         |       / |
      //         |      /  |
      // y_top   -------   |
      //         |    /    |
      //         |   /     |
      // y_bottom----      |
      //         | /       |
      //         |/        |
      //         ************
      //            x_left
      //               x_right
      //
      double x_right = (Hx/Hy) * (y_top    - y0) + x0;
      double x_left  = (Hx/Hy) * (y_bottom - y0) + x0;

      const double tg = tan(_angle_rad_abs); // for brevity
      _Y[0] = y_bottom + (x_left - x0) * tg;
      _Y[1] = y_bottom - (x1 - x_left) * tg;
      _Y[2] = y_top + (x_right - x0) * tg;
      _Y[3] = y_top - (x1 - x_right) * tg;
    }
  }

  // now when we know all real coordinates of the layer's vertices,
  // we can calculate the coefficients of equations describing
  // layer's sloped sides (remember that it has 2 sides parallel to OY axis).
  // bottom side connects 0 and 1 vertices
  _a_bottom = (_Y[1] - _Y[0]) / Hx;
  _b_bottom = _Y[0] - _a_bottom * _X[0];
  // top side connects 2 and 3 vertices
  _a_top = (_Y[3] - _Y[2]) / Hx;
  _b_top = _Y[2] - _a_top * _X[2];
}



bool Layer::contains_element(const fem::Rectangle &cell, const std::vector<fem::Point> &points) const
{
  // we think that a cell belongs to a layer if a center of the cell belongs to the layer
  double xc = 0., yc = 0.; // center of the cell

  for (unsigned int i = 0; i < cell.n_vertices; ++i)
  {
    xc += points[cell.vertex(i)].coord(0);
    yc += points[cell.vertex(i)].coord(1);
  }
  xc /= cell.n_vertices;
  yc /= cell.n_vertices;

  expect(xc > _min_point.coord(0), "X-center of the cell is less than left limit. That's wrong");
  expect(xc < _max_point.coord(0), "X-center of the cell is more than right limit. That's wrong");
//  expect(yc > _min_point.coord(1), "Y-center of the cell is less than bottom limit. That's wrong");
//  expect(yc < _max_point.coord(1), "Y-center of the cell is more than top limit. That's wrong");

  // calculate the intersection of y=xc line with bottom and top slope side of the layer
  const double y_bottom = _a_bottom * xc + _b_bottom;
  const double y_top = _a_top * xc + _b_top;

  // the cell belong to the layer if cell's center is within received bounds,
  // and it also has to be inside domain with layers
  if (yc >= _min_point.coord(1) && yc <= _max_point.coord(1) &&
      yc >= y_bottom && yc <= y_top)
    return true;

  return false; // the layer doesn't contain this cell
}



double Layer::thickness() const
{
  return _thickness;
}
