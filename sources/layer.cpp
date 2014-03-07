#include "layer.h"
#include "fem/rectangle.h"
#include "fem/math_functions.h"
#include "fem/auxiliary_functions.h"


Layer::Layer() { }

Layer::Layer(const Layer &layer)
  : Quadrangle(layer),
    _number(layer._number)
{ }

Layer& Layer::operator =(const Layer &layer)
{
  Quadrangle::operator =(layer);
  _number = layer._number;
}

Layer::~Layer() { }



void Layer::init(unsigned int number, const std::vector<double> &thickness_percent,
                 const Point &min_point, const Point &max_point,
                 double angle)
{
  const double right_angle = 90.; // right angle in degrees
  require(fabs(angle) < right_angle, "Angle doesn't belong to correct range: (-90, 90), its value : " + d2s(angle));
  require(fabs(fabs(angle) - right_angle) > FLOAT_NUMBERS_EQUALITY_TOLERANCE,
          "Angle is equal to right angle (90), what is prohibited");

  _number = number;
  _thickness_percent = thickness_percent[number];
  _min_point = min_point;
  _max_point = max_point;

  _angle = angle; // the angle between horizontal (OX) axis
  _angle_rad_abs = fabs(_angle * PI / 180.); // abs angle in radians

  // the numberation of vertices is the same as in case of rectangle
  // 2 --- 3
  // |     |
  // 0 --- 1
  // and since the layers are distributed (nearly) horizontally
  // the x-components of the (0-th and 2-nd) and (1-st and 3-rd) points coincide
  // and they coincide to x-limits of the domain
  _X[0] = _X[2] = _min_point.coord(0);
  _X[1] = _X[3] = _max_point.coord(0);

  // define y-coordinates of the layer vertices
  const double Hy = _max_point.coord(1) - _min_point.coord(1); // the length (y-direction) of the domain
  const double Hx = _max_point.coord(0) - _min_point.coord(0); // the length (x-direction) of the domain
  const unsigned int n_layers = thickness_percent.size(); // the total amount of layers
  double y_prev, y_cur = _min_point.coord(1);
  for (int i = 0; i < _number + 1; ++i)
  {
    y_prev = y_cur;
    y_cur = (i == n_layers - 1 ? _max_point.coord(1) : y_prev + 0.01*thickness_percent[i]*Hy);
  }
  //_y_prev = y_prev;
  //_y_cur  = y_cur;

  if (fabs(_angle) < FLOAT_NUMBERS_EQUALITY_TOLERANCE) // if angle is zero, there is no stretching, and transformation is not required
  {
    _Y[0] = _Y[1] = y_prev;
    _Y[2] = _Y[3] = y_cur;

   // _x_pos_left = _x_pos_right = 0.; // they don't matter in this case
    //for (int i = 0; i < n_vertices; ++i)
    //  _transform_y[i] = 1.;
  }
  else // there is an angle
  {
    const double x0 = _min_point.coord(0);
    const double x1 = _max_point.coord(0);
    const double y0 = _min_point.coord(1);

    if (_angle > 0) // if the angle is positive
    {
      double x_prev = (Hx/Hy) * (-y_prev + y0) + x1;
      double x_cur  = (Hx/Hy) * (-y_cur  + y0) + x1;

      const double tg = tan(_angle_rad_abs);
      _Y[0] = y_prev - (x_prev - x0) * tg;
      _Y[1] = y_prev + (x1 - x_prev) * tg;
      _Y[2] = y_cur - (x_cur - x0) * tg;
      _Y[3] = y_cur + (x1 - x_cur) * tg;

      //_transform_y[0] = (fabs(y_prev - _Y[0]) < FLOAT_NUMBERS_EQUALITY_TOLERANCE ? 1. : y_prev / _Y[0]);
      //_transform_y[1] = (fabs(y_prev - _Y[1]) < FLOAT_NUMBERS_EQUALITY_TOLERANCE ? 1. : y_prev / _Y[1]);
      //_transform_y[2] = (fabs(y_cur  - _Y[2]) < FLOAT_NUMBERS_EQUALITY_TOLERANCE ? 1. : y_cur  / _Y[2]);
      //_transform_y[3] = (fabs(y_cur  - _Y[3]) < FLOAT_NUMBERS_EQUALITY_TOLERANCE ? 1. : y_cur  / _Y[3]);

      //_x_pos_left = x_cur;
      //_x_pos_right = x_prev;

      _ab = (_Y[1] - _Y[0]) / Hx;
      _at = (_Y[3] - _Y[2]) / Hx;
      _bb = _Y[0] - _ab * _X[0];
      _bt = _Y[2] - _at * _X[2];

    }
    else // if the angle is negative
    {
      require(false, "Not implemented");
    }

    //calc_transformation();
  }


//  double prev_layers_thickness = 0.; // total thickness of all previous layers
//  for (int i = 0; i < _number; ++i)
//    prev_layers_thickness += thicknesses[i];
//  const double beg_point_y = _min_point.coord(1) + prev_layers_thickness;
//  double beg_point_x = _min_point.coord(0); // in case of negative or 0 angle
//  if (_angle > 0)
//    beg_point_x = _max_point.coord(0); // in case of positive angle

//  _beg_point = Point(beg_point_x, beg_point_y); // beginning point

}



//void Layer::calc_transformation()
//{
//  // the length of the domain. the layer is supposed to lye over the whole
//  // (but in case of non-zero angle for some layers it is not true)
//  const double Hx = _max_point.coord(0) - _min_point.coord(0);

//  // the stretching in y-direction leads to this extra length in respect to normal height (depth)
//  const double dy = Hx * tan(fabs(_angle_rad));


//  // if dy is not zero, then angle is not zero, and there is a stretching, so we need a transformation

//  // the transformation matrix is different when the angle is either positive or negative
//  if (_angle > 0)
//  {
//    _transform_y[0] = _beg_point.coord(1) / (_beg_point.coord(1) - dy);
//    _transform_y[1] = 1.;
//    _transform_y[2] = 1.;
//    _transform_y[3] = (_beg_point.coord(1) + _thickness) / (_beg_point.coord(1) + _thickness + dy);
//  }
//  else // if angle is negative or 0
//  {
//    _transform_y[0] = 1.;
//    _transform_y[1] = _beg_point.coord(1) / (_beg_point.coord(1) - dy);
//    _transform_y[2] = (_beg_point.coord(1) + _thickness) / (_beg_point.coord(1) + _thickness + dy);
//    _transform_y[3] = 1.;
//  }
//}



bool Layer::contains_element(const Rectangle &cell, const std::vector<Point> &points) const
{
  // we think that a cell belongs to a layer if a center of the cell belongs to the layer
  double xc = 0., yc = 0.; // center of the cell

  //const double tg = tan(_angle_rad_abs); // precomputed tangent

  //double x[n_vertices], y[n_vertices];
  //double ynew[n_vertices];

  for (int i = 0; i < cell.n_vertices; ++i)
  {
    const double x = points[cell.vertex(i)].coord(0); // coordinates of the vertex
    const double y = points[cell.vertex(i)].coord(1);

    xc += x;
    yc += y;
  }

//  x[0] = x[2] = 0;
//  x[1] = x[3] = 1;
//  y[0] = 0;
//  y[1] = 1;
//  y[2] = 1;
//  y[3] = 2;


//  ynew[0] = y[0] - x[0]*tg;//(_x_pos_right - x[0]) * tg;
//  ynew[1] = y[1] - x[1]*tg;//(_x_pos_right - x[1]) * tg;
//  ynew[2] = y[2] - x[2]*tg;//(_x_pos_left  - x[2]) * tg;
//  ynew[3] = y[3] - x[3]*tg;//(_x_pos_left  - x[3]) * tg;

//  for (int i = 0; i < cell.n_vertices; ++i)
//  {
//    xc += x[i];
//    yc += ynew[i];
//  }
  xc /= cell.n_vertices;
  yc /= cell.n_vertices;

//  expect(xc > _min_point.coord(0), "X-center of the cell is less than left limit. That's wrong");
//  expect(xc < _max_point.coord(0), "X-center of the cell is more than right limit. That's wrong");
//  expect(yc > _min_point.coord(1), "Y-center of the cell is less than bottom limit. That's wrong");
//  expect(yc < _max_point.coord(1), "Y-center of the cell is more than top limit. That's wrong");

//  if (yc >= _y_prev &&
//      yc <= _y_cur)
//    return true; // the layer contains this cell


  const double yb = _ab * xc + _bb;
  const double yt = _at * xc + _bt;

  if (yc >= yb && yc <= yt)
    return true;

  return false; // the layer doesn't contain this cell
}
