#ifndef INCLUDED_GEOMETRIC_SUNA
#define INCLUDED_GEOMETRIC_SUNA


float xlinear(float x1,float y1,float x2,float y2,float y)
{
  return ((x1-x2)/(y1-y2)*(y-y1)+x1);
}

float ylinear(float x1,float y1,float x2,float y2,float x)
{
  return (((y1-y2)/(x1-x2))*(x-x1)+y1);
}

float Trapezoid(float a,float b,float c)
{
  return(b+a)*(c)*0.5;
}

#endif // GEOMETRIC_H
