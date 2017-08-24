//
// Created by dodo on 23/08/17.
//

#ifndef NUMERICAL_ANALYSIS_VOLUMOUSOBJECT_HPP
#define NUMERICAL_ANALYSIS_VOLUMOUSOBJECT_HPP

#include <string>

class Point3D {
 private:
  double x, y, z;
 public:
  Point3D() {
    x = y = z = 0;
  }

  Point3D(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  double getX() const {
    return x;
  }

  void setX(double x) {
    this->x = x;
  }

  double getY() const {
    return y;
  }

  void setY(double y) {
    this->y = y;
  }

  double getZ() const {
    return z;
  }

  void setZ(double z) {
    this->z = z;
  }

  std::string toString() {
    return "\tX: " + std::to_string(x) + "\n\tY: " + std::to_string(y) + "\n\tZ: " + std::to_string(z);
  }
};

class VolumousObject {
 private:
  double volume, weight;
  Point3D centerOfMass;
 public:

  VolumousObject() {
    volume = weight = 0;
  }

  double getVolume() const {
    return volume;
  }

  void setVolume(double volume) {
    VolumousObject::volume = volume;
  }

  double getWeight() const {
    return weight;
  }

  void setWeight(double weight) {
    VolumousObject::weight = weight;
  }

  Point3D &getCenterOfMass() {
    return centerOfMass;
  }

  void setCenterOfMass(const Point3D &centerOfMass) {
    VolumousObject::centerOfMass = centerOfMass;
  }

  std::string toString() {
    return "Object details:\n\tVolume: " + std::to_string(volume) + "\n\tWeight:" + std::to_string(weight)
        + "\nCenter of mass:\n" + centerOfMass.toString();
  }
};

#endif //NUMERICAL_ANALYSIS_VOLUMOUSOBJECT_HPP
