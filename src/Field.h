#ifndef FIELD_H
#define FIELD_H

#include <iostream>
#include <SFML/Graphics.hpp>
#include "evector.h"
#include "matrix.h"

using namespace sf;
using namespace std;

class Field : public Drawable
{
	public:
		//ctor & dtor
		Field(const int& rows, const int& columns, const int& height, const int& width, const int hMap[], const int& dsStep = 1);
		virtual ~Field();
		//operator
		friend ostream& operator << (ostream& os, const Field &F);

		//public method
		void setHeightMap(const int hMap[]);
		void setScreenSize(const int& height, const int& width);
		void diamondSquareStep();
		void setObsPoint(const Vector& watchingPoint, const Vector& sightPoint);
		void translation(const double& a, const double& b, const double& c);
		void rotation(const Vector& vect, const double& angle);
		void render();

	protected:
		void draw(RenderTarget &target, RenderStates states) const;
	private:
		//attribute
		Vector* surf;
		Vertex* vertices;
		int x;
		int y;
		int screenHeight;
		int screenWidth;
		int fieldHeight;
		int fieldWidth;
		Vector2u dimensions;
		Matrix4 trans;

		//method
		void setColors();
		void setVerticesPos();
		void perspective();
};

#endif // FIELD_H
