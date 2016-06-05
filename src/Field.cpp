#include "Field.h"
#include <time.h>

//-----------------------------------------------------------------------------------------Ctor & Dtor
Field::Field(const int& rows, const int& columns, const int& width, const int& height, const int hMap[], const int& dsStep)
{
	x = columns;
	y = rows;

	fieldHeight = height;
	fieldWidth = width;

	surf = new Vector[height*width];
	vertices = new Vertex[height*width];

	for(int i = 0; i < y; i++)
	{
		for (int j = 0; j < x; j++)
		{
			surf[i*x+j] = Vector(((double)i/(double)(x-1))*width, ((double)j/(double)(y-1))*height, hMap[i*x+j]);
		}
	}

	//init rand
	srand(time(NULL));

	int i;
	for(i = 0; i < dsStep; i++)
	{
		diamondSquareStep();
//		cout<<"Diamond-Square n"<<i+1<<endl;
	}
//	cout<<"Field generated!"<<endl;

	setColors();

	trans = Matrix4::Identity;
}

Field::~Field()
{
	delete[] vertices;
	delete[] surf;
}


//-----------------------------------------------------------------------------------------Operator
ostream& operator<<(ostream& os, const Field &F)
{
	int i,j;
	for(i = 0; i < F.x; i++)
	{
		for (j = 0; j < F.y; j++)
		{
			os<<F.surf[i*F.x + j]<<endl;;
		}
		os<<endl;
	}

	return os;
}


//-----------------------------------------------------------------------------------------Public methods
void Field::setScreenSize(const int& height, const int& width)
{
	screenHeight = height;
	screenWidth = width;
	setVerticesPos();
}

//a little randomized
void Field::diamondSquareStep()
{
	//new field size definition
	int newx = 2*x - 1;
	int newy = 2*y - 1;
	Vector* nRes = new Vector[newx*newy];

	//Copying unchanged points
//	cout<<"Copying unchanged points"<<endl;
	int i,j;
	for(i = 0; i < newx; i+=2)
	{
		for(j = 0; j < newy; j+=2)
		{
			nRes[i*newx + j] = surf[(i/2)*x + j/2];
		}
	}
//	cout<<"OK"<<endl;

	//diamond
//	cout<<"Diamond pass"<<endl;
	Vector average(0.0);
	for(i = 1; i < newx; i+=2)
	{
		for(j = 1; j < newx; j+=2)
		{
			average = nRes[(i-1)*newx + (j-1)] + nRes[(i-1)*newx + (j+1)] + nRes[(i+1)*newx + (j-1)] + nRes[(i+1)*newx + (j+1)];
			average+=Vector(0,0,rand() % 100 - 10);
			average/=4;
			nRes[i*newx + j] = average;
		}
	}
//	cout<<"OK"<<endl;
	//square
//	cout<<"Square pass"<<endl;
	int neighCount;
	for(i = 0; i < newx; i++)
	{
		for(j = (i+1)%2; j < newy; j+=2)
		{
			average = Vector(0.0);
			neighCount = 0;

			//first row has no top neighbour
			if(i != 0)
			{
				average += nRes[(i-1)*newx + j];
				neighCount++;
			}
			//last row has no bottom neighbour
			if(i != newx-1)
			{
				average += nRes[(i+1)*newx + j];
				neighCount++;
			}
			//first column has no left neighbour
			if(j != 0)
			{
				average += nRes[i*newx + (j-1)];
				neighCount++;
			}
			//last column has not right neighbour
			if(j != newy-1)
			{
				average += nRes[i*newx + (j+1)];
				neighCount++;
			}

			average+=Vector(0,0,rand() % 20 - 10);
			average/=neighCount;

			if(i == 0)
			{
				average = Vector(0, average[1], average[2]);
			}
			else if (i == newx - 1)
			{
				average = Vector(fieldWidth, average[1], average[2]);
			}
			if (j == 0)
			{
				average = Vector(average[0], 0, average[2]);
			}
			else if (j == newy - 1)
			{
				average = Vector(average[0], fieldHeight, average[2]);
			}

			nRes[i*newx + j] = average;
		}
	}
//	cout<<"OK"<<endl;

	//Saving new datas
//	cout<<"Saving new datas"<<endl;
	x = newx;
	y = newy;
	delete[] vertices;
	vertices = new Vertex[x*y];
	delete[] surf;
	surf = nRes;
//	cout<<"OK"<<endl;
}

void Field::draw(RenderTarget &target, RenderStates states) const
{
	int i,j;
	Vertex v[3];
	//half the trinagles
	for(i = 0; i < x-1; i++)
	{
		for(j = 0; j < y-1; j++)
		{
			v[0] = vertices[i*x + j];
			v[1] = vertices[i*x + j+1];
			v[2] = vertices[(i+1)*x + j+1];
			target.draw(v, 3, PrimitiveType::Triangles);
			//This is sufficient to draw the field,

			//This outlines each drawn trianlge (debug purpose)
			v[0].color = Color::Black;
			v[1].color = Color::Black;
			v[2].color = Color::Black;
			Vertex l1[2] = {v[0], v[1]};
			Vertex l2[2] = {v[0], v[2]};
			Vertex l3[2] = {v[1], v[2]};
			target.draw(l1, 2, PrimitiveType::Lines);
			target.draw(l2, 2, PrimitiveType::Lines);
			target.draw(l3, 2, PrimitiveType::Lines);
		}
	}
	//the other half
	for(i = x-1; i > 0; i--)
	{
		for(j = y-1; j > 0; j--)
		{
			v[0] = vertices[i*x + j];
			v[1] = vertices[i*x + j-1];
			v[2] = vertices[(i-1)*x + j-1];
			target.draw(v, 3, PrimitiveType::Triangles);
		}
	}
}

void Field::setObsPoint(const Vector& watchingPoint, const Vector& sightPoint)
{
	Vector diff = sightPoint - watchingPoint;
	translation(diff[0], diff[1], diff[2]);

	Normalize(diff);
	Vector normal = Orthogonal(diff);

	double crossProduct = Cosine(Vector(0,0,1), diff);
	double angle = acos(crossProduct);

	rotation(normal, angle);
}

void Field::translation(const double& a, const double& b, const double& c)
{
	Matrix4 mat = Matrix4::Translate(Vector(a, b, c));

	trans = trans*mat;
}

void Field::rotation(const Vector& axis, const double& angle)
{
	Matrix4 mat = Matrix4::Rotate(axis, angle);

	trans = trans*mat;
}

void Field::render()
{

	for(int i = 0; i < x*y; i++)
	{
		//Perspective
		Matrix4 mat (1,0,0,0,
					 0,1,0,0,
					 0,0,1,0,
					 0,0,1/Norm(surf[i]),1);

		Matrix4 tmp = trans * mat;

		surf[i] = tmp * surf[i];
	}
	setVerticesPos();
}
//-----------------------------------------------------------------------------------------Private methods
void Field::setColors()
{
	int i;
	for(i = 0; i < x*y; i++)
	{
		if(surf[i][2] < 800)
			vertices[i].color = Color(242,121,0);
		else
			vertices[i].color = Color::White;
		if(surf[i][2] < 500)
			vertices[i].color = Color::Green;
		if(surf[i][2] < 300)
		{
			vertices[i].color = Color::Blue;
			surf[i].setz(300);
		}

	}
}

void Field::setVerticesPos()
{
	int i,j;
	for(i = 0; i < x; i++)
	{
		for(j = 0; j < y; j++)
		{
			vertices[i*x + j].position = Vector2f(surf[i*x+j][0], surf[i*x+j][1]);
		}
	}
}

void Field::perspective()
{

}
