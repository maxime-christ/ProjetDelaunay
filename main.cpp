#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <iostream>
#include "Field.h"
#include "matrix.h"

using namespace sf;
using namespace std;

int main()
{
    RenderWindow window (VideoMode(1200, 768), "SFML works!");

	int hMap[25] = {700, 700, 700, 600, 600,
					800, 1000, 750, 700, 650,
					750, 800, 700, 450, 400,
					600, 700, 500, 200, 300,
					900, 850, 250, 250, 275};

    Field field(5,5, 1200, 768, hMap, 3);
    field.setScreenSize(1200, 768);

	field.translation(150, 150, 0);
	field.rotation(Vector(1,0,0), -45);
	field.translation(0,750,0);

//	field.setObsPoint(Vector(500, 500, 1500), Vector(0,0,700));
	field.render();

	window.draw(field);

    window.display();

    while (window.isOpen())
    {
        Event event;
        while (window.pollEvent(event))
        {
			switch (event.type)
			{
				case Event::Closed :
					window.close();
					break;

				case Event::KeyPressed :
					if(event.key.code == Keyboard::Escape)
						window.close();
					break;

				case Event::MouseButtonPressed :

					break;

				default :
					break;
			}
        }
    }

//    free(vertices);

    return 0;
}
